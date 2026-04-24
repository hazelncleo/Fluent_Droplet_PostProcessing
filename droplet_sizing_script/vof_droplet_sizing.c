#include "vof_droplet_sizing.h"


/*
 *   UDF Functions
 */


DEFINE_EXECUTE_ON_LOADING(on_loading, libname) {

    if (udm_offset == UDM_UNRESERVED) {
        udm_offset = Reserve_User_Memory_Vars(NUM_UDM);
    }

    if (udm_offset == UDM_UNRESERVED) {
        Message("\nThe library: \"%s\" requires %d UDMs to be defined before loading.\n", libname, NUM_UDM);
    } else {
        #if !RP_NODE
        Message("\n%d UDMs have been reserved for the library: \"%s\", with offset: \"%d\".\n", NUM_UDM, libname, udm_offset);
        Set_User_Memory_Name(udm_offset, "droplet-id");
        #endif
    }
}


DEFINE_EXECUTE_AT_END(multithreaded_droplet_sizes_runtime) {
    multithreaded_calculation();
}


DEFINE_EXECUTE_AT_END(singlethreaded_droplet_sizes_runtime) {
    singlethreaded_calculation();
}


DEFINE_ON_DEMAND(droplet_sizes_ondemand) {

    if (MULTIPLE_COMPUTE_NODE_P) {
        multithreaded_calculation();
    } else {
        singlethreaded_calculation();
    }
}


/*
 *   Main control functions
 */


void multithreaded_calculation() {

    if (N_TIME > SKIP_TIMESTEP){

        #if !RP_HOST

            cell_stack cells_to_reexplore;
            initialize_stack(&cells_to_reexplore);

            init_udm();

            compute_droplet_data(&cells_to_reexplore);

            if I_AM_NODE_ZERO_P {
                node_zero_send_data();
            }



            assemble_droplets(&cells_to_reexplore);

            if I_AM_NODE_ZERO_P {
                node_zero_send_droplet_connections();
            }

        #endif

        #if !RP_NODE

            double time_spent;

            clock_t begin, end;
            begin = clock();

            Message("\n#################################################################\n");
            host_process_multithreaded();

            end = clock();
            time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

            Message("Completed successfully in %.3f seconds\n", time_spent);
            Message("\n#################################################################\n");
        #endif
    }
}


void singlethreaded_calculation() {

    if (N_TIME > SKIP_TIMESTEP){

        #if !RP_HOST

            init_udm();

            compute_droplet_data_singlethreaded();

        #endif

        #if !RP_NODE

            double time_spent;

            clock_t begin, end;
            begin = clock();

            Message("\n#################################################################\n");
            host_process_singlethreaded();

            end = clock();
            time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

            Message("Completed successfully in %.3f seconds\n", time_spent);
            Message("\n#################################################################\n");
        #endif
    }
}


/*
 *   Compute node functions
 */


/* Initialize the values of the UDM to 0 */
void init_udm() {

    Domain *mixture_domain = Get_Domain(1);

    Thread *cell_thread;
    Thread **pt;

    cell_t cell;

    if (udm_offset != UDM_UNRESERVED) {

        mp_thread_loop_c(cell_thread, mixture_domain, pt) {

            if (FLUID_THREAD_P(cell_thread)) {

                begin_c_loop_int(cell, cell_thread) {

                    /* Set UDM to 0 */
                    C_UDMI(cell, cell_thread, udm_offset) = 0;

                } end_c_loop_int(cell, cell_thread)
            }
        }

    } else {
        Error("UDM has not been reserved\n");
    }

    PRF_GSYNC();

    EXCHANGE_SVAR_MESSAGE(mixture_domain, (SV_UDM_I, SV_NULL));

}


/* Calculate droplet sizes, positions and velocities, for singlethreaded analysis */
void compute_droplet_data_singlethreaded() {

    Domain *mixture_domain = Get_Domain(1);

    Thread *cell_thread, *phase_thread;
    Thread **pt;

    cell_t cell;

    int cell_explored, message, droplet_id;
    droplet_id = 1;
    message = -1;

    real vof;

    mp_thread_loop_c(cell_thread, mixture_domain, pt) {
        if (FLUID_THREAD_P(cell_thread)) {

            phase_thread = pt[DROPLET_PHASE];

            begin_c_loop_int(cell, cell_thread) {

                vof = C_VOF(cell, phase_thread);
                cell_explored = C_UDMI(cell, cell_thread, udm_offset);

                /* If cell is in new droplet, perform floodfill algorithm */
                if (cell_explored == 0 && vof > 0.5) {

                    found_new_droplet_singlethreaded(cell, cell_thread, phase_thread, droplet_id);
                    ++droplet_id;

                } else if (cell_explored == 0) {

                    /* If cell is not in droplet, set to explored */
                    C_UDMI(cell, cell_thread, udm_offset) = -1;
                }
            }
            end_c_loop_int(cell,cell_thread)
        }
    }

    /* Send finished computing message to host node */
    PRF_CSEND_INT(node_host, &message, 1, droplet_id);
}


/* Performs floodfill on a droplet, for singlethreaded analysis */
void found_new_droplet_singlethreaded(cell_t first_cell, Thread *cell_thread, Thread *phase_thread, int droplet_id) {

    cell_stack candidates_stack;
    initialize_stack(&candidates_stack);

    Thread *cell_face_thread;

    cell_t current_cell, adjacent_cell;

    face_t cell_face;

    int cell_explored, local_face_id, message;
    message = 0;

    real vof;
    real *droplet_values;
    vof = C_VOF(first_cell, phase_thread);
    droplet_values = (real*)calloc(8, sizeof(real)); /* 0 = vol, 1 = mass, 2-4 = centroid, 5-7 = velocity */

    /* Set the droplet id of the first cell and push it to the stack */
    C_UDMI(first_cell, cell_thread, udm_offset) = droplet_id;

    push_to_stack(&candidates_stack, first_cell);

    first_value_update(droplet_values, vof, first_cell, cell_thread);

    /* Perform floodfill until the stack is empty */
    while (!stack_is_empty(&candidates_stack)) {

        /* Get next cell on the stack */
        current_cell = pop_from_stack(&candidates_stack);

        /* Loop over faces of current cell */
        c_face_loop(current_cell, cell_thread, local_face_id) {

            /* Get global face & face thread values */
            cell_face = C_FACE(current_cell, cell_thread, local_face_id);
            cell_face_thread = C_FACE_THREAD(current_cell, cell_thread, local_face_id);

            /* If the face is not on the boundary */
            if (!BOUNDARY_FACE_THREAD_P(cell_face_thread)) {

                /* Get the adjacent cell */
                adjacent_cell = F_C0(cell_face, cell_face_thread);
                if (adjacent_cell == current_cell) {
                    adjacent_cell = F_C1(cell_face, cell_face_thread);
                }

                /* Get adjacent cells vof & UDM value */
                vof = C_VOF(adjacent_cell, phase_thread);
                cell_explored = C_UDMI(adjacent_cell, cell_thread, udm_offset);


                /* If adjacent cell has not been explored and is in droplet */
                if (cell_explored == 0 && vof > 0.5) {

                    /* Push cell to the stack, set its droplet id and add values to the droplet */
                    push_to_stack(&candidates_stack, adjacent_cell);
                    C_UDMI(adjacent_cell, cell_thread, udm_offset) = droplet_id;
                    subsequent_value_update(droplet_values, vof, adjacent_cell, cell_thread);

                } else if (cell_explored == 0) {

                    /* Cell not in droplet */
                    C_UDMI(adjacent_cell, cell_thread, udm_offset) = -1;
                }
            }
        }
    }

    /* Divide centroid & velocity by mass */
    final_value_update(droplet_values);

    /* Send message to the host & send final droplet values */
    PRF_CSEND_INT(node_host, &message, 1, droplet_id);
    PRF_CSEND_REAL(node_host, droplet_values, 8, droplet_id);
}


/* Calculate droplet sizes, positions and velocities, for multithreaded analysis */
void compute_droplet_data(cell_stack *cells_to_reexplore) {

    Domain *mixture_domain = Get_Domain(1);

    Thread *cell_thread, *phase_thread;
    Thread **pt;

    cell_t cell;

    int cell_explored, message, droplet_id, receiving_node;
    droplet_id = myid + 1;
    message = -1;
    receiving_node = (I_AM_NODE_ZERO_P ? node_host : node_zero); /* If node zero send to host, otherwise send to node zero */

    real vof;

    mp_thread_loop_c(cell_thread, mixture_domain, pt) {
        if (FLUID_THREAD_P(cell_thread)) {

            phase_thread = pt[DROPLET_PHASE];

            begin_c_loop_int(cell, cell_thread) {

                vof = C_VOF(cell, phase_thread);
                cell_explored = C_UDMI(cell, cell_thread, udm_offset);

                /* If cell is in new droplet, perform floodfill algorithm */
                if (cell_explored == 0 && vof > 0.5) {

                    found_new_droplet(cell, cell_thread, phase_thread, droplet_id, cells_to_reexplore);
                    droplet_id += compute_node_count;

                } else if (cell_explored == 0) {

                    /* If cell is not in droplet, set to explored */
                    C_UDMI(cell, cell_thread, udm_offset) = -1;
                }
            }
            end_c_loop_int(cell,cell_thread)
        }
    }

    /* Send finished computing message to receiving node & sync exterior cell UDM values */
    PRF_CSEND_INT(receiving_node, &message, 1, droplet_id);

    PRF_GSYNC();

    EXCHANGE_SVAR_MESSAGE(mixture_domain, (SV_UDM_I, SV_NULL));
}


/* Performs floodfill on a droplet, for multithreaded analysis */
void found_new_droplet(cell_t first_cell, Thread *cell_thread, Thread *phase_thread, int droplet_id, cell_stack *cells_to_reexplore) {

    cell_stack candidates_stack;
    initialize_stack(&candidates_stack);

    Thread *cell_face_thread;

    cell_t current_cell, adjacent_cell;

    face_t cell_face;

    int cell_explored, local_face_id, adjacent_cell_node_id, message, receiving_node, droplet_outside_node;
    receiving_node = (I_AM_NODE_ZERO_P ? node_host : node_zero);
    droplet_outside_node = 0;

    real vof;
    real *droplet_values;
    vof = C_VOF(first_cell, phase_thread);
    droplet_values = (real*)calloc(8, sizeof(real)); /* 0 = vol, 1 = mass, 2-4 = centroid, 5-7 = velocity */

    /* Set the droplet id of the first cell and push it to the stack */
    C_UDMI(first_cell, cell_thread, udm_offset) = droplet_id;
    push_to_stack(&candidates_stack, first_cell);

    first_value_update(droplet_values, vof, first_cell, cell_thread);

    /* Perform floodfill until the stack is empty */
    while (!stack_is_empty(&candidates_stack)) {

        /* Get next cell on the stack */
        current_cell = pop_from_stack(&candidates_stack);

        /* Loop over faces of current cell */
        c_face_loop(current_cell, cell_thread, local_face_id) {

            /* Get global face & face thread values */
            cell_face = C_FACE(current_cell, cell_thread, local_face_id);
            cell_face_thread = C_FACE_THREAD(current_cell, cell_thread, local_face_id);

            /* If the face is not on the boundary */
            if (!BOUNDARY_FACE_THREAD_P(cell_face_thread)) {

                /* Get the adjacent cell */
                adjacent_cell = F_C0(cell_face, cell_face_thread);
                if (adjacent_cell == current_cell) {
                    adjacent_cell = F_C1(cell_face, cell_face_thread);
                }

                /* Get adjacent cells vof & UDM value */
                vof = C_VOF(adjacent_cell, phase_thread);
                cell_explored = C_UDMI(adjacent_cell, cell_thread, udm_offset);

                /* If adjacent cell has not been explored and is in droplet */
                if (cell_explored == 0 && vof > 0.5) {

                    /* Get ID of the adjacent cells node */
                    adjacent_cell_node_id = C_PART(adjacent_cell, cell_thread);

                    /* Check if adjacent cell is in the current node */
                    if (adjacent_cell_node_id != myid) {

                        /* Mark the node for combining and save cell id */
                        droplet_outside_node = 1;
                        push_to_stack(cells_to_reexplore, adjacent_cell);

                    } else {

                        /* Push cell to the stack, set its droplet id and add values to the droplet */
                        push_to_stack(&candidates_stack, adjacent_cell);
                        C_UDMI(adjacent_cell, cell_thread, udm_offset) = droplet_id;
                        subsequent_value_update(droplet_values, vof, adjacent_cell, cell_thread);
                    }

                } else if (cell_explored == 0) {

                    /* Cell not in droplet */
                    C_UDMI(adjacent_cell, cell_thread, udm_offset) = -1;
                }
            }
        }
    }

    /* Check if droplet crosses node boundaries */
    if (!droplet_outside_node) {

        /* Divide centroid & velocity by mass */
        final_value_update(droplet_values);

        /* Send message to the host & send final droplet values */
        message = 0;
        PRF_CSEND_INT(receiving_node, &message, 1, droplet_id);
        PRF_CSEND_REAL(receiving_node, droplet_values, 8, droplet_id);

    } else {

        /* Send message to the host & send final droplet values */
        message = 1;
        PRF_CSEND_INT(receiving_node, &message, 1, droplet_id);
        PRF_CSEND_REAL(receiving_node, droplet_values, 8, droplet_id);

        /* Push the first cell of the droplet to be rexplored & add a gap in the stack */
        push_to_stack(cells_to_reexplore, first_cell);
        add_gap_to_stack(cells_to_reexplore);
    }
}


/* Find droplet ids of cells connected to boundary droplets */
void assemble_droplets(cell_stack *cells_to_reexplore) {

    Domain *mixture_domain = Get_Domain(1);

    Thread *cell_thread;
    Thread **pt;

    cell_t cell;

    int *attached_droplets, n_droplets, new_droplet_id, message, receiving_node;
    n_droplets = -1;
    message = -1;
    attached_droplets = (int*)calloc(MAX_COMBINE_DROPLETS, sizeof(int));
    receiving_node = (I_AM_NODE_ZERO_P ? node_host : node_zero);

    mp_thread_loop_c(cell_thread, mixture_domain, pt) {
        if (FLUID_THREAD_P(cell_thread)) {

            while(!stack_is_empty(cells_to_reexplore)) {

                /* Check if the top of the stack is a gap */
                if (stack_is_gap(cells_to_reexplore)) {

                    /* If number of droplets to be combined in greater than 0 send data */
                    if (n_droplets >= 0) {
                        send_droplet_message(n_droplets, attached_droplets, receiving_node);
                    }

                    /* Remove the gap and start exploring the next droplets connections */
                    rem_gap_from_stack(cells_to_reexplore);
                    cell = pop_from_stack(cells_to_reexplore);
                    attached_droplets[0] = C_UDMI(cell, cell_thread, udm_offset);
                    n_droplets = 1;

                } else {

                    /* Get droplet ID of a neighbouring nodes cell */
                    cell = pop_from_stack(cells_to_reexplore);
                    new_droplet_id = C_UDMI(cell, cell_thread, udm_offset);

                    /* Check if the droplet ID is not already connected */
                    if (droplet_in_array(attached_droplets, n_droplets, new_droplet_id)) {
                        attached_droplets[n_droplets++] = new_droplet_id;
                    }
                }
            }

            /* Send final droplet information */
            if (n_droplets > 0) {
                send_droplet_message(n_droplets, attached_droplets, receiving_node);
            }

            /* Send finished computing message to receiving node */
            PRF_CSEND_INT(receiving_node, &message, 1, myid);
        }
    }
}


/* Check if a droplet id is already in an array of droplet ids */
int droplet_in_array(int *attached_droplets, int n_droplets, int new_droplet_id) {

    for (int droplet = 0; droplet < n_droplets; droplet++) {
        if (attached_droplets[droplet] == new_droplet_id) {
            return 0;
        }
    }
    return 1;
}


/* Send array of droplet ids to the receiving node */
void send_droplet_message(int array_size, int *attached_droplets, int receiving_node) {

    /* Allocate int array and assign values */
    int *droplet_ids = (int*)malloc(array_size * sizeof(int));

    for (int droplet = 0; droplet < array_size; ++droplet) {
        droplet_ids[droplet] = attached_droplets[droplet];
    }

    /* Send size of array and then the array to the receiving node */
    PRF_CSEND_INT(receiving_node, &array_size, 1, myid);
    PRF_CSEND_INT(receiving_node, droplet_ids, array_size, myid);
}


/*
 *   Stack functions
 */


/* Initialize the stack */
void initialize_stack(cell_stack *stack_obj) {
    stack_obj->top = -1;
    stack_obj->n_gaps = 0;
}


/* Check if the stack is empty */
int stack_is_empty(cell_stack *stack_obj) {
    if (stack_obj->top == -1) {
        return 1;
    } else {
        return 0;
    }
}


/* Check if the stack is full */
int stack_is_full(cell_stack *stack_obj) {
    if (stack_obj->top >= MAX_STACK_SIZE - 1) {
        return 1;
    } else {
        return 0;
    }
}


/* Check if a gap is on top of the stack */
int stack_is_gap(cell_stack *stack_obj) {
    if (stack_obj->n_gaps > 0 && stack_obj->gaps[stack_obj->n_gaps-1] == stack_obj->top) {
        return 1;
    } else {
        return 0;
    }
}


/* Push a cell onto the stack */
void push_to_stack(cell_stack *stack_obj, cell_t cell) {
    if (stack_is_full(stack_obj)) {
        Error("The stack on Node %d is full\n", myid);
    }
    stack_obj->arr[++stack_obj->top] = cell;
}


/* Add a gap onto the stack to seperate values */
void add_gap_to_stack(cell_stack *stack_obj) {
    stack_obj->gaps[stack_obj->n_gaps++] = ++stack_obj->top;
}


/* Remove a gap from the top of the stack */
void rem_gap_from_stack(cell_stack *stack_obj) {
    if (stack_is_gap(stack_obj)) {
        stack_obj->n_gaps--;
        stack_obj->top--;
    } else {
        Error("Tried to remove gap from stack on node %d, but top is not a gap\n", myid);
    }
}


/* Pop the top value off the stack */
cell_t pop_from_stack(cell_stack *stack_obj) {

    cell_t popped;

    if (stack_is_empty(stack_obj)) {
        Error("Tried to pop value off an empty stack on node %d\n", myid);
    } else if (stack_is_gap(stack_obj)) {
        Error("Tried to pop value off stack, but was a gap on node %d\n", myid);
    }

    popped = stack_obj->arr[stack_obj->top];
    stack_obj->top--;
    return popped;
}


/*
 *   Data storage functions
 */


/* Initialize the data storage values to 0 */
void initialize_datastorage(datastorage *droplets_datastorage) {

    droplets_datastorage->n_droplets = 0;                                                    /* Total number of droplet values that have been saved */
    droplets_datastorage->n_to_combine = 0;                                                  /* Total number of droplets after combining */
    memset(droplets_datastorage->combination_ids, 0, MAX_COMBINE_DROPLETS * sizeof(int));    /* Array of temporary IDs to help with combining */
}


/* Add a droplets values and ID to the storage object */
void add_droplet_values_to_datastorage(datastorage *droplets_datastorage, real *droplet_values, int droplet_id) {

    /* Save droplet values */
    for (int droplet_value = 0; droplet_value < 8; ++droplet_value) {
        droplets_datastorage->droplet_values[droplets_datastorage->n_droplets][droplet_value] = droplet_values[droplet_value];
    }

    /* Save droplet id & Increment number of droplets saved */
    droplets_datastorage->droplet_ids[droplets_datastorage->n_droplets] = droplet_id;
    droplets_datastorage->n_droplets++;
}


/* Get the index of a given droplet id */
int get_droplet_index_from_datastorage(datastorage *droplets_datastorage, int droplet_id) {

    for (int index = 0; index < droplets_datastorage->n_droplets; ++index) {
        if (droplets_datastorage->droplet_ids[index] == droplet_id) {
            return index;
        }
    }
    Error("Error finding droplet");
    return -1;
}


/*
 *  Check if any of droplets in array "droplets" have already been assigned
 *  If they have then return that id
 *  otherwise return an unassigned id
 */
void check_droplets_already_assigned(datastorage *droplets_datastorage, int n_to_assign, int *droplets, int *reassign_combinations) {

    int index;

    reassign_combinations[0] = 0;
    reassign_combinations[1] = droplets_datastorage->n_to_combine + 1; /* Set to unassigned id */

    for (int current_droplet = 0; current_droplet < n_to_assign; ++current_droplet) {

        index = get_droplet_index_from_datastorage(droplets_datastorage, droplets[current_droplet]); /* Get index of current droplet in datastorage */

        /* If droplet is assigned and its id is less than the current id reassign it */
        if (droplets_datastorage->combination_ids[index] > 0) {
            if (reassign_combinations[1] == droplets_datastorage->n_to_combine + 1) {

                /* If first connected id then reduce */
                reassign_combinations[1] = droplets_datastorage->combination_ids[index];

            } else if (droplets_datastorage->combination_ids[index] < reassign_combinations[1]) {

                /* Add previous combination id to list of ids to reduce */
                reassign_combinations[2 + reassign_combinations[0]] = reassign_combinations[1];
                reassign_combinations[1] = droplets_datastorage->combination_ids[index];
                reassign_combinations[0]++;

            } else if (droplets_datastorage->combination_ids[index] > reassign_combinations[1]) {

                /* Add found combination id to list of ids to reduce */
                reassign_combinations[2 + reassign_combinations[0]] = droplets_datastorage->combination_ids[index];
                reassign_combinations[0]++;

            }
        }
    }

    /* If droplet is not attached to any current droplets increment number of combined droplets */
    if (reassign_combinations[1] == droplets_datastorage->n_to_combine + 1) {
        droplets_datastorage->n_to_combine++;
    }
}

/* Assign each droplet in array a combined id */
/* TODO MIGHT NEED TO CATCH SOME EDGE CASES */
void assign_droplets_to_combine(datastorage *droplets_datastorage, int n_to_assign, int *droplets) {

    int *reassign_combinations;
    reassign_combinations = (int*)malloc(MAX_COMBINE_DROPLETS * sizeof(int));

    check_droplets_already_assigned(droplets_datastorage, n_to_assign, droplets, reassign_combinations);

    for (int droplet = 0; droplet < n_to_assign; ++droplet) {
        droplets_datastorage->combination_ids[get_droplet_index_from_datastorage(droplets_datastorage, droplets[droplet])] = reassign_combinations[1];
    }

    for (int reassign = 2; reassign < (2 + reassign_combinations[0]); ++reassign) {
        for (int combination = 0; combination < droplets_datastorage->n_to_combine; ++combination) {
            if (droplets_datastorage->combination_ids[combination] == reassign_combinations[reassign]) {
                droplets_datastorage->combination_ids[combination] = reassign_combinations[1];
            }
        }
    }
}


/* Combine the values for a combination droplet */
int get_values_from_datastorage(datastorage *droplets_datastorage, real *droplet_values, int combination_id) {

    int droplet_id = 0;

    for (int current_droplet = 0; current_droplet < droplets_datastorage->n_droplets; ++current_droplet) {

        /* If the droplet is part of the combination droplet add its values */
        if (droplets_datastorage->combination_ids[current_droplet] == combination_id) {

            add_vectors(droplet_values, droplets_datastorage->droplet_values[current_droplet]);

            /* Check if new id is less than current id */
            if ((droplet_id == 0) || (droplets_datastorage->droplet_ids[current_droplet] < droplet_id)) {
                droplet_id = droplets_datastorage->droplet_ids[current_droplet];
            }
        }
    }

    /* No droplets of this combination id found */
    if (droplet_id == 0) {
        return -1;
    }

    /* Divide by the mass */
    final_value_update(droplet_values);

    return droplet_id;
}


/* Add one droplet values vector to another */
void add_vectors(real *droplet_values, real *values_to_add) {
    for (int value = 0; value < 8; ++value) {
        droplet_values[value] += values_to_add[value];
    }
}


/* Calculate initial values for first cell in droplet */
void first_value_update(real *droplet_values, real vof, cell_t cell, Thread *cell_thread) {

    real temp[3];
    C_CENTROID(temp, cell, cell_thread);

    droplet_values[0] = C_VOLUME(cell, cell_thread) * vof;      /* Calculate volume */
    droplet_values[1] = droplet_values[0] * SECONDARY_DENSITY;  /* Calculate mass */

    /* Calculate centroid */
    droplet_values[2] = temp[0] * droplet_values[1];
    droplet_values[3] = temp[1] * droplet_values[1];
    droplet_values[4] = temp[2] * droplet_values[1];

    /* Calculate velocity */
    droplet_values[5] = C_U(cell, cell_thread) * droplet_values[1];
    droplet_values[6] = C_V(cell, cell_thread) * droplet_values[1];
    droplet_values[7] = C_W(cell, cell_thread) * droplet_values[1];
}


/* Update the values of a droplet with a new cells values */
void subsequent_value_update(real *droplet_values, real vof, cell_t cell, Thread *cell_thread) {

    real temp_value;
    temp_value = C_VOLUME(cell, cell_thread) * vof; /* Calculate volume of current cell */
    real temp_vector[3];
    C_CENTROID(temp_vector, cell, cell_thread);

    droplet_values[0] += temp_value;                /* Calculate total droplet volume */
    temp_value *= SECONDARY_DENSITY;                /* Calculate current cell mass */
    droplet_values[1] += temp_value;                /* Calculate total droplet mass */

    /* Update centroid position */
    droplet_values[2] += temp_vector[0] * temp_value;
    droplet_values[3] += temp_vector[1] * temp_value;
    droplet_values[4] += temp_vector[2] * temp_value;

    /* Update velocity */
    droplet_values[5] += C_U(cell, cell_thread) * temp_value;
    droplet_values[6] += C_V(cell, cell_thread) * temp_value;
    droplet_values[7] += C_W(cell, cell_thread) * temp_value;
}


/* Calculate final centroid and velocities by dividing by mass */
void final_value_update(real *droplet_values) {

    /* Final centroid values */
    droplet_values[2] /= droplet_values[1];
    droplet_values[3] /= droplet_values[1];
    droplet_values[4] /= droplet_values[1];

    /* Final velocity values */
    droplet_values[5] /= droplet_values[1];
    droplet_values[6] /= droplet_values[1];
    droplet_values[7] /= droplet_values[1];
}


/*
 *   Node zero functions
 */


/* Receive data from compute nodes and send to host */
void node_zero_send_data() {

    int message = 0, droplet_id = 2, all_nodes_completed = 0;
    int* nodes_completed;
    nodes_completed = (int*)calloc((compute_node_count - 1), sizeof(int));

    real droplet_values[8] = {0.};

    while (!all_nodes_completed) {

        all_nodes_completed = 1;

        for (int current_node = 1; current_node < compute_node_count; ++current_node) {
            if (!nodes_completed[current_node - 1]) {

                /*
                 * message = -1 for no more droplets
                 * message =  0 for droplet within single domain
                 * message =  1 for droplet in more domains
                 */
                PRF_CRECV_INT(current_node, &message, 1, droplet_id);
                PRF_CSEND_INT(node_host, &message, 1, droplet_id);

                if (message == -1) {
                    nodes_completed[current_node - 1] = 1;    /* Node calculations completed */
                } else if (message >= 0) {

                    /* Receive values for droplet */
                    PRF_CRECV_REAL(current_node, droplet_values, 8, droplet_id);
                    PRF_CSEND_REAL(node_host, droplet_values, 8, droplet_id);

                    all_nodes_completed = 0;    /* Keep looping */
                }
            }
            droplet_id++;   /* Increment to next node */
        }
        droplet_id++;   /* Skip node zero */
    }
}


/* Receive connection data from compute nodes and send to host */
void node_zero_send_droplet_connections() {

    int message = 0, all_nodes_completed = 0;
    int *nodes_completed;
    nodes_completed = (int*)calloc((compute_node_count - 1), sizeof(int));

    while (!all_nodes_completed) {

        all_nodes_completed = 1;

        for (int current_node = 1; current_node < compute_node_count; ++current_node) {
            if (!nodes_completed[current_node - 1]) {

                /*
                 * message = -1 for no more connections to receive
                 * message >  0 message is equal to size of array to receive
                 */
                PRF_CRECV_INT(current_node, &message, 1, current_node);
                PRF_CSEND_INT(node_host, &message, 1, myid);

                if (message == -1) {
                    nodes_completed[current_node-1] = 1;
                } else {

                    int *droplet_ids = (int*)malloc(message * sizeof(int));

                    /* Receive connections for droplet and send to host */
                    PRF_CRECV_INT(current_node, droplet_ids, message, current_node);
                    PRF_CSEND_INT(node_host, droplet_ids, message, myid);

                    all_nodes_completed = 0;    /* Keep looping */
                }
            }
        }
    }
}


/*
 *   Host node functions
 */


/* Gets file name, creates file and writes header line */
FILE *file_handler() {

    char fname[60];

    FILE *fptr = NULL;

    int max_len, timestep, sim_id, n_cycles, x_grid_position, y_grid_position;
    max_len = sizeof(fname);
    timestep = N_TIME;

    real vibration_amplitude, vibration_frequency, noise_amplitude, noise_frequency;

    /* Set filename to droplets_<timestep_number>.csv */
    if (snprintf(fname, max_len, "droplets_data/droplets_%d.csv", timestep) >= max_len) {
        Error("filename is larger than allocated buffer\n");
    }

    fptr = fopen(fname, "w");

    if (fptr == NULL) {
        Error("Error writing to file.\n");
    }

    Message("File: \"%s\" opened on host process.\n", fname);

    /* Prints header to the csv file */
    if (PRINT_SIMDATA) {

        sim_id = RP_Get_Integer("user/sim_id");
        n_cycles = RP_Get_Integer("user/n_cycles");
        vibration_amplitude = RP_Get_Real("user/vibration_amplitude");
        vibration_frequency = RP_Get_Real("user/vibration_frequency");

        /* 1 = simple vibration, 2 = simple vibration + noise, 3 = solid coupled */
        if (sim_id == 1) {

            fprintf(fptr, "#################################################################\n");
            fprintf(fptr, "    Simulation type: Simple vibration, timestep: %d, time: %e\n", timestep, CURRENT_TIME);
            fprintf(fptr, "-----------------------------------------------------------------\n");
            fprintf(fptr, "    Parameter Values:\n");
            fprintf(fptr, "        Number of cycles = %d\n", n_cycles);
            fprintf(fptr, "        Vibration amplitude = %e\n", vibration_amplitude);
            fprintf(fptr, "        Vibration frequency = %e\n", vibration_frequency);
            fprintf(fptr, "#################################################################\n");

        } else if (sim_id == 2) {

            noise_amplitude = RP_Get_Real("user/noise_amplitude");
            noise_frequency = RP_Get_Real("user/noise_frequency");

            fprintf(fptr, "#################################################################\n");
            fprintf(fptr, "    Simulation type: Simple vibration + noise, timestep: %d, time: %e\n", timestep, CURRENT_TIME);
            fprintf(fptr, "-----------------------------------------------------------------\n");
            fprintf(fptr, "    Parameter Values:\n");
            fprintf(fptr, "        Number of cycles = %d\n", n_cycles);
            fprintf(fptr, "        Vibration amplitude = %e\n", vibration_amplitude);
            fprintf(fptr, "        Vibration frequency = %e\n", vibration_frequency);
            fprintf(fptr, "        Noise amplitude = %e\n", noise_amplitude);
            fprintf(fptr, "        Noise frequency = %e\n", noise_frequency);
            fprintf(fptr, "#################################################################\n");

        } else if (sim_id == 3) {

            x_grid_position = RP_Get_Integer("user/x_grid_position");
            y_grid_position = RP_Get_Integer("user/y_grid_position");

            fprintf(fptr, "#################################################################\n");
            fprintf(fptr, "    Simulation type: Solid coupled, timestep: %d, time: %e\n", timestep, CURRENT_TIME);
            fprintf(fptr, "-----------------------------------------------------------------\n");
            fprintf(fptr, "    Parameter Values:\n");
            fprintf(fptr, "        Number of cycles = %d\n", n_cycles);
            fprintf(fptr, "        Vibration amplitude = %e\n", vibration_amplitude);
            fprintf(fptr, "        Vibration frequency = %e\n", vibration_frequency);
            fprintf(fptr, "        X grid position = %d\n", x_grid_position);
            fprintf(fptr, "        Y grid position = %d\n", y_grid_position);
            fprintf(fptr, "#################################################################\n");

        }
    }

    /* droplet_id, volume of droplet, mass of droplet, x coord, y coord, z coord, u velocity, v velocity, w velocity */
    fprintf(fptr, "droplet_id,volume,mass,x,y,z,u,v,w\n");

    return fptr;
}


/* Saves droplet data to the file */
void save_line_to_file(FILE *fptr, real *droplet_values, int droplet_id) {
    fprintf(
        fptr,
        "%d,%e,%e,%e,%e,%e,%e,%e,%e\n",
        droplet_id,        /* droplet id     */
        droplet_values[0], /* droplet volume */
        droplet_values[1], /* mass           */
        droplet_values[2], /* x centroid     */
        droplet_values[3], /* y centroid     */
        droplet_values[4], /* z centroid     */
        droplet_values[5], /* x velocity     */
        droplet_values[6], /* y velocity     */
        droplet_values[7]  /* z velocity     */
    );
}


/* Receive data from single compute node */
void receive_node_data_singlethreaded(FILE *fptr) {

    int message, droplet_id;
    droplet_id = 1;

    real droplet_values[8] = {0.};

    while (1) {

        /*
         * Get message from node zero
         * -1 complete
         * 0 droplet inside node
         * 1 droplet outside node
         */
        PRF_CRECV_INT(node_zero, &message, 1, droplet_id);

        if (message == -1) {
            break;
        } else if (message == 0) {
            /* Receive values for droplet & save to file */
            PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
            save_line_to_file(fptr, droplet_values, droplet_id);
        }
        droplet_id++;
    }
}


/* Receive data from node zero first pass */
void receive_node_zero_data(FILE *fptr, datastorage *droplets_datastorage) {

    int message, droplet_id = 1;
    real droplet_values[8] = {0.};

    while (1) {

        /*
         * Get message from node zero
         * -1 complete
         *  0 droplet inside node
         *  1 droplet outside node
         */
        PRF_CRECV_INT(node_zero, &message, 1, droplet_id);

        if (message == -1) {
            break;
        } else if (message == 0) {

            /* Receive values for droplet & save to file */
            PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
            save_line_to_file(fptr, droplet_values, droplet_id);

        } else if (message == 1) {

            /* Receive values for droplet & store */
            PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
            add_droplet_values_to_datastorage(droplets_datastorage, droplet_values, droplet_id);

        }
        droplet_id += compute_node_count;
    }
}


/* Receive data from each compute node through node zero (excluding node zero) */
void receive_compute_node_data(FILE *fptr, datastorage *droplets_datastorage) {

    int message = 0, droplet_id = 2, all_nodes_completed = 0;
    int* nodes_completed;
    nodes_completed = (int*)calloc((compute_node_count - 1), sizeof(int));

    real droplet_values[8] = {0.};

    while (!all_nodes_completed) {

        all_nodes_completed = 1;

        for (int current_node = 1; current_node < compute_node_count; ++current_node) {
            if (!nodes_completed[current_node - 1]) {

                /*
                 * Get message from node zero
                 * -1 complete
                 *  0 droplet inside node
                 *  1 droplet outside node
                 */
                PRF_CRECV_INT(node_zero, &message, 1, droplet_id);

                if (message == -1) {

                    nodes_completed[current_node - 1] = 1; /* Node calculations completed */

                } else if (message == 0) {

                    /* Receive values for droplet and save to file */
                    PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
                    save_line_to_file(fptr, droplet_values, droplet_id);

                    all_nodes_completed = 0;    /* Keep looping */

                } else if (message == 1) {

                    /* Receive values for droplet and & store */
                    PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
                    add_droplet_values_to_datastorage(droplets_datastorage, droplet_values, droplet_id);

                    all_nodes_completed = 0;    /* Keep looping */

                }
            }
            droplet_id++;   /* Increment to next node */
        }
        droplet_id++;   /* Skip node zero */
    }
}


/* Receive node zero droplet connection data */
void receive_node_zero_connections(datastorage *droplets_datastorage) {

    int message, *droplet_ids;

    while (1) {

        /*
         * message = -1 for no more connections to receive
         * message >  0 message is equal to size of array to receive
         */
        PRF_CRECV_INT(node_zero, &message, 1, node_zero);

        if (message == -1) {
            break;
        } else {

            /* Received node connections and assign droplets */
            droplet_ids = (int*)malloc(message * sizeof(int));
            PRF_CRECV_INT(node_zero, droplet_ids, message, node_zero);
            assign_droplets_to_combine(droplets_datastorage, message, droplet_ids);

        }
    }
}


/* Receive other compute nodes droplet connection data through node zero (excluding node zero) */
void receive_compute_node_connections(datastorage *droplets_datastorage) {

    int message, all_nodes_completed, *nodes_completed, *droplet_ids;
    all_nodes_completed = 0;
    nodes_completed = (int*)calloc((compute_node_count - 1), sizeof(int));

    while (!all_nodes_completed) {

        all_nodes_completed = 1;

        for (int current_node = 1; current_node < compute_node_count; ++current_node) {
            if (nodes_completed[current_node - 1] == 0) {

                /*
                 * Receive message from node zero
                 * -1 = no more messages
                 * >0 = size of array to receive
                 */
                PRF_CRECV_INT(node_zero, &message, 1, node_zero);

                if (message == -1) {
                    nodes_completed[current_node - 1] = 1; /* current node has no more messages to send */
                } else {

                    /* Allocate an array and assign the droplets to be combined */
                    droplet_ids = (int*)malloc(message * sizeof(int));
                    PRF_CRECV_INT(node_zero, droplet_ids, message, node_zero);
                    assign_droplets_to_combine(droplets_datastorage, message, droplet_ids);
                    all_nodes_completed = 0;

                }
            }
        }
    }
}


/* Combine connected droplets on node boundaries */
void combine_boundary_droplets(FILE *fptr, datastorage *droplets_datastorage) {

    int droplet_id;

    real droplet_values[8];

    for (int combination_id = 1; combination_id < droplets_datastorage->n_to_combine; ++combination_id) {

        memset(droplet_values, (real) 0., 8 * sizeof(real)); /* Reset values to 0 */

        droplet_id = get_values_from_datastorage(droplets_datastorage, droplet_values, combination_id);

        if (droplet_id > 0) {
            save_line_to_file(fptr, droplet_values, droplet_id);
        }
    }
}


/* Host node main compute process singlethreaded */
void host_process_singlethreaded() {

    FILE *fptr = file_handler();

    receive_node_data_singlethreaded(fptr);

    fclose(fptr);
}


/* Host node main compute process multithreaded */
void host_process_multithreaded() {

    FILE *fptr = file_handler();                    /* Create file and open for saving */

    datastorage droplets_datastorage;
    initialize_datastorage(&droplets_datastorage);             /* Initialize datastorage struct */

    receive_node_zero_data(fptr, &droplets_datastorage);     /* Write node 0 first */

    receive_compute_node_data(fptr, &droplets_datastorage);  /* Write other compute nodes */

    receive_node_zero_connections(&droplets_datastorage);    /* Receive connection data from node zero */

    receive_compute_node_connections(&droplets_datastorage); /* Receive connection data from other compute nodes */

    combine_boundary_droplets(fptr, &droplets_datastorage);  /* Combine droplets on boundaries */

    fclose(fptr);
}