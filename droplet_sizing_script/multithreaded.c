#include "multithreaded.h"

// Calculate droplet sizes function
DEFINE_ON_DEMAND(calculate_droplet_sizes){
    clock_t begin = clock();
    #if !RP_HOST

        Stack cells_to_reexplore;
        initialize(&cells_to_reexplore);

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
        Message("\n");
        host_process();
    #endif

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    #if !RP_NODE
        Message("Host node, completed successfully in %.3f seconds\n", time_spent);
    #else
        Message("Node %d, completed successfully in %.3f seconds\n", myid, time_spent);
    #endif
    
}


void found_new_droplet(cell_t first_cell, Thread *cell_thread, Thread *water_thread, int droplet_id, Stack *cells_to_reexplore){

    real vof= C_VOF(first_cell, water_thread), temp_value, droplet_values[8] = {0.}; // 0 = vol, 1 = mass, 2-4 = centroid, 5-7 = velocity
    int cell_explored, local_face_id, adjacent_cell_node_id, message;
    int receiving_node = (I_AM_NODE_ZERO_P ? node_host : node_zero);
    bool droplet_outside_cell = false;
    Thread *cell_face_thread;
    cell_t current_cell, adjacent_cell;
    face_t cell_face;
    
    Stack stack;
    initialize(&stack);
    C_UDMI(first_cell, cell_thread, 0) = droplet_id; // Cell explored
    push(&stack, first_cell); // Push to the stack

    first_value_update(droplet_values, vof, first_cell, cell_thread);

    while (!isEmpty(&stack)) {

        current_cell = pop(&stack);

        // Loop over faces of current cell
        c_face_loop(current_cell, cell_thread, local_face_id){

            cell_face = C_FACE(current_cell, cell_thread, local_face_id);
            cell_face_thread = C_FACE_THREAD(current_cell, cell_thread, local_face_id);

            // If the face is not on the boundary
            if (!BOUNDARY_FACE_THREAD_P(cell_face_thread)){

                // Get the new adjacent cell
                adjacent_cell = F_C0(cell_face, cell_face_thread);
                if (adjacent_cell == current_cell) {
                    adjacent_cell = F_C1(cell_face, cell_face_thread);
                }
                
                vof = C_VOF(adjacent_cell, water_thread);
                cell_explored = C_UDMI(adjacent_cell, cell_thread, 0);
                
                // If new cell has not been explored and is droplet
                if (cell_explored == 0 && vof > 0.5) {

                    adjacent_cell_node_id = C_PART(adjacent_cell, cell_thread);

                    if (adjacent_cell_node_id == myid) {

                        push(&stack, adjacent_cell);
                        C_UDMI(adjacent_cell, cell_thread, 0) = droplet_id;
                        subsequent_value_update(droplet_values, vof, adjacent_cell, cell_thread);

                    } else {
                        droplet_outside_cell = true;
                        push(cells_to_reexplore, adjacent_cell);
                    }

                } else if (cell_explored == 0) {
                    C_UDMI(adjacent_cell, cell_thread, 0) = -1;
                }
            } 
        }
    }

    if (!droplet_outside_cell){
        
        final_value_update(droplet_values);

        message = 0;
        PRF_CSEND_INT(receiving_node, &message, 1, droplet_id);
        PRF_CSEND_REAL(receiving_node, droplet_values, 8, droplet_id);
        
    } else {

        message = 1;
        PRF_CSEND_INT(receiving_node, &message, 1, droplet_id);
        PRF_CSEND_REAL(receiving_node, droplet_values, 8, droplet_id);

        push(cells_to_reexplore, first_cell);
        addGap(cells_to_reexplore);
    }
}


void compute_droplet_data(Stack *cells_to_reexplore){

    Message("Node %d initializing\n", myid);
    Domain *mixture_domain = Get_Domain(1);
    Domain *water_domain = DOMAIN_SUB_DOMAIN(mixture_domain, S_PHASE);
    Thread *cell_thread, *water_thread, *adjacent_cell_thread, *cell_face_thread;
    Thread **pt;
    face_t cell_face;
    cell_t cell, current_cell, adjacent_cell;
    bool droplet_outside_cell = false;
    int cell_explored, message, droplet_id = myid + 1;
    int receiving_node = (I_AM_NODE_ZERO_P ? node_host : node_zero);
    real vof;

    Message("Node %d initialized\n", myid);

    mp_thread_loop_c(cell_thread, mixture_domain, pt){
        if (FLUID_THREAD_P(cell_thread)){

            Message("Node %d starting calculation\n", myid);

            water_thread = pt[S_PHASE];

            begin_c_loop_int(cell, cell_thread){

                vof = C_VOF(cell, water_thread);
                cell_explored = C_UDMI(cell, cell_thread, 0);
                
                // Cell in new droplet
                if (cell_explored == 0 && vof > 0.5) { 

                    found_new_droplet(cell, cell_thread, water_thread, droplet_id, cells_to_reexplore);
                    droplet_id += compute_node_count;

                } else if (cell_explored == 0) { // Cell not in droplet
                    C_UDMI(cell, cell_thread, 0) = -1; 
                }
            }
            end_c_loop_int(cell,cell_thread)
        }
    }
    message = -1;
    PRF_CSEND_INT(receiving_node, &message, 1, droplet_id);
    EXCHANGE_SVAR_MESSAGE(mixture_domain, (SV_UDM_I, SV_NULL));
}


void init_udm(){
    Domain *mixture_domain = Get_Domain(1);
    Thread *cell_thread, **pt;
    cell_t cell;
    mp_thread_loop_c(cell_thread, mixture_domain, pt){

        if (FLUID_THREAD_P(cell_thread)){
            begin_c_loop_int(cell, cell_thread){
                C_UDMI(cell, cell_thread, 0) = 0;
            } end_c_loop_int(cell, cell_thread)
        }
    }
}


bool unique(int *attached_droplets, int n_droplets, int new_droplet_id){

    for (int droplet = 0; droplet < n_droplets; droplet++){
        if (attached_droplets[droplet] == new_droplet_id) {
            return false;
        }
    }
    return true;
}


void send_droplet_message(int message, int *attached_droplets, int receiving_node){

    int *droplet_ids = (int*)malloc(message * sizeof(int));
    
    for (int i = 0; i < message; i++){
        droplet_ids[i] = attached_droplets[i];
    }
    
    PRF_CSEND_INT(receiving_node, &message, 1, myid);
    PRF_CSEND_INT(receiving_node, droplet_ids, message, myid);
}


void assemble_droplets(Stack *cells_to_reexplore){
    Domain *mixture_domain = Get_Domain(1);
    Thread *cell_thread, **pt;
    cell_t cell;
    int current_droplet_id, attached_droplets[MAX_DROPLET_COMBINES] = {0}, n_droplets = -1, new_droplet_id, message;
    int receiving_node = (I_AM_NODE_ZERO_P ? node_host : node_zero);
    mp_thread_loop_c(cell_thread, mixture_domain, pt){

        if (FLUID_THREAD_P(cell_thread)){

            while(!isEmpty(cells_to_reexplore)) {

                if (isGap(cells_to_reexplore)){

                    if (n_droplets >= 0){
                        send_droplet_message(n_droplets, attached_droplets, receiving_node);  
                    }

                    remGap(cells_to_reexplore);
                    cell = pop(cells_to_reexplore);
                    attached_droplets[0] = C_UDMI(cell, cell_thread, 0);
                    n_droplets = 1;

                } else {
                    cell = pop(cells_to_reexplore);
                    new_droplet_id = C_UDMI(cell, cell_thread, 0);

                    if (unique(attached_droplets, n_droplets, new_droplet_id)){
                        attached_droplets[n_droplets++] = new_droplet_id;
                    }
                }
            }

            // Send final droplet information
            if (n_droplets > 0){
                send_droplet_message(n_droplets, attached_droplets, receiving_node);
            }

            message = -1;
            PRF_CSEND_INT(receiving_node, &message, 1, myid);
        }
    }
}   