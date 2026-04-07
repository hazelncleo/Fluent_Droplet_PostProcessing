
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdbool.h>
#include "udf.h"
#include "stack.h"
#include "data_storage.h"
#include "multithreaded.h"


// Calculate droplet sizes function
DEFINE_ON_DEMAND(calculate_droplet_sizes){

    // Open file & write first line
    #if !RP_NODE
        FILE *fptr = file_handler();
    #endif

    // Sync all nodes
    PRF_GSYNC();

    calculate_droplets();

    fclose(fptr);
}


void calculate_droplets(){
    #if !RP_HOST
    /*Domain *mixture_domain = Get_Domain(1);
    Domain *water_domain = DOMAIN_SUB_DOMAIN(mixture_domain, S_PHASE);
    Thread *cell_thread, *water_thread, *adjacent_cell_thread, *cell_face_thread;
    Thread **pt;
    face_t cell_face;
    cell_t cell, current_cell, adjacent_cell;    
    int local_face_id;

    Stack stack;
    initialize(&stack);

    DataStorage datastorage;
    initialize_storage(&datastorage);

    mp_thread_loop_c(cell_thread, mixture_domain, pt){
        if (FLUID_THREAD_P(cell_thread)){

            // Initially all cells unexplored (UDM=0)
            begin_c_loop(cell, cell_thread){
                C_UDMI(cell, cell_thread, 0) = 0;
            }
            end_c_loop(cell, cell_thread)

            water_thread = pt[S_PHASE];

            begin_c_loop(cell, cell_thread){

                datastorage.vof_water = C_VOF(cell, water_thread);
                datastorage.explored = C_UDMI(cell, cell_thread, 0);
                
                // Cell in new droplet
                if (datastorage.explored == 0 && datastorage.vof_water > 0.5) { 

                    C_UDMI(cell, cell_thread, 0) = datastorage.droplet_id; // Cell explored
                    push(&stack, cell); // Push to the stack
                    

                    first_value_update(&datastorage, cell, cell_thread);

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

                                datastorage.vof_water = C_VOF(adjacent_cell, water_thread);
                                datastorage.explored = C_UDMI(adjacent_cell, cell_thread, 0);
                                
                                // If new cell has not been explored and is droplet
                                if (datastorage.explored == 0 && datastorage.vof_water > 0.5) {

                                    push(&stack, adjacent_cell);
                                    C_UDMI(adjacent_cell, cell_thread, 0) = datastorage.droplet_id;

                                    subsequent_value_update(&datastorage, adjacent_cell, cell_thread);

                                } else if (datastorage.explored == 0) {
                                    C_UDMI(adjacent_cell, cell_thread, 0) = -1;
                                }
                            } 
                        }
                    }

                    final_value_update(&datastorage);

                    if I_AM_NODE_ZERO_P {
                        write_droplet_data_to_file(&datastorage, fptr);
                    }
                    datastorage.droplet_id++;

                } else if (datastorage.explored == 0) { // Cell not in droplet
                    C_UDMI(cell, cell_thread, 0) = -1; 
                }
            }
            end_c_loop(cell,cell_thread)
        }
    }
    if I_AM_NODE_ZERO_P {
        Message("\n------------------\nCalculation completed! (number of droplets = %d)\n------------------\n", datastorage.droplet_id);
    }*/
    #endif

    #if !RP_NODE

        int nodes_completed[compute_node_count] = {0}, message[2];
        bool all_nodes_completed = false;

        while (!all_nodes_completed) {

            // Receive status from node_zero
            // message[0] = id of node message from
            // message[1] = -1 for no more droplets
            // message[1] = 0 for droplet within single domain
            // message[1] = 1 for droplet in more domains
            PRF_CRECV_INT(node_zero, message, 2, node_zero);

        }





    #endif
}


FILE *file_handler(){

    char fname[50];
    int max_len = sizeof fname, timestep = N_TIME;

    if (snprintf(fname, max_len, "droplets_%d.csv", timestep) >= max_len){
        Error("filename is larger than allocated buffer\n");
    } 

    FILE *fptr = NULL;

    fptr = fopen(fname, "w");
    if (fptr == NULL){
        Error("Error writing to file.");
    }

    // droplet_id,volume,u_x,u_y,u_z,v_x,v_y,v_z
    fprintf(fptr, "droplet_id,volume,ux,uy,uz,vx,vy,vz\n");

    return fptr;
}


void combine_droplets(){

}