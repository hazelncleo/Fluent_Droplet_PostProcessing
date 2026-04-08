
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdbool.h>
#include "udf.h"
#include "stack.h"
#include "data_storage.h"
#include "multithreaded.h"


// Calculate droplet sizes function
DEFINE_ON_DEMAND(calculate_droplet_sizes){
    calculate_droplets();
}


void calculate_droplets(){

    // Node zero first moves other nodes data to host
    #if !RP_HOST

        compute_droplet_data();

        if I_AM_NODE_ZERO_P {
            node_zero_send_data();
        }
        
        
    #endif

    #if !RP_NODE
        host_write_to_file();
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
    fprintf(fptr, "droplet_id,volume,ux,uy,uz,vx,vy,vz,mass\n");
    return fptr;
}


void node_zero_send_data(){
    
    Message("Node zero!\n");
    int message, i, droplet_id, j = 0, nodes_completed[compute_node_count];
    bool all_nodes_completed = false;
    real droplet_values[8] = {0.};
    memset(nodes_completed,0,sizeof(nodes_completed));

    while (!all_nodes_completed) {

        all_nodes_completed = true;

        for (i = 1; i < compute_node_count; i++){
            if (nodes_completed[i] == 0){

                // message = -1 for no more droplets
                // message = 0 for droplet within single domain
                // message = 1 for droplet in more domains
                droplet_id = (i + 1) + (j * compute_node_count);
                PRF_CRECV_INT(i, &message, 1, droplet_id);
                PRF_CSEND_INT(node_host, &message, 1, droplet_id);
                
                if (message == -1){

                    // Node calculations completed
                    nodes_completed[i] = 1;

                } else if (message == 0){
                    
                    // Receive values for droplet
                    PRF_CRECV_REAL(i, droplet_values, 8, droplet_id);
                    PRF_CSEND_REAL(node_host, droplet_values, 8, droplet_id);
                    
                    // Keep looping
                    all_nodes_completed = false;

                } else if (message == -1) {
                    // TODO
                    all_nodes_completed = false;
                }
            }
        }
        j++;
    }
}


void compute_droplet_data(){

    Message("\nNode %d initializing\n", myid);

    Domain *mixture_domain = Get_Domain(1);
    Domain *water_domain = DOMAIN_SUB_DOMAIN(mixture_domain, S_PHASE);
    Thread *cell_thread, *water_thread, *adjacent_cell_thread, *cell_face_thread;
    Thread **pt;
    face_t cell_face;
    cell_t cell, current_cell, adjacent_cell;    
    int local_face_id, cell_node_id, droplet_id[2] = {0}; // 0 = id, 1 = explored
    real droplet_values[8] = {0.}; // 0 = vol, 1-3 = centroid, 4-6 = velocity, 7 = mass
    real calc_values[2] = {0.}; // 0 = vof, 1 = temp value
    bool droplet_outside_cell = false;
    droplet_id[0] = myid + 1;
    
    Message("\nNode %d initialized\n", myid);

    Stack stack;
    initialize(&stack);

    mp_thread_loop_c(cell_thread, mixture_domain, pt){
        if (FLUID_THREAD_P(cell_thread)){

            Message("\nNode %d starting calculation\n", myid);
            // Initially all cells unexplored (UDM=0)
            begin_c_loop(cell, cell_thread){
                C_UDMI(cell, cell_thread, 0) = 0;
            }
            end_c_loop(cell, cell_thread)

            water_thread = pt[S_PHASE];

            begin_c_loop(cell, cell_thread){

                calc_values[0] = C_VOF(cell, water_thread);
                droplet_id[1] = C_UDMI(cell, cell_thread, 0);
                
                // Cell in new droplet
                if (droplet_id[1] == 0 && calc_values[0] > 0.5) { 

                    C_UDMI(cell, cell_thread, 0) = droplet_id[0]; // Cell explored
                    push(&stack, cell); // Push to the stack

                    first_value_update(droplet_values, calc_values, cell, cell_thread);

                    cell_node_id = C_PART(cell, cell_thread);
                    if (cell_node_id != myid) {
                        droplet_outside_cell = true;
                    }

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
                                
                                calc_values[0] = C_VOF(adjacent_cell, water_thread);
                                droplet_id[1] = C_UDMI(adjacent_cell, cell_thread, 0);
                                
                                // If new cell has not been explored and is droplet
                                if (droplet_id[1] == 0 && calc_values[0] > 0.5) {

                                    push(&stack, adjacent_cell);
                                    C_UDMI(adjacent_cell, cell_thread, 0) = droplet_id[0];

                                    subsequent_value_update(droplet_values, calc_values, adjacent_cell, cell_thread);

                                    cell_node_id = C_PART(adjacent_cell, cell_thread);
                                    if (cell_node_id != myid) {
                                        droplet_outside_cell = true;
                                    }

                                } else if (droplet_id[1] == 0) {
                                    C_UDMI(adjacent_cell, cell_thread, 0) = -1;
                                }
                            } 
                        }
                    }

                    final_value_update(droplet_values, calc_values);

                    int temp = 0;
                    if I_AM_NODE_ZERO_P {
                        PRF_CSEND_INT(node_host, &temp, 1, droplet_id[0]);
                        PRF_CSEND_REAL(node_host, droplet_values, 8, droplet_id[0]);
                    } else {
                        PRF_CSEND_INT(node_zero, &temp, 1, droplet_id[0]);
                        PRF_CSEND_REAL(node_zero, droplet_values, 8, droplet_id[0]);
                    }
                    

                    droplet_id[0] += compute_node_count;
                    reset(droplet_values, calc_values);

                } else if (droplet_id[1] == 0) { // Cell not in droplet
                    C_UDMI(cell, cell_thread, 0) = -1; 
                }
            }
            end_c_loop(cell,cell_thread)
        }
    }
    int temp = -1;
    if I_AM_NODE_ZERO_P {
        PRF_CSEND_INT(node_host, &temp, 1, droplet_id[0]);
    } else {
        PRF_CSEND_INT(node_zero, &temp, 1, droplet_id[0]);
    }
    Message("\nNode %d completed!\n", myid);
}


void host_write_to_file(){
    FILE *fptr = file_handler();
    Message("\nFile opened on host process\n");

    int message, i, droplet_id, j = 0, nodes_completed[compute_node_count];
    bool all_nodes_completed = false;
    real droplet_values[8] = {0.};
    memset(nodes_completed,0,sizeof(nodes_completed));
    Message("\nHost process initialized\n");

    // Write node 0 first
    droplet_id = 1;
    while (nodes_completed[0] == 0) {
        PRF_CRECV_INT(node_zero, &message, 1, droplet_id);

        if (message == -1){

            // Node calculations completed
            break;

        } else if (message == 0){
            
            // Receive values for droplet
            PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
            fprintf(
                fptr, 
                "%d,%e,%e,%e,%e,%e,%e,%e,%e\n",
                droplet_id,
                droplet_values[0], // droplet volume
                droplet_values[1], // x centroid
                droplet_values[2], // y centroid
                droplet_values[3], // z centroid
                droplet_values[4], // x velocity
                droplet_values[5], // y velocity
                droplet_values[6], // z velocity
                droplet_values[7]  // mass
            );

        } else if (message == -1) {
        }
        droplet_id += compute_node_count;
    }

    // Write all other nodes
    while (!all_nodes_completed) {

        all_nodes_completed = true;

        for (i = 1; i < compute_node_count; i++){
            if (nodes_completed[i] == 0){

                // message = -1 for no more droplets
                // message = 0 for droplet within single domain
                // message = 1 for droplet in more domains
                droplet_id = (i + 1) + (j * compute_node_count);
                PRF_CRECV_INT(node_zero, &message, 1, droplet_id);
                
                if (message == -1){

                    // Node calculations completed
                    nodes_completed[i] = 1;

                } else if (message == 0){
                    
                    // Receive values for droplet
                    PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
                    fprintf(
                        fptr, 
                        "%d,%e,%e,%e,%e,%e,%e,%e,%e\n",
                        droplet_id,
                        droplet_values[0], // droplet volume
                        droplet_values[1], // x centroid
                        droplet_values[2], // y centroid
                        droplet_values[3], // z centroid
                        droplet_values[4], // x velocity
                        droplet_values[5], // y velocity
                        droplet_values[6], // z velocity
                        droplet_values[7]  // mass
                    );
                    
                    // Keep looping
                    all_nodes_completed = false;

                } else if (message == -1) {
                    // TODO
                    all_nodes_completed = false;
                }
            }
        }
        
        j++;
    }

    fclose(fptr);
}