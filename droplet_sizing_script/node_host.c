#include "node_host.h"



FILE *file_handler(){

    char fname[50];
    int max_len = sizeof fname, timestep = N_TIME;

    if (snprintf(fname, max_len, "droplets_%d.csv", timestep) >= max_len){
        Error("filename is larger than allocated buffer\n");
    } 

    FILE *fptr = NULL;

    fptr = fopen(fname, "w");
    if (fptr == NULL){
        Error("Error writing to file.\n");
    }

    Message("File: \"%s\" opened on host process.\n", fname);

    // droplet_id,volume,mass,u_x,u_y,u_z,v_x,v_y,v_z
    fprintf(fptr, "droplet_id,volume,mass,ux,uy,uz,vx,vy,vz\n");
    return fptr;
}


void save_line_to_file(FILE *fptr, real *droplet_values, int droplet_id){
    fprintf(
        fptr, 
        "%d,%e,%e,%e,%e,%e,%e,%e,%e\n",
        droplet_id,
        droplet_values[0], // droplet volume
        droplet_values[1], // mass
        droplet_values[2], // x centroid
        droplet_values[3], // y centroid
        droplet_values[4], // z centroid
        droplet_values[5], // x velocity
        droplet_values[6], // y velocity
        droplet_values[7]  // z velocity
    );
}


void receive_node_zero_data(FILE *fptr, Datastorage *datastorage){
    
    int message, cell_to_explore_from, droplet_id = 1;
    real droplet_values[8] = {0.};

    while (true) {

        // Get message from node zero
        // -1 complete
        //  0 droplet inside node
        //  1 droplet outside node
        PRF_CRECV_INT(node_zero, &message, 1, droplet_id);

        if (message == -1){
            break;
        } else if (message == 0){
            // Receive values for droplet & save to file
            PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
            save_line_to_file(fptr, droplet_values, droplet_id);

        } else if (message == 1) {
            // Receive values for droplet & save values & id
            PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
            addValues(datastorage, droplet_values, droplet_id);
        }
        droplet_id += compute_node_count;
    }
}


void receive_compute_node_data(FILE *fptr, Datastorage *datastorage){

    int message, cell_to_explore_from, droplet_id = 2;
    int *nodes_completed = (int*)calloc((compute_node_count - 1), sizeof(int));
    bool all_nodes_completed = false;
    real droplet_values[8] = {0.};

    // Write all other nodes
    while (!all_nodes_completed) {

        all_nodes_completed = true;

        for (int current_node = 1; current_node < compute_node_count; current_node++){
            if (nodes_completed[current_node - 1] == 0){

                // message = -1 for no more droplets
                // message = 0 for droplet within single domain
                // message = 1 for droplet in more domains
                PRF_CRECV_INT(node_zero, &message, 1, droplet_id);
                
                if (message == -1){

                    // Node calculations completed
                    nodes_completed[current_node - 1] = 1;

                } else if (message == 0){
                    
                    // Receive values for droplet
                    PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
                    save_line_to_file(fptr, droplet_values, droplet_id);
                    
                    // Keep looping
                    all_nodes_completed = false;

                } else if (message == 1) {
                    PRF_CRECV_REAL(node_zero, droplet_values, 8, droplet_id);
                    addValues(datastorage, droplet_values, droplet_id);
                    all_nodes_completed = false;
                }
            }
            droplet_id++;
        }
        droplet_id++;
    }
}


void combine_boundary_droplets(FILE *fptr, Datastorage *datastorage){
    int message;
    int *nodes_completed = (int*)calloc((compute_node_count - 1), sizeof(int));
    bool all_nodes_completed = false;

    // Receive node zero droplet data
    while (true) {

        PRF_CRECV_INT(node_zero, &message, 1, node_zero);

        if (message == -1){
            break;
        } else {
            int *droplet_ids = (int*)malloc(message * sizeof(int));
            PRF_CRECV_INT(node_zero, droplet_ids, message, node_zero);
            assignDroplets(datastorage, message, droplet_ids);
        }
    }

    // Receive other compute nodes droplet data
    while (!all_nodes_completed) {

        all_nodes_completed = true;

        for (int current_node = 1; current_node < compute_node_count; current_node++){
            if (nodes_completed[current_node - 1] == 0){

                // Receive message, -1 = done, else is size of int array to receive
                PRF_CRECV_INT(node_zero, &message, 1, node_zero);

                if (message == -1){
                    nodes_completed[current_node-1] = 1;
                } else {
                    int *droplet_ids = (int*)malloc(message * sizeof(int));
                    PRF_CRECV_INT(node_zero, droplet_ids, message, node_zero);
                    assignDroplets(datastorage, message, droplet_ids);
                    all_nodes_completed = false;
                }
            }
        }

    }

    real droplet_values[8];
    int droplet_id;
    Message("N Droplets: %d\n", datastorage->n_droplets);

    for (int combination_id = 1; combination_id < datastorage->n_to_combine; combination_id++){

        memset(droplet_values, (real) 0., 8 * sizeof(real));
        droplet_id = getValues(datastorage, droplet_values, combination_id);

        if (droplet_id > 0){
            save_line_to_file(fptr, droplet_values, droplet_id);
        }
        
    }
}


void host_process(){

    FILE *fptr = file_handler();

    Datastorage datastorage; 

    // Write node 0 first
    receive_node_zero_data(fptr, &datastorage);

    // Write other compute nodes
    receive_compute_node_data(fptr, &datastorage);

    // Combine droplets on boundaries
    combine_boundary_droplets(fptr, &datastorage);

    fclose(fptr);
}