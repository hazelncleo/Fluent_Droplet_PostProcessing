#include "node_zero.h"



void node_zero_send_data(){
    
    int message, cell_to_explore_from, droplet_id = 2;
    int *nodes_completed = (int*)calloc((compute_node_count - 1), sizeof(int));
    bool all_nodes_completed = false;
    real droplet_values[8] = {0.};

    while (!all_nodes_completed) {

        all_nodes_completed = true;

        for (int current_node = 1; current_node < compute_node_count; current_node++){
            if (nodes_completed[current_node] == 0){

                // message = -1 for no more droplets
                // message = 0 for droplet within single domain
                // message = 1 for droplet in more domains
                PRF_CRECV_INT(current_node, &message, 1, droplet_id);
                PRF_CSEND_INT(node_host, &message, 1, droplet_id);
                
                if (message == -1){

                    // Node calculations completed
                    nodes_completed[current_node] = 1;

                } else if (message == 0){
                    
                    // Receive values for droplet
                    PRF_CRECV_REAL(current_node, droplet_values, 8, droplet_id);
                    PRF_CSEND_REAL(node_host, droplet_values, 8, droplet_id);
                    
                    // Keep looping
                    all_nodes_completed = false;

                } else if (message == 1) {
                    
                    PRF_CRECV_REAL(current_node, droplet_values, 8, droplet_id);
                    PRF_CSEND_REAL(node_host, droplet_values, 8, droplet_id);

                    all_nodes_completed = false;
                }
            }
            droplet_id++;
        }
        droplet_id++;
    }
}