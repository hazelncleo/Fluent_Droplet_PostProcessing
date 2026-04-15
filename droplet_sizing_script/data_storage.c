#include "udf.h"
#include "data_storage.h"

#define SECONDARY_DENSITY 997. // water density kg/m^3


void initializeDatastorage(Datastorage *datastorage){
    datastorage->n_droplets = 0;
    datastorage->n_to_combine = 0;
    memset(datastorage->combination_ids, 0, MAX_COMBINE_DROPLETS * sizeof(int));
}


void addValues(Datastorage *datastorage, real *droplet_values, int droplet_id){

    for (int i = 0; i < 8; i++){
        datastorage->droplet_values[datastorage->n_droplets][i] = droplet_values[i];
    }
    datastorage->droplet_ids[datastorage->n_droplets] = droplet_id;
    datastorage->n_droplets++;
}


int getIndex(Datastorage *datastorage, int droplet_id){

    for (int i = 0; i < datastorage->n_droplets; i++){
        if (datastorage->droplet_ids[i] == droplet_id){
            return i;
        }
    }
    Error("Error finding droplet");
    return -1;
}


int checkDropletsAssigned(Datastorage *datastorage, int n_to_assign, int *droplets){

    int index, combination_id = datastorage->n_to_combine + 1;

    for (int i = 0; i < n_to_assign; i++){

        index = getIndex(datastorage, droplets[i]);

        if (datastorage->combination_ids[index] > 0 && datastorage->combination_ids[index] < combination_id){
            combination_id = datastorage->combination_ids[index];
        }
    }

    if (combination_id == datastorage->n_to_combine + 1){
        datastorage->n_to_combine++;
    }

    return combination_id;
}


void assignDroplets(Datastorage *datastorage, int n_to_assign, int *droplets){

    int combination_id = checkDropletsAssigned(datastorage, n_to_assign, droplets);

    for (int i = 0; i < n_to_assign; i++){
        datastorage->combination_ids[getIndex(datastorage, droplets[i])] = combination_id;
    }
}


int getValues(Datastorage *datastorage, real *droplet_values, int combination_id){

    int droplet_id = 999999, inc=0;

    for (int i = 0; i < datastorage->n_droplets; i++){
        if (datastorage->combination_ids[i] == combination_id) {

            add_vec(droplet_values, datastorage->droplet_values[i]);
            inc++;
            if (datastorage->droplet_ids[i] < droplet_id){
                droplet_id = datastorage->droplet_ids[i];
            }
        }
    }

    if (droplet_id == 999999){
        return -1;
    }

    final_value_update(droplet_values);

    return droplet_id;
}


void add_vec(real *droplet_values, real *values_to_add){
    for (int i = 0; i < 8; i++){
        droplet_values[i] += values_to_add[i];
    }
}


// Calculate initial values for first cell in droplet
void first_value_update(real *droplet_values, real vof, cell_t cell, Thread *cell_thread){

    droplet_values[0] = C_VOLUME(cell, cell_thread) * vof;
    droplet_values[1] = droplet_values[0] * SECONDARY_DENSITY;

    real temp[3];
    C_CENTROID(temp, cell, cell_thread);
    droplet_values[2] = temp[0] * droplet_values[1];
    droplet_values[3] = temp[1] * droplet_values[1];
    droplet_values[4] = temp[2] * droplet_values[1];

    droplet_values[5] = C_U(cell, cell_thread) * droplet_values[1];
    droplet_values[6] = C_V(cell, cell_thread) * droplet_values[1];
    droplet_values[7] = C_W(cell, cell_thread) * droplet_values[1];
}


// Calculate values for each subsequent cell in droplet
void subsequent_value_update(real *droplet_values, real vof, cell_t adjacent_cell, Thread *cell_thread){

    real temp_value = C_VOLUME(adjacent_cell, cell_thread) * vof; 

    droplet_values[0] += temp_value; // Add volume

    temp_value *= SECONDARY_DENSITY;

    droplet_values[1] += temp_value; // Add mass

    real temp_vector[3];
    C_CENTROID(temp_vector, adjacent_cell, cell_thread);
    droplet_values[2] += temp_vector[0] * temp_value;
    droplet_values[3] += temp_vector[1] * temp_value;
    droplet_values[4] += temp_vector[2] * temp_value;

    droplet_values[5] += C_U(adjacent_cell, cell_thread) * temp_value;
    droplet_values[6] += C_V(adjacent_cell, cell_thread) * temp_value;
    droplet_values[7] += C_W(adjacent_cell, cell_thread) * temp_value;
}


// Calculate final values
void final_value_update(real *droplet_values){

    droplet_values[2] /= droplet_values[1];
    droplet_values[3] /= droplet_values[1];
    droplet_values[4] /= droplet_values[1];
    
    droplet_values[5] /= droplet_values[1];
    droplet_values[6] /= droplet_values[1];
    droplet_values[7] /= droplet_values[1];
}