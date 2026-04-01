#include "udf.h"
#include <stdbool.h>
#include <stdio.h>

#define SECONDARY_DENSITY 997. // water density kg/m^3

typedef struct {
    int droplet_id;
    real vof_water;
    int explored;
    real temp;
    real droplet_volume;
    real centroid[ND_ND];
    real temp_vector[ND_ND];
    real velocity[ND_ND];
    real mass;
} DataStorage;

// Initialise values
void initialize_storage(DataStorage *datastorage){
    datastorage->droplet_id = 1;
    datastorage->vof_water  = 0.;
    datastorage->explored   = 0;
    datastorage->temp       = 0.;
    datastorage->droplet_volume = 0.;
    datastorage->centroid[0] = 0.;
    datastorage->centroid[1] = 0.;
    datastorage->centroid[2] = 0.;
    datastorage->temp_vector[0] = 0.;
    datastorage->temp_vector[1] = 0.;
    datastorage->temp_vector[2] = 0.;
    datastorage->velocity[0] = 0.;
    datastorage->velocity[1] = 0.;
    datastorage->velocity[2] = 0.;
    datastorage->mass = 0.;
}

// Reset values
void reset(DataStorage *datastorage){
    datastorage->vof_water  = 0.;
    datastorage->explored   = 0;
    datastorage->temp       = 0.;
    datastorage->droplet_volume = 0.;
    datastorage->centroid[0] = 0.;
    datastorage->centroid[1] = 0.;
    datastorage->centroid[2] = 0.;
    datastorage->temp_vector[0] = 0.;
    datastorage->temp_vector[1] = 0.;
    datastorage->temp_vector[2] = 0.;
    datastorage->velocity[0] = 0.;
    datastorage->velocity[1] = 0.;
    datastorage->velocity[2] = 0.;
    datastorage->mass = 0.;
}

// Calculate initial values for first cell in droplet
void first_value_update(DataStorage *datastorage, cell_t cell, Thread *cell_thread){

    datastorage->droplet_volume = C_VOLUME(cell, cell_thread)*datastorage->vof_water;
    datastorage->mass = datastorage->droplet_volume * SECONDARY_DENSITY;

    C_CENTROID(datastorage->centroid, cell, cell_thread);
    datastorage->centroid[0] *= datastorage->mass;
    datastorage->centroid[1] *= datastorage->mass;
    datastorage->centroid[2] *= datastorage->mass;

    datastorage->velocity[0] = C_U(cell, cell_thread) * datastorage->mass;
    datastorage->velocity[1] = C_V(cell, cell_thread) * datastorage->mass;
    datastorage->velocity[2] = C_W(cell, cell_thread) * datastorage->mass;
}

// Calculate values for each subsequent cell in droplet
void subsequent_value_update(DataStorage *datastorage, cell_t adjacent_cell, Thread *cell_thread){
    datastorage->temp = C_VOLUME(adjacent_cell, cell_thread) * datastorage->vof_water;

    datastorage->droplet_volume += datastorage->temp;

    datastorage->temp *= SECONDARY_DENSITY;

    datastorage->mass += datastorage->temp;

    C_CENTROID(datastorage->temp_vector, adjacent_cell, cell_thread);

    datastorage->centroid[0] += datastorage->temp_vector[0] * datastorage->temp;
    datastorage->centroid[1] += datastorage->temp_vector[1] * datastorage->temp;
    datastorage->centroid[2] += datastorage->temp_vector[2] * datastorage->temp;

    datastorage->velocity[0] += C_U(adjacent_cell, cell_thread) * datastorage->temp;
    datastorage->velocity[1] += C_V(adjacent_cell, cell_thread) * datastorage->temp;
    datastorage->velocity[2] += C_W(adjacent_cell, cell_thread) * datastorage->temp;
}

// Calculate final values
void final_value_update(DataStorage *datastorage){

    datastorage->centroid[0] /= datastorage->mass;
    datastorage->centroid[1] /= datastorage->mass;
    datastorage->centroid[2] /= datastorage->mass;

    datastorage->velocity[0] /= datastorage->mass;
    datastorage->velocity[1] /= datastorage->mass;
    datastorage->velocity[2] /= datastorage->mass;
}

// Write droplet data to the file
void write_droplet_data_to_file(DataStorage *datastorage, FILE *fptr){
    fprintf(
        fptr, 
        "%d,%e,%e,%e,%e,%e,%e,%e\n",
        datastorage->droplet_id,
        datastorage->droplet_volume,
        datastorage->centroid[0],
        datastorage->centroid[1],
        datastorage->centroid[2],
        datastorage->velocity[0],
        datastorage->velocity[1],
        datastorage->velocity[2]
    );
}