#ifndef DATA_STORAGE_H_   /* Include guard */
#define DATA_STORAGE_H_
#include "udf.h"
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
void initialize_storage(DataStorage *datastorage);

// Reset values
void reset(DataStorage *datastorage);

// Calculate initial values for first cell in droplet
void first_value_update(DataStorage *datastorage, cell_t cell, Thread *cell_thread);

// Calculate values for each subsequent cell in droplet
void subsequent_value_update(DataStorage *datastorage, cell_t adjacent_cell, Thread *cell_thread);

// Calculate final values
void final_value_update(DataStorage *datastorage);

// Write droplet data to the file
void write_droplet_data_to_file(DataStorage *datastorage, FILE *fptr);

#endif