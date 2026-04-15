#ifndef DATA_STORAGE_H_   /* Include guard */
#define DATA_STORAGE_H_

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#include <stdio.h>
#include <stdbool.h>
#include "udf.h"
#include "config_parameters.h"

typedef struct {
    int droplet_ids[MAX_COMBINE_DROPLETS];
    int combination_ids[MAX_COMBINE_DROPLETS];
    int n_droplets;
    int n_to_combine;
    real droplet_values[MAX_COMBINE_DROPLETS][8];
} Datastorage;

void initializeDatastorage(Datastorage *datastorage);

void addValues(Datastorage *datastorage, real *droplet_values, int droplet_id);

int getIndex(Datastorage *datastorage, int droplet_id);

int checkDropletsAssigned(Datastorage *datastorage, int n_to_assign, int *droplets);

void assignDroplets(Datastorage *datastorage, int n_to_assign, int *droplets);

int getValues(Datastorage *datastorage, real *droplet_values, int combination_id);

void add_vec(real *droplet_values, real *values_to_add);

// Calculate initial values for first cell in droplet
void first_value_update(real *droplet_values, real vof, cell_t cell, Thread *cell_thread);

// Calculate values for each subsequent cell in droplet
void subsequent_value_update(real *droplet_values, real vof, cell_t adjacent_cell, Thread *cell_thread);

// Calculate final values
void final_value_update(real *droplet_values);

#endif