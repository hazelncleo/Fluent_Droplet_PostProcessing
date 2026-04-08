#ifndef DATA_STORAGE_H_   /* Include guard */
#define DATA_STORAGE_H_
#include "udf.h"
#include <stdio.h>

#define SECONDARY_DENSITY 997. // water density kg/m^3

// Reset values
void reset(real *droplet_values, real *calc_values);

// Calculate initial values for first cell in droplet
void first_value_update(real *droplet_values, real *calc_values, cell_t cell, Thread *cell_thread);

// Calculate values for each subsequent cell in droplet
void subsequent_value_update(real *droplet_values, real *calc_values, cell_t adjacent_cell, Thread *cell_thread);

// Calculate final values
void final_value_update(real *droplet_values, real *calc_values);

#endif