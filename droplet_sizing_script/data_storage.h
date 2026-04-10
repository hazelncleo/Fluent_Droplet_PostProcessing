#ifndef DATA_STORAGE_H_   /* Include guard */
#define DATA_STORAGE_H_

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#include "udf.h"
#include <stdio.h>

#define SECONDARY_DENSITY 997. // water density kg/m^3

// Calculate initial values for first cell in droplet
void first_value_update(real *droplet_values, real vof, cell_t cell, Thread *cell_thread);

// Calculate values for each subsequent cell in droplet
void subsequent_value_update(real *droplet_values, real vof, cell_t adjacent_cell, Thread *cell_thread);

// Calculate final values
void final_value_update(real *droplet_values);

#endif