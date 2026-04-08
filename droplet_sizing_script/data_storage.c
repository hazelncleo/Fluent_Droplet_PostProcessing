#include "udf.h"
#include "data_storage.h"

#define SECONDARY_DENSITY 997. // water density kg/m^3

// Reset values
void reset(real *droplet_values, real *calc_values){
    int i;
    for (i = 0; i < 7; i++){
        droplet_values[i] = 0.;
    }
    for (i = 0; i < 3; i++){
        calc_values[i] = 0.;
    }
}

// Calculate initial values for first cell in droplet
void first_value_update(real *droplet_values, real *calc_values, cell_t cell, Thread *cell_thread){

    droplet_values[0] = C_VOLUME(cell, cell_thread) * calc_values[0];
    calc_values[1] = droplet_values[0] * SECONDARY_DENSITY;

    real temp[3];
    C_CENTROID(temp, cell, cell_thread);
    droplet_values[1] = temp[0] * calc_values[1];
    droplet_values[2] = temp[1] * calc_values[1];
    droplet_values[3] = temp[2] * calc_values[1];

    droplet_values[4] = C_U(cell, cell_thread) * calc_values[1];
    droplet_values[5] = C_V(cell, cell_thread) * calc_values[1];
    droplet_values[6] = C_W(cell, cell_thread) * calc_values[1];
}

// Calculate values for each subsequent cell in droplet
void subsequent_value_update(real *droplet_values, real *calc_values, cell_t adjacent_cell, Thread *cell_thread){

    calc_values[2] = C_VOLUME(adjacent_cell, cell_thread) * calc_values[0];

    droplet_values[0] += calc_values[2];

    calc_values[2] *= SECONDARY_DENSITY;

    calc_values[1] += calc_values[2];

    real temp[3];
    C_CENTROID(temp, adjacent_cell, cell_thread);
    droplet_values[1] += temp[0] * calc_values[2];
    droplet_values[2] += temp[1] * calc_values[2];
    droplet_values[3] += temp[2] * calc_values[2];

    droplet_values[4] += C_U(adjacent_cell, cell_thread) * calc_values[2];
    droplet_values[5] += C_V(adjacent_cell, cell_thread) * calc_values[2];
    droplet_values[6] += C_W(adjacent_cell, cell_thread) * calc_values[2];
}

// Calculate final values
void final_value_update(real *droplet_values, real *calc_values){

    droplet_values[1] /= calc_values[1];
    droplet_values[2] /= calc_values[1];
    droplet_values[3] /= calc_values[1];

    droplet_values[4] /= calc_values[1];
    droplet_values[5] /= calc_values[1];
    droplet_values[6] /= calc_values[1];
}