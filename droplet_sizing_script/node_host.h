#ifndef NODE_HOST_H_   /* Include guard */
#define NODE_HOST_H_

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#include <stdio.h>
#include <stdbool.h>
#include "udf.h"

FILE *file_handler();

void save_line_to_file(FILE *fptr, real *droplet_values, int droplet_id);

int receive_node_zero_data(FILE *fptr);

int receive_compute_node_data(FILE *fptr, int n_droplets_outside);

void combine_boundary_droplets(FILE *fptr, int n_droplets_outside);

void host_process();

#endif