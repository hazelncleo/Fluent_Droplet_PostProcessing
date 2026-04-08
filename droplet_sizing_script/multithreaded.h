#ifndef MULTITHREADED_H_   /* Include guard */
#define MULTITHREADED_H_

#include <stdio.h>
#include "udf.h"

void calculate_droplets();

FILE *file_handler();

void node_zero_send_data();

void compute_droplet_data();

void host_write_to_file();

#endif
