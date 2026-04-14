#ifndef MULTITHREADED_H_   /* Include guard */
#define MULTITHREADED_H_

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include "udf.h"
#include "stack.h"
#include "data_storage.h"
#include "node_host.h"
#include "node_zero.h"

#define MAX_DROPLET_COMBINES 250


void found_new_droplet(cell_t first_cell, Thread *cell_thread, Thread *water_thread, int droplet_id, Stack *cells_to_reexplore);

void compute_droplet_data(Stack *cells_to_reexplore);

void init_udm();

bool unique(int *attached_droplets, int n_droplets, int new_droplet_id);

void send_droplet_message(int message, int *attached_droplets, int receiving_node);

void assemble_droplets(Stack *cells_to_reexplore);

#endif
