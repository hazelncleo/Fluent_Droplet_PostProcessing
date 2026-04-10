#ifndef MULTITHREADED_H_   /* Include guard */
#define MULTITHREADED_H_

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#include <stdio.h>
#include <stdbool.h>
#include "udf.h"
#include "stack.h"
#include "data_storage.h"
#include "node_host.h"
#include "node_zero.h"

void found_new_droplet(cell_t first_cell, Thread *cell_thread, Thread *water_thread, int droplet_id);

void compute_droplet_data();

#endif
