#ifndef VOF_DROPLET_SIZING_H_   /* Include guard */
#define VOF_DROPLET_SIZING_H_

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS

#include <stdio.h>
#include <time.h>
#include "udf.h"

#define SECONDARY_DENSITY 997. /* kg/m^3 */
#define MAX_COMBINE_DROPLETS 20000
#define MAX_STACK_SIZE 2000000
#define MAX_N_GAPS 10000
#define DROPLET_PHASE 1 /* 0 = primary, 1 = secondary, and so on... */
#define PRINT_SIMDATA 1
#define SKIP_TIMESTEP 50

#define NUM_UDM 1
static int udm_offset = UDM_UNRESERVED;


/*
 *   Structs
 */


typedef struct {
    cell_t arr[MAX_STACK_SIZE];
    int top;
    int gaps[MAX_N_GAPS];
    int n_gaps;
} cell_stack;


typedef struct {
    int droplet_ids[MAX_COMBINE_DROPLETS];
    int combination_ids[MAX_COMBINE_DROPLETS];
    int n_droplets;
    int n_to_combine;
    real droplet_values[MAX_COMBINE_DROPLETS][8];
} datastorage;


/*
 *   Main control functions
 */

void multithreaded_calculation();

void singlethreaded_calculation();

/*
 *   Compute node functions
 */

void init_udm();

void compute_droplet_data_singlethreaded();

void found_new_droplet_singlethreaded(cell_t first_cell, Thread *cell_thread, Thread *phase_thread, int droplet_id);

void compute_droplet_data(cell_stack *cells_to_reexplore);

void found_new_droplet(cell_t first_cell, Thread *cell_thread, Thread *water_thread, int droplet_id, cell_stack *cells_to_reexplore);

void assemble_droplets(cell_stack *cells_to_reexplore);

int droplet_in_array(int *attached_droplets, int n_droplets, int new_droplet_id);

void send_droplet_message(int message, int *attached_droplets, int receiving_node);

/*
 *   Stack functions
 */

void initialize_stack(cell_stack *stack_obj);

int stack_is_empty(cell_stack *stack_obj);

int stack_is_full(cell_stack *stack_obj);

int stack_is_gap(cell_stack *stack_obj);

void push_to_stack(cell_stack *stack_obj, cell_t cell);

void add_gap_to_stack(cell_stack *stack_obj);

void rem_gap_from_stack(cell_stack *stack_obj);

cell_t pop_from_stack(cell_stack *stack_obj);

/*
 *   Data storage functions
 */

void initialize_datastorage(datastorage *droplets_datastorage);

void add_droplet_values_to_datastorage(datastorage *droplets_datastorage, real *droplet_values, int droplet_id);

int get_droplet_index_from_datastorage(datastorage *droplets_datastorage, int droplet_id);

void check_droplets_already_assigned(datastorage *droplets_datastorage, int n_to_assign, int *droplets, int *reassign_combinations);

void assign_droplets_to_combine(datastorage *droplets_datastorage, int n_to_assign, int *droplets);

int get_values_from_datastorage(datastorage *droplets_datastorage, real *droplet_values, int combination_id);

void add_vectors(real *droplet_values, real *values_to_add);

void first_value_update(real *droplet_values, real vof, cell_t cell, Thread *cell_thread);

void subsequent_value_update(real *droplet_values, real vof, cell_t cell, Thread *cell_thread);

void final_value_update(real *droplet_values);

/*
 *   Node zero functions
 */

void node_zero_send_data();

void node_zero_send_droplet_connections();

/*
 *   Host node functions
 */

FILE *file_handler();

void save_line_to_file(FILE *fptr, real *droplet_values, int droplet_id);

void receive_node_data_singlethreaded(FILE *fptr);

void receive_node_zero_data(FILE *fptr, datastorage *droplets_datastorage);

void receive_compute_node_data(FILE *fptr, datastorage *droplets_datastorage);

void receive_node_zero_connections(datastorage *droplets_datastorage);

void receive_compute_node_connections(datastorage *droplets_datastorage);

void combine_boundary_droplets(FILE *fptr, datastorage *droplets_datastorage);

void host_process_singlethreaded();

void host_process_multithreaded();

#endif