
#include "udf.h"
#include <stdbool.h>
#include <stdio.h>

// Modifiable parameters
#define MAX_STACK_SIZE 1000000
#define SECONDARY_DENSITY 997. // water density kg/m^3




// Define a structure for the stack
typedef struct {
    cell_t arr[MAX_STACK_SIZE];  
    int top;        
} Stack;

// Function to initialize the stack
void initialize(Stack *stack) {
    stack->top = -1;  
}

// Function to check if the stack is empty
bool isEmpty(Stack *stack) {
    return stack->top == -1;  
}

// Function to check if the stack is full
bool isFull(Stack *stack) {
    return stack->top >= MAX_STACK_SIZE - 1;  
}

// Function to push an element onto the stack
void push(Stack *stack, cell_t value) {
    if (isFull(stack)) {
        return;
    }
    stack->arr[++stack->top] = value;
}

// Function to pop an element from the stack
cell_t pop(Stack *stack) {
    if (isEmpty(stack)) {
        return -1;
    }

    cell_t popped = stack->arr[stack->top];
    stack->top--;
    return popped;
}

// Function to peek the top element of the stack
cell_t peek(Stack *stack) {
    if (isEmpty(stack)) {
        return -1;
    }
    return stack->arr[stack->top];
}






// Calculate droplet sizes function
DEFINE_ON_DEMAND(calculate_droplet_sizes){
    
    Domain *mixture_domain = Get_Domain(1);
    Domain *water_domain = DOMAIN_SUB_DOMAIN(mixture_domain,S_PHASE);
    Thread *cell_thread, *water_thread, *adjacent_cell_thread, *cell_face_thread;
    Thread **pt;
    real vof_water, explored, temp_volume, droplet_volume, centroid[ND_ND], centroid_temp[ND_ND], velocity[ND_ND], temp_mass, mass;
    int local_face_id, droplet_id = 1;
    face_t cell_face;
    cell_t cell, current_cell, adjacent_cell;

    FILE *fptr;

    fptr = fopen("droplets.csv", "w");

    if I_AM_NODE_ZERO_P {
        fprintf(fptr, "droplet_id,volume,ux,uy,uz,vx,vy,vz\n");// droplet_id,volume,diameter,u_x,u_y,u_z,v_x,v_y,v_z
    }

    Stack stack;
    initialize(&stack);

    mp_thread_loop_c(cell_thread, mixture_domain, pt){
        if (FLUID_THREAD_P(cell_thread)){

            // Initially all cells unexplored (UDM=0)
            begin_c_loop(cell, cell_thread){
                C_UDMI(cell, cell_thread, 0) = 0;
            }
            end_c_loop(cell, cell_thread)

            water_thread = pt[S_PHASE];

            begin_c_loop(cell, cell_thread){

                vof_water = C_VOF(cell, water_thread);
                explored = C_UDMI(cell, cell_thread, 0);

                if (explored == 0 && vof_water > 0.5) { // Cell in new droplet

                    C_UDMI(cell, cell_thread, 0) = droplet_id; // Cell explored
                    push(&stack, cell); // Push to the stack
                    
                    droplet_volume = C_VOLUME(cell, cell_thread)*vof_water;
                    mass = droplet_volume*SECONDARY_DENSITY; 
                    C_CENTROID(centroid, cell, cell_thread);
                    centroid[0] *= mass;
                    centroid[1] *= mass;
                    centroid[2] *= mass;
                    velocity[0] = C_U(cell, cell_thread) * mass;
                    velocity[1] = C_V(cell, cell_thread) * mass;
                    velocity[2] = C_W(cell, cell_thread) * mass;

                    while (!isEmpty(&stack)) {

                        current_cell = pop(&stack);

                        // Loop over faces of current cell
                        c_face_loop(current_cell, cell_thread, local_face_id){

                            cell_face = C_FACE(current_cell, cell_thread, local_face_id);
                            cell_face_thread = C_FACE_THREAD(current_cell, cell_thread, local_face_id);

                            // If the face is not on the boundary
                            if (!BOUNDARY_FACE_THREAD_P(cell_face_thread)){

                                // Get the new adjacent cell
                                adjacent_cell = F_C0(cell_face, cell_face_thread);
                                if (adjacent_cell == current_cell) {
                                    adjacent_cell = F_C1(cell_face, cell_face_thread);
                                }

                                vof_water = C_VOF(adjacent_cell, water_thread);
                                explored = C_UDMI(adjacent_cell, cell_thread, 0);
                                
                                // If new cell has not been explored and is droplet
                                if (explored == 0 && vof_water > 0.5) {

                                    push(&stack, adjacent_cell);
                                    C_UDMI(adjacent_cell, cell_thread, 0) = droplet_id;

                                    temp_volume = C_VOLUME(adjacent_cell, cell_thread) * vof_water;
                                    droplet_volume += temp_volume;
                                    temp_mass = temp_volume * SECONDARY_DENSITY;
                                    mass += temp_mass;
                                    C_CENTROID(centroid_temp, adjacent_cell, cell_thread);
                                    centroid[0] += centroid_temp[0]*temp_mass;
                                    centroid[1] += centroid_temp[1]*temp_mass;
                                    centroid[2] += centroid_temp[2]*temp_mass;
                                    velocity[0] += C_U(adjacent_cell, cell_thread) * temp_mass;
                                    velocity[1] += C_V(adjacent_cell, cell_thread) * temp_mass;
                                    velocity[2] += C_W(adjacent_cell, cell_thread) * temp_mass;

                                } else if (explored == 0) {
                                    C_UDMI(adjacent_cell, cell_thread, 0) = -1;
                                }
                            } 
                        }
                    }

                    centroid[0] /= mass;
                    centroid[1] /= mass;
                    centroid[2] /= mass;
                    velocity[0] /= mass;
                    velocity[1] /= mass;
                    velocity[2] /= mass;

                    if I_AM_NODE_ZERO_P {
                        fprintf(fptr, "%d,%e,%e,%e,%e,%e,%e,%e\n",droplet_id,droplet_volume,centroid[0],centroid[1],centroid[2],velocity[0],velocity[1],velocity[2]);
                    }
                    droplet_id++;

                } else if (explored == 0) { // Cell not in droplet
                    C_UDMI(cell, cell_thread, 0) = -1; 
                }
            }
            end_c_loop(cell,cell_thread)
        }
    }
    if I_AM_NODE_ZERO_P {
        Message("\n------------------\nCalculation completed! (number of droplets = %d)\n------------------\n", droplet_id);
    }
}


