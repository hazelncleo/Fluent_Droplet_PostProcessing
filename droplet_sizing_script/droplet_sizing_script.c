
#include "udf.h"
#include <stdbool.h>
#include <stdio.h>

// Modifiable parameters
#define DOMAIN_ID_CALCULATE 4
#define MAX_STACK_SIZE 2000000



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
    
    int droplet_id = 1;
    Domain *mixture_domain = Get_Domain(1);
    Domain *water_domain = DOMAIN_SUB_DOMAIN(mixture_domain,S_PHASE);
    Thread *cell_thread;
    Thread **pt;
    cell_t cell;
    real vof_water;
    real explored;
    Thread *water_thread;

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

                if (explored == 0 && vof_water > 0.5) { // Cell in droplet

                    C_UDMI(cell, cell_thread, 0) = droplet_id; 
                    push(&stack, cell);

                    int local_face_id;
                    face_t cell_face;
                    cell_t current_cell;
                    cell_t adjacent_cell_0;
                    cell_t adjacent_cell_1;

                    while (!isEmpty(&stack)){

                        current_cell = pop(&stack);

                        c_face_loop(current_cell, cell_thread, local_face_id){

                            cell_face = C_FACE(cell, cell_thread, local_face_id);
                            //adjacent_cell_0 = F_C0(cell_face, cell_thread);
                            
                            Message("Here! %d\n", cell_face);

                        }
                        
                    }

                } else { // Cell not in droplet
                    C_UDMI(cell, cell_thread, 0) = -1; 
                }
            }
            end_c_loop(cell,cell_thread)
        }
    }
}


