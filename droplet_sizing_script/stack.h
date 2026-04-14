#ifndef STACK_H_   /* Include guard */
#define STACK_H_

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#include <stdbool.h>
#include <stdio.h>
#include "udf.h"

#define MAX_STACK_SIZE 1000000
#define MAX_N_GAPS 1000

// Define a structure for the stack
typedef struct {
    cell_t arr[MAX_STACK_SIZE];
    int top;
    int gaps[MAX_N_GAPS];
    int n_gaps;
} Stack;

// Function to initialize the stack
void initialize(Stack *stack);

// Function to check if the stack is empty
bool isEmpty(Stack *stack);

// Function to check if the stack is full
bool isFull(Stack *stack);

// Function to check if current value is a gap
bool isGap(Stack *stack);

// Function to push an element onto the stack
void push(Stack *stack, cell_t value);

// Function to add a gap to the 
void addGap(Stack *stack);

int remGap(Stack *stack);

// Function to pop an element from the stack
cell_t pop(Stack *stack);

// Function to peek the top element of the stack
cell_t peek(Stack *stack);

#endif