#ifndef STACK_H_   /* Include guard */
#define STACK_H_

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#include "udf.h"
#include <stdbool.h>

#define MAX_STACK_SIZE 1000000

// Define a structure for the stack
typedef struct {
    cell_t arr[MAX_STACK_SIZE];  
    int top;        
} Stack;

// Function to initialize the stack
void initialize(Stack *stack);

// Function to check if the stack is empty
bool isEmpty(Stack *stack);

// Function to check if the stack is full
bool isFull(Stack *stack);

// Function to push an element onto the stack
void push(Stack *stack, cell_t value);

// Function to pop an element from the stack
cell_t pop(Stack *stack);

// Function to peek the top element of the stack
cell_t peek(Stack *stack);

#endif