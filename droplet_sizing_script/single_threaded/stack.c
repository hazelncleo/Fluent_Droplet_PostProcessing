#include <stdbool.h>
#include <stdio.h>
#include "udf.h"
#include "stack.h"

#define MAX_STACK_SIZE 1000000

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


