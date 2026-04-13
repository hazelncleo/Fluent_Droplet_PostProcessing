#include "stack.h"

#define MAX_STACK_SIZE 1000000

// Function to initialize the stack
void initialize(Stack *stack) {
    stack->top = -1;
    stack->n_gaps = 0;  
}

// Function to check if the stack is empty
bool isEmpty(Stack *stack) {
    return stack->top == -1;  
}

// Function to check if the stack is full
bool isFull(Stack *stack) {
    return stack->top >= MAX_STACK_SIZE - 1;  
}

bool isGap(Stack *stack) {
    if (stack->n_gaps > 0){
        return (stack->gaps[stack->n_gaps-1] == stack->top);
    } else {
        return false;
    }
}

// Function to push an element onto the stack
void push(Stack *stack, cell_t value) {
    if (isFull(stack)) {
        return;
    }
    stack->arr[++stack->top] = value;
}

void addGap(Stack *stack) {
    stack->gaps[stack->n_gaps++] = ++stack->top;
}

int remGap(Stack *stack) {
    if (isGap(stack)){
        stack->n_gaps--;
        stack->top--;
    } else {
        return -1;
    }
}

// Function to pop an element from the stack
cell_t pop(Stack *stack) {
    if (isEmpty(stack)) {
        return -1;
    } else if (isGap(stack)) {
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


