//
// FILE: pulse_node.h
//
// DESCRIPTION: 
//   Function definitions for pulse_node.c
// 

#ifndef PULSE_NODE_H
#define PULSE_NODE_H

typedef struct node_tag {

    struct node_tag *succ; // next node
    struct node_tag *pred; // previous node
    double time;           // pulse location for individual pulses
    double theta[2];       // 0: individual pulse mass; 1: individual pulse variance
    double *mean_contrib;  // vector of secretion contribution for this pulse at
                           // each timepoint
    double eta[2];         // scaling parameter distr. gamma(4/2, 4/2), for
                           // making the truncated t-distribution prior on
                           // individual pulse mass/width
} Node_type;


Node_type *initialize_node(void);

Node_type *sequential_search(Node_type *list, double time);

void insert_node(Node_type *new_node, Node_type *list);

void delete_node(Node_type *node, Node_type *list);

void print_list(Node_type *list);

void destroy_list(Node_type *list);

#endif 


//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

