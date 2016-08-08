///-----------------------------------------------------------------------------
///
/// FILE: hash.h
///
/// DESCRIPTION: 
///   Function definitions for hash.c
/// 
///-----------------------------------------------------------------------------

#ifndef HASH_H
#define HASH_H

// Include this here since defn of structures used in arguments/returns
#include "deconvolution_main.h"

Node_type *initialize_node(void);

Node_type *sequential_search(Node_type *list, double time);

void insert_node(Node_type *new_node, Node_type *list);

void delete_node(Node_type *node, Node_type *list);

void print_list(Node_type *list);

void destroy_list(Node_type *list);

#endif // HASH_H

//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

