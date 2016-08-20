//
// pulse_node.c
//
// this module implements the linked list.  in this implementation, the first 
// node in the linked list is a dummy node.  one must initialize the list with
// the dummy node as the first.                                               
//
// Subroutines that exist in this program
//      *initialize_node
//      *sequential_search
//      insert_node
//      delete_node
//      print_list
//      destroy_list
//

#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include "pulse_node.h"

//#include "r_interface.h"
extern double fitstart;


//
// initialize_node:
//      this allocates memory for a new pulse and gives it some initial
//      parameters
//
//      ARGUMENTS: 
//          None
//
//      RETURNS: 
//          p    - the created pulse
//
// Variable definitions
//      i  - generic counter
//      *p - initialized pulse
//
Node_type *initialize_node(void) {

  int i;
  Node_type *p;

  /* initialize memory for a node */
  if ((p = (Node_type *)malloc(sizeof(Node_type))) == NULL) {
    Rprintf("Initialize_node: out of memory, can't allocate node\n");
    exit(1);
  }

  /* all nodes initialized with a pointer to the NULL pointer and time = 0,
   * the parameters all initialized to zero it is up to the programmer to
   * insure that the correct values are inserted into the node at the
   * appropriate time */

  p->succ = p->pred = NULL;
  p->mean_contrib   = NULL;
  p->time           = fitstart;
  for (i = 0; i < 2; i++) {
    p->theta[i] = 0;
  }

  return p;

}




//
// sequential_search: 
//      finds where a pulse's location fits within an established linked list
//
//      ARGUMENTS: 
//          Node_type *list - this is the current list of pulses that exist
//          double time     - this is the location of a new pulse, note that we
//                            do not actually input a new pulse, we just input
//                            its location
//      RETURNS: 
//          *loc            - the newly created pulse's prececessor
//
// Variable definitions
//      loc   - the newly created pulse's predecessor
//      locm1 - the newly created pulse's successor
//
Node_type *sequential_search(Node_type *list, double time) {

  Node_type *loc, *locm1;

  /* searches through the linked list "list" for the first occassion of
   * "list->time" that is less than "time".  The new node will be inserted
   * between "locm1" and "loc".  "locm1"  preceeds "loc" in "list" by one
   * position */

  locm1 = list;
  for (loc=list;loc;loc=loc->succ){
    if (loc->time >= time) {
      return locm1;
    } else {
      locm1 = loc;
    }
  }

  return locm1;

}




//
// insert_node: 
//      integrates a newly created pulse into the linked list
//
//      ARGUMENTS: 
//          Node_type *new_node - the newly created pulse
//          Node_type *list     - this is the current list of pulses that exist
//
//      RETURNS: 
//          None                - all updates are made internally
//
//
// Variable definitions
//      node              - new_node's predecessor
// 
// Subroutines used
//      sequential_search - found in this file; identifies inputted pulse's
//                          predecessor
//
void insert_node(Node_type *new_node, Node_type *list) {

  Node_type *node;

  /* insert the node "new_node" in the linked list "list" find the position in
   * which to insert "new_node" */

  node = sequential_search(list, new_node->time);
  /********************************************************/

  /* reassign pointers so that everyone is pointing to the correct node */

  new_node->pred = node;
  new_node->succ = node->succ;
  if (node->succ != NULL) {
    (node->succ)->pred = new_node;
  }
  node->succ = new_node;

}




//
// delete_node: 
//      frees memory associated with inputted pulse and reassigns pointers of
//      remaining pulses
//
//      ARGUMENTS: 
//          Node_type *node - the pulse to be deleted
//          Node_type *list - this is the current list of pulses that exist
//
//      RETURNS: 
//          None            - all updates are made internally
//
// Variable definitions
//      None
// 
// Subroutines used
//      None
//
void delete_node(Node_type *node, Node_type *list) {

  free(node->mean_contrib);

  if (node->succ == NULL) {
    (node->pred)->succ = node->succ;
  } else {
    (node->pred)->succ = node->succ;
    (node->succ)->pred = node->pred;
  }

  free(node);

}




//
// print_list: 
//      prints information about all pulses in the linked list
//
//      ARGUMENTS: 
//          Node_type *list - this is the current list of pulses that exist;
//
//      RETURNS: 
//          None            - all updates are made internally
//
// Variable definitions
//      i     - generic counter
//      *node - counter through pulses
// 
// Subroutines used
//      None
//
void print_list(Node_type *list) {

  int i;
  Node_type *node;

  /* traverses the linked list and prints the contents of each node */
  i = 1;
  node = list->succ;
  Rprintf("Pulse No. Time  Mass  Width\n");
  while (node != NULL) {
    Rprintf("%2d %8.2lf %8.2lf %8.2lf \n", i, node->time, node->theta[0],
           node->theta[1]);
    node = node->succ;
    i++;
  }
}




// 
// destroy_list: 
//      frees memory throughout the linked list
//
//      ARGUMENTS: 
//          Node_type *list - this is the current list of pulses that exist;
//          
//      RETURNS: 
//          None            - all updates are made internally
//
// Variable definitions
//      loc  - counter through pulses
// 
void destroy_list(Node_type *list) {

  Node_type *loc;

  while (list->succ != NULL) {
    loc = list->succ;
    free(list->mean_contrib);
    free(list);
    list = loc;
  }

  free(list->mean_contrib);
  free(list);

}



//******************************************************************************
// END OF FILE
//******************************************************************************

