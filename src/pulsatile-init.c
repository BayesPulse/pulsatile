/*
 * Creates entry points for .Call routines
 *   - New suggestion/request as of R 3.4
 *   - Adds layer of protection + bit of speed
 */
#include "r_interface.h"
#include <R_ext/Rdynload.h>


// Define methods called with .Call
static const
R_CallMethodDef callMethods[] = {
  {"decon_r_interface", (DL_FUNC) &decon_r_interface, 40},
  {NULL,                NULL,                         0}
};


// Initialize DLL
void 
R_init_pulsatile(DllInfo *info) {

	// Register routines (only using .Call)
  R_registerRoutines(info,
                     NULL,
                     callMethods,
                     NULL,
                     NULL);

	// Restrict to only explicitly defined names
  R_useDynamicSymbols(info, FALSE);

	// Make call available to other packages
  R_RegisterCCallable("pulsatile", "decon_r_interface", 
											(DL_FUNC) &decon_r_interface);

}


