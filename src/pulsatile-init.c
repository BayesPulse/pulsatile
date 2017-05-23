/*
 * Creates entry points for .Call routines
 *   - New suggestion/request as of R 3.4
 *   - Adds layer of protection + bit of speed
 */
#include "r_interface.h"
#include <R_ext/Rdynload.h>

// Hadley's R packages approach
void R_init_pulsatile(DllInfo *info) {
    R_RegisterCCallable(info, "decon_r_interface", 
                        (DL_FUNC) &decon_r_interface);
}

