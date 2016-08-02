#-------------------------------------------------------------------------------
# decon.R - Wrapper for fitting decon
#-------------------------------------------------------------------------------

#decon <- function(data, priors, starting.vals, proposal.vars, seeds) {
#
#  .Call(decon_, data, priors, starting.vals, proposal.vars, seeds)
#
#}


#' test_inout
#'
#' testing reading in data, converting to c data types and converting back.
#'   test list/dataframe to 
#' 
#' @param data named vector of parameters
#' @keywords pulse simulation
#' @export
#' @useDynLib pulsatile testc
#' @useDynLib pulsatile testspec
#' @examples none currently
#' test_inout()
test_inout <- function(x, ...) {
  
  if ("pulse_spec" %in% class(x)) { 
    .Call(testspec, x, 
          PACKAGE = "pulsatile")

  } else if (is.data.frame(x)) {
    .Call(testc, x, 
          PACKAGE = "pulsatile")
  }

}


#' show_args
#'
#' Example via R source code: http://git.io/v0zFx
#' 
#' @useDynLib pulsatile showArgs1
#' @param data named vector of parameters
#' @keywords pulse simulation
#' @export
#' @examples none currently
#' show_args()
show_args <- function(.data) {
  
  .Call(showArgs1, .data, PACKAGE = "pulsatile")

}

