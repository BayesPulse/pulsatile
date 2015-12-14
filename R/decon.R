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
#' Testing reading in data, converting to C data types and converting back.
#'   Test List/Dataframe to 
#' 
#' @param data Named vector of parameters
#' @keywords pulse simulation
#' @export
#' @examples none currently
#' @useDynLib pulsatile test_inout_
#' test_inout()
test_inout <- function(data) {
  
  .Call(test_inout_, data, PACKAGE = "pulsatile")

}
