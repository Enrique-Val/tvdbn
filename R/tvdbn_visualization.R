library(dbnR)

#' @title Plots a dynamic Bayesian network in a hierarchical way
#'
#' @description Plot an Time-Vaying Dynamic Bayesian Network in an ordered manner.
#' The nodes are horizontally ordered in time slices, where the leftmost slice is
#' the oldest and the rightmost is the present. The nodes are ordered alphabetically
#' in each time slice, top to down.
#' @details This is function is a copy-paste of the function `dbnR::plot_dynamic_network()`,
#' from the package `dbnR` with minor modifications. The original author,
#'  *David Quesada*, deserves all the credit for this function.
#' @param structure the structure or fit of the network.
#' @param offset the blank space between time slices
#' @return the visualization of the DBN
#' @import dbnR
#' @export
plot_time_varying_network <- function(structure, offset = 200){
  structure = remove_zeros(structure)
  if (class(structure)[1] == "tvdbn.fit") {
    class(structure)[1] = "dbn.fit"
  }
  else if (class(structure)[1] == "tvdbn") {
    class(structure)[1] = "dbn"
  }
  ret = dbnR::plot_dynamic_network(structure = structure, offset = offset)
  ret$x$nodes["x"] = rev(unlist(ret$x$nodes["x"]))
  return(ret)
}

