# In this file, we put all the eXplaniable Artificial Intelligence (XAI) related methods.

#' @title  The hamming distances accross the transition networks.
#'
#' @description  Get the difference between all the adyacent transition networks
#' of a given Time-varying Dynamic Bayesian Network
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`.
#' @return A vector with the Hamming distance between adyadcents transition networks.
#' Each position i represents the Hamming distance between the `i-th` and `i+1-th`
#' transition network.
#'
#' @export
hamming_changes <- function(tvdbn.fit) {
  tn_list = tvdbn::all_transition_network_graph(tvdbn.fit)
  diff = c()
  for (i in 1:(length(tn_list)-1)) {
    diff = c(diff, bnlearn::hamming(tn_list[[i]],tn_list[[i+1]]))
  }

  return(diff)

}
