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
hamming_changes <- function(tvdbn.fit, time = NULL) {
  tn_list = tvdbn::all_transition_network_graph(tvdbn.fit)
  diff = c()
  if (is.null(time)) {
    for (i in 1:(length(tn_list)-1)) {
      diff = c(diff, bnlearn::hamming(tn_list[[i]],tn_list[[i+1]]))
    }
  }
  else {
    for (i in 1:(length(tn_list))) {
      diff = c(diff, bnlearn::hamming(tn_list[[time]],tn_list[[i]]))
    }
  }
  return(diff)

}



#' @title  The global graph of two Bayesian networks
#'
#' @description  Get the global graph of two Bayesian networks.
#'
#' The global graph contains just the arcs that are present in both networks. Such
#' arcs are referred as persistent arcs or persistent edges
#' @param bn1 An input Bayesian network
#' @param bn2 An input Bayesian network with the same set of vertex than `bn1`
#' @return The global graph of the two input Bayesian networks
#'
#' @export
global_graph_2g <- function(bn1, bn2) {
  # Convert to bn if the input object is type bn.fit
  if (any(class(bn1) == "bn.fit")) {
    bn1 = tvdbn::graph(bn1)
  }
  if (any(class(bn2) == "bn.fit")) {
    bn2 = tvdbn::graph(bn2)
  }
  global_graph = bn1
  for (i in 1:length(bn1$arcs[,1])) {
    #Check if the arcs is in the other network bn2
    if (! any(apply(bn2$arcs, 1, function(x,want) isTRUE(all.equal(x,want)),   bn1$arcs[i,]))) {
      # If it is not, delete that arc from the global grapj
      global_graph = bnlearn::drop.arc(global_graph, bn1$arcs[i,1],bn1$arcs[i,2])
    }
  }
  return(global_graph)
}

#' @title  The global graph of a TV-DBN or a list of BNs
#'
#' @description  Get the global graph of a time-varying dynamic Bayesian network or of a list of Bayesian networks
#'
#' The global graph contains just the arcs that are present in all the networks. Such
#' arcs are referred as persistent arcs or persistent edges
#' @param tvdbn A  time-varying DBN of type `tvdbn.fit` or `tvdbn` OR a list of Bayesian networks
#' with the same set of nodes
#' @return A (transition) Bayesian network, the global graph of the input networks
#'
#' @export
global_graph <- function(tvdbn, t_0 = NULL, t_f = NULL) {
  bn_list = tvdbn
  # If the input is a TV-DBN, we have to separate the transition networks
  if (any(class(tvdbn) == "tvdbn") || any(class(tvdbn) == "tvdbn.fit")) {
    if (is.null(t_0)) {
      t_0 = 1
    }
    if (is.null(t_f)) {
      t_f = tvdbn::get_time_points(tvdbn)
    }
    bn_list = tvdbn::all_transition_network_graph(tvdbn)[t_0:(t_f-1)]
  }
  global_graph = bn_list[[1]]
  for (i in 2:length(bn_list)) {
    global_graph = tvdbn::global_graph_2g(global_graph, bn_list[[i]])
    if (is_empty(global_graph$arcs)) {
      print(paste("Warning: From intial time ",t_0," to time ",i+1,", no persistent arc was found. Specify a lower time interval",sep=""))
      break
    }
  }
  return(global_graph)
}



#' @title  The Hamming graph of two Bayesian networks
#'
#' @description  Get the Hamming graph of two Bayesian networks.
#'
#' The Hamming graph contains the arcs of both networks that are absent in the other. Such
#' arcs are referred as contingent arcs or contingent edges
#' @param bn1 An input Bayesian network
#' @param bn2 An input Bayesian network with the same set of vertex than `bn1`
#' @return The Hamming graph of the two input Bayesian networks
#'
#' @export
hamming_graph_2g <- function(bn1, bn2) {
  # Convert to bn if the input object is type bn.fit
  if (any(class(bn1) == "bn.fit")) {
    bn1 = tvdbn::graph(bn1)
  }
  if (any(class(bn2) == "bn.fit")) {
    bn2 = tvdbn::graph(bn2)
  }
  len = length(bn1$arcs[,1])
  bn1_arcs = bn1$arcs
  for (i in 1:length(bn1_arcs[,1])) {
    #Check if the arcs is in the other network bn2
    if (any(apply(bn2$arcs, 1, function(x,want) isTRUE(all.equal(x,want)),   bn1_arcs[i,]))) {
      # If it is, delete that arc from both graph
      bn1 = bnlearn::drop.arc(bn1, bn1_arcs[i,1],bn1_arcs[i,2])
      bn2 = bnlearn::drop.arc(bn2, bn1_arcs[i,1],bn1_arcs[i,2])
    }
  }
  # Merge the two networks
  for (i in 1:length(bn2$arcs[,1])) {
    bn1 = bnlearn::set.arc(bn1,bn2$arcs[i,1],bn2$arcs[i,2])
  }
  return(bn1)
}






