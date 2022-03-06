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

  if (length(bn1$arcs) == 0) {
    return(bn1)
  }
  if (length(bn2$arcs) == 0) {
    return(bn2)
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
#' @param t_0 Start computing the global graph at the instant t_0. (Use only if the
#' type of `tvdbn` is `tvdbn` or `tvdbn.fit`. With list, the paramter becomes useless)
#' @param t_f Finish computing the global graph at the instant t_f. (Use only if the
#' type of `tvdbn` is `tvdbn` or `tvdbn.fit`. With list, the paramter becomes useless)
#' @param threshold A threshold of tolerance. This allows for adding arcs to the global arc
#' that were missing on a percentage of the transition networks. The input should be a real number between
#' 0 and 1. 0 means that all arcs of the transition networks will be added to the global_graph (
#' very relaxed threshold) and 1 means that only the arcs present in every single transition network
#' will be added to the global graph (very strict threshold).
#' @return A (transition) Bayesian network, the global graph of the input networks
#'
#' @export
global_graph <- function(tvdbn, t_0 = NULL, t_f = NULL, threshold = 1) {
  assertthat::assert_that(0 <= threshold & 1 >= threshold)
  t0 = Sys.time()
  bn_list = tvdbn
  arc_matrix = NULL
  n_trans_network = NULL
  vars = NULL
  # If the input is a TV-DBN
  if (any(class(tvdbn) == "tvdbn") || any(class(tvdbn) == "tvdbn.fit")) {
    if (is.null(t_0)) {
      t_0 = 0
    }
    if (is.null(t_f)) {
      t_f = tvdbn::get_time_points(tvdbn)-1
    }
    # Correct t_0 and t_f
    t_0 = max(0,t_0)
    t_f = min(tvdbn::get_time_points(tvdbn)-1, t_f)
    vars = get_variables(tvdbn)
    n_vars = length(vars)
    n_trans_network = t_f-t_0
    arcs = arcs(tvdbn)
    arc_matrix = matrix(0,nrow = n_vars, ncol = n_vars, dimnames = list(vars,vars) )
    for (i in 1:length(arcs(tvdbn)[,1])) {
      t = as.integer(remove_time_name(arcs[i,1])[2])
      if (t >= t_0 & t < t_f) {
        from = remove_time_name(arcs[i,1])[1]
        to = remove_time_name(arcs[i,2])[1]
        arc_matrix[from,to] = arc_matrix[from,to] + 1
      }
    }

    tn = Sys.time()
    print(tn-t0)


  }
  else {
    # If the input is a transition network list
    vars = get_variables(bn_list[[1]])
    n_vars = length(vars)
    n_trans_network = length(bn_list)
    # A matrix in which we will count the arcs that appears in every transition network.
    # The position (i,j) means an arc from the variable i to the variable j
    arc_matrix = matrix(0,nrow = n_vars, ncol = n_vars, dimnames = list(vars,vars) )
    for (i in 1:length(bn_list)) {
      arcs_i = arcs(bn_list[[i]])
      if (length(arcs_i) > 0) {
        for (j in 1:nrow(arcs_i)) {
          from = remove_time_name(arcs_i[j,1])[1]
          to = remove_time_name(arcs_i[j,2])[1]
          arc_matrix[from,to] = arc_matrix[from,to] + 1
        }
      }
    }
    tn = Sys.time()
    print(tn-t0)
  }

  # Now we have the arc_matrix and we are going to build the global_graph with said matrix
  # First, we threshold the matrix so that it becomes and adjacency matrix
  arc_matrix = arc_matrix >= threshold*n_trans_network
  rownames(arc_matrix) = sapply(rownames(arc_matrix), time_name,0)
  colnames(arc_matrix) = sapply(colnames(arc_matrix), time_name,1)
  global_graph = bnlearn::empty.graph(c(rownames(arc_matrix), colnames(arc_matrix)))
  arc_set = matrix(ncol=2, dimnames = list(NULL, c("from", "to")))[-1,]
  for (i in rownames(arc_matrix)) {
    for (j in colnames(arc_matrix)) {
      if (arc_matrix[i,j])
        arc_set = rbind(arc_set, c(i,j))
    }
  }
  arcs(global_graph) = arc_set
  class(global_graph) = c("tvdbn",class(global_graph))
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





# summary_network <- function(tvdbn, frequency) {
#   t0 = Sys.time()
#   assertthat::assert_that("tvdbn" %in% class(tvdbn) || "tvdbn.fit" %in% class(tvdbn))
#   assertthat::assert_that(frequency > 1)
#   assertthat::assert_that(get_time_points(tvdbn) > frequency)
#   if ("tvdbn.fit" %in% class(tvdbn)) {
#     tvdbn = graph(tvdbn)
#   }
#   gg = list()
#   tns = all_transition_network_graph(tvdbn)
#   for (i in 1:(length(tns)%/%frequency)) {
#     gg[[i]] = global_graph(tns[(1+(i-1)*frequency):(i*frequency)])
#   }
#   if (length(tns) %% frequency != 0) {
#     gg[[round(length(tns)/frequency)]] =global_graph(tns[(1+length(tns)-(length(tns) %% frequency)):length(tns)])
#   }
#   summary = gg[[1]]
#   for (i in 2:length(gg)) {
#     summary = append_networks(summary, gg[[i]])
#   }
#   tn = Sys.time()
#   print(tn-t0)
#   return(summary)
# }

summary_network <- function(tvdbn, frequency, threshold=1) {
  t0 = Sys.time()
  assertthat::assert_that("tvdbn" %in% class(tvdbn) || "tvdbn.fit" %in% class(tvdbn))
  assertthat::assert_that(frequency > 1)
  assertthat::assert_that(get_time_points(tvdbn)+1 > frequency)

  tn_number = get_time_points(tvdbn)-1
  gg = list()
  for (i in 1:round(tn_number/frequency)) {
    gg[[i]] = global_graph(tvdbn,(i-1)*frequency,i*frequency, threshold = threshold)
  }
  summary = gg[[1]]
  for (i in 2:length(gg)) {
    summary = append_networks(summary, gg[[i]])
  }
  tn = Sys.time()
  print(tn-t0)
  return(summary)
}


#' @title  Graph that contains every arc of the transition networks of a TV-DBN or a list of BNs
#'
#' @description  Get a graph that contains every arc of the transition networks of a time-varying
#'  dynamic Bayesian network or of a list of Bayesian networks with the same set of nodes
#'
#' @param tvdbn A  time-varying DBN of type `tvdbn.fit` or `tvdbn` OR a list of Bayesian networks
#' with the same set of nodes
#' @param t_0 Start computing the global graph at the instant t_0. (Use only if the
#' type of `tvdbn` is `tvdbn` or `tvdbn.fit`. With list, the paramter becomes useless)
#' @param t_f Finish computing the global graph at the instant t_f. (Use only if the
#' type of `tvdbn` is `tvdbn` or `tvdbn.fit`. With list, the paramter becomes useless)
#' @return A (transition) Bayesian network, the all-arc graph of the input networks
#'
#' @export
all_arc_graph <- function(tvdbn, t_0 = NULL, t_f = NULL) {
  return(tvdbn::global_graph(tvdbn = tvdbn, t_0 = t_0, t_f = t_f, threshold = 0))
}

