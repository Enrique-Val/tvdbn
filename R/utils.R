#' @title  Construct the time name of a variable
#'
#' @description  Construct the time name of a variable (i.e. the name of the variable with the time point correctly appended)
#'
#' @param name Name of the variable
#' @param t Time instant (integer)
#' @return The name of the variable with the time point correctly appended.
#'
#' @export
time_name <- function(name, t, len = NULL) {
  if (! is.null(len)) {
    # Type checking
    assertthat::assert_that(is.numeric(len))
    assertthat::assert_that(len > t)
    assertthat::assert_that(t >= 0)
    # Compute the number of zeros to append. be careful if t=0
    zeros = 0
    if (t == 0) {
      zeros = trunc(log10(len))
    }
    else {
      zeros = trunc(log10(len)) - trunc(log10(t))
    }
    if (zeros > 0) {
      for (i in 1:zeros) {
        t = paste(0,t,sep="")
      }
    }
  }
  return(paste(name,t,sep="_t_"))
}



#' @title  Get the name of a variable and the time point from a time name
#'
#' @description  Get the name of a variable and the time point from a time name.
#' @param time_name The name of a variable with a time point correctly appended.
#' Be careful to supply a string with the correct format
#' @return A vector with two elements:
#'
#'    1) The name of the variable (without the time appended)
#'
#'    2) The time point (be careful, since it the type is `character`)
#'
#' @export
remove_time_name <- function(time_name) {
  split = strsplit(time_name, "_")[[1]]
  var_name = paste(head(split,length(split)-2), collapse="_")
  time = split[length(split)]
  return(c(var_name,time))
}

#' @title  Nodes of a Time-varying DBN
#'
#' @description  Get the nodes of a time-varying DBN
#' @param tvdbn A TV-DBN, fitted or just the structure
#' @return The nodes of the network
#'
#' @export
nodes <- function(tvdbn) {
  class(tvdbn) = class(tvdbn)[-1]
  return(bnlearn::nodes(tvdbn))
}




#' @title  Structure of a node subset of a Time-varying DBN
#'
#' @description  Get the structure of a Time-varying DBN given a fitted Time-varying DBN,
#'  considering only a subset of nodes of the network
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`
#' @param nodes The nodes whose structure will be obtained.
#' @return The structure of the TV-DBN, type = `tvdbn`
#'
#' @export
subgraph <- function(tvdbn.fit, nodes) {
  dag = bnlearn::subgraph(x = tvdbn.fit, nodes = nodes)
  return(bn_to_tvdbn(dag))
}


#' @title  Structure of a Time-varying DBN
#'
#' @description  Get the structure of a Time-varying DBN given a fitted Time-varying DBN
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`
#' @return The structure of the TV-DBN, type = `tvdbn`
#'
#' @export
graph <- function(tvdbn.fit) {
  return(tvdbn::subgraph(tvdbn.fit = tvdbn.fit, nodes = nodes(tvdbn.fit)))
}


#' @title  Transition network in a time point (TODO)
#'
#' @description  Get the transition network of a certain time instant.
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`
#' @param time The transition network to obtain
#' @return A transition network with the marginal and conditional probabilities for each node
#' (at this moment in the implementation, it only return the conditional probabilities)
#' @import purrr
#'
#' @export
transition_network <- function(tvdbn.fit, time) {
  variables = get_variables(tvdbn.fit)
  # First, we store the conditional probabilities of nodes of the transition network
  prob_list = list()
  for (i in variables) {
    prob_list[[time_name(i,time, get_time_points(tvdbn.fit))]] = tvdbn.fit[[time_name(i,time,get_time_points(tvdbn.fit))]]
  }
  # Next, we compute the marginal probabilities of the parent nodes
  # TODO
  # class(prob_list) = c("tvdbn.fit","bn.fit","bn.fit.gnet")
  return(prob_list)
}


#' @title  Structure of the transition network in a time point
#'
#' @description  Get the transition network of a certain time instant.
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`.
#' @param time The transition network to obtain.
#' @param normalize A flag that indicates if we normalize the time points obtained (starting
#' from zero) or we leave them as they are.
#' @return The structure of the transition network of the TV-DBN in
#' the given time instant.
#' @import purrr
#'
#' @export
transition_network_graph <- function(tvdbn.fit, time, normalize = TRUE) {
  variables = get_variables(tvdbn.fit)
  # We obtain the variables of the time of our interest
  var_list = unlist(map(variables, time_name, time, get_time_points(tvdbn.fit)))
  # It said time is higher than 0, then we also consider the previous time slice
  # (just the previous if we assume Markov order 1. In the future, we will expand for
  # higher Markovian orders)
  if (time > 0) {
    var_list = c(var_list,unlist(map(variables, time_name, time-1, get_time_points(tvdbn.fit))))
  }
  # return the subgraph with the nodes of interest (selected and previous time slices)
  to_ret = tvdbn::subgraph(tvdbn.fit = tvdbn.fit, nodes = sort(var_list))
  # We might want to normalize the result
  if (normalize) {
    to_ret = tvdbn::normalize_time(to_ret)
  }
  return(to_ret)
}


#' @title  Markovian order of a transition network
#'
#' @description  Get the Markovian order of a transition network. In
#' other words, the number of past instant from which the transition
#' network model is dependent.
#' @param trans_network A Time-varying DBN structure, type `tvdbn`.
#' @return The Markovian order of the transition network.
#'
#' @export
transition_network_markovian_order <- function(trans_network) {
  # Get the Markovian Order of the Transition Network
  # For the time being, it does exactly the same thing as the function "get_time_points", as we
  # are only considering first Markov order relations
  order = length(trans_network$nodes)/length(get_variables(trans_network)) -1
  return(order)
}


#' @title  Normalize the name of the time variables of a TN
#'
#' @description  Normalize the name of the time variables of a Transition Network.
#' This allows comparing them using the Hamming distance
#' @param trans_network A fitted Time-varying DBN of type `tvdbn`.
#' @return The transition network with the naming normalized
#'
#' @export
normalize_time <- function(trans_network) {
  return(tvdbn::change_time(trans_network, ini = 0))
}

#' @title  List of the structure of all the transition networks
#'
#' @description  Get all the transition network structures of the Time-varying DBN
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`.
#' @param normalize A flag that indicates if we should normalize the time instants in
#' the obtained transition networks.
#' @return A list `tvdbn` containing the structure (type `tvdbn`) of all the
#'  transition networks
#'
#' @export
all_transition_network_graph <- function(tvdbn.fit, normalize = TRUE) {
  tvdbn = tvdbn::graph(tvdbn.fit)
  tn_list = list()
  time_instants = tvdbn::get_time_points(tvdbn.fit)
  if (normalize) {
    for (i in 1:(time_instants-1)) {
      tn_list[[i]] = tvdbn::normalize_time(tvdbn::transition_network_graph(tvdbn, i))
    }
  }
  else {
    for (i in 1:(time_instants-1)) {
      tn_list[[i]] = tvdbn::transition_network_graph(tvdbn, i)
    }
  }

  return(tn_list)
}


#' @title  List of the structure of all the transition networks
#'
#' @description  Get all the transition network structures of the Time-varying DBN
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`.
#' @param normalize A flag that indicates if we should normalize the time instants in
#' the obtained transition networks.
#' @return A list `tvdbn` containing the structure (type `tvdbn`) of all the
#'  transition networks
#'
#' @export
append_networks <- function(tvdbn1, tvdbn2) {
  assertthat::assert_that("tvdbn" %in% class(tvdbn1) && "tvdbn" %in% class(tvdbn2))
  # First, we normalize the time names
  tvdbn1 = normalize_time(tvdbn1)
  tvdbn2 = change_time(tvdbn2, ini = get_time_points(tvdbn1)-1)

  tvdbn2_nodes_toadd = nodes(tvdbn2)
  to_del = c()

  # Eliminate from the list of nodes of tvdbn2 all the nodes of the first time point
  for (i in 1:length(get_variables(tvdbn2))-1) {
    to_del = c(to_del,nodes(tvdbn2)[1+i*get_time_points(tvdbn2)])
  }
  tvdbn2_nodes_toadd = tvdbn2_nodes_toadd[ ! tvdbn2_nodes_toadd %in% to_del]

  ## Add all the nodes and arcs of tvdbn2 to tvdbn1
  for (i in tvdbn2_nodes_toadd) {
    tvdbn1 = add.node(tvdbn1, i)
  }

  if (length(tvdbn2$arcs[,1]) > 0) {
    for (i in 1:length(tvdbn2$arcs[,1])) {
      tvdbn1 = bnlearn::set.arc(tvdbn1, tvdbn2$arcs[i,1], tvdbn2$arcs[i,2])
    }
  }

  return(tvdbn1)

}



#' @title Change the initial time of a TV-DBN
#'
#' @description  Change the initial time of a TV-DBN. Be EXTREMELY careful, as the function may interfere
#' with other functionality, such as visualization. The best choice is to always have a normalized TV-DBN
#' (starting at time 0)
#' @param trans_network A Time-varying DBN structure of type `tvdbn`.
#' @return The TV-DBN with the time points shifted
#'
#' @export
change_time <- function(trans_network, ini = 0) {
  # Get the number of time instants of the network
  n_time_points = get_time_points(trans_network)
  # Time points is a list with all the time instant of the transition. For instance, it may range from 3 to 5
  time_points = c()
  sorted_list = sort(names(trans_network$nodes))
  for (i in sorted_list[1:(n_time_points)]) {
    time_points = c(time_points,as.numeric(remove_time_name(i)[2]))
  }
  time_points = sort(time_points)
  # bnlearn has a bug and does not let us use the function nodes with inherited class
  # We cast to class "bn" and later we will revert this change
  class(trans_network) = "bn"


  # Number of variables
  n_vars = length(bnlearn::nodes(trans_network))/(n_time_points)
  time_nodes = bnlearn::nodes(trans_network)
  new_nodes = c()

  for (i in time_nodes) {
    # Get the old time instant
    decomposed = remove_time_name(i)
    old_tp = as.numeric(decomposed[2])
    # Build the new time instant
    new_tp = match(old_tp, time_points)-1
    var_name = decomposed[1]
    new_name = time_name(var_name,new_tp+ini, len = n_time_points+ini)
    new_nodes = c(new_nodes, new_name)
  }
  trans_network = bnlearn::rename.nodes(trans_network, new_nodes)

  class(trans_network) = c("tvdbn",class(trans_network))
  return(trans_network)

}

remove_zeros <- function(tvdbn) {
  tmp = class(tvdbn)
  class(tvdbn) = "bn"
  for (i in bnlearn::nodes(tvdbn)) {
    name = remove_time_name(i)
    bnlearn::nodes(tvdbn)[match(i,bnlearn::nodes(tvdbn))] = time_name(name[1],as.integer(name[2]))
    print(time_name(name[1],as.integer(name[2])))
  }
  class(tvdbn) = tmp
  return(tvdbn)
}

remove.node <- function(x,node) {
  assertthat::assert_that("tvdbn" %in% class(x) || "tvdbn.fit" %in% class(x))
  tmp = NULL
  if ("tvdbn" %in% class(x)) {
    assertthat::assert_that("bn" %in% class(x))
    tmp = class(x)
  }
  else if ("tvdbn.fit" %in% class(x)) {
    assertthat::assert_that("bn.fit" %in% class(x))
    tmp = class(x)
  }
  x = bnlearn::remove.node(x,node)
  class(x) = tmp
  return(x)
}

sort_nodes <- function(nodes, order = "variable") {

}
