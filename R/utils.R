#' @title  Construct the time name of a variable
#'
#' @description  Construct the time name of a variable (i.e. the name of the variable with the time point correctly appended)
#'
#' @param name Name of the variable
#' @param t Time instant (integer)
#' @param len The length of the dataset or of the TVDBN (number of slices). This parameter is optional and will append as many left
#' zeros as necessary so that all the time stamps have the same length
#' @return The name of the variable with the time point correctly appended.
#'
#' @export
time_name <- function(name, t, len = NULL) {
  if (! is.null(len)) {
    # Type checking
    assertthat::assert_that(is.numeric(len))
    assertthat::assert_that(len > 0)
    assertthat::assert_that(t >= 0)
    assertthat::assert_that(len > t)
    # Compute the number of zeros to append. be careful if t=0
    zeros = 0
    if (t == 0) {
      zeros = trunc(log10(len-1))
    }
    else {
      zeros = trunc(log10(len-1)) - trunc(log10(t))
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
  var_list = NULL
  network_length = NULL
  if (is_padded(tvdbn.fit)) {
    network_length = get_time_points(tvdbn.fit)
  }

  var_list = unlist(map(variables, time_name, time, network_length))

  # It said time is higher than 0, then we also consider the previous time slice
  # (just the previous if we assume Markov order 1. In the future, we will expand for
  # higher Markovian orders)
  if (time > 0) {
    var_list = c(var_list,unlist(map(variables, time_name, time-1, network_length)))
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
  for (i in 1:(time_instants-1)) {
    tn_list[[i]] = tvdbn::normalize_time(tvdbn::transition_network_graph(tvdbn, i, normalize = normalize))
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

  tvdbn2_nodes_toadd = nodes(tvdbn2)[as.integer(sapply(nodes(tvdbn2),remove_time_name)[2,]) != (get_time_points(tvdbn1)-1)]

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
  # Check if the time points are padded
  new_length = NULL
  if (is_padded(trans_network)) {
    new_length = get_time_points(trans_network)+ini
  }

  # Get the old initial time
  old_ini = get_initial_time(trans_network)

  # Get the number of time instants of the network
  n_time_points = get_time_points(trans_network)

  # Set of nodes
  time_nodes = nodes(trans_network)
  new_nodes = c()

  for (i in time_nodes) {
    # Get the old time instant
    decomposed = remove_time_name(i)
    var_name = decomposed[1]
    time_point = as.numeric(decomposed[2])
    new_name = time_name(var_name,time_point+ini-old_ini, len = new_length)
    new_nodes = c(new_nodes, new_name)
  }
  trans_network = rename.nodes(trans_network, new_nodes)
  return(trans_network)

}

remove_zeros <- function(tvdbn) {
  # Set of nodes
  time_nodes = nodes(tvdbn)
  new_nodes = c()

  for (i in time_nodes) {
    # Get the old time instant
    decomposed = remove_time_name(i)
    var_name = decomposed[1]
    time_point = as.numeric(decomposed[2])
    new_name = time_name(var_name,time_point)
    new_nodes = c(new_nodes, new_name)
  }
  tvdbn = rename.nodes(tvdbn, new_nodes)
  return(tvdbn)
}

add_zeros <- function(tvdbn) {
  network_length = get_time_points(tvdbn)

  # Set of nodes
  time_nodes = nodes(tvdbn)
  new_nodes = c()

  for (i in time_nodes) {
    # Get the old time instant
    decomposed = remove_time_name(i)
    var_name = decomposed[1]
    time_point = as.numeric(decomposed[2])
    new_name = time_name(var_name,time_point, len = network_length)
    new_nodes = c(new_nodes, new_name)
  }
  tvdbn = rename.nodes(tvdbn, new_nodes)
  return(tvdbn)
}

sort_nodes <- function(nodes, order = "name") {
  if(order == "name") {
    nodes = nodes[order(as.integer(sapply(nodes, remove_time_name)[2,]))]
    nodes = nodes[order(sapply(nodes, remove_time_name)[1,])]
  }
  else if (order == "time") {
    nodes = nodes[order(sapply(nodes, remove_time_name)[1,])]
    nodes = nodes[order(as.integer(sapply(nodes, remove_time_name)[2,]))]
  }
  else {
    stop("\"order\" parameter should be \"name\" or \"time\"\n")
  }
  return(nodes)
}


