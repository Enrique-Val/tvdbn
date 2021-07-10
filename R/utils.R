#' @title  Construct the time name of a variable
#'
#' @description  Construct the time name of a variable (i.e. the name of the variable with the time point correctly appended)
#'
#' @param name Name of the variable
#' @param t Time instant (integer)
#' @return The name of the variable with the time point correctly appended.
#'
#' @export
time_name <- function(name, t) {
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
  return(subgraph(tvdbn.fit = tvdbn.fit, nodes = bnlearn::nodes(tvdbn.fit)))
}


#' @title  Transition network in a time point
#'
#' @description  Get the transition network of a certain time instant.
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`
#' @return The structure of the TV-DBN, type = `tvdbn`
#' @import purrr
#'
#' @export
transition_network <- function(tvdbn.fit, time) {
  variables = get_variables(tvdbn.fit)
  # First, we store the conditional probabilities of nodes of the transition network
  prob_list = list()
  for (i in variables) {
    prob_list[[time_name(i,time)]] = tvdbn.fit[[time_name(i,time)]]
  }
  return(prob_list)
}


#' @title  Structure of the transition network in a time point
#'
#' @description  Get the transition network of a certain time instant.
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`.
#' @param time The transition network to obtain.
#' @return The structure of the transition network of the TV-DBN in
#' the given time instant.
#' @import purrr
#'
#' @export
transition_network_graph <- function(tvdbn.fit, time) {
  variables = get_variables(tvdbn.fit)
  # First, we store the conditional probabilities of nodes of the transition network
  var_list = unlist(map(variables, time_name, time))
  if (time > 0) {
    var_list = c(var_list,unlist(map(variables, time_name, time-1)))
  }
  return(tvdbn::subgraph(tvdbn.fit = tvdbn.fit, nodes = var_list))
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
  order = length(trans_network$nodes)/length(get_variables(trans_network)) -1
  return(order)
}


#' @title  Structure of the transition network in a time point
#'
#' @description  Get the transition network of a certain time instant.
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`.
#' @param time The transition network to obtain.
#' @return The structure of the transition network of the TV-DBN in
#' the given time instant.
#'
#' @import purrr
#' @export
transition_network_normalize_name <- function(trans_network) {
  # Get the Markovian Order of the Transition Network
  order = transition_network_markovian_order(trans_network)
  min = Inf
  sorted_list = sort(names(trans_network$nodes))
  for (i in sorted_list[0:order+1]) {
    if (as.numeric(remove_time_name(i)[2]) < min) {
      min = as.numeric(remove_time_name(i)[2])
    }
  }
  mstring = modelstring(trans_network)

  for (i in 0:order) {
    for (separator in c("]","\\|",":")){
      old = paste("_t_",as.character(min+i),separator,sep="")
      new = paste("_t_",as.character(i),separator,sep="")
      mstring = gsub(old, new, mstring)
    }
  }

  return(model2network(mstring))

}
