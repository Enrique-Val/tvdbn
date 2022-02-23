# This file define the class tvdbn and some basic functionality of it

# We selected the S3 class convention instead of S4 since bnlearn using S3. We highly recommend adopting
# the S4 in future projects, but in this project it was unfeasible.


#' @title Cast from a regular Bayesian Network to a Time-varying DBN
#'
#' @description  Method to cast from a regular Bayesian Network to a Time-varying Dynamic Bayesian
#' Network. This method should be used withgreat care or outright ignored in favour of learn_tvdbn
#'  or similars. If used, be sure to follow the next conditions:
#'
#'    (1) All the variables should have their corresponding time point appended (e.g variable_t_0)
#'     and all the variables must be present with any time point
#'
#'    (2) Only allow arcs from a time instant t to a time instant t+1 (intraslices arcs are also
#'     forbidden in the current release).
#'
#' @details We selected the S3 class convention instead of S4 since bnlearn using S3. We highly
#'  recommend adopting the S4 in future projects, but in this project was unfeasible.
#'
#' The only attributes that can be of our interest are the name of the variables and number of
#'time points, which can be calculated
#'
#' @param bn The Bayesian Network that we desire to convert to a TV-DBN. The type can be `bn` or `bn.fit`
#' @return A fitted Time-varying Dynamic Bayesian Network based on the Bayesian Network provided.
#'  The type of the return object is `tvdbn.fit`, that inherits from `bn.fit` and from `bn.fit.dnet`
#'  or `bn.fit.gnet`, depending on the type (discrete or gaussian) of the nodes. Alternatively, if
#'  a graph (untrained network) of type `bn` is provided, a object of type `tvdbn` that inherits
#'  from `bn` is returned
#' @import bnlearn
#' @export
bn_to_tvdbn <- function(bn) {
  if (class(bn)[1] == "bn.fit") {
    class(bn) <- c("tvdbn.fit",class(bn))
  }
  if (class(bn)[1] == "bn") {
    class(bn) <- c("tvdbn",class(bn))
  }
  return(bn)
}


#######################
# GETTERS AND SETTERS #
#######################

#' @title Get the number of time points of a Time-Varying DBN
#'
#' @description  Get the number of time points of time instants of the Time-Varying DBN.
#'
#' @param tvdbn The Time-varying DBN. Can be an object of type `tvdbn` or `tvdbn.fit`
#' @return The number of time instants of the Time-Varying DBN.
#' @import bnlearn
#' @export
get_time_points <- function(tvdbn) {
  assertthat::assert_that(any(class(tvdbn) == "tvdbn") || any(class(tvdbn) == "tvdbn.fit"))

  nodes_tvdbn = sort(nodes(tvdbn))
  first_var = remove_time_name(nodes_tvdbn[1])[1]
  time_points = 0

  for (i in nodes_tvdbn) {
    if (first_var == remove_time_name(i)[1]) {
      time_points = time_points+1
    }
    else {
      return(time_points)
    }
  }
}


#' @title Get the template variables
#'
#' @description  Get the template variables of our TV-DBN, i.e. the variables not
#' instantiated in each time point.
#'
#' @param tvdbn The Time-varying DBN. Can be an object of type `tvdbn` or `tvdbn.fit`
#' @return The number of template variables of the Time-Varying DBN.
#' @import bnlearn
#' @export
get_variables <- function(tvdbn) {
  time_instants = get_time_points(tvdbn)
  node_vector = sort(nodes(tvdbn))

  variable_vector = c()
  i = 1
  while(i <= length(node_vector)) {
    variable_vector = c(variable_vector, tvdbn::remove_time_name(node_vector[i])[1])
    i = i + time_instants
  }

  return(variable_vector)
}


get_initial_time <- function(tvdbn) {
  return(as.integer(remove_time_name(sort_nodes(nodes(tvdbn))[1])[2]))
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

is_padded <- function(tvdbn) {
  assertthat::assert_that("tvdbn" %in% class(tvdbn) || "tvdbn.fit" %in% class(tvdbn))
  first_node = sort_nodes(nodes(tvdbn))[1]
  if (nchar(remove_time_name(first_node)[2]) > 1) {
    return(TRUE)
  }
  return(FALSE)
}





###################
## FUNCTIONALITY ##
###################

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
