# This file define the class tvdbn and some basic functionality of it

# We selected the S3 class convention instead of S4 since bnlearn using S3. We highly recommend adopting
# the S4 in future projects, but in this project was unfeasible.


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



#' @title Get the number of time points of a Time-Varying DBN
#'
#' @description  Get the number of time points of time instants of the Time-Varying DBN.
#'
#' @param tvdbn The Time-varying DBN. Can be an object of type `tvdbn` or `tvdbn.fit`
#' @return The number of time points of time instants of the Time-Varying DBN.
#' @import bnlearn
#' @export
get_time_points <- function(tvdbn) {
  if (class(tvdbn)[1] == "tvdbn.fit") {
    first_var = remove_time_name(nodes(tvdbn)[1])[1]
    time_points = 0
    for (i in nodes(tvdbn)) {
      if (first_var == remove_time_name(i)[1]) {
        time_points = time_points+1
      }
      else {
        return(time_points)
      }
    }
  }

  else if (class(tvdbn)[1] == "tvdbn") {
    first_var = remove_time_name(names(tvdbn$nodes)[1])[1]
    time_points = 0
    for (i in names(tvdbn$nodes)) {
      if (first_var == remove_time_name(i)[1]) {
        time_points = time_points+1
      }
      else {
        return(time_points)
      }
    }
  }
}
