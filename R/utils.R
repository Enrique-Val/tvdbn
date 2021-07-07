#' @title  Construct the time name of a variable
#'
#' @description  Construct the time name of a variable (i.e. the name of the variable with the time point correctly appended)
#'
#' @param name Name of the variable
#' @param t Time instant (integer)
#' @return The name of the variable with the time point correctly appended.
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
#'    2) The time point
remove_time_name <- function(time_name) {
  split = strsplit(time_name, "_")[[1]]
  return(c(paste(head(split,length(split)-2), collapse="_"),split[length(split)]))
}

#' @title  Structure of a Time-varying DBN
#'
#' @description  Get the structure of a Time-varying DBN given a fitted Time-varying DBN
#' @param tvdbn.fit A fitted Time-varying DBN of type `tvdbn.fit`
#' @return The structure of the TV-DBN, type = `tvdbn`
graph <- function(tvdbn.fit) {
  dag = bnlearn::subgraph(x = tvdbn.fit, nodes = bnlearn::nodes(tvdbn.fit))
  return(bn_to_tvdbn(dag))
}
