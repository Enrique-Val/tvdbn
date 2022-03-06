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
plot_time_varying_dbn <- function(structure, offset=100) {
  net_nodes = sort_nodes(nodes(structure), order = "time")
  initial_time = get_initial_time(structure)
  x = c()
  y = c()
  prior_t = -1
  for (i in 1:length(net_nodes)) {
    node = net_nodes[i]
    t = as.integer(remove_time_name(node)[2]) - initial_time
    x = c(x, 100 + (offset+100) * t)
    if (t != prior_t) {
      y = c(y,100)
      prior_t = t
    }
    else {
      y = c(y, y[length(y)]+100)
    }
  }

  time_points = get_time_points(structure)
  color = NULL
  if(utils::packageVersion("grDevices") < "3.6.0")
    color <- grDevices::heat.colors(time_points)
  else
    color <- grDevices::hcl.colors(time_points, palette = "RdYlBu")

  color <- grDevices::col2rgb(color) / 255
  color <- apply(color, 2, function(x){do.call(grDevices::rgb,as.list(x))})

  all_colors = c()
  for (node in net_nodes) {
    t = as.integer(remove_time_name(node)[2]) - initial_time
    all_colors = c(all_colors, color[t+1])
  }


  nodes <- data.frame(id = net_nodes,
                      label = net_nodes,
                      x = x,
                      y = y,
                      #level = node_levels(structure, nodes_uniq)[2,],
                      color.background = all_colors,
                      color.border = "black",
                      borderWidth = 2,
                      #shape = "star",
                      shadow = FALSE,
                      physics = FALSE)

  edges <- data.frame(from = bnlearn::arcs(structure)[,1],
                      to = bnlearn::arcs(structure)[,2],
                      arrows = "to",
                      #color = "orange",
                      #smooth = TRUE,
                      shadow = FALSE,
                      color = "black")

  ret <- visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visOptions(highlightNearest = list(enabled = T, hover = T),
                           nodesIdSelection = T)

  return(ret)

}


plot_spatial_tvdbn <- function(structure, offset = 300) {
  net_nodes = sort_nodes(nodes(structure), order = "time_only") # Sorting is needed is the net is manipulated
  initial_time = get_initial_time(structure)
  x = c()
  y = c()
  prior_t = -1
  prior_lat = -6000
  internal_shift = 100
  real_offset = 100+offset+length(tvdbn::get_dimensions(get_variables(structure))[[1]])*100
  for (i in 1:length(net_nodes)) {
    node = net_nodes[i]
    t = as.integer(remove_time_name(node)[2]) - initial_time

    if (t != prior_t) {
      prior_t = t
      prior_lat = get_coordinates(node)[1]
      internal_shift = 100
      y = c(y,100)
    }
    else if (get_coordinates(node)[1] != prior_lat) {
      prior_lat = get_coordinates(node)[1]
      internal_shift = 100
      y = c(y, y[length(y)]+100)
    }
    else {
      y = c(y, y[length(y)])
    }
    internal_shift = internal_shift+100

    x = c(x, internal_shift + (real_offset+100) * t)
  }

  time_points = get_time_points(structure)
  color = NULL
  if(utils::packageVersion("grDevices") < "3.6.0")
    color <- grDevices::heat.colors(time_points)
  else
    color <- grDevices::hcl.colors(time_points, palette = "RdYlBu")

  color <- grDevices::col2rgb(color) / 255
  color <- apply(color, 2, function(x){do.call(grDevices::rgb,as.list(x))})

  all_colors = c()
  for (node in net_nodes) {
    t = as.integer(remove_time_name(node)[2]) - initial_time
    all_colors = c(all_colors, color[t+1])
  }


  nodes <- data.frame(id = net_nodes,
                      label = net_nodes,
                      x = x,
                      y = y,
                      #level = node_levels(structure, nodes_uniq)[2,],
                      color.background = all_colors,
                      color.border = "black",
                      borderWidth = 2,
                      #shape = "star",
                      shadow = FALSE,
                      physics = FALSE)

  edges <- data.frame(from = bnlearn::arcs(structure)[,1],
                      to = bnlearn::arcs(structure)[,2],
                      arrows = "to",
                      #color = "orange",
                      #smooth = TRUE,
                      shadow = FALSE,
                      color = "black")

  ret <- visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visOptions(highlightNearest = list(enabled = T, hover = T),
                           nodesIdSelection = T)

  return(ret)

}
