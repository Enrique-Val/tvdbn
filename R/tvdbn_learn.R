#' Imports

library(visNetwork)
library(bnlearn)
library(stringr)
library(dbnR)
library(glmnet)
library(purrr)
library(data.table)

#######################
# Auxiliary functions #
#######################

# Kernel of the function
kernel_function <- function(value, kernel_bandwidth) {
  return(exp(-(value^2)/kernel_bandwidth))
}

causal_kernel_function <- function(value) {
  theta = 2
  K = 2
  eta = -2
  to_ret = ((value-eta)^(K-1)*exp(-(value-eta)/theta)) / (theta^K*factorial(K-1))
  return( to_ret )
}


weight_time_point <- function(t_star,t, dataset_length,kernel_bandwidth = NULL) {
  if (is.null(kernel_bandwidth)) {
    kernel_bandwidth = (dataset_length^2)/49
  }
  tmp <- 0
  for (t_alt in 1:dataset_length) {
    tmp = tmp+kernel_function(t_alt-t_star,kernel_bandwidth)
  }
  return(kernel_function(t-t_star,kernel_bandwidth)/tmp)
}

weight_time_series <- function(t_star, dataset_length,kernel_bandwidth = NULL) {
  if (is.null(kernel_bandwidth)) {
    kernel_bandwidth = (dataset_length^2)/49
  }
  tmp = 1:dataset_length
  tmp = kernel_function(tmp-t_star,kernel_bandwidth)
  tmp = tmp/sum(tmp)
  return(tmp)
}

weight_causal_time_series <- function(t_star, dataset_length) {
  if (t_star == 1) {
    tmp = 1:dataset_length
    tmp = causal_kernel_function(tmp-1)

  }
  else {
    tmp = 1:(dataset_length-t_star+2)
    tmp = c(rep(0,t_star-2), causal_kernel_function(tmp-2))
  }
  tmp = tmp/sum(tmp)
  return(tmp)
}

learn_tvdbn_coefficients <- function(x, type = "relaxed", kernel_bandwidth = NULL, blacklist = list(), whitelist = list(), max_parents = n_variables, padding = FALSE, spatial_penalty = NULL) {
  start.time <- Sys.time()
  A = list()
  intercept = list()
  sd = list()
  if (length(dim(x)) == 2) {
    # A single time series was input
    x = array(x,dim = c(nrow(x),ncol(x),1), dimnames = list(NULL, colnames(x), NULL))
  }
  n_time_series = length(x[1,1,])
  ts_length = nrow(x[,,1]) # Number of time points in the series
  variables = colnames(x[,,1]) # List of the variables of the time series
  n_variables = length(variables) # Number of variables of the series
  previous_lambda = list()
  sd_cv = c()
  sd_b = c()

  # Compute the spatial penalty (in case that we are dealing with spatial data)
  spatial_sigma = NULL
  latitudes = NULL
  longitudes = NULL
  n_lats = NULL
  n_lons = NULL
  if (!is.null(spatial_penalty)) {
    assertthat::assert_that(is.numeric(spatial_penalty))
    # Create a covariance matrix
    spatial_sigma = matrix(c(spatial_penalty,0,0,spatial_penalty),2,2)
    tmp = tvdbn::get_dimensions(colnames(x))
    latitudes = tmp[[1]]
    longitudes = tmp[[2]]
    n_lats = length(latitudes)
    n_lons = length(longitudes)
  }

  # Find the marginal distributions of the elements of the first time point
  weights_t0 = c()
  for (i in 1:n_time_series){
    weights_t0 = c(weights_t0,weight_time_series(1,ts_length, kernel_bandwidth))
  }
  mean_t0 = c()
  sd_t0 = c()
  x_2d = aperm(x,c(1,3,2))
  dim(x_2d)<- c(ts_length*n_time_series, n_variables)
  for (i in 1:n_variables) {
    mean_t0_i = sum(weights_t0*x_2d[,i])/sum(weights_t0)
    sd_t0_i = sum(abs(mean_t0_i-x_2d[,i])*weights_t0)/sum(weights_t0)
    mean_t0 = c(mean_t0, mean_t0_i)
    sd_t0 = c(sd_t0, sd_t0_i)
  }

  intercept[[1]] = mean_t0
  sd[[1]] = sd_t0
  for (t_star in 2:ts_length) {
    A[[t_star-1]] = matrix(nrow = n_variables, ncol = n_variables,
                           dimnames = list(map(variables,time_name,t_star-1, if (padding) ts_length else NULL), map(variables,time_name,t_star-2, if (padding) ts_length else NULL)))
    intercept[[t_star]] = matrix(nrow = n_variables, dimnames = list(map(variables,time_name,t_star-1, if (padding) ts_length else NULL)))
    sd[[t_star]] = matrix(nrow = n_variables, dimnames = list(map(variables,time_name,t_star-1, if (padding) ts_length else NULL)))
  }

  # Process the blacklist
  blacklist_processed = vector(mode = "list", length = n_variables)

  for (forbidden_arc in blacklist) {
    i = match(forbidden_arc[1], variables)
    j = match(forbidden_arc[2], variables)
    blacklist_processed[[j]] = c(blacklist_processed[[j]],i)
  }

  #Process the whitelist
  # Process the blacklist
  whitelist_processed = matrix(data=1, nrow = n_variables, ncol = n_variables, dimnames=list(variables,variables))



  for (mandatory_arc in whitelist) {
    i = match(mandatory_arc[1], variables)
    j = match(mandatory_arc[2], variables)
    whitelist_processed[j,i] = 0
  }

  # Check for other arcs restrictions
  # In the original article, there cannot be arcs between the same variable in different time points.
  if (type == "original") {
    for (name_i in variables) {
      i = match(name_i, variables)
      if (!(i %in% blacklist_processed[[i]])) {
        blacklist_processed[[i]] = c(blacklist_processed[[i]],i)
      }
    }
  }

  # If autoregressive is selected, the model forces x_t-1 as the parent of x_t (opposite of "original")
  else if (type == "autoregressive") {
    for (i in variables) {
      whitelist_processed[i,i] = 0
    }
  }


  #Iterate for all times and variables
  lambda=list()
  for (i in 1:n_variables) {
    # Compute the spatial penalty. The same variable regardless of the time point will have the same penalty.
    # If the net is not spatial (or no spatial penalty was specified), we give every feature the same penalty (if whitelisted, no penalty)
    penalty_i = whitelist_processed[i,]
    if (!is.null(spatial_sigma)) {
      # Compute the mu (that is, the coordinates of the variable)
      coords = tvdbn::get_coordinates(variables[i])
      mu = c(match(coords[1], latitudes), match(coords[2],longitudes))
      penalty_i =penalty_i* (1/mvtnorm::dmvnorm(matrix(c(rep(1:n_lats,each = n_lons),rep(1:n_lons,n_lats)),ncol=2), mean = mu, sigma = spatial_sigma))
    }

    for (t_star in 2:ts_length) {
      cvfit = NULL
      index = NULL
      weights = NULL
      x_i_t = NULL
      y_i_t = NULL
      if (type == "causal") {
        # Reduce time series to the causal boundary (previous + current + following instants)
        # If the current instant is the first or the second, there is nothing to do
        if (t_star == 2) {
          weights = weight_causal_time_series(2, ts_length)
          weights = weights[2:ts_length]
          y_i_t = x[2:ts_length,i]
          x_i_t = x[1:(ts_length-1),]
        }
        else {
          weights = weight_causal_time_series(2, ts_length-t_star+2)
          y_i_t = x[(t_star-1):ts_length,i]
          x_i_t = x[(t_star-2):(ts_length-1),]
        }
      }
      else {
        #Reweight time series
        for (ts in 1:n_time_series) {
          weights_ts = weight_time_series(t_star, ts_length, kernel_bandwidth)
          weights = c(weights, weights_ts[2:ts_length])
          y_i_t = c(y_i_t, x[2:ts_length,i,ts])
          x_i_t = rbind(x_i_t, x[1:(ts_length-1),,ts])
        }
        #y_i_t = as.matrix(y_i_t)
        #x_i_t = as.matrix(x_i_t)

      }
      # L1-regression
      cvfit = NULL
      if (t_star>=ts_length-4 && type == "causal") {
        A[[t_star-1]][i,] = A[[t_star-2]][i,]
        intercept[[t_star]][i] = intercept[[t_star-1]][i]
        sd[[t_star]][i] = sd[[t_star-1]][i]
      }
      else {
        set.seed(0)
        cvfit = glmnet::cv.glmnet(x_i_t,y_i_t,weights= weights,
                          exclude = blacklist_processed[[i]], penalty.factor = penalty_i,
                          nfolds = 3, dfmax = max_parents)

        if (i == 1) {
          lambda[[t_star-1]]=cvfit$lambda
        }
        # Index of the of the fitting that yielded the best results
        index = cvfit$index[2]
        #print(cvfit$glmnet.fit$beta[,index])
        #print(cvfit$glmnet.fit$beta
        # Store in a triple: (coefficients, intercept, sd).
        A[[t_star-1]][i,] =  as.matrix(cvfit$glmnet.fit$beta)[,index]
        intercept[[t_star]][i] = cvfit$glmnet.fit$a0[index]
        sd[[t_star]][i] = sqrt(cvfit$cvm[index])
        previous_lambda[[i]] = cvfit$lambda

        # Compare crossvalidation with biased evaluation
        asse = glmnet::assess.glmnet(cvfit,x_i_t, y_i_t,weights)
        sd_cv = c(sd_cv, cvfit$cvm[index])
        sd_b = c(sd_b, asse$mse[1])
      }
    }
  }
  end.time <- Sys.time()
  print(end.time-start.time)

  #Compare evaluation with biased evaluation
  if (FALSE) {
  print("cv mean and sd")
  print(mean(sd_cv))
  print(sd(sd_cv))
  print("biased mean and sd")
  print(mean(sd_b))
  print(sd(sd_b))
  print("diff mean and sd")
  print(mean(abs(sd_cv-sd_b)))
  print(sd(abs(sd_cv-sd_b)))
  print("diff std mean and sd")
  mean_comp = mean(c(sd_cv,sd_b))
  sd_comp = sd(c(sd_cv,sd_b))
  sd_cv = (sd_cv-mean_comp)/sd_comp
  sd_b = (sd_b-mean_comp)/sd_comp
  diff_list = abs(sd_cv-sd_b)
  print(mean(diff_list))
  print(sd(diff_list))
  }
  return(list("A" = A, "intercept" = intercept, "sd" = sd, "lambda" = lambda))

}


learn_tvdbn_structure <- function(A) {
  dag = NULL
  node_set = dimnames(A[[1]])[[2]]
  for (i in 1:length(A)) {
    node_set = c(node_set, dimnames(A[[i]])[[1]])
  }
  dag = bnlearn::empty.graph(nodes = node_set)

  arc_set = matrix(ncol = 2, dimnames = list(NULL, c("from", "to")))[-1,]
  for (t in 1:(length(A))) {
    A.t = A[[t]]
    for (i in dimnames(A.t)[[1]]) {
      for (j in dimnames(A.t)[[2]]) {
        if (A.t[i,j] != 0) {
          arc_set = rbind(arc_set, c(j, i))
        }
      }
    }
  }

  arcs(dag) = arc_set
  return(tvdbn::bn_to_tvdbn(dag))
}



learn_tvdbn_parameters <- function(dag, A, intercept, sd) {
  distribution_list = list()
  for (i in 1:length(dimnames(A[[1]])[[2]])) {
    distribution_list[[dimnames(A[[1]])[[2]][i]]] = list(coef = c("(Intercept)"=intercept[[1]][i]), sd = sd[[1]][i])
  }

  for (t in 1:(length(A))) {
    A.t = A[[t]]
    for (i in 1:length(dimnames(A.t)[[1]])) {
      i.name = dimnames(A.t)[[1]][i]
      non_zero_entries = A.t[i,]!=0
      coefficient_list = A.t[i,non_zero_entries]
      names(coefficient_list) = names(A.t[i,])[non_zero_entries]
      # Append the intercept to the list.
      coefficient_list = c(coefficient_list, "(Intercept)" = intercept[[t+1]][i])
      # In t.i.dist we store the normal conditioned distribution of the variable i in the time t
      t.i.dist = list(coef = coefficient_list, sd = sd[[t+1]][i])
      distribution_list[[i.name]] = t.i.dist
    }
  }
  bayesian_network = custom.fit(dag, distribution_list)

  tvd_bayesian_network = tvdbn::bn_to_tvdbn(bayesian_network)
  #tvd_bayesian_network = list("bn" = bayesian_network, "n" = length(A)+1, "variables" = variables)
  #attr(tvd_bayesian_network, "class") = c("tvdbn.fit",class(bayesian_network))
  return (tvd_bayesian_network)
}

learn_tvdbn <- function(x, type = "relaxed", kernel_bandwidth = NULL, blacklist = list(), whitelist = list(), max_parents = ncol(x), padding = FALSE, spatial_penalty = NULL) {
  ret = learn_tvdbn_coefficients(x, type = type, kernel_bandwidth = kernel_bandwidth, blacklist = blacklist,
                                 whitelist = whitelist, max_parents = max_parents, padding = padding, spatial_penalty = spatial_penalty)
  A = ret[["A"]]
  intercept = ret[["intercept"]]
  sd = ret[["sd"]]
  dag = learn_tvdbn_structure(A)
  tvdbn = learn_tvdbn_parameters(dag, A, intercept, sd)
  return(tvdbn)
}


