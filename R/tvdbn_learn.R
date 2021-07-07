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

learn_tvdbn_coefficients <- function(x, type = "relaxed", blacklist = list(), whitelist = list()) {
  print(length(x[,1]))
  print(nrow(x))
  start.time <- Sys.time()
  A = list()
  intercept = list()
  variance = list()
  previous_lambda = list()

  # Find the marginal distributions of the elements of the first time point
  t1_weights = weight_time_series(1,nrow(x))
  mean = c()
  sd = c()
  for (i in 1:ncol(x)) {
    mean_i = sum(t1_weights*x[,i])/sum(t1_weights)
    sd_i = sum((mean_i-x[,i])^2*t1_weights)/sum(t1_weights)
    mean = c(mean, mean_i)
    sd = c(sd, sd_i)
  }

  for (t_star in 2:nrow(x)) {
    A[[t_star-1]] = matrix(nrow = length(x[1,]), ncol = length(x[1,]), dimnames = list(map(dimnames(x)[[2]],time_name,t_star-1), map(dimnames(x)[[2]],time_name,t_star-2)))
    intercept[[t_star-1]] = matrix(nrow = length(x[1,]), dimnames = list(map(dimnames(x)[[2]],time_name,t_star-1)))
    variance[[t_star-1]] = matrix(nrow = length(x[1,]), dimnames = list(map(dimnames(x)[[2]],time_name,t_star-1)))
  }

  # Process the blacklist
  blacklist_processed = vector(mode = "list", length = ncol(x))


  print(blacklist_processed)

  for (forbidden_arc in blacklist) {
    i = match(forbidden_arc[1], colnames(x))
    j = match(forbidden_arc[2], colnames(x))
    print(i)
    print(j)
    blacklist_processed[[j]] = c(blacklist_processed[[j]],i)
  }

  #Process the whitelist
  # Process the blacklist
  whitelist_processed = matrix(data=1, nrow = ncol(x), ncol = ncol(x), dimnames=list(dimnames(x)[[2]],dimnames(x)[[2]]))


  print(whitelist_processed)

  for (mandatory_arc in whitelist) {
    i = match(mandatory_arc[1], colnames(x))
    j = match(mandatory_arc[2], colnames(x))
    print(i)
    print(j)
    whitelist_processed[j,i] = 0
  }

  #Iterate for all times and variables
  lambda=list()
  for (i in 1:ncol(x)) {
    for (t_star in 2:nrow(x)) {
      cvfit = NULL
      index = NULL
      weights = NULL
      x_i_t = NULL
      y_i_t = NULL
      ts_length = nrow(x)
      if (type != "causal1" && type != "causal2" ) {
        #Reweight time series
        weights = weight_time_series(t_star, ts_length)
        weights = weights[2:ts_length]
        y_i_t = x[2:length(x[,1]),i]
        x_i_t = x[1:(length(x[,1])-1),]
      }
      else if (type == "causal1") {
        #Reweight time series using the causal function
        weights = weight_causal_time_series(t_star, ts_length)
        weights = weights[2:ts_length]
        y_i_t = x[2:length(x[,1]),i]
        x_i_t = x[1:(length(x[,1])-1),]
        #print(weights)
      }
      else {
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
      # L1-regression
      cvfit = NULL
      exclude_list = blacklist_processed[[i]]
      # In the original article, there cannot be arcs between the same variable in different time points.
      if (type == "original" && !(i %in% exclude_list)) {
        exclude_list = append(exclude_list,i)
      }
      if (t_star>=nrow(x)-4 && (type == "causal1" || type == "causal2")) {
        A[[t_star-1]][i,] = A[[t_star-2]][i,]
        intercept[[t_star-1]][i] = intercept[[t_star-2]][i]
        variance[[t_star-1]][i] = variance[[t_star-2]][i]
      }
      else {
        if (TRUE) {
          cvfit = glmnet::cv.glmnet(x_i_t,y_i_t,weights= weights,
                            exclude = exclude_list, penalty.factor = whitelist_processed[i,],
                            nfolds = 3)
        }
        else {
          cvfit = glmnetPlus::cv.glmnet(x_i_t,y_i_t,weights= weights,
                                    exclude = exclude_list, penalty.factor = whitelist_processed[i,],
                                     lambda = previous_lambda[[i]], beta0 =A[[t_star-2]][i,],
                                    type.gaussian = "naive")
        }
        if (i == 1) {
          lambda[[t_star-1]]=cvfit$lambda
        }
        # Index of the of the fitting that yielded the best results
        index = cvfit$index[2]
        print(cvfit$glmnet.fit$beta[,index])
        #print(cvfit$glmnet.fit$beta
        # Store in a triple: (coefficients, intercept, variance).
        A[[t_star-1]][i,] =  as.matrix(cvfit$glmnet.fit$beta)[,index]
        intercept[[t_star-1]][i] = cvfit$glmnet.fit$a0[index]
        variance[[t_star-1]][i] = sqrt(cvfit$cvm[index])
        previous_lambda[[i]] = cvfit$lambda
      }
    }
  }
  end.time <- Sys.time()
  print(end.time-start.time)
  return(list("A" = A, "intercept" = intercept, "variance" = variance, "lambda" = lambda))

}


learn_tvdbn_structure <- function(A) {
  string_model = ""
  # Variables in the first time instant
  for (i in dimnames(x)[[2]]) {
    string_model = paste(string_model,"[",time_name(i,0),"]",sep="")
  }
  for (t in 1:(length(A))) {
    A.t = A[[t]]
    string_model.t = ""
    for (i in dimnames(A.t)[[1]]) {
      string_model = paste(string_model,"[",i,"|",sep="")

      for (j in dimnames(A.t)[[2]]) {
        if (A.t[i,j] != 0) {
          string_model = paste(string_model,j,":",sep="")
        }
      }
      string_model = str_sub(string_model,1,nchar(string_model)-1)
      string_model = paste(string_model,"]",sep="")
    }
    string_model = paste(string_model,string_model.t,sep="")
  }
  print(string_model)
  dag <-model2network(string_model)
  return(bn_to_tvdbn(dag))
}

learn_tvdbn_parameters <- function(dag, A, intercept, variance) {
  distribution_list = list()
  for (i in dimnames(x)[[2]]) {
    distribution_list[[time_name(i,0)]] = list(coef = c("(Intercept)"=0), sd = 1)
  }

  for (t in 1:(length(A))) {
    A.t = A[[t]]
    for (i in 1:length(dimnames(A.t)[[1]])) {
      i.name = dimnames(A.t)[[1]][i]
      non_zero_entries = A.t[i,]!=0
      coefficient_list = A.t[i,non_zero_entries]
      names(coefficient_list) = names(A.t[i,])[non_zero_entries]
      # Append the intercept to the list.
      coefficient_list = c(coefficient_list, "(Intercept)" = intercept[[t]][i])
      # In t.i.dist we store the normal conditioned distribution of the variable i in the time t
      t.i.dist = list(coef = coefficient_list, sd = variance[[t]][i])
      distribution_list[[i.name]] = t.i.dist
    }
  }
  bayesian_network = custom.fit(dag, distribution_list)
  variables = c()
  for (i in dimnames(A[[1]])[[1]]) {
    variables = c(variables, substring(i,1,str_length(i)-4))

  }

  tvd_bayesian_network = bn_to_tvdbn(bayesian_network)
  #tvd_bayesian_network = list("bn" = bayesian_network, "n" = length(A)+1, "variables" = variables)
  #attr(tvd_bayesian_network, "class") = c("tvdbn.fit",class(bayesian_network))
  return (tvd_bayesian_network)
}

learn_tvdbn <- function(x, type = "relaxed", blacklist = list(), whitelist = list()) {
  ret = learn_tvdbn_coefficients(x, type = type, blacklist = blacklist, whitelist = whitelist)
  A = ret[["A"]]
  intercept = ret[["intercept"]]
  variance = ret[["variance"]]
  dag = learn_tvdbn_structure(A)
  tvdbn = learn_tvdbn_parameters(dag, A, intercept, variance)
  return(tvdbn)
}


