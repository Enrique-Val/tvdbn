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

learn_tvdbn_coefficients <- function(x, type = "relaxed", blacklist = list(), whitelist = list(), max_parents = ncol(x)) {
  print(length(x[,1]))
  print(nrow(x))
  start.time <- Sys.time()
  A = list()
  intercept = list()
  sd = list()
  previous_lambda = list()
  sd_cv = c()
  sd_b = c()

  # Find the marginal distributions of the elements of the first time point
  weights_t0 = weight_time_series(1,nrow(x))
  mean_t0 = c()
  sd_t0 = c()
  for (i in 1:ncol(x)) {
    mean_t0_i = sum(weights_t0*x[,i])/sum(weights_t0)
    sd_t0_i = sum(abs(mean_t0_i-x[,i])*weights_t0)/sum(weights_t0)
    mean_t0 = c(mean_t0, mean_t0_i)
    sd_t0 = c(sd_t0, sd_t0_i)
  }

  intercept[[1]] = mean_t0
  sd[[1]] = sd_t0
  for (t_star in 2:nrow(x)) {
    A[[t_star-1]] = matrix(nrow = length(x[1,]), ncol = length(x[1,]), dimnames = list(map(dimnames(x)[[2]],time_name,t_star-1), map(dimnames(x)[[2]],time_name,t_star-2)))
    intercept[[t_star]] = matrix(nrow = length(x[1,]), dimnames = list(map(dimnames(x)[[2]],time_name,t_star-1)))
    sd[[t_star]] = matrix(nrow = length(x[1,]), dimnames = list(map(dimnames(x)[[2]],time_name,t_star-1)))
  }

  # Process the blacklist
  blacklist_processed = vector(mode = "list", length = ncol(x))


  print(blacklist_processed)

  for (forbidden_arc in blacklist) {
    i = match(forbidden_arc[1], colnames(x))
    j = match(forbidden_arc[2], colnames(x))
    blacklist_processed[[j]] = c(blacklist_processed[[j]],i)
  }

  #Process the whitelist
  # Process the blacklist
  whitelist_processed = matrix(data=1, nrow = ncol(x), ncol = ncol(x), dimnames=list(dimnames(x)[[2]],dimnames(x)[[2]]))


  print(whitelist_processed)

  for (mandatory_arc in whitelist) {
    i = match(mandatory_arc[1], colnames(x))
    j = match(mandatory_arc[2], colnames(x))
    whitelist_processed[j,i] = 0
  }

  # Check for other arcs restrictions
  # In the original article, there cannot be arcs between the same variable in different time points.
  if (type == "original") {
    for (name_i in colnames(x)) {
      i = match(name_i, colnames(x))
      if (!(i %in% blacklist_processed[[i]])) {
        blacklist_processed[[i]] = c(blacklist_processed[[i]],i)
      }
    }
  }

  # If autoregressive is selected, the model forces x_t-1 as the parent of x_t (opposite of "original")
  else if (type == "autoregressive") {
    for (i in colnames(x)) {
      whitelist_processed[i,i] = 0
    }
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
        weights = weight_time_series(t_star, ts_length)
        weights = weights[2:ts_length]
        y_i_t = x[2:length(x[,1]),i]
        x_i_t = x[1:(length(x[,1])-1),]
      }
      # L1-regression
      cvfit = NULL
      if (t_star>=nrow(x)-4 && type == "causal") {
        A[[t_star-1]][i,] = A[[t_star-2]][i,]
        intercept[[t_star]][i] = intercept[[t_star-1]][i]
        sd[[t_star]][i] = sd[[t_star-1]][i]
      }
      else {
        cvfit = glmnet::cv.glmnet(x_i_t,y_i_t,weights= weights,
                          exclude = blacklist_processed[[i]], penalty.factor = whitelist_processed[i,],
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
  return(list("A" = A, "intercept" = intercept, "sd" = sd, "lambda" = lambda))

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



learn_tvdbn_parameters <- function(dag, A, intercept, sd) {
  distribution_list = list()
  for (i in 1:length(dimnames(x)[[2]])) {
    distribution_list[[time_name(dimnames(x)[[2]][i],0)]] = list(coef = c("(Intercept)"=intercept[[1]][i]), sd = sd[[1]][i])
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
  variables = c()
  for (i in dimnames(A[[1]])[[1]]) {
    variables = c(variables, substring(i,1,str_length(i)-4))

  }

  tvd_bayesian_network = tvdbn::bn_to_tvdbn(bayesian_network)
  #tvd_bayesian_network = list("bn" = bayesian_network, "n" = length(A)+1, "variables" = variables)
  #attr(tvd_bayesian_network, "class") = c("tvdbn.fit",class(bayesian_network))
  return (tvd_bayesian_network)
}

learn_tvdbn <- function(x, type = "relaxed", blacklist = list(), whitelist = list(), max_parents = ncol(x)) {
  ret = learn_tvdbn_coefficients(x, type = type, blacklist = blacklist,
                                 whitelist = whitelist, max_parents = max_parents)
  A = ret[["A"]]
  intercept = ret[["intercept"]]
  sd = ret[["sd"]]
  dag = learn_tvdbn_structure(A)
  tvdbn = learn_tvdbn_parameters(dag, A, intercept, sd)
  return(tvdbn)
}


