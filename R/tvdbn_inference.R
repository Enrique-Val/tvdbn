inference <- function(tvdbn.fit, folded_dt, ini, end) {
  network_length = NULL
  if (is_padded(tvdbn.fit)) {
    network_length = get_time_points(tvdbn.fit)
  }

  mean_results = matrix(nrow = end-ini+1, ncol = ncol(folded_dt), dimnames = list(NULL,colnames(folded_dt)))
  sd_results = matrix(nrow = end-ini+1, ncol = ncol(folded_dt), dimnames = list(NULL,colnames(folded_dt)))
  error_results = matrix(nrow = end-ini+1, ncol = ncol(folded_dt), dimnames = list(NULL,colnames(folded_dt)))

  mean_results[1,] = folded_dt[ini+1,]
  sd_results[1,] = rep(0,ncol(folded_dt))
  error_results[1,] = rep(0,ncol(folded_dt))

  for (row in 2:(end-ini+1)) {
    for (col in colnames(folded_dt)) {
      variable = time_name(col, row-1+ini, network_length)
      coefs = tvdbn.fit[[variable]][["coefficients"]]
      name_coefs = names(tvdbn.fit[[variable]][["coefficients"]])
      sum = coefs[1]
      for (k in name_coefs[-1]) {
        sum = sum + coefs[k]*mean_results[row-1, remove_time_name(k)[1]]
      }
      mean_results[row, col] = sum

      sd_i = tvdbn.fit[[variable]][["sd"]]^2
      for (k in name_coefs[-1]) {
        sd_i = sd_i + (coefs[k]*sd_results[row-1, remove_time_name(k)[1]])^2
      }
      sd_results[row,col] = sqrt(sd_i)

      error_results[row,col] =(sum - folded_dt[row+ini,col])^2
    }
  }
  mse = mean(error_results[-1,])
  if (FALSE) {
    rmse = c()
    for (i in 1:ncol(error_results)) {
      rmse = c(rmse, mean(error_results[-1,i])/sd(folded_dt[,i]))
    }
    rmse = mean(rmse)
  }
  return(list(mean_results, sd_results, error_results, mse))
}

