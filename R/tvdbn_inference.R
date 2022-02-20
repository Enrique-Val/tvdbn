inference <- function(tvdbn.fit, folded_dt, ini, end) {
  mean_results = matrix(nrow = end-ini+1, ncol = ncol(folded_dt), dimnames = list(NULL,colnames(folded_dt)))
  sd_results = matrix(nrow = end-ini+1, ncol = ncol(folded_dt), dimnames = list(NULL,colnames(folded_dt)))
  mse = c()

  mean_results[1,] = folded_dt[ini+1,]

  for (row in 2:(end-ini+1)) {
    for (col in colnames(folded_dt)) {
      variable = time_name(col, row-1+ini, get_time_points(tvdbn.fit))
      coefs = tvdbn.fit[[variable]][["coefficients"]]
      name_coefs = names(tvdbn.fit[[variable]][["coefficients"]])
      sum = coefs[1]
      for (k in name_coefs[-1]) {
        sum = sum + coefs[k]*mean_results[row-1, remove_time_name(k)[1]]
      }
      mean_results[row, col] = sum
      mse = c(mse, (sum - folded_dt[row+ini,col])^2)
    }
  }
  return(list(mean_results, sd_results, mean(mse)))
}
