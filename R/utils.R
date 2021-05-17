# The most frequent result for one categorical variable
# Used in combine_boot with method = 'factor'
Mode_cat <- function(x) {
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# Calculate the Wilcox's VarNC for a vector of categorical values
VA_fact <- function(x) {
  lev_miss <- length(levels(x)) - length(unique(x))
  freq <- table(x) / length(x)
  return(VA(c(freq, rep(0, lev_miss))))
}


# Calculate the missing proportion in MAR3, used in 'generate_miss'
# Solve (1-x)^p + (1-m)*p*x -1 = 0 in (0, 1)
# where m = miss_perc, p = num_co
monot_quantil <- function(miss_perc, num_col) {
  m <- miss_perc
  p <- num_col
  tmp_result <- c()
  tempt <- linspace(0.01, 1, n = 1000)
  i <- 1
  for (x in tempt) {
    tmp_result[i] <- abs((1 - x)^p + (1 - m) * p * x - 1)
    i <- i + 1
  }
  return(tempt[which.min(tmp_result)])
}
