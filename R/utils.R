#' ordinal_encode
#' @description Ordinal encode the categorical columns in a dataframe
#' @param df Dataframe
#' @param idx_col_cat Categorical columns index
#' @export
#' @return Dataframe after ordinal encoding the categorical columns
ordinal_encode <- function(df, idx_col_cat) {
  for (j in idx_col_cat) {
    df[, j] <- as.numeric(factor(df[, j], levels = levels(df[, j])))
  }
  return(df)
}

#' factor_encode
#' @description Factor encode the categorical columns in a dataframe
#' @param df Dataframe
#' @param idx_col_cat Categorical columns index
#' @return Dataframe after factor encoding the categorical columns
#' @export
factor_encode <- function(df, idx_col_cat) {
  for (j in idx_col_cat) {
    df[, j] <- factor(df[, j])
  }
  return(df)
}



#' factor_ordinal_encode
#' @description Change the levels of categorical columns into integer numbers
#' @param df Dataframe
#' @param idx_col_cat Categorical columns index
#' @export
#' @return Dataframe after encoding the categorical columns
factor_ordinal_encode <- function(df, idx_col_cat) {
  for (j in idx_col_cat) {
    levels(df[, j]) <- c(1:nlevels(df[, j]))
  }
  return(df)
}

#' normalize_num
#' @description Scale the numerical columns in a dataframe
#' @param df Dataframe
#' @param idx_col_num Numerical columns index
#' @return Dataframe after scaling the numerical columns
normalize_num <- function(df, idx_col_num) {
  df[, idx_col_num] <- data.frame(apply(df[, idx_col_num], 2, scale))
  return(df)
}


#' VA_fact
#' @description Calculate the Wilcox's VarNC for a vector of categorical values
#' @param x A vector of categorical variable
#' @return Wilcox's VarNC
VA_fact <- function(x) {
  lev_miss <- nlevels(x) - length(unique(x))
  freq <- table(x) / length(x)
  return(qualvar::VA(c(freq, rep(0, lev_miss))))
}


#' dict_onehot
#' @description Create the dictionary for categorical columns. For example, {Y7:Y7_1,Y7_2,...}.
#' @param df Dataframe.
#' @param col_cat categorical column index
#' @export
#' @return A dictionary that connects onehot form and factor form of the categorical columns.
dict_onehot <- function(df, col_cat) {
  dict_name <- list()
  ls_col_name <- colnames(df)
  for (i in col_cat) {
    col_cat_name <- ls_col_name[i]
    dict_name[[col_cat_name]] <- paste(col_cat_name, levels(df[[col_cat_name]]), sep = "_")
  }
  return(dict_name)
}

#' Mode_cat
#' @description Find the most frequent result for one categorical variable.
#' This function is used in combine_boot with method = 'factor'
#' @param x Vector of categorical results. For example c('A','B','A').
#' @return Most frequent category in vector \code{x}.
Mode_cat <- function(x) {
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


#' dict_level
#' @description Create the dictionary for categorical columns levels. For example, {Y7:A,B,C,...}.
#' @param df Dataframe.
#' @param col_cat categorical column index
#' @export
#' @return A dictionary that connects categorical columns with their original levels.
dict_level <- function(df, col_cat) {
  ls_col_cat_name <- colnames(df)[col_cat]
  res <- list()
  for (col_name in ls_col_cat_name) {
    res[[col_name]] <- c(levels(df[[col_name]]))
  }
  return(res)
}


#' which.max.random
#' @description Determines the location, i.e., index of the maximum of a numeric (or logical) vector.
#' In case of tie, the index is selected randomly.
#' @param vec Numeric (logical, integer or double) vector.
#' @export
#' @return Index of the maximum value.
which.max.random <- function(vec) {
  maxima <- which(vec == max(vec))
  if (length(maxima) > 1) {
    maxima <- sample(maxima, 1)
  }
  return(maxima)
}


which_max_cat <- function(x, name, dict_cat) {
  return(dict_cat[[name]][which.max.random(x)])
}
