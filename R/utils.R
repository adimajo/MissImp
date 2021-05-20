#' ordinal_encode
#' @description Ordinal encode the categorical columns in a dataframe
#' @param df Dataframe
#' @param idx_col_cat Categorical columns index
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
factor_encode <- function(df, idx_col_cat) {
  for (j in idx_col_cat) {
    df[, j] <- factor(df[, j])
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
  lev_miss <- length(levels(x)) - length(unique(x))
  freq <- table(x) / length(x)
  return(qualvar::VA(c(freq, rep(0, lev_miss))))
}


#Create a dictionary {Y7:Y7_1,Y7_2,...}, {Y_7: Y_7_1, Y_7_2,...}
dict_onehot = function(df, col_cat){
  dict_name = list()
  ls_col_name = colnames(df)
  for(i in col_cat){
    col_cat_name = ls_col_name[i]
    dict_name[[col_cat_name]]=paste(col_cat_name, levels(df[[col_cat_name]]), sep = "_")
  }
  return(dict_name)
}
