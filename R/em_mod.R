#' prepare_df_for_em
#' @description Change the position of categorical columns to adapt to the em.mix input form.
#' @param df Original dataframe
#' @param col_cat Categorical columns index in \code{df}
#' @return Dataframe where the categorical columns are the first columns
prepare_df_for_em <- function(df, col_cat) {
  allcol <- c(1:length(colnames(df)))
  col_rest <- allcol[!allcol %in% col_cat]
  return(df[, c(col_cat, col_rest)])
}


#' prob_vector_cat
#' @description Extract the probability vector from a tensor of probability. For example,
#' if there are only 2 categorical variables \code{(Y1, Y2)} with \code{nlevels(Y1)=2, nlevels(Y2)=3},
#'  then the tensor is a matrix of probability. \code{tensor[i,j]} is the probability of \code{Y1=i, Y2=j}.
#'  If \code{obs=(NA,2)}, then \code{prob_vector_cat} returns the probability vector
#'  of \code{Y1} when \code{Y2=2} and the onehot probability vector of \code{Y2=2}.
#'  If \code{obs=(NA,NA)}, then \code{prob_vector_cat} returns the marginal probability vector
#'  of \code{Y1} and that of \code{Y2=2}.
#' @param obs One observation of all categorical variables.
#' @param tensor A tensor (multidimensional array) of probability, with dimension \code{nlevels(Y1) x nlevels(Y2) x...x nlevels(Yk)} if there are k categorical columns.
#' @param dict_name_cat A dictionary of categorical variable names and their corresponding onehot name. For example, \code{Y7: Y7_1, Y7_2,...}.
#' @export
#' @return Probability vector for all categorical columns.
prob_vector_cat <- function(obs, tensor, dict_name_cat) {
  cat_name <- unlist(dict_name_cat)
  obs.disj <- data.frame(matrix(rep(0, length(cat_name)), nrow = 1))
  colnames(obs.disj) <- cat_name
  for (name in names(dict_name_cat)) {
    if (!is.na(obs[[name]])) {
      cat <- paste0(name, "_", obs[[name]])
      obs.disj[[cat]] <- 1
      obs[[name]] <- cat
    }
  }
  if (any(is.na(obs))) { # if there are NA(s)
    obs2 <- as.list(obs)
    obs2[is.na(obs)] <- TRUE
    prob_matrix <- data.frame(do.call(`[`, c(list(tensor), obs2)))
    dimt <- dim(prob_matrix)
    dimt <- dimt[dimt != 1]
    vec <- c()
    for (t in c(1:length(dimt))) {
      vec <- c(vec, apply(prob_matrix, t, sum))
    }
    obs.disj[names(vec)] <- vec
  }
  return(unlist(obs.disj))
}

#' em_mod: modified EM Imputation with probability vector
#'
#' @description \code{em_mod} is a em imputation function that returns categorical columns results both in factor and in onehot probability vector form.
#' Please find the detailed documentation of \code{em.mix} and \code{imp.mix} in the 'mix' package. Only the modifications are explained on this page.
#' After the estimation of parameter \code{pi} in \code{em.mix}, we change it into a tensor (multidimensional array) and
#' extract the probability vector from this tensor with the help of function \code{prob_vector_cat}.
#' @param df Data matrix with missing values.
#' @param col_cat Categorical columns index
#' @export
#' @return \code{ximp} imputed data matrix.
#' @return \code{ximp.disj} imputed data matrix of same type as 'ximp' for the numeric columns.
#'  For the categorical columns, the prediction of probability for each category is shown in form of onehot probability vector.
em_mod <- function(df, col_cat) {
  exist_cat <- !all(c(0, col_cat) == c(0))
  if (!exist_cat) { # if there are only numerical columns
    s <- mix::prelim.mix(df, length(col_cat)) # do preliminary manipulations
    thetahat <- mix::em.mix(s) # ML estimate for unrestricted model
    mix::rngseed(43) # set random number generator seed
    ximp <- mix::imp.mix(s, thetahat, df) # impute under newtheta (one draw)
    return(list(ximp = ximp, ximp.disj = ximp))
  }

  dict_name_cat <- dict_onehot(df, col_cat)
  df.cat <- df[, col_cat]
  # prepare for em and imputation
  # The categorical columns must be ordinal encoded and they must be the first columns of the dataframe
  df_for_em <- prepare_df_for_em(df, col_cat)
  s <- mix::prelim.mix(df_for_em, length(col_cat)) # do preliminary manipulations
  thetahat <- mix::em.mix(s) # ML estimate for unrestricted model
  mix::rngseed(43) # set random number generator seed
  ximp <- mix::imp.mix(s, thetahat, df_for_em) # impute under newtheta (one draw)
  ximp <- factor_encode(data.frame(ximp), c(1:length(col_cat)))
  # rearrange columns
  cols <- colnames(df)
  ximp <- ximp[cols]

  # create disjonctive table
  dims <- c()
  dimnames <- list()
  dummy <- dummyVars(" ~ .", data = ximp, sep = "_")
  ximp.disj <- data.frame(predict(dummy, newdata = ximp))

  # construction of tensor pi
  i <- 1
  for (name in names(dict_name_cat)) {
    # ximp.disj[[name]] <- df_with_mv[[name]]
    lev <- dict_name_cat[[name]]
    dims <- c(dims, length(lev))
    dimnames[[i]] <- lev
    i <- i + 1
  }
  pi <- thetahat$pi
  tensor_pi <- array(pi, dims, dimnames)

  # ximp.disj onehot probability
  cat_name <- unlist(dict_name_cat)
  ximp.disj[cat_name] <- t(apply(as.matrix(df.cat), 1,
    FUN = function(x) prob_vector_cat(x, tensor_pi, dict_name_cat)
  ))
  return(list(ximp = data.frame(ximp), ximp.disj = data.frame(ximp.disj)))
}
