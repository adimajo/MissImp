#' MI_EM_amelia: Multiple imputation with EM
#'
#' @description \code{MI_EM_amelia} is a multiple imputation function with EM algorithm.
#' This function returns both multiple imputation results of the original incomplete data and
#' the final imputed data (derived from multiple imputation results by Robin's Rule).
#' As for categorical columns, both factor form and onehot probability vector form are returned.
#'
#' More details about the EM implementation could be found in the documentation of function \code{amelia} from
#'  'Amelia' package.
#' @param df_with_mv Data matrix with missing values.
#' @param col_cat Categorical columns index
#' @param col_num Numerical columns index
#' @param num_imp Number of multiple imputations
#' @export
#' @return \code{ximp} Final imputed data matrix.
#' @return \code{ximp.disj} Final imputed data matrix of same type as 'ximp' for the numeric columns.
#'  For the categorical columns, the prediction of probability for each category is shown in form of onehot probability vector.
#' @return \code{ls_ximp} List imputed data matrix from multiple imputation procedure.
#' @return \code{ls_ximp.disj} List imputed data matrix from multiple imputation procedure, with categorical columns
#' in onehot probability vector form.
#' @references
#' Honaker, J., King, G., Blackwell, M. (2011).
#' Amelia II: A Program for Missing Data.
#' \emph{Journal of Statistical Software}, \bold{45(7)}, 1--47.
#' \url{http://www.jstatsoft.org/v45/i07/}.
#'
#' Statistical Analysis with Missing Data, by Little and Rubin, 2002
MI_EM_amelia <- function(df_with_mv, col_num, col_cat = NULL, num_imp = 5) {
  exist_cat <- !all(c(0, col_cat) == c(0))
  if (exist_cat) {
    name_cat <- colnames(df_with_mv)[col_cat]
    # Deal with the problem that nlevels(df[[col]]) > length(unique(df[[col]]))
    for (col in name_cat) {
      df_with_mv[[col]] <- factor(as.character(df_with_mv[[col]]))
    }
    # remember the levels for each categorical column
    dict_lev <- dict_level(df_with_mv, col_cat)
    # preserve colnames for ximp.disj
    dummy <- dummyVars(" ~ .", data = df_with_mv, sep = "_")
    col_names.disj <- colnames(data.frame(predict(dummy, newdata = df_with_mv)))
    # represent the factor columns with their ordinal levels
    df_with_mv <- factor_ordinal_encode(df_with_mv, col_cat)
    # Create dict_cat with categroical columns
    dict_name_cat <- dict_onehot(df_with_mv, col_cat)


    # imputation
    imp_amelia <- Amelia::amelia(df_with_mv, m = num_imp, p2s = 0, noms = col_cat, boot.type = "none")
    # Extract the imputations with categorical columns in onehot form
    imp_amelia_disj <- list()
    i <- 1
    for (imp in imp_amelia$imputations) {
      imp$index <- as.numeric(row.names(imp))
      imp_num <- imp[, col_num]
      dummy <- dummyVars(" ~ .", data = imp[, col_cat, drop = FALSE], sep = "_")
      imp_cat <- data.frame(predict(dummy, newdata = imp[, col_cat, drop = FALSE]))
      imp_num$index <- as.numeric(row.names(imp_num))
      imp_cat$index <- as.numeric(row.names(imp_cat))
      imp_amelia_disj[[i]] <- merge(x = imp_num, y = imp_cat, by = "index", all = TRUE)
      i <- i + 1
    }
    # Combine multiple imputations
    imp_merge_disj <- Reduce(function(dtf1, dtf2) {
      rbind(dtf1, dtf2)
    }, imp_amelia_disj)

    ximp.disj <- stats::aggregate(. ~ index, data = imp_merge_disj, mean)
    ximp.all <- ximp.disj
    # Derive the final categorical results from the probability vectors
    names_cat <- names(dict_name_cat)
    for (name in names_cat) {
      ximp.all[[name]] <- apply(ximp.all[dict_name_cat[[name]]], 1, which_max_cat, name, dict_name_cat)
      ximp.all[[name]] <- unlist(ximp.all[[name]])
      ximp.all[[name]] <- factor(ximp.all[[name]])
    }
    ximp <- ximp.all[colnames(df_with_mv)]
    ximp.disj <- ximp.disj[-c(1)] # remove "index" column
  }
  else {
    imp_amelia <- Amelia::amelia(df_with_mv, m = num_imp, p2s = 0, boot.type = "none")
    imp_amelia_disj <- list()
    i <- 1
    for (imp in imp_amelia$imputations) {
      imp$index <- as.numeric(row.names(imp))
      imp_amelia_disj[[i]] <- imp
      i <- i + 1
    }

    imp_merge <- Reduce(function(dtf1, dtf2) {
      rbind(dtf1, dtf2)
    }, imp_amelia_disj)
    ximp <- stats::aggregate(. ~ index, data = imp_merge, mean)
    ximp.disj <- ximp[-c(1)]
    ximp <- ximp.disj
  }
  if (exist_cat) {
    for (col in name_cat) {
      levels(ximp[[col]]) <- dict_lev[[col]]
    }
    colnames(ximp.disj) <- col_names.disj
    for (i in seq(length(imp_amelia$imputations))) {
      for (col in name_cat) {
        levels(imp_amelia$imputations[[i]][[col]]) <- dict_lev[[col]]
      }
      colnames(imp_amelia_disj[[i]]) <- col_names.disj
    }
  }

  return(list(ximp = ximp, ximp.disj = ximp.disj, ls_ximp = imp_amelia$imputations, ls_ximp.disj = imp_amelia_disj))
}
