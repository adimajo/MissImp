
#' Imputation of Missing Values with variance calculation
#'
#' @description \code{MissImp} is a function that could impute the mixed-type missing values with a chosen imputation method
#' (single imputation or multiple imputation). It can be used to impute continuous, integer and/or categorical data.
#'  Bootstrap or Jackknife resampling method could also be applied to estimate the variance of prediction.
#'  This function returns not only the final imputation result \code{df_result}, but also the estimated variance for each imputed value \code{df_result_var}.
#'  For the categorical variables, the factor formed result is returned in \code{df_result}, while the probability vector (onehot) formed result
#'  is returned in \code{df_result_disj}. If the original complete dataset is also given, a MSE calculated on numerical columns is returned,
#'  as well as a F1-score calculated on the categorical columns.
#'
#' @param df A \code{data.frame} with missing values to impute. For categorical variables, most of the preprocessing is included in the package,
#' but no special character ('+', '$', '-' for example) is allowed in the levels of each categorical variable.
#' @param imp_method Imputation method to be chosen.
#' \itemize{
#' \item \strong{missRanger}: Random-forest based imputation method. This is a fast version of 'missForest'
#' \item \strong{missForest}: Random-forest based imputation method.
#' \item \strong{kNN}: k-Nearest Neighbors based imputation method.
#' \item \strong{PCA}: Principle Component Analysis based imputation method.
#' Factor Analysis of Mixed Data (FAMD) is used for mixed-type data. With large proportion of missing data, this method may generate some convergence error.
#' \item \strong{EM}: Expectation-Maximization imputation method. This method can not deal with large number of categorical variables.
#' \item \strong{MI_EM}: Multiple imputation with EM method.
#' \item \strong{MI_PCA}: Multiple imputation with PCA-based method.
#' \item \strong{MICE}: Multiple Imputation by Chained Equations (MICE). This is a bayesian based method.
#' \item \strong{MI_Ranger}: Multiple imputation with Random Forest based method.
#' }
#' @param resample_method Resampling method to be chosen.
#' \itemize{
#' \item \strong{bootstrap}: Generate \code{n_resample} incomplete datasets by draws with replacement.
#' Each incomplete dataset has the same number of rows as the original one.
#' \item \strong{jackknife}: Generate \code{n_resample} incomplete datasets by removing nrows/\code{n_resample} rows of the original dataset.
#' \item \strong{none}: No resampling.
#' }
#' @param n_resample Number of datasets created by resampling.
#' @param col_cat Index of categorical columns.
#' @param col_dis Index of discrete columns with integer value.
#' @param maxiter_tree Max number of iterations with tree-based methods ("missRanger", "missForest", "MI_Ranger", "MI_Ranger_bis").
#' @param maxiter_pca Max number of iterations with PCA-based methods ("PCA", "MI_PCA").
#' @param maxiter_mice Max number of iterations with MICE imputation method.
#' @param ncp_pca Number of component for PCA-based methods ("PCA", "MI_PCA").
#' @param learn_ncp If \code{ncp_pca} is learned. This could lead to long execution time.
#' @param cat_combine_by Combine method for the categorical part of imputed resampled datasets.
#' \itemize{
#' \item \strong{factor} The mode value of several predictions is taken as the final result.
#' \item \strong{onehot} The average of the probability vectors is taken before choosing the category that maximise the probability as the final result.
#' }
#' @param var_cat Method used to caculate the 'variance' on a prediction of categorical variable. 'wilcox_va' means the Wilcox's VarNC index and 'unalike' means unalikeability.
#' @param df_complete Complete dataset without missing values (if known). With this dataset, the performance of each imputation method could be estimated.
#' @param num_mi Number of multiple imputation.
#' @return \code{df_result} The final imputation result.
#' @return \code{df_result_disj} The imputation result with probability vector for the categorical columns.
#' @return \code{df_result_var} The variance for each imputed value. These variances are calculated based on imputation results of the incomplete datasets after resampling.
#' @return \code{df_result_var_disj} The variance for each imputed value and for the probability vector of the categorical columns.
#' @return \code{MSE} MSE calculated with the numerical imputations.
#' @return \code{F1} F1-score calculated with the categorical imputations.
#' @export
#' @importFrom stats var reformulate terms.formula predict setNames
#'
#' @export
MissImp <- function(df, imp_method = "missRanger", resample_method = "bootstrap",
                    n_resample = 2 * round(log(nrow(df))), col_cat = c(), col_dis = c(),
                    maxiter_tree = 10, maxiter_pca = 100, maxiter_mice = 10, ncp_pca = ncol(df) / 2,
                    learn_ncp = FALSE, cat_combine_by = "factor", var_cat = "wilcox_va",
                    df_complete = NULL, num_mi = NULL) {
  if (all(!is.na(df))) {
    stop("The input dataframe is complete. Imputation is not needed.\n")
  }
  imp_method <- match.arg(imp_method, c("missRanger", "kNN", "missForest", "PCA", "EM", "MI_EM", "MI_PCA", "MICE", "MI_Ranger", "MI_Ranger_bis"))
  resample_method <- match.arg(resample_method, c("bootstrap", "jackknife", "none"))
  cat_combine_by <- match.arg(cat_combine_by, c("factor", "onehot"))
  var_cat <- match.arg(var_cat, c("wilcox_va", "unalike"))
  if (imp_method == "MI_EM" && is.null(num_mi)) {
    stop("For multiple imputation method, 'num_mi' is needed.\n")
  }
  exist_cat <- !all(c(0, col_cat) == c(0))
  num_col <- ncol(df)
  num_row <- nrow(df)
  col_con <- c(1:num_col)
  col_con <- col_con[!col_con %in% c(col_cat, col_dis)]




  ## 0. Preparation
  if (exist_cat) {
    dict_lev <- dict_level(df, col_cat)
    dict_name_cat <- dict_onehot(df, col_cat)
  }


  ## 1. Create several datasets. (Resampling)
  if (resample_method == "bootstrap") {
    ls_df <- bootsample(df, n_resample)
  } else if (resample_method == "jackknife") {
    ls_df <- jacksample(df, n_resample)
  } else { # if resample_method=='none', there will be no resampling
    ls_df <- list(df)
  }


  ## 2. Single Imputation
  ## Input: ls_df list of dataframes after resampling
  ## Output: ls.imp.onehot, ls.imp.fact two list of imputed dataframe, the only difference is the form of the categorical columns. If there is not any categorical column, then these two lists are identical.
  ls.imp.onehot <- list()
  ls.imp.fact <- list()
  i <- 1
  for (dfi in ls_df) {
    dummy <- dummyVars(" ~ .", data = dfi, sep = "_")
    tmp.disj <- data.frame(predict(dummy, newdata = dfi))
    if (imp_method == "missRanger") {
      res <- missRanger_mod(dfi, col_cat = col_cat, maxiter = maxiter_tree)
      imp.onehot_i <- res$ximp.disj
      imp.fact_i <- res$ximp
    } else if (imp_method == "missForest") {
      res <- missForest_mod(xmis = dfi, maxiter = maxiter_tree, col_cat = col_cat)
      imp.onehot_i <- res$ximp.disj
      imp.fact_i <- res$ximp
    } else if (imp_method == "kNN") {
      res <- suppressWarnings((kNN_mod(dfi, col_cat = col_cat, weightDist = TRUE)))
      imp.onehot_i <- res$ximp.disj
      imp.fact_i <- res$ximp
    } else if (imp_method == "EM") {
      res <- em_mod(dfi, col_cat = col_cat) # Error: matrix not sinugular, the categorical variable with too many categories may ends in a sparse matrix
      imp.onehot_i <- res$ximp.disj
      imp.fact_i <- res$ximp
    } else if (imp_method == "PCA") {
      if (learn_ncp) {
        ncp_pca <- estim_ncpFAMD_mod(dfi, method.cv = "Kfold", verbose = F, maxiter = maxiter_pca)$ncp
      }
      res <- imputeFAMD_mod(dfi, ncp = ncp_pca, maxiter = maxiter_pca)
      imp.onehot_i <- res$tab.disj
      imp.fact_i <- res$completeObs
    } else if (imp_method == "MI_EM") {
      res <- MI_EM_amelia(dfi, col_num = c(col_con, col_dis), col_cat = col_cat, num_imp = num_mi)
      imp.onehot_i <- res$ximp.disj
      imp.fact_i <- res$ximp
    } else if (imp_method == "MI_PCA") {
      res <- MIFAMD_mod(dfi, ncp = ncp_pca, maxiter = maxiter_pca, nboot = num_mi)
      imp.onehot_i <- res$ximp.disj
      imp.fact_i <- res$ximp
    } else if (imp_method == "MICE") {
      col_cat <- which(!sapply(dfi, is.numeric))
      exist_cat <- !all(c(0, col_cat) == c(0))
      if (exist_cat) {
        name_cat <- colnames(dfi)[col_cat]
        # Deal with the problem that nlevels(df[[col]]) > length(unique(df[[col]]))
        for (col in name_cat) {
          dfi[[col]] <- factor(as.character(dfi[[col]]))
        }
        # remember the levels for each categorical column
        dict_lev <- dict_level(dfi, col_cat)
        # preserve colnames for ximp.disj
        dummy <- dummyVars(" ~ .", data = dfi, sep = "_")
        col_names.disj <- colnames(data.frame(predict(dummy, newdata = dfi)))
        # represent the factor columns with their ordinal levels
        dfi <- factor_ordinal_encode(dfi, col_cat)
      }
      res0 <- mice::mice(dfi, m = num_mi, maxit = maxiter_mice)
      res <- result_mice(res0, impnum = num_mi, col_cat = col_cat)
      # change back to original levels
      if (exist_cat) {
        for (col in name_cat) {
          levels(res$ximp[[col]]) <- dict_lev[[col]]
        }
        colnames(res$ximp.disj) <- col_names.disj
      }
      imp.onehot_i <- res$ximp.disj
      imp.fact_i <- res$ximp
    } else if (imp_method == "MI_Ranger") {
      res <- MI_missRanger(data.frame(dfi), col_cat = col_cat, num_mi = num_mi, maxiter = maxiter_tree)
      imp.onehot_i <- res$ximp.disj
      imp.fact_i <- res$ximp
    }
    # } else if (imp_method == "MI_Ranger_bis") {
    #   res <- MI_missRanger_bis(data.frame(dfi), col_cat = col_cat, num_mi = num_mi, maxiter = maxiter_tree)
    #   imp.onehot_i <- res$ximp.disj
    #   imp.fact_i <- res$ximp
    # }
    # in case there is one category that doesn't appear in dfi
    tmp.disj[colnames(imp.onehot_i)] <- imp.onehot_i
    tmp.disj[is.na(tmp.disj)] <- 0
    ls.imp.onehot[[i]] <- data.frame(tmp.disj)
    ls.imp.fact[[i]] <- data.frame(imp.fact_i)
    i <- i + 1
  }

  ## 3. Final result
  if (resample_method == "bootstrap") {
    if (cat_combine_by == "factor") {
      ls.imp.tmp <- ls.imp.fact
      col_cat_boot <- col_cat
      col_con_boot <- col_con
      col_dis_boot <- col_dis
    } else {
      ls.imp.tmp <- ls.imp.onehot
      # With onehot form, the column index have changed. Y7 -> Y7_1,...Y7_6
      num_name <- colnames(df[, c(col_con, col_dis), drop = FALSE])
      con_name <- colnames(df[, col_con, drop = FALSE])
      dis_name <- colnames(df[, col_dis, drop = FALSE])
      col_cat_boot <- which(!colnames(ls.imp.onehot[[1]]) %in% num_name)
      col_con_boot <- which(colnames(ls.imp.onehot[[1]]) %in% con_name)
      col_dis_boot <- which(colnames(ls.imp.onehot[[1]]) %in% dis_name)
    }
    res <- combine_boot(
      ls_df = ls.imp.tmp,
      col_con = col_con_boot, col_dis = col_dis_boot, col_cat = col_cat_boot,
      num_row_origin = num_row, method = cat_combine_by, dict_cat = dict_name_cat,
      var_cat = var_cat
    )
  } else if (resample_method == "jackknife") {
    if (cat_combine_by == "factor") {
      stop("Please choose cat_combine_by='onehot' when resample_method=='jackknife'.")
    }
    # With onehot form, the column index have changed. Y7 -> Y7_1,...Y7_6
    num_name <- colnames(df[, c(col_con, col_dis), drop = FALSE])
    con_name <- colnames(df[, col_con, drop = FALSE])
    dis_name <- colnames(df[, col_dis, drop = FALSE])
    col_cat_jack <- which(!colnames(ls.imp.onehot[[1]]) %in% num_name)
    col_con_jack <- which(colnames(ls.imp.onehot[[1]]) %in% con_name)
    col_dis_jack <- which(colnames(ls.imp.onehot[[1]]) %in% dis_name)
    res <- combine_jack(ls.imp.onehot,
      col_con = col_con_jack,
      col_dis = col_dis_jack, col_cat = col_cat_jack, method = cat_combine_by,
      dict_cat = dict_name_cat, var_cat = var_cat
    )
  } else { # if resample_method=='none', there will be no resampling
    res <- list()
    res[["df_result"]] <- ls.imp.fact
    res[["df_result_disj"]] <- ls.imp.onehot
    res[["df_result_var_disj"]] <- NA
    res[["df_result_var"]] <- NA
  }


  # Change back the categorical variable levels
  if (exist_cat) {
    name_cat <- names(dict_lev)
    for (name in name_cat) {
      res$imp[[name]] <- apply(as.array(as.character(res$imp[[name]])), 1, function(x) {
        unlist(strsplit(x, "_"))[-1]
      })
    }
  }




  ## 4. Evaluation matrix
  if (!is.null(df_complete)) { # original complete dataset is provided
    mask <- data.frame(is.na(df))
    colnames(mask) <- colnames(df)
    MSE_imp <- ls_MSE(df_complete, ls.imp.fact,
      mask = mask, col_num_comp = c(col_con, col_dis),
      resample_method = resample_method
    )
    if (exist_cat && cat_combine_by == "factor") {
      F1_imp <- ls_F1(df_complete, ls.imp.fact,
        mask = mask, col_cat_comp = col_cat, col_cat_imp = col_cat,
        resample_method = resample_method, combine_method = cat_combine_by
      )
    } else if (exist_cat && cat_combine_by == "onehot") {
      # With onehot form, the column index have changed. Y7 -> Y7_1,...Y7_6
      num_name <- colnames(df[, c(col_con, col_dis)])
      col_cat_F1 <- which(!colnames(ls.imp.onehot[[1]]) %in% num_name)
      F1_imp <- ls_F1(df_complete, ls.imp.onehot,
        mask = mask, col_cat_comp = col_cat, col_cat_imp = col_cat_F1,
        resample_method = resample_method, combine_method = cat_combine_by,
        dict_cat = dict_name_cat
      )
    } else {
      F1_imp <- NA
    }
    res[["MSE"]] <- MSE_imp
    res[["F1"]] <- F1_imp
  }

  return(res)
}



# Version jackknife estimate using estimation on the full incomplete dataset
# MissImp_bis <- function(df, imp_method = "missRanger", resample_method = "bootstrap",
#                     n_resample = 2 * round(log(nrow(df))), col_cat = c(), col_dis = c(),
#                     maxiter_tree = 10, maxiter_pca = 100, maxiter_mice = 10, ncp_pca = ncol(df) / 2,
#                     learn_ncp = FALSE, cat_combine_by = "factor", var_cat = "wilcox_va",
#                     df_complete = NULL, num_mi = NULL) {
#   if (all(!is.na(df))) {
#     stop("The input dataframe is complete. Imputation is not needed.\n")
#   }
#   imp_method <- match.arg(imp_method, c("missRanger", "kNN", "missForest", "PCA", "EM", "MI_EM", "MI_PCA", "MICE", "MI_Ranger", "MI_Ranger_bis"))
#   resample_method <- match.arg(resample_method, c("bootstrap", "jackknife", "none"))
#   cat_combine_by <- match.arg(cat_combine_by, c("factor", "onehot"))
#   var_cat <- match.arg(var_cat, c("wilcox_va", "unalike"))
#   if (imp_method == "MI_EM" && is.null(num_mi)) {
#     stop("For multiple imputation method, 'num_mi' is needed.\n")
#   }
#   exist_cat <- !all(c(0, col_cat) == c(0))
#   num_col <- ncol(df)
#   num_row <- nrow(df)
#   col_con <- c(1:num_col)
#   col_con <- col_con[!col_con %in% c(col_cat, col_dis)]
#
#
#
#
#   ## 0. Preparation
#   if (exist_cat) {
#     dict_lev <- dict_level(df, col_cat)
#     dict_name_cat <- dict_onehot(df, col_cat)
#   }
#
#
#   ## 1. Create several datasets. (Resampling)
#   if (resample_method == "bootstrap") {
#     ls_df <- bootsample(df, n_resample)
#   } else if (resample_method == "jackknife") {
#     ls_df <- jacksample(df, n_resample)
#   } else { # if resample_method=='none', there will be no resampling
#     ls_df <- list(df)
#   }
#
#
#   ## 2. Single Imputation
#   ## Input: ls_df list of dataframes after resampling
#   ## Output: ls.imp.onehot, ls.imp.fact two list of imputed dataframe, the only difference is the form of the categorical columns. If there is not any categorical column, then these two lists are identical.
#   ls.imp.onehot <- list()
#   ls.imp.fact <- list()
#   i <- 1
#   if(resample_method == "jackknife"){
#     ls_df[[n_resample+1]] = df
#   }
#   for (dfi in ls_df) {
#     dummy <- dummyVars(" ~ .", data = dfi, sep = "_")
#     tmp.disj <- data.frame(predict(dummy, newdata = dfi))
#     if (imp_method == "missRanger") {
#       res <- missRanger_mod(dfi, col_cat = col_cat, maxiter = maxiter_tree)
#       imp.onehot_i <- res$ximp.disj
#       imp.fact_i <- res$ximp
#     } else if (imp_method == "missForest") {
#       res <- missForest_mod(xmis = dfi, maxiter = maxiter_tree, col_cat = col_cat)
#       imp.onehot_i <- res$ximp.disj
#       imp.fact_i <- res$ximp
#     } else if (imp_method == "kNN") {
#       res <- suppressWarnings((kNN_mod(dfi, col_cat = col_cat, weightDist = TRUE)))
#       imp.onehot_i <- res$ximp.disj
#       imp.fact_i <- res$ximp
#     } else if (imp_method == "EM") {
#       res <- em_mod(dfi, col_cat = col_cat) # Error: matrix not sinugular, the categorical variable with too many categories may ends in a sparse matrix
#       imp.onehot_i <- res$ximp.disj
#       imp.fact_i <- res$ximp
#     } else if (imp_method == "PCA") {
#       if (learn_ncp) {
#         ncp_pca <- estim_ncpFAMD_mod(dfi, method.cv = "Kfold", verbose = F, maxiter = maxiter_pca)$ncp
#       }
#       res <- imputeFAMD_mod(dfi, ncp = ncp_pca, maxiter = maxiter_pca)
#       imp.onehot_i <- res$tab.disj
#       imp.fact_i <- res$completeObs
#     } else if (imp_method == "MI_EM") {
#       res <- MI_EM_amelia(dfi, col_num = c(col_con, col_dis), col_cat = col_cat, num_imp = num_mi)
#       imp.onehot_i <- res$ximp.disj
#       imp.fact_i <- res$ximp
#     } else if (imp_method == "MI_PCA") {
#       res <- MIFAMD_mod(dfi, ncp = ncp_pca, maxiter = maxiter_pca, nboot = num_mi)
#       imp.onehot_i <- res$ximp.disj
#       imp.fact_i <- res$ximp
#     } else if (imp_method == "MICE") {
#       col_cat <- which(!sapply(dfi, is.numeric))
#       exist_cat <- !all(c(0, col_cat) == c(0))
#       if (exist_cat) {
#         name_cat <- colnames(dfi)[col_cat]
#         # Deal with the problem that nlevels(df[[col]]) > length(unique(df[[col]]))
#         for (col in name_cat) {
#           dfi[[col]] <- factor(as.character(dfi[[col]]))
#         }
#         # remember the levels for each categorical column
#         dict_lev <- dict_level(dfi, col_cat)
#         # preserve colnames for ximp.disj
#         dummy <- dummyVars(" ~ .", data = dfi, sep = "_")
#         col_names.disj <- colnames(data.frame(predict(dummy, newdata = dfi)))
#         # represent the factor columns with their ordinal levels
#         dfi <- factor_ordinal_encode(dfi, col_cat)
#       }
#       res0 <- mice::mice(dfi, m = num_mi, maxit = maxiter_mice)
#       res <- result_mice(res0, impnum = num_mi, col_cat = col_cat)
#       # change back to original levels
#       if (exist_cat) {
#         for (col in name_cat) {
#           levels(res$ximp[[col]]) <- dict_lev[[col]]
#         }
#         colnames(res$ximp.disj) <- col_names.disj
#       }
#       imp.onehot_i <- res$ximp.disj
#       imp.fact_i <- res$ximp
#     } else if (imp_method == "MI_Ranger") {
#       res <- MI_missRanger(data.frame(dfi), col_cat = col_cat, num_mi = num_mi, maxiter = maxiter_tree)
#       imp.onehot_i <- res$ximp.disj
#       imp.fact_i <- res$ximp
#     }
#     # } else if (imp_method == "MI_Ranger_bis") {
#     #   res <- MI_missRanger_bis(data.frame(dfi), col_cat = col_cat, num_mi = num_mi, maxiter = maxiter_tree)
#     #   imp.onehot_i <- res$ximp.disj
#     #   imp.fact_i <- res$ximp
#     # }
#     # in case there is one category that doesn't appear in dfi
#     tmp.disj[colnames(imp.onehot_i)] <- imp.onehot_i
#     tmp.disj[is.na(tmp.disj)] <- 0
#     ls.imp.onehot[[i]] <- data.frame(tmp.disj)
#     ls.imp.fact[[i]] <- data.frame(imp.fact_i)
#     i <- i + 1
#   }
#
#   ## 3. Final result
#   if (resample_method == "bootstrap") {
#     if (cat_combine_by == "factor") {
#       ls.imp.tmp <- ls.imp.fact
#       col_cat_boot <- col_cat
#       col_con_boot <- col_con
#       col_dis_boot <- col_dis
#     } else {
#       ls.imp.tmp <- ls.imp.onehot
#       # With onehot form, the column index have changed. Y7 -> Y7_1,...Y7_6
#       num_name <- colnames(df[, c(col_con, col_dis), drop = FALSE])
#       con_name <- colnames(df[, col_con, drop = FALSE])
#       dis_name <- colnames(df[, col_dis, drop = FALSE])
#       col_cat_boot <- which(!colnames(ls.imp.onehot[[1]]) %in% num_name)
#       col_con_boot <- which(colnames(ls.imp.onehot[[1]]) %in% con_name)
#       col_dis_boot <- which(colnames(ls.imp.onehot[[1]]) %in% dis_name)
#     }
#     res <- combine_boot(ls.imp.tmp,
#                         col_con = col_con_boot, col_dis = col_dis_boot, col_cat = col_cat_boot,
#                         num_row_origin = num_row, method = cat_combine_by, dict_cat = dict_name_cat,
#                         var_cat = var_cat
#     )
#   } else if (resample_method == "jackknife") {
#     if (cat_combine_by == "factor") {
#       stop("Please choose cat_combine_by='onehot' when resample_method=='jackknife'.")
#     }
#     # With onehot form, the column index have changed. Y7 -> Y7_1,...Y7_6
#     num_name <- colnames(df[, c(col_con, col_dis), drop = FALSE])
#     con_name <- colnames(df[, col_con, drop = FALSE])
#     dis_name <- colnames(df[, col_dis, drop = FALSE])
#     col_cat_jack <- which(!colnames(ls.imp.onehot[[1]]) %in% num_name)
#     col_con_jack <- which(colnames(ls.imp.onehot[[1]]) %in% con_name)
#     col_dis_jack <- which(colnames(ls.imp.onehot[[1]]) %in% dis_name)
#     ls.df.imp <- ls.imp.onehot[c(1:n_resample)]
#     df_full <- ls.imp.onehot[[n_resample+1]]
#     res <- combine_jack_bis(ls.df.imp,
#                             df_full,
#                         col_con = col_con_jack,
#                         col_dis = col_dis_jack, col_cat = col_cat_jack, method = cat_combine_by,
#                         dict_cat = dict_name_cat, var_cat = var_cat
#     )
#   } else { # if resample_method=='none', there will be no resampling
#     res <- list()
#     res[["df_result"]] <- ls.imp.fact
#     res[["df_result_disj"]] <- ls.imp.onehot
#     res[["df_result_var_disj"]] <- NA
#     res[["df_result_var"]] <- NA
#   }
#
#
#   # Change back the categorical variable levels
#   if (exist_cat) {
#     name_cat <- names(dict_lev)
#     for (name in name_cat) {
#       res$imp[[name]] <- apply(as.array(as.character(res$imp[[name]])), 1, function(x) {
#         unlist(strsplit(x, "_"))[-1]
#       })
#     }
#   }
#
#
#
#
#   ## 4. Evaluation matrix
#   if (!is.null(df_complete)) { # original complete dataset is provided
#     mask <- data.frame(is.na(df))
#     colnames(mask) <- colnames(df)
#     ls.imp <- ls.imp.fact[c(1:n_resample)]
#     df_full <- ls.imp.fact[[n_resample+1]]
#     MSE_imp <- ls_MSE_bis(df_complete, ls.imp.fact, df_full,
#                       mask = mask, col_num_comp = c(col_con, col_dis),
#                       resample_method = resample_method
#     )
#     if (exist_cat && cat_combine_by == "factor") {
#       F1_imp <- ls_F1(df_complete, ls.imp.fact,
#                       mask = mask, col_cat_comp = col_cat, col_cat_imp = col_cat,
#                       resample_method = resample_method, combine_method = cat_combine_by
#       )
#     } else if (exist_cat && cat_combine_by == "onehot") {
#       # With onehot form, the column index have changed. Y7 -> Y7_1,...Y7_6
#       num_name <- colnames(df[, c(col_con, col_dis)])
#       col_cat_F1 <- which(!colnames(ls.imp.onehot[[1]]) %in% num_name)
#       F1_imp <- ls_F1(df_complete, ls.imp.onehot,
#                       mask = mask, col_cat_comp = col_cat, col_cat_imp = col_cat_F1,
#                       resample_method = resample_method, combine_method = cat_combine_by,
#                       dict_cat = dict_name_cat
#       )
#     } else {
#       F1_imp <- NA
#     }
#     res[["MSE"]] <- MSE_imp
#     res[["F1"]] <- F1_imp
#   }
#
#   return(res)
# }
