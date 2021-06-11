single_imp <- function(df, imp_method = "missRanger", resample_method = "bootstrap", 
                       n_resample = 2 * round(log(nrow(df))), col_cat = c(), col_dis = c(), 
                       maxiter_tree = 10, maxiter_pca = 100, ncp_pca = ncol(df) / 2, 
                       learn_ncp = TRUE, cat_combine_by = "factor", var_cat = "wilcox_va",
                       df_complete=NULL) {
  if (all(!is.na(df))) {
    stop("The input dataframe is complete. Imputation is not needed.")
  }
  imp_method <- match.arg(imp_method, c("missRanger", "kNN", "missForest", "PCA", "EM"))
  resample_method <- match.arg(resample_method, c("bootstrap", "jackknife", "none"))
  cat_combine_by <- match.arg(cat_combine_by, c("factor", "onehot"))
  var_cat <- match.arg(var_cat, c("wilcox_va", "unalike"))
  exist_cat <- !all(c(0, col_cat) == c(0))
  num_col <- ncol(df)
  num_row <- nrow(df)
  col_con <- c(1:num_col)
  col_con <- col_con[!col_con %in% c(col_cat, col_dis)]

  ## 0. Preparation
  dict_name_cat <- list()
  if (exist_cat) {
    # remember the levels for each categorical column
    dict_lev <- dict_level(df, col_cat)
    # represent the factor columns with their ordinal levels
    df <- ordinal_encode(df, col_cat)
    df <- factor_encode(df, col_cat)
    # create dictionary for the onehot columns
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
    if (imp_method == "missRanger") {
      res <- missRanger_mod(dfi, col_cat = col_cat, maxiter = maxiter_tree)
      ls.imp.onehot[[i]] <- data.frame(res$ximp.disj)
      ls.imp.fact[[i]] <- data.frame(res$ximp)
    }
    else if (imp_method == "missForest") {
      res <- missForest_mod(xmis = dfi, maxiter = maxiter_tree, col_cat = col_cat)
      ls.imp.onehot[[i]] <- data.frame(res$ximp.disj)
      ls.imp.fact[[i]] <- data.frame(res$ximp)
    }
    else if (imp_method == "kNN") {
      res <- kNN_mod(dfi, col_cat = col_cat, weightDist=TRUE)
      ls.imp.onehot[[i]] <- data.frame(res$ximp.disj)
      ls.imp.fact[[i]] <- data.frame(res$ximp)
    }
    else if (imp_method == "EM") {
      res <- em_mod(dfi, col_cat = col_cat)
      ls.imp.onehot[[i]] <- data.frame(res$ximp.disj)
      ls.imp.fact[[i]] <- data.frame(res$ximp)
    }
    else if (imp_method == "PCA") {
      if (learn_ncp) {
        ncp_pca <- estim_ncpFAMD_mod(dfi, method.cv = "Kfold", verbose = F, maxiter = maxiter_pca)$ncp
      }
      res <- imputeFAMD_mod(dfi, ncp = ncp_pca, maxiter = maxiter_pca)
      ls.imp.onehot[[i]] <- data.frame(res$tab.disj)
      ls.imp.fact[[i]] <- data.frame(res$completeObs)
    }
    i <- i + 1
  }

  ## In case of Jackknife, one imputation for the original incomplete dataset is needed.
  if (resample_method == "jackknife") {
    if (imp_method == "missRanger") {
      res <- missRanger_mod(df, col_cat = col_cat, maxiter = maxiter_tree)
      imp.full.onehot <- data.frame(res$ximp.disj)
    }
    else if (imp_method == "missForest") {
      res <- missForest_mod(xmis = df, maxiter = maxiter_tree, col_cat = col_cat)
      imp.full.onehot <- data.frame(res$ximp.disj)
    }
    else if (imp_method == "kNN") {
      res <- kNN_mod(df, col_cat = col_cat, weightDist=TRUE)
      imp.full.onehot <- data.frame(res$ximp.disj)
    }
    else if (imp_method == "EM") {
      res <- em_mod(df, col_cat = col_cat)
      imp.full.onehot <- data.frame(res$ximp.disj)
    }
    else if (imp_method == "PCA") {
      if (learn_ncp) {
        ncp_pca <- estim_ncpFAMD_mod(df, method.cv = "Kfold", verbose = F, maxiter = maxiter_pca)$ncp
      }
      res <- imputeFAMD_mod(df, ncp = ncp_pca, maxiter = maxiter_pca)
      imp.full.onehot <- data.frame(res$tab.disj)
    }
  }
  


  
  ## 3. Final result
  if (resample_method == "bootstrap") {
    if (cat_combine_by == "factor") {
      ls.imp.tmp <- ls.imp.fact
      col_cat_boot <- col_cat
    }
    else {
      ls.imp.tmp <- ls.imp.onehot
      # With onehot form, the categorical column index has changed. Y7 -> Y7_1,...Y7_6
      col_cat_boot <- c(1:ncol(ls.imp.onehot[[1]]))
      col_cat_boot <- col_cat_boot[!col_cat_boot %in% c(col_con, col_dis)]
    }
    res <- combine_boot(ls.imp.tmp,
      col_con = col_con, col_dis = col_dis, col_cat = col_cat_boot,
      num_row_origin = num_row, method = cat_combine_by, dict_cat = dict_name_cat,
      var_cat = var_cat
    )
  }
  else if (resample_method == "jackknife") {
    if (cat_combine_by == "factor") {
      stop("Please choose cat_combine_by='onehot' when resample_method=='jackknife'.")
    }
    col_cat_jack <- c(1:ncol(ls.imp.onehot[[1]]))
    col_cat_jack <- col_cat_jack[!col_cat_jack %in% c(col_con, col_dis)]
    res <- combine_jack(ls.imp.onehot, imp.full.onehot,
      col_con = col_con,
      col_dis = col_dis, col_cat = col_cat_jack, method = cat_combine_by,
      dict_cat = dict_name_cat, var_cat = var_cat
    )
  }
  else { # if resample_method=='none', there will be no resampling
    res <- list()
    res[[df_result]] <- ls.imp.fact
    res[[df_result_disj]] <- ls.imp.onehot
    res[[df_result_var_disj]] <- NA
    res[[df_result_var]] <- NA
  }
  # Change back the categorical variable levels
  if(exist_cat){
    name_cat = names(dict_lev)
    for(name in name_cat){
      levels(res$imp[[name]]) <- dict_lev[[name]]
    }
  }
  
  
  ## 4. Evaluation matrix
  if(!is.null(df_complete)){ #original complete dataset is provided
    mask <- data.frame(is.na(df))
    colnames(mask) <- colnames(mask)
    df_imp_full <- NULL
    if(resample_method=='jackknife'){
      df_imp_full <- imp.full.onehot
    }
    MSE_imp <- ls_MSE(df_complete, ls.imp.fact, mask = mask, col_num = c(col_con,col_dis), 
                  resample_method = resample_method, df_imp_full=df_imp_full)
    if(exist_cat && cat_combine_by=='factor'){
      F1 <- ls_F1(df_complete, ls.imp.fact, mask = mask, col_cat = col_cat, dict_lev = dict_lev,
                  resample_method = resample_method, combine_method = cat_combine_by)
    }else if (exist_cat && cat_combine_by=='onehot'){
      col_cat_oh <- c(1:length(colnames(ls.imp.onehot[[1]])))
      col_cat_oh <- col_cat_oh[!col_cat_oh %in% c(col_con, col_dis)]
      F1_imp <-ls_F1(df_complete, ls.imp.onehot, mask = mask, col_cat = col_cat_oh, dict_lev = dict_lev,
                 resample_method = resample_method, combine_method = cat_combine_by, 
                 dict_cat = dict_name_cat, df_imp_full=df_imp_full)
    }else{
      F1_imp <- NA
    }
    res[["MSE"]] <- MSE_imp
    res[["F1"]] <- F1_imp
  }
  
  return(res)
}
