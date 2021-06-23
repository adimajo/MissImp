# There are density comparing functions for some of the imputation method, but here we try to define the function that could be used for every method
# Input: two datasets with same column names
# Output: the comparison of the density distribution for each column
dens_comp <- function(df_comp, df_imp) {
  ls_col_name <- colnames(df_comp)
  ls_p <- list()
  for (ycol in ls_col_name) {
    dat <- data.frame(
      y = c(df_comp[[ycol]], df_imp[[ycol]]),
      lines = rep(c("complete", "imputed"), each = 100)
    )
    p <- ggplot(dat, aes(x = y, fill = lines)) +
      geom_density(alpha = 0.3) +
      labs(x = ycol)
    show(p)
  }
}

# MSE (average over missing number)
# As we did a scale operation on the complete dataset, we may need to rescale the output imputed data set then calculate the MSE
# Input: Original dataset, list of imputed dataset, mask of missingness, column index for numerical columns
# Output: A list of MSE values for each imputation result, a average MSE value, variance of MSE
ls_MSE <- function(df_comp, ls_df_imp, mask, col_num, resample_method = "bootstrap", df_imp_full = NULL) {
  resample_method <- match.arg(resample_method, c("bootstrap", "jackknife", "none", "mi"))
  if (resample_method == "jackknife" && is.null(df_imp_full)) {
    stop("With jackknife resampling method, df_imp_full is required.\n")
  }
  ls_mse_result <- c()
  mask_num <- mask[, col_num] * 1
  df_comp_num <- df_comp[, col_num]
  if (resample_method == "jackknife") {
    df_imp_full_num <- df_imp_full[, col_num]
  }
  col_name_num <- colnames(df_comp_num)
  n_sample <- length(ls_df_imp)
  i <- 1
  for (df_imp in ls_df_imp) {
    df_imp_num <- df_imp[, col_num]
    df_imp_num$index <- as.numeric(row.names(df_imp_num))
    if (resample_method == "bootstrap") {
      df_imp_num$index <- floor(df_imp_num$index)
      df_num <- stats::aggregate(. ~ index, data = df_imp_num[c("index", col_name_num)], mean)
      df_imp_i <- df_num[col_name_num]
      df_comp_i <- df_comp_num[df_num$index, ]
      mask_num_i <- mask_num[df_num$index, ]
    }
    else if (resample_method == "jackknife") {
      # df_imp_full_num$index <- as.numeric(row.names(df_imp_full_num))
      # df_minus <- n_sample * df_imp_full_num[df_imp_num$index, ]
      # -(n_sample - 1) * df_imp_num
      df_minus <- df_imp_num
      df_minus$index <- df_imp_num$index
      df_imp_i <- df_minus[col_name_num]
      df_comp_i <- df_comp_num[df_minus$index, ]
      mask_num_i <- mask_num[df_minus$index, ]
    }
    else {
      df_imp_i <- df_imp_num
      df_comp_i <- df_comp_num
      mask_num_i <- mask_num
    }

    mse_result <- (sqrt(sum((as.matrix(df_comp_i) * mask_num_i - as.matrix(df_imp_i) * mask_num_i)^2) / sum(mask_num_i)))
    ls_mse_result[i] <- mse_result
    i <- i + 1
  }
  return(list(
    "list_MSE" = ls_mse_result,
    "Mean_MSE" = mean(ls_mse_result),
    "Variance_MSE" = var(ls_mse_result)
  ))
}


# F1
# Input: Original dataset, list of imputed dataset, mask of missingness, column index for categorical columns
# Output: A list of F1-scores for each imputation result, a average F1-score, variance of F1-score
ls_F1 <- function(df_comp, ls_df_imp, mask, col_cat, dict_lev, resample_method = "bootstrap", combine_method = "onehot", df_imp_full = NULL, dict_cat = NULL) {
  resample_method <- match.arg(resample_method, c("bootstrap", "jackknife", "none", "mi"))
  combine_method <- match.arg(combine_method, c("onehot", "factor"))
  if (resample_method == "jackknife" && is.null(df_imp_full)) {
    stop("With jackknife resampling method, df_imp_full is required.\n")
  }
  if (resample_method == "jackknife" && combine_method == "factor") {
    stop("With jackknife resampling method, combine_method could only be 'onehot'.\n")
  }
  if (combine_method == "onehot" && is.null(dict_cat)) {
    stop("With onehot combining method, dict_cat is needed.\n")
  }

  ls_f1_result <- c()
  i <- 1

  n_sample <- length(ls_df_imp)

  for (df_imp in ls_df_imp) {
    # Take only the categorical part
    df_imp_cat <- df_imp[, col_cat]
    col_name_cat <- colnames(df_imp_cat)
    df_imp_cat$index <- as.numeric(row.names(df_imp_cat))
    if (resample_method == "bootstrap") {
      df_imp_cat$index <- floor(df_imp_cat$index)
      if (combine_method == "onehot") { # df_imp_cat in form of onehot
        # take average of the onehot probabilities for the doublons
        df_cat <- stats::aggregate(. ~ index, data = df_imp_cat[c("index", col_name_cat)], mean)
        # convert onehot to factor form
        names_cat <- names(dict_cat)
        which_max_cat <- function(x, name) {
          return(dict_cat[[name]][which.max.random(x)])
        }
        for (name in names_cat) {
          df_cat[[name]] <- apply(df_cat[dict_cat[[name]]], 1, which_max_cat, name)
          df_cat[[name]] <- unlist(df_cat[[name]])
          df_cat[[name]] <- factor(df_cat[[name]])
          levels(df_cat[[name]]) <- dict_lev[[name]]
        }
        df_imp_i <- df_cat[names_cat]
        df_comp_i <- df_comp[df_cat$index, ][names_cat]
        mask_cat_i <- mask[df_cat$index, ][names_cat] * 1
      }
      else { # df_imp_cat in form of factor
        df_cat <- stats::aggregate(. ~ index, data = df_imp_cat[c("index", col_name_cat)], Mode_cat)
        names_cat <- names(dict_lev)
        for (name in names_cat) {
          levels(df_cat[[name]]) <- dict_lev[[name]]
        }
        df_imp_i <- df_cat[col_name_cat]
        df_comp_i <- df_comp[df_cat$index, col_cat]
        mask_cat_i <- mask[df_cat$index, col_cat] * 1
      }
    }
    else if (resample_method == "jackknife") {
      df_imp_full_cat <- df_imp_full[col_name_cat]
      df_imp_full_cat$index <- as.numeric(row.names(df_imp_full_cat))
      df_minus <- n_sample * df_imp_full_cat[which(df_imp_full_cat$index %in% df_imp_cat$index), ]
      -(n_sample - 1) * df_imp_cat
      df_minus$index <- df_imp_cat$index
      # convert onehot to factor form
      names_cat <- names(dict_cat)
      which_max_cat <- function(x, name) {
        return(dict_cat[[name]][which.max.random(x)])
      }
      for (name in names_cat) {
        df_minus[[name]] <- apply(df_minus[dict_cat[[name]]], 1, which_max_cat, name)
        df_minus[[name]] <- unlist(df_minus[[name]])
        df_minus[[name]] <- factor(df_minus[[name]])
        levels(df_minus[[name]]) <- dict_lev[[name]]
      }
      df_imp_i <- df_minus[names_cat]
      df_comp_i <- df_comp[df_minus$index, ][names_cat]
      mask_cat_i <- mask[df_minus$index, ][names_cat] * 1
    }
    else {
      df_imp_i <- df_imp_cat
      df_comp_i <- df_comp[, col_cat]
      mask_cat_i <- mask_cat[, col_cat]
    }

    # Change to one vector
    mask_concat <- as.vector(as.matrix(mask_cat_i))
    df_comp_concat <- as.vector(as.matrix(df_comp_i))
    df_comp_true <- df_comp_concat[mask_concat == 1]
    df_imp_concat <- as.vector(as.matrix(df_imp_i))
    df_imp_predict <- df_imp_concat[mask_concat == 1]
    # micro f1 score : https://sebastianraschka.com/faq/docs/multiclass-metric.html
    f1_result <- F1_Score_micro(factor(df_comp_true), factor(df_imp_predict))
    ls_f1_result[i] <- mean(f1_result)
    i <- i + 1
  }
  return(list(
    "list_F1" = ls_f1_result,
    "Mean_F1" = mean(ls_f1_result),
    "Variance_F1" = var(ls_f1_result)
  ))
}
