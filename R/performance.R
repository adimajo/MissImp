# dens_comp
# @description There are density comparing functions for some of the imputation
# method, but here we try to define the function that could be used for every
# method.
# @param df_comp One dataset.
# @param df_imp Another dataset with the same columns as \code{df_comp}.
# @return The comparison of the density distribution for each column.
# dens_comp <- function(df_comp, df_imp) {
#   ls_col_name <- colnames(df_comp)
#   ls_p <- list()
#   for (ycol in ls_col_name) {
#     dat <- data.frame(
#       "y" = c(df_comp[[ycol]], df_imp[[ycol]]),
#       "lines" = rep(c("complete", "imputed"), each = nrow(df_comp))
#     )
#     p <- ggplot2::ggplot(dat, ggplot2::aes(x = y, fill = lines)) +
#       ggplot2::geom_density(alpha = 0.3) +
#       ggplot2::labs(x = ycol)
#     show(p)
#   }
# }

#' List of MSE
#' @description \code{ls_MSE} is a function that returns a list of MSE
#' corresponding to the given list of imputed datasets.
#' \code{resample_method} is needed because with 'bootstrap' method, we could
#' have repeated lines in the imputed datasets, and with both 'jackknife' and
#' 'bootstrap', the imputed datasets could not cover all the lines.
#' With the purpose of giving every variable the same weight, we scale each
#' variable with the mean and variance calculated from the complete dataset.
#'
#' If the complete and imputed datasets are mix-typed, then only the numerical
#' parts are taken into account.
#' @param df_comp The original complete dataset.
#' @param ls_df_imp List of imputed dataset.
#' @param mask Mask of missingness (1 means missing value and 0 means
#' observed value)
#' @param resample_method Default value is 'bootstrap', could also be
#' 'jackknife' or 'none'.
#' @param col_num_comp Indices of numerical columns in the complete dataset
#' @export
#' @return \code{list_MSE} List of MSE corresponding to the given list of
#' imputed datasets.
#' @return \code{Mean_MSE} Mean value of MSE.
#' @return \code{Variance_MSE} Variance of MSE.
#' @return  \code{list_MSE_scale} List of scaled MSE corresponding to the given
#' list of imputed datasets. Before performing the calculation of MSE,
#' the imputed data set and complete dataset are both scaled with Min-Max scale
#' using the parameter from complete dataset.
#' @return \code{Mean_MSE_scale} Mean value of scaled MSE.
#' @return \code{Variance_MSE} Variance of scaled MSE.
ls_MSE <- function(df_comp,
                   ls_df_imp,
                   mask,
                   col_num_comp,
                   resample_method = "bootstrap") {
  resample_method <- match.arg(resample_method, c(
    "bootstrap",
    "jackknife", "none"
  ))
  ls_mse_result <- c()
  ls_mse_result_scale <- c()
  mask_num <- mask[, col_num_comp] * 1
  df_comp_num <- df_comp[, col_num_comp]
  col_name_num <- colnames(df_comp_num)
  n_sample <- length(ls_df_imp)
  i <- 1
  for (df_imp in ls_df_imp) {
    df_imp_num <- df_imp[col_name_num]
    if (resample_method == "bootstrap") {
      df_imp_num[["index"]] <- as.numeric(row.names(df_imp_num))
      df_imp_num$index <- floor(df_imp_num$index)
      df_num <- stats::aggregate(. ~ index, data = df_imp_num[
        c("index", col_name_num)
      ], mean)
      df_imp_i <- df_num[col_name_num]
      df_comp_i <- df_comp_num[df_num$index, ]
      mask_num_i <- mask_num[df_num$index, ]
    } else if (resample_method == "jackknife") {
      df_imp_num[["index"]] <- as.numeric(row.names(df_imp_num))
      df_imp_i <- df_imp_num[col_name_num]
      df_comp_i <- df_comp_num[df_imp_num$index, ]
      mask_num_i <- mask_num[df_imp_num$index, ]
    } else {
      df_imp_i <- df_imp_num
      df_comp_i <- df_comp_num
      mask_num_i <- mask_num
    }

    mse_result <- sum((df_comp_i - df_imp_i)^2) / sum(mask_num_i)
    ls_mse_result[i] <- mse_result

    # MinMax Scale each variable based on estimated parameters from
    # corresponding complete dataset
    min_comp <- apply(df_comp_i, 2, min)
    max_min_comp <- apply(df_comp_i, 2, function(x) {
      max(x) - min(x)
    })
    df_comp_i_scale <- as.matrix(df_comp_i - min_comp) / t(
      matrix(rep(max_min_comp, nrow(df_comp_i)), nrow = length(max_min_comp))
    )
    df_imp_i_scale <- as.matrix(df_imp_i - min_comp) / t(
      matrix(rep(max_min_comp, nrow(df_comp_i)), nrow = length(max_min_comp))
    )
    mse_result_scale <- sum((df_comp_i_scale - df_imp_i_scale)^2) / sum(
      mask_num_i
    )
    ls_mse_result_scale[i] <- mse_result_scale

    i <- i + 1
  }
  return(list(
    list_MSE = ls_mse_result,
    Mean_MSE = mean(ls_mse_result),
    Variance_MSE = var(ls_mse_result),
    list_MSE_scale = ls_mse_result_scale,
    Mean_MSE_scale = mean(ls_mse_result_scale),
    Variance_MSE_scale = var(ls_mse_result_scale)
  ))
}

# F1
# Input: Original dataset, list of imputed dataset, mask of missingness, column
# index for categorical columns
# Output: A list of F1-scores for each imputation result, a average F1-score,
# variance of F1-score
#' List of MSE
#' @description \code{ls_F1} is a function that returns a list of F1-Score
#' corresponding to the given list of imputed datasets.
#' \code{resample_method} is needed because with 'bootstrap' method, we could
#' have repeated lines in the imputed datasets,
#' and with both 'jackknife' and 'bootstrap', the imputed datasets could not
#' cover all the lines.
#' @param df_comp The original complete dataset.
#' @param ls_df_imp List of imputed dataset.
#' @param mask Mask of missingness (1 means missing value and 0 means observed
#' value).
#' @param col_cat_comp Indices of categorical columns in the complete dataset.
#' @param col_cat_imp Indices of categorical columns in the imputed dataset.
#' @param resample_method Default value is 'bootstrap', could also be
#' 'jackknife' or 'none'.
#' @param combine_method When \code{resample_method} = 'bootstrap',
#' \code{combine_method} could be 'factor' or 'onehot'.
#' When \code{method} = 'onehot', \code{ls_F1} takes the average of the one-hot
#' probability vector for each observation, then choose the position of maximum
#' probability as the predicted category. When \code{method} = 'factor', or
#' each observation, \code{ls_F1} chooses the mode value over the imputed
#' dataframes as the predicted category.
#' @param dict_cat The dictionary of categorical columns names if "onehot"
#' method is applied. For example, it could be list("Y7"=c("Y7_1","Y7_2"),
#' "Y8"=c("Y8_1","Y8_2","Y8_3")).
#' @export
#' @return \code{list_F1} List of F1 corresponding to the given list of
#' imputed datasets.
#' @return \code{Mean_F1} Mean value of F1.
#' @return \code{Variance_F1} Variance of F1.
ls_F1 <- function(df_comp,
                  ls_df_imp,
                  mask,
                  col_cat_comp,
                  col_cat_imp,
                  resample_method = "bootstrap",
                  combine_method = "onehot",
                  dict_cat = NULL) {
  resample_method <- match.arg(resample_method, c(
    "bootstrap",
    "jackknife", "none"
  ))
  combine_method <- match.arg(combine_method, c("onehot", "factor"))
  # if (resample_method == "jackknife" && is.null(df_imp_full)) {
  #   stop("With jackknife resampling method, df_imp_full is required.\n")
  # }
  if (resample_method == "jackknife" && combine_method == "factor") {
    stop("With jackknife resampling method, combine_method could only be
         'onehot'.\n")
  }
  if (resample_method != "none" && combine_method == "onehot" && is.null(
    dict_cat
  )) {
    stop("With onehot combining method, dict_cat is needed.\n")
  }
  dict_lev <- dict_level(df_comp, col_cat_comp)

  ls_f1_result <- c()
  i <- 1

  n_sample <- length(ls_df_imp)

  for (df_imp in ls_df_imp) {
    # Take only the categorical part
    df_imp_cat <- df_imp[, col_cat_imp]
    col_name_cat <- colnames(df_imp_cat)

    if (resample_method == "bootstrap") {
      df_imp_cat[["index"]] <- as.numeric(row.names(df_imp_cat))
      df_imp_cat$index <- floor(df_imp_cat$index)
      if (combine_method == "onehot") { # df_imp_cat in form of onehot
        # take average of the onehot probabilities for the doublons
        df_cat <- stats::aggregate(. ~ index, data = df_imp_cat[
          c("index", col_name_cat)
        ], mean)
        # convert onehot to factor form
        names_cat <- names(dict_cat)
        for (name in names_cat) {
          df_cat[[name]] <- apply(
            df_cat[dict_cat[[name]]], 1,
            which_max_cat, name, dict_cat
          )
          df_cat[[name]] <- unlist(df_cat[[name]])
          df_cat[[name]] <- apply(as.array(df_cat[[name]]), 1, function(x) {
            unlist(strsplit(x, "_"))[-1]
          })
          df_cat[[name]] <- factor(df_cat[[name]])
        }
        df_imp_i <- df_cat[names_cat]
        df_comp_i <- df_comp[df_cat$index, ][names_cat]
        mask_cat_i <- mask[df_cat$index, ][names_cat] * 1
      } else { # df_imp_cat in form of factor
        df_cat <- df_imp_cat[c("index", col_name_cat)] %>%
          dplyr::group_by(.data$index) %>%
          dplyr::summarise(dplyr::across(dplyr::all_of(col_name_cat), Mode_cat))
        names_cat <- names(dict_lev)
        for (name in names_cat) {
          df_cat[[name]] <- apply(as.array(df_cat[[name]]), 1, function(x) {
            unlist(strsplit(x, "_"))[-1]
          })
          df_cat[[name]] <- factor(df_cat[[name]])
        }
        df_imp_i <- df_cat[col_name_cat]
        df_comp_i <- df_comp[df_cat$index, col_cat_comp]
        mask_cat_i <- mask[df_cat$index, col_cat_comp] * 1
      }
    } else if (resample_method == "jackknife") {
      # convert onehot to factor form
      df_imp_cat[["index"]] <- as.numeric(row.names(df_imp_cat))
      names_cat <- names(dict_cat)
      for (name in names_cat) {
        df_imp_cat[[name]] <- apply(
          df_imp_cat[dict_cat[[name]]], 1,
          which_max_cat, name, dict_cat
        )
        df_imp_cat[[name]] <- unlist(df_imp_cat[[name]])
        df_imp_cat[[name]] <- apply(
          as.array(df_imp_cat[[name]]), 1,
          function(x) {
            unlist(strsplit(x, "_"))[-1]
          }
        )
        df_imp_cat[[name]] <- factor(df_imp_cat[[name]])
      }
      df_imp_i <- df_imp_cat[names_cat]
      df_comp_i <- df_comp[df_imp_cat$index, ][names_cat]
      mask_cat_i <- mask[df_imp_cat$index, ][names_cat] * 1
    } else {
      df_imp_i <- df_imp_cat
      df_comp_i <- df_comp[, col_cat_comp]
      mask_cat_i <- data.frame(mask[, col_cat_comp] * 1)
    }

    # Change to one vector
    ls_f1 <- c()
    for (col in colnames(df_comp_i)) {
      mask_concat <- as.vector(mask_cat_i[[col]])
      df_comp_concat <- as.vector(df_comp_i[[col]])
      df_comp_true <- df_comp_concat[mask_concat == 1]
      df_imp_concat <- as.vector(df_imp_i[[col]])
      df_imp_predict <- df_imp_concat[mask_concat == 1]
      ls_f1 <- c(ls_f1, F1_Score_micro(df_comp_true, df_imp_predict))
    }
    # micro f1 score:
    # https://sebastianraschka.com/faq/docs/multiclass-metric.html
    ls_f1_result[i] <- mean(ls_f1)
    i <- i + 1
  }
  return(list(
    list_F1 = ls_f1_result,
    Mean_F1 = mean(ls_f1_result),
    Variance_F1 = var(ls_f1_result)
  ))
}

# Version jackknife estimate using estimation on the full incomplete dataset
# ls_MSE_bis <- function(df_comp,
#                    ls_df_imp,
#                    df_full,
#                    mask,
#                    col_num_comp,
#                    resample_method = "bootstrap") {
#   resample_method <- match.arg(resample_method, c("bootstrap",
# "jackknife", "none"))
#   ls_mse_result <- c()
#   ls_mse_result_scale <- c()
#   mask_num <- mask[, col_num_comp] * 1
#   col_names <- colnames(df_comp)
#   df_comp_num <- df_comp[, col_num_comp]
#   col_name_num <- colnames(df_comp_num)
#   df_full_num <- df_full[col_name_num]
#   n_sample <- length(ls_df_imp)
#   i <- 1
#   for (df_imp in ls_df_imp) {
#     df_imp_num <- df_imp[col_name_num]
#     df_imp_num[["index"]] <- as.numeric(row.names(df_imp_num))
#     if (resample_method == "bootstrap") {
#       df_imp_num$index <- floor(df_imp_num$index)
#       df_num <- stats::aggregate(. ~ index, data = df_imp_num[c("index",
# col_name_num)], mean)
#       df_imp_i <- df_num[col_name_num]
#       df_comp_i <- df_comp_num[df_num$index, ]
#       mask_num_i <- mask_num[df_num$index, ]
#     } else if (resample_method == "jackknife") {
#       df_full_num$index <- as.numeric(row.names(df_full_num))
#       df_imp_i <- n_sample * df_full_num[df_imp_num$index, ] - (
# n_sample - 1) * df_imp_num
#       df_imp_i <- df_imp_i[col_name_num]
#       df_comp_i <- df_comp_num[df_imp_num$index, ]
#       mask_num_i <- mask_num[df_imp_num$index, ]
#     } else {
#       df_imp_i <- df_imp_num
#       df_comp_i <- df_comp_num
#       mask_num_i <- mask_num
#     }
#
#     mse_result <- sum((df_comp_i - df_imp_i)^2) / sum(mask_num_i)
#     ls_mse_result[i] <- mse_result
#
#     # MinMax Scale each variable based on estimated parameters from
#     # corresponding complete dataset
#     min_comp <- apply(df_comp_i, 2, min)
#     max_min_comp <- apply(df_comp_i, 2, function(x) {
#       max(x) - min(x)
#     })
#     df_comp_i_scale <- as.matrix(df_comp_i - min_comp) / t(
# matrix(rep(max_min_comp, nrow(df_comp_i)), nrow = length(max_min_comp)))
#     df_imp_i_scale <- as.matrix(df_imp_i - min_comp) / t(
# matrix(rep(max_min_comp, nrow(df_comp_i)), nrow = length(max_min_comp)))
#     mse_result_scale <- sum((df_comp_i_scale - df_imp_i_scale)^2) / sum
# (mask_num_i)
#     ls_mse_result_scale[i] <- mse_result_scale
#
#     i <- i + 1
#   }
#   return(list(
#     list_MSE = ls_mse_result,
#     Mean_MSE = mean(ls_mse_result),
#     Variance_MSE = var(ls_mse_result),
#     list_MSE_scale = ls_mse_result_scale,
#     Mean_MSE_scale = mean(ls_mse_result_scale),
#     Variance_MSE_scale = var(ls_mse_result_scale)
#   ))
# }
