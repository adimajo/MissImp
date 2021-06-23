#' jacksample: create several dataframes by jackknife resampling
#'
#' @description
#' \code{jacksample} function creates \code{num_sample} samples without the group of observation k (k=1,...,\code{num_sample}).
#' We regroup the original data into \code{num_sample} groups by their index.
#' @param df A complete or incomplete dataframe
#' @param num_sample Number of jackknifed samples.
#' @return A list of jackknifed dataframes.
#' @export
#' @references Statistical Analysis with Missing Data, by Little and Rubin, 2002
#' @examples
#' n <- 10000
#' mu.X <- c(1, 2, 3)
#' Sigma.X <- matrix(c(9, 3, 2, 3, 4, 0, 2, 0, 1), nrow = 3)
#' X.complete.cont <- MASS::mvrnorm(n, mu.X, Sigma.X)
#' rs <- generate_miss(X.complete.cont, 0.5, mechanism = "MNAR2")
#' ls_boot <- jacksample(rs$X.incomp, 4)
jacksample <- function(df, num_sample) {
  num_row <- nrow(df)
  row_tranch <- num_row %/% num_sample
  rest <- num_row %% num_sample
  ls_df_new <- list()
  i <- 1
  start <- 1
  while (i <= num_sample - rest) {

    # num_sample-rest samples with row_tranch rows
    end <- start + row_tranch - 1
    if (start == 1) {
      df_new <- df[(end + 1):num_row, ]
    }
    else if (end == num_row) {
      df_new <- df[1:(start - 1), ]
    }
    else {
      df_new <- df[c(1:(start - 1), (end + 1):num_row), ]
    }
    # df_new[["old_idx"]] = c(start:end)
    ls_df_new[[i]] <- df_new
    i <- i + 1
    start <- end + 1
  }

  while (i <= num_sample) {
    # rest samples with row_tranch+1 rows
    end <- start + (row_tranch + 1) - 1
    if (end >= num_row) {
      df_new <- df[c(1:(start - 1)), ]
    }
    else {
      df_new <- df[c(1:(start - 1), (end + 1):num_row), ]
    }

    # df_new[["old_idx"]] = c(start:end)
    ls_df_new[[i]] <- df_new
    i <- i + 1
    start <- end + 1
  }
  return(ls_df_new)
}





#' combine_jack: combine imputed jackknife datasets
#'
#' @description
#' \code{combine_jack} function combines several imputed jackknifed dataframes
#' into the final imputed dataframe and provide the variance for each imputed value.
#' @param ls_df A list of imputed jackknifed dataframes.
#' @param df_full Imputation result on the original dataframe before jackknife resampling.
#' @param col_con Continous columns index.
#' @param col_dis Discret columns index.
#' @param col_cat Categorical columns index.
#' @param method The encoded method of categorical columns in the imputed dataframes.
#' This function is only coded for "onehot" situation.
#' @param dict_cat The dictionary of categorical columns names if "onehot" method is applied.
#' For example, it could be list("Y7"=c("Y7_1","Y7_2"), "Y8"=c("Y8_1","Y8_2","Y8_3")).
#' @param var_cat The method of variance calculation for the categorical columns.
#' "unalike" will lead to the calculation of unalikeability, while "wilcox_va" will lead to the calculation of Wilcox index: VarNC.
#' @return \code{df_result_disj} The final imputed dataframe with the categorical columns in onehot form.
#' @return \code{df_result_var_disj} The variance matrix for the final imputation dataframe with the categorical columns in onehot form.
#' @return \code{df_result} The final imputed dataframe with the categorical columns in factor form.
#' @return \code{df_result_var} The variance matrix for the final imputation dataframe with the categorical columns in factor form.
#' @export
#' @references Statistical Analysis with Missing Data, by Little and Rubin, 2002

combine_jack <- function(ls_df,
                         df_full,
                         col_con,
                         col_dis = c(),
                         col_cat = c(),
                         method = "onehot",
                         dict_cat = NULL,
                         var_cat = "unalike") {
  ls_df_minus <- list()
  ls_col_name <- colnames(ls_df[[1]])
  exist_dis <- !all(c(0, col_dis) == c(0))
  if (exist_dis) {
    col_name_dis <- ls_col_name[col_dis]
  }
  exist_cat <- !all(c(0, col_cat) == c(0))
  if (exist_cat) {
    col_name_cat <- ls_col_name[col_cat]
    var_cat <- match.arg(var_cat, c("unalike", "wilcox_va"))
    is_unalike <- (var_cat == "unalike")
    if (is.null(dict_cat)) {
      stop("dict_cat is needed when there are onehot categorical columns")
    }
  }
  i <- 1
  n_sample <- length(ls_df)
  df_full$index <- as.numeric(row.names(df_full))
  while (i <= n_sample) {
    ls_df[[i]]$index <- as.numeric(row.names(ls_df[[i]]))
    # only for numeric columns
    ls_df_minus[[i]] <- n_sample * df_full[ls_df[[i]]$index, ] - (n_sample -
      1) * ls_df[[i]]
    if (exist_cat) {
      ls_df_minus[[i]][col_name_cat][ls_df_minus[[i]][col_name_cat] < 0] <- 0
    }
    ls_df_minus[[i]]$index <- ls_df[[i]]$index
    i <- i + 1
  }
  # Put all imputed datasets together
  df_new_merge <- Reduce(function(dtf1, dtf2) {
    rbind(dtf1, dtf2)
  }, ls_df_minus)
  # Add the combined result for categorical variables deducted from the onehot result
  if (exist_cat) {
    names_cat <- names(dict_cat)
    which_max_cat <- function(x, name) {
      return(dict_cat[[name]][which.max.random(x)])
    }
    for (name in names_cat) {
      df_new_merge[[name]] <- apply(df_new_merge[dict_cat[[name]]], 1, which_max_cat, name)
      df_new_merge[[name]] <- unlist(df_new_merge[[name]])
      df_new_merge[[name]] <- factor(df_new_merge[[name]])
    }
  }



  # By Jackknife combining rules for disjunctive part
  df_new <- stats::aggregate(. ~ index, data = df_new_merge[c("index", ls_col_name)], mean) # Here for each indice, there are n_sample-1 values
  if (exist_dis) {
    df_new[col_name_dis] <- round(df_new[col_name_dis]) # round for the discret variables
  }
  df_new_var_disj <- stats::aggregate(. ~ index, data = df_new_merge[c("index", ls_col_name)], var)
  df_new_var_disj[c(ls_col_name)] <- df_new_var_disj[c(ls_col_name)] / (n_sample -
    1)


  # Final result according to the mean of onehot probability
  if (exist_cat) {
    for (name in names_cat) {
      df_new[[name]] <- apply(df_new[dict_cat[[name]]], 1, which_max_cat, name)
      df_new[[name]] <- unlist(df_new[[name]])
      df_new[[name]] <- factor(df_new[[name]])
    }
    df_result <- df_new[-c(1, col_cat + 1)] # remove the index and the onehot columns
    if (is_unalike) {
      df_new_cat_var <- stats::aggregate(. ~ index, data = df_new_merge[c("index", names_cat)], uwo4419::unalike)
    }
    else {
      df_new_cat_var <- df_new_merge[c("index", names_cat)] %>%
        dplyr::group_by(index) %>%
        dplyr::summarise(across(all_of(names_cat), VA_fact))
    }

    df_new_var <- merge(df_new_var_disj, df_new_cat_var, by = "index")
    df_result_var <- df_new_var[-c(1, col_cat + 1)] # remove the index and the onehot columns


    # remove index column and the combined categorical columns
    df_result_disj <- df_new[, !(names(df_new) %in% c(names_cat, "index"))]
    df_result_var_disj <- df_new_var[, !(names(df_new) %in% c(names_cat, "index"))]
    return(
      list(
        "imp.disj" = df_result_disj,
        "uncertainty.disj" = df_result_var_disj,
        "imp" = df_result,
        "uncertainty" = df_result_var
      )
    )
  }
  else {
    return(list("imp" = df_new[-c(1)], "uncertainty" = df_new_var_disj[-c(1)]))
  }
}
