#' bootsample: create several dataframes by bootstrap resampling
#'
#' @description
#' \code{bootsample} function generates several dataframes by drawing samples of size n by simple random sampling with replacement.
#'
#' In average, if we want to cover all the n rows in the original dataframe,
#' the number of samples should be greater than log(n).
#' @param df A complete or incomplete dataframe
#' @param num_sample Number of bootstrapped samples. We suggest that num_sample > Log(n), where n is the number of rows.
#' @return A list of bootstrapped dataframes.
#' @export
#' @references Statistical Analysis with Missing Data, by Little and Rubin, 2002
#' @examples
#' n <- 10000
#' mu.X <- c(1, 2, 3)
#' Sigma.X <- matrix(c(9, 3, 2, 3, 4, 0, 2, 0, 1), nrow = 3)
#' X.complete.cont <- MASS::mvrnorm(n, mu.X, Sigma.X)
#' rs <- generate_miss(X.complete.cont, 0.5, mechanism = "MNAR2")
#' ls_boot <- bootsample(rs$X.incomp, 4)
bootsample <- function(df, num_sample) {
  num_row <- nrow(df)
  if (num_sample < log(num_row)) {
    warning(
      "Warning: We suggest that the bootstrap number > log(number of rows), if not, all the rows may not be covered in the bootstrap samples.\n"
    )
  }
  ls_df_new <- list()
  i <- 1
  ls_idx <- c()
  while (i <= num_sample - 1) {
    chosen_idx <- gdata::resample(c(1:num_row), num_row, replace = TRUE)
    ls_idx <- c(ls_idx, chosen_idx)
    ls_df_new[[i]] <- df[chosen_idx, ]
    i <- i + 1
  }
  # last sample
  ls_idx_miss <- c(1:num_row)
  ls_idx_miss <- ls_idx_miss[!ls_idx_miss %in% unique(ls_idx)]
  chosen_idx <- gdata::resample(c(1:num_row), num_row - length(ls_idx_miss), replace = TRUE)
  chosen_idx <- c(chosen_idx, ls_idx_miss)
  ls_df_new[[i]] <- df[chosen_idx, ]

  return(ls_df_new)
}


#' combine_boot: combine imputed bootstrap datasets
#'
#' @description
#' \code{combine_boot} function combines several imputed bootstrapped dataframes
#' into the final imputed dataframe and provide the variance for each imputed value.
#' @param ls_df A list of imputed bootstrapped dataframes.
#' @param col_con Continous columns index.
#' @param col_dis Discret columns index.
#' @param col_cat Categorical columns index.
#' @param num_row_origin Number of rows in the original incomplete dataframe before bootstrapping.
#' @param method The encoded method of categorical columns in the imputed dataframes.
#' This function is coded for both "onehot" and "factor" encoded situations.
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
combine_boot <- function(ls_df,
                         col_con,
                         col_dis = c(),
                         col_cat = c(),
                         num_row_origin,
                         method = "onehot",
                         dict_cat = NULL,
                         var_cat = "unalike") {
  method <- match.arg(method, c("onehot", "factor"))
  var_cat <- match.arg(var_cat, c("unalike", "wilcox_va"))
  is_unalike <- (var_cat == "unalike")
  is_onehot <- (method == "onehot")
  if (is_onehot & is.null(dict_cat)) {
    stop("dict_cat is needed when the method is onehot.")
  }
  ls_df_new <- list()

  ls_col_name <- colnames(ls_df[[1]])
  col_num <- c(col_con, col_dis)
  exist_dis <- !all(c(0, col_dis) == c(0))
  if (exist_dis) {
    col_name_dis <- ls_col_name[col_dis]
  }
  exist_cat <- !all(c(0, col_cat) == c(0))
  if (exist_cat) {
    col_name_cat <- ls_col_name[col_cat]
    if (is_onehot) {
      col_name_cat <- ls_col_name[col_cat]
      # We treat everything as numeric if the categorical variables are encoded by onehot probability
      col_num <- c(col_num, col_cat)
    }
  }
  col_name_num <- ls_col_name[col_num]
  # Deal with doublets in each imputed dataset
  i <- 1
  for (df in ls_df) {
    df$index <- as.numeric(row.names(df))
    df$index <- floor(df$index)
    df_num <- stats::aggregate(. ~ index, data = df[c("index", col_name_num)], mean)
    if (exist_cat && !is_onehot) {
      df_cat <- stats::aggregate(. ~ index, data = df[c("index", col_name_cat)], Mode_cat)
      ls_df_new[[i]] <- merge(df_num, df_cat, by = "index")
    }
    else {
      ls_df_new[[i]] <- df_num
    }

    # the expected number of rows is n(1-(n-1)^k/n^k) where k = n
    # expectation of unique rows for sampling with replacement
    # We need a varaince matrix here to pass ? No, see https://stats.stackexchange.com/questions/399382/bootstrap-rubins-rules-and-uncertainty-of-sub-estimates

    i <- i + 1
  }

  # Combine all imputed datasets together
  df_new_merge <- Reduce(function(dtf1, dtf2) {
    rbind(dtf1, dtf2)
  }, ls_df_new)
  # Add the combined result for categorical variables deducted from the onehot result
  if (exist_cat && is_onehot) {
    names_cat <- names(dict_cat)
    which_max_cat <- function(x, name) {
      return(dict_cat[[name]][which.max(x)])
    }
    for (name in names_cat) {
      df_new_merge[[name]] <- apply(df_new_merge[dict_cat[[name]]], 1, which_max_cat, name)
      df_new_merge[[name]] <- unlist(df_new_merge[[name]])
      df_new_merge[[name]] <- factor(df_new_merge[[name]])
    }
  }

  # By Bootstrap combining rules
  df_new_mean_num <- stats::aggregate(. ~ index, data = df_new_merge[c("index", col_name_num)], mean)
  df_new_num_var <- stats::aggregate(. ~ index, data = df_new_merge[c("index", col_name_num)], var)
  if (exist_dis) {
    df_new_mean_num[col_name_dis] <- round(df_new_mean_num[col_name_dis]) # for the discret variables
  }

  if (exist_cat && !is_onehot) {
    df_new_merge <- factor_encode(df_new_merge, col_cat + 1)
    df_new_mode_cat <- stats::aggregate(. ~ index, data = df_new_merge[c("index", col_name_cat)], Mode_cat)
    df_new <- merge(df_new_mean_num, df_new_mode_cat, by = "index")
    df_new <- factor_encode(df_new, col_cat + 1) # +1 is for the index column
    if (is_unalike) {
      df_new_cat_var <- stats::aggregate(. ~ index, data = df_new_merge[c("index", col_name_cat)], uwo4419::unalike)
    }
    else {
      df_new_cat_var <- df_new_merge[c("index", col_name_cat)] %>%
        dplyr::group_by(index) %>%
        dplyr::summarise(across(col_name_cat, VA_fact))
      # df_new_cat_var = aggregate(.~index , data =df_new_merge[c("index",col_name_cat)], VA_fact)
    }
    df_new_var <- merge(df_new_num_var, df_new_cat_var, by = "index")
    # We need to use unalikeability to mesure the uncertainty of the categorical variables
  }
  else {
    df_new <- df_new_mean_num
    df_new_var <- df_new_num_var
  }
  # Add the combined result for categorical variables deducted from the onehot result to df_new
  if (is_onehot && exist_cat) {
    for (name in names_cat) {
      df_new[[name]] <- apply(df_new[dict_cat[[name]]], 1, which_max_cat, name)
      df_new[[name]] <- unlist(df_new[[name]])
    }
    if (is_unalike) {
      df_new_cat_var <- stats::aggregate(. ~ index, data = df_new_merge[c("index", names_cat)], uwo4419::unalike)
    }
    else {
      df_new_cat_var <- df_new_merge[c("index", names_cat)] %>%
        dplyr::group_by(index) %>%
        dplyr::summarise(across(names_cat, VA_fact))
      # df_new_cat_var = aggregate(.~index , simplify=FALSE, data =df_new_merge[c("index",names_cat)], VA_fact)
    }
    df_new_var <- merge(df_new_num_var, df_new_cat_var, by = "index")
  }

  # Add the uncovered rows, fill in with NA
  num_row_new <- nrow(df_new)
  if (num_row_new < num_row_origin) {
    print(paste0("Covered number of rows: ", num_row_new))
    print(paste0("Original number of rows: ", num_row_origin))
    warning(
      "Warning: All the rows in the original dataset is not selected in the bootstrap samples. An increase on the number of bootstrap samples is suggested.\n"
    )
  }
  df_idx <- data.frame("index" = c(1:num_row_origin))
  df_result_with_idx <- merge(
    x = df_idx,
    y = df_new,
    by = "index",
    all = TRUE
  )
  df_result_var_with_idx <- merge(
    x = df_idx,
    y = df_new_var,
    by = "index",
    all = TRUE
  )

  # remove index column and the combined categorical columns
  if (is_onehot && exist_cat) {
    df_result_disj <- df_result_with_idx[, !(names(df_result_with_idx) %in% c(names_cat, "index"))]
    df_result_var_disj <- df_result_var_with_idx[, !(names(df_result_with_idx) %in% c(names_cat, "index"))]
  }
  else {
    df_result_disj <- df_result_with_idx[-c(1)]
    df_result_var_disj <- df_result_var_with_idx[-c(1)]
  }



  # Final result for the categorical variables
  if (is_onehot) {
    df_result <- df_result_with_idx[-c(1, col_cat + 1)] # remove index and onehot columns, +1 is for the adjustment
    df_result_var <- df_result_var_with_idx[-c(1, col_cat + 1)]
  }
  else {
    df_result <- df_result_disj
    df_result_var <- df_result_var_disj
  }



  return(
    list(
      "imp.disj" = df_result_disj,
      "uncertainty.disj" = df_result_var_disj,
      "imp" = df_result,
      "uncertainty" = df_result_var
    )
  )
}
