


# Generate missing values with different mechanisms
# Input: A complete dataset, the mechanism and the proportion of missingness
# Output: A dataset with missing values according to the mechanism and the proportion of missingness
generate_miss <- function(df,
                          miss_perc,
                          mechanism = "MCAR",
                          # c("MCAR", "MAR1", "MAR2","MAR3","MNAR1","MNAR2")
                          # We could add here the parameters for produce_NA function in MAR1 and MNAR1
                          mar2.col.ctrl = 1) {
  mechanism <-
    match.arg(mechanism, c("MCAR", "MAR1", "MAR2", "MAR3", "MNAR1", "MNAR2"))
  ls_col_name <- colnames(df)
  num_col <- length(ls_col_name)
  if (mechanism == "MCAR") {
    # Bernouilli
    mcar <- produce_NA(df, mechanism = "MCAR", perc.missing = miss_perc)
    X.mcar <- mcar$data.incomp
    R.mcar <- data.frame(mcar$idx_newNA)
    real_miss_perc <- sum(R.mcar * 1) / prod(dim(R.mcar * 1))
    return(list(
      "X.incomp" = X.mcar,
      "R.mask" = R.mcar,
      "real_miss_perc" = real_miss_perc
    ))
  }
  else if (mechanism == "MAR1") {
    # Logistic regression to determinate the missingness
    # The options in this produce_NA function could be added to the main function
    mar1 <- produce_NA(
      df,
      mechanism = "MAR",
      perc.missing = miss_perc,
      by.patterns = F,
      logit.model = "MID"
    )
    X.mar1 <- mar1$data.incomp
    R.mar1 <- data.frame(mar1$idx_newNA)
    real_miss_perc <- sum(R.mar1 * 1) / prod(dim(R.mar1 * 1))
    return(list(
      "X.incomp" = X.mar1,
      "R.mask" = R.mar1,
      "real_miss_perc" = real_miss_perc
    ))
  }
  else if (mechanism == "MAR2") {
    # Censoring algorithm. Everything depends on the quantile of one specified complete column.
    # For example, the Y2[i] will be removed if Y1[i]<q(30%) of Y1)
    # For the categorical variable cat, the quantile is taken on the levels(cat)
    idx_ctrl <- mar2.col.ctrl
    X.mar2 <- df
    for (coll in ls_col_name[1:num_col]) {
      if (coll == ls_col_name[idx_ctrl]) {
        next
      }
      X.mar2[, coll] <- delete_MAR_censoring(X.mar2, miss_perc,
        coll,
        cols_ctrl = ls_col_name[idx_ctrl]
      )[, coll]
    }
    R.mar2 <- data.frame(is.na(X.mar2))
    real_miss_perc <- sum(R.mar2 * 1) / prod(dim(R.mar2 * 1))
    return(list(
      "X.incomp" = X.mar2,
      "R.mask" = R.mar2,
      "real_miss_perc" = real_miss_perc
    ))
  }
  else if (mechanism == "MAR3") {
    # Monotone with censoring mechanism, each column depends on the quantile of the observed data of the column before
    X.mar3 <- df
    if (miss_perc * num_col >= (num_col - 1)) {
      stop("Error: MAR3 mechanism cannot work with this miss_perc")
    }
    perc <- monot_quantil(miss_perc = miss_perc, num_col = num_col)
    i <- 1
    while (i < num_col) {
      ls_row <- which(!is.na(X.mar3[, ls_col_name[i]]))
      # if(i != num_col-1){
      for (coll in ls_col_name[(i + 1):num_col]) {
        X.mar3[ls_row, coll] <- delete_MAR_censoring(X.mar3[ls_row, ], perc, coll,
          cols_ctrl = ls_col_name[i]
        )[ls_row, coll]
      }
      # }
      # else{ # the last column adjust the missing percentage to approach the target missing percentage
      #   R = data.frame(is.na(X.mar3))
      #   p_adjust = (prod(dim(R*1)) * miss_perc - sum(R*1))/length(ls_row)
      #   X.mar3[ls_row,ls_col_name[num_col]] = delete_MAR_censoring(X.mar3[ls_row,], p_adjust, ls_col_name[num_col],
      #                                              cols_ctrl = ls_col_name[i])[ls_row,ls_col_name[num_col]]
      # }

      i <- i + 1
    }
    R.mar3 <- data.frame(is.na(X.mar3))
    real_miss_perc <- sum(R.mar3 * 1) / prod(dim(R.mar3 * 1))
    return(list(
      "X.incomp" = X.mar3,
      "R.mask" = R.mar3,
      "real_miss_perc" = real_miss_perc
    ))
  }
  else if (mechanism == "MNAR1") {
    # logistic regression to determinate the missingness, with num_patt_mnar random patterns
    mnar1 <- produce_NA(
      df,
      mechanism = "MNAR",
      perc.missing = miss_perc,
      by.patterns = F,
      logit.model = "LEFT"
    )
    X.mnar1 <- mnar1$data.incomp
    R.mnar1 <- data.frame(mnar1$idx_newNA)
    real_miss_perc <- sum(R.mnar1 * 1) / prod(dim(R.mnar1 * 1))
    return(list(
      "X.incomp" = X.mnar1,
      "R.mask" = R.mnar1,
      "real_miss_perc" = real_miss_perc
    ))
  }
  else if (mechanism == "MNAR2") {
    # A missing value in "X", if the x-value is below the miss_perc % quantile of "the first column "X"
    X.mnar2 <- df
    for (coll in ls_col_name) {
      X.mnar2[, coll] <- delete_MNAR_censoring(X.mnar2, miss_perc, coll)[, coll]
    }
    R.mnar2 <- data.frame(is.na(X.mnar2))
    real_miss_perc <- sum(R.mnar2 * 1) / prod(dim(R.mnar2 * 1))
    return(list(
      "X.incomp" = X.mnar2,
      "R.mask" = R.mnar2,
      "real_miss_perc" = real_miss_perc
    ))
  }
}



# Generate a list of incomplete data with all the mechanisms that are in generate_miss function
generate_miss_ls <- function(df, miss_perc) {
  return(
    list(
      "mcar" = generate_miss(df, miss_perc, mechanism = "MCAR"),
      "mar1" = generate_miss(df, miss_perc, mechanism = "MAR1"),
      "mar2" = generate_miss(df, miss_perc, mechanism = "MAR2"),
      "mar3" = generate_miss(df, miss_perc, mechanism = "MAR3"),
      "mnar1" = generate_miss(df, miss_perc, mechanism = "MNAR1"),
      "mnar2" = generate_miss(df, miss_perc, mechanism = "MNAR2")
    )
  )
}


# Encoding functions
ordinal_encode <- function(df, idx_col_cat) {
  for (j in idx_col_cat) {
    df[, j] <- as.numeric(factor(df[, j], levels = levels(df[, j])))
  }
  return(df)
}

factor_encode <- function(df, idx_col_cat) {
  for (j in idx_col_cat) {
    df[, j] <- factor(df[, j])
  }
  return(df)
}

normalize_num <- function(df, idx_col_num) {
  df[, idx_col_num] <- data.frame(apply(df[, idx_col_num], 2, scale))
  return(df)
}
