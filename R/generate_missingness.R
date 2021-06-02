#' generate_miss: Generate missing values with different mechanisms
#' 
#' @description
#' \code{generate_miss} function generates missing values in a complete 
#' dataframe with the chosen mechanism and the chosen proportion of missingness.
#' @param df Complete dataframe.
#' @param miss_perc Desired percentage of missing values. Note that all the mechanisms could only 
#' approach this proportion. The real proportion of missingness will be returned with the incomplete dataframe.
#' When the number of columns is small, some mechanisms could be incompatible with large missing percentage. 
#' @param mechanism Desired missing mechanism.
#' \itemize{
#'  \item \strong{MCAR} generates missing values by missing completely at random mechanism with Bernouilli distribution. 
#'  \item \strong{MAR1} generates missing values by missing at random mechanism with logistic regression on the observed data.
#'  \item \strong{MAR2} generates missing values by missing at random mechanism with censoring algorithm. 
#' The missingness in every other column depends on the quantile of one specified complete column \code{mar2.col.ctrl}.
#' For example, on row i, Y2[i] will be removed if Y1[i]<q(30\%) of Y1.
#'  \item \strong{MAR3} generates missing values by missing at random mechanism with monotone censoring mechanism.
#' The missingness in every column depends on the quantile of the observed data of the column before.
#' For example, on row i, Y2[i],Y3[i],...,Yn[i] will be removed if Y1[i]<q(30\%) of Y1.
#' And then on row s, Y3[s],...,Yn[s] will be removed if Y2[s]<q(30\%) of observed Y2 (Those who are not removed in step 1).
#'  \item \strong{MNAR1} generates missing values by missing not at random mechanism with logistic regression on the observed and missing part of data.
#'  \item \strong{MNAR2} generates missing values by missing not at random mechanism with censoring mechanism.
#'  For example, for each column j, on row i, Yj[i] will be removed if Yj[i]<q(30\%) of Y1.
#' }
#' @param mar2.col.ctrl Control column in mechanism MAR2
#' 
#' @return \code{X.incomp} Generated incomplete dataframe.
#' @return \code{R.mask} Mask of missingness. R[2,3]=1 means that df[2,3] is missing.
#' @return \code{real_miss_perc} Real proportion of missingness.
#' @export
#' @references
#' \itemize{
#' \item \url{https://rmisstastic.netlify.app/how-to/generate/misssimul}
#' \item \url{https://cran.r-project.org/web/packages/missMethods/vignettes/Generating-missing-values.html}
#' \item Santos, M. S., R. C. Pereira, A. F. Costa, J. P. Soares, J. Santos, and P. H. Abreu. 2019. Generating Synthetic Missing Data: A Review by Missing Mechanism. IEEE Access 7: 11651–67. \url{https://doi.org/10.1109/ACCESS.2019.2891360}.
#' }
#' @examples
#' n = 10000
#' mu.X = c(1, 2, 3)
#' Sigma.X = matrix(c(9, 3, 2, 3, 4, 0, 2, 0, 1), nrow = 3)
#' X.complete.cont = MASS::mvrnorm(n, mu.X, Sigma.X)
#' rs = generate_miss(X.complete.cont, 0.5, mechanism = "MNAR2")
#' rs$X.incomp
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
    colnames(R.mcar) <- colnames(X.mcar)
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
    colnames(R.mar1) <- colnames(X.mar1)
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
      X.mar2[, coll] <- missMethods::delete_MAR_censoring(X.mar2, miss_perc,
        coll,
        cols_ctrl = ls_col_name[idx_ctrl]
      )[, coll]
    }
    R.mar2 <- data.frame(is.na(X.mar2))
    real_miss_perc <- sum(R.mar2 * 1) / prod(dim(R.mar2 * 1))
    colnames(R.mar2) <- colnames(X.mar2)
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
        X.mar3[ls_row, coll] <-  missMethods::delete_MAR_censoring(X.mar3[ls_row, ], perc, coll,
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
    colnames(R.mar3) <- colnames(X.mar3)
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
    colnames(R.mnar1) <- colnames(X.mnar1)
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
      X.mnar2[, coll] <-  missMethods::delete_MNAR_censoring(X.mnar2, miss_perc, coll)[, coll]
    }
    R.mnar2 <- data.frame(is.na(X.mnar2))
    real_miss_perc <- sum(R.mnar2 * 1) / prod(dim(R.mnar2 * 1))
    colnames(R.mnar2) <- colnames(X.mnar2)
    return(list(
      "X.incomp" = X.mnar2,
      "R.mask" = R.mnar2,
      "real_miss_perc" = real_miss_perc
    ))
  }
}


#' generate_miss_ls: Generate a list of incomplete dataframes with different missing mechanisms
#' 
#' @description
#' \code{generate_miss_ls} function generates a list of incomplete dataframes. Missing values are generated in the given complete 
#' dataframe with the all mechanisms in \code{generate_miss} function and the chosen proportion of missingness.
#' @param df Complete dataframe.
#' @param miss_perc Desired percentage of missing values. Note that all the mechanisms could only 
#' approach this proportion. The real proportion of missingness will be returned with the incomplete dataframe.
#' When the number of columns is small, some mechanisms could be incompatible with large missing percentage. 
#' 
#' @return \code{mcar} Incomplete dataframe object with MCAR mechanism. The detailed description could be 
#' found in \code{generate_miss} documentation. Each object contains three attributes: 
#' \code{X.incomp}, \code{R.mask}, and \code{real_miss_perc}
#' @return \code{mar1} Incomplete dataframe object with MAR1 mechanism.
#' @return \code{mar2} Incomplete dataframe object with MAR2 mechanism.
#' @return \code{mar3} Incomplete dataframe object with MAR3 mechanism.
#' @return \code{mnar1} Incomplete dataframe object with MNAR1 mechanism.
#' @return \code{mnar2} Incomplete dataframe object with MNAR2 mechanism.
#' @export
#' @references
#' @examples
#' n = 10000
#' mu.X = c(1, 2, 3)
#' Sigma.X = matrix(c(9, 3, 2, 3, 4, 0, 2, 0, 1), nrow = 3)
#' X.complete.cont = MASS::mvrnorm(n, mu.X, Sigma.X)
#' rs = generate_miss_ls(X.complete.cont, 0.4)
#' rs$mcar$X.incomp
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


#' monot_quantil
#' @description Calculate the missing proportion in MAR3, used in 'generate_miss'. 
#' Solve the function (1-x)^p + (1-m)*p*x -1 = 0 in (0, 1), where m = miss_perc, p = num_co
#' @param miss_perc Desired proportion of missingness
#' @param num_col Number of columns
#' @return The required percentage for the quantil in MAR3. Solution to the function (1-x)^p + (1-m)*p*x -1 = 0 in (0, 1)
monot_quantil <- function(miss_perc, num_col) {
  m <- miss_perc
  p <- num_col
  tmp_result <- c()
  tempt <- pracma::linspace(0.01, 1, n = 1000)
  i <- 1
  for (x in tempt) {
    tmp_result[i] <- abs((1 - x)^p + (1 - m) * p * x - 1)
    i <- i + 1
  }
  return(tempt[which.min(tmp_result)])
}

