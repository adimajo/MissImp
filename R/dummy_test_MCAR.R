#' dummy_test_matrix: Create the matrix of p-value for dummy t-chi-test
#'
#' @description
#' \code{dummy_test_matrix} generates a matrix of p-values for dummy t-chi-test. The null hypothesis(H0) is that the missing mechanism is MCAR.
#' The position [i,j] of this matrix shows the p-value of the test that the missigness in Yi does not depend on the value of Yj.
#'
#' We note Yj_1 as the part of Yj where Yi is missing, and Yj_0 as the part of Yj where Yi is observed.
#' Mj_1 and Mj_0 correspond to the mask of missingness where Yi is missing or observed. Mi is the mask of missingness for Yi.
#' For example, if Yi[3] is missing and Yj[3] is observed, then Mj_1[3]=0, Mi[3]=1.
#'
#' There are four situations:
#' \itemize{
#'  \item Yj is completely missing. In this case, no test will be done.
#'  \item Yj is partially observed, but Yj_1 (or Yj_0) is completely missing(or with only one observed data).
#'  In this case, a t-test is performed to test if the mean of Mj_0 (or Mj_1) is 1.
#'  \item Yj is numerical, Yj_1 and Yj_0 are both partially observed.
#'  In this case, a paired t-test is performed to test if Yj_1 and Yj_0 have the same mean.
#'  \item Yj is categorical, Yj_1 and Yj_0 are both partially observed.
#'  In this case, a chi-squared test is performed to test if Yj and Mi are independent.
#' }
#' @param df An incomplete dataframe.
#' @param col_cat The categorical columns index.
#' @return A matrix of p-value, where the position [i,j] shows the p-value of the test that the missigness in Yi does not depend on the value of Yj.
#' @export
#' @references
#' Missing value analysis & Data imputation, G. David Garson, 2015.
#' @examples
#' n <- 10000
#' mu.X <- c(1, 2, 3)
#' Sigma.X <- matrix(c(9, 3, 2, 3, 4, 0, 2, 0, 1), nrow = 3)
#' X.complete.cont <- MASS::mvrnorm(n, mu.X, Sigma.X)
#' rs <- generate_miss(X.complete.cont, 0.5, mechanism = "MNAR2")
#' dummy_test_matrix(rs$X.incomp, c())
#'
#' # dummy_test_matrix(airquality, col_cat = c(5:6))
dummy_test_matrix <- function(df, col_cat = c()) {
  test_result_dummy <- data.frame()
  df <- factor_encode(df, col_cat) # in case the categorical column is not a factor
  R_df <- data.frame(is.na(df)) * 1.0
  ls_col_name <- colnames(df)
  ls_row_name <- c()
  num_col <- length(ls_col_name)
  i <- 1
  while (i <= num_col) {
    col_ctr <- ls_col_name[i] # control column
    df_1 <- df[R_df[[col_ctr]] == 1, ]
    df_0 <- df[R_df[[col_ctr]] == 0, ]
    R_1 <- R_df[R_df[[col_ctr]] == 1, ] # this matrix is needed for the t-test in a special case
    R_0 <- R_df[R_df[[col_ctr]] == 0, ]
    row1 <- length(row.names(df_1))
    row0 <- length(row.names(df_1))
    ls_row_name[i] <- col_ctr
    if (row1 == 0 |
      row0 == 0) {
      # if one column contains only NA or doesn't contain any NA
      test_result_dummy[i, ] <- NA
      i <- i + 1
      next
    }
    j <- 1
    while (j <= num_col) {
      if (i != j) {
        col_test <- ls_col_name[j]
        # if the test column is all NA, the test is not done on this column
        if (all(is.na(df[[col_test]]))) {
          test_result_dummy[i, col_test] <- NA
          j <- j + 1
          next
        }

        # if when the control column is NA, the test column is all NA(or with only one value)
        # or when the control column is observed , the test column is all NA(or with only one value)
        # the paired t-test is done on the indicator matrix R_1 an R_0
        critc_1 <- all(is.na(df_1[[col_test]])) || (sum(!is.na(df_1[[col_test]])) == 1)
        critc_0 <- all(is.na(df_0[[col_test]])) || (sum(!is.na(df_0[[col_test]])) == 1)
        if (critc_1 ||
          critc_0) {
          # the situation of critic_1 && critic_0 is discussed before
          if (critc_1 && !critc_0) {
            test_result_dummy[i, col_test] <- stats::t.test(R_0[[col_test]] * 1, mu = 1)$p.value
          } else {
            test_result_dummy[i, col_test] <- stats::t.test(R_1[[col_test]] * 1, mu = 1)$p.value
          }

          j <- j + 1
          next
        }

        # if the test column is categorical, we use chi 2 test
        if (j %in% col_cat) {
          test_table <- data.frame(col_test = df[[col_test]], R = R_df[[col_ctr]])
          test_result_dummy[i, col_test] <- suppressWarnings(with(test_table, stats::chisq.test(col_test, R))$p.value)
          # warnings about the approximation of normal distribution if the sample is small
          j <- j + 1
          next
        }
        test_result_dummy[i, col_test] <- stats::t.test(df_1[[col_test]], df_0[[col_test]])$p.value
      }
      j <- j + 1
    }
    i <- i + 1
  }
  row.names(test_result_dummy) <- ls_row_name
  return(test_result_dummy)
}





#' dummy_test: dummy t-chi-test for MCAR
#'
#' @description
#' \code{dummy_test} combines the p-values from \code{dummy_test_matrix} by Fisher's method.
#'
#' \code{dummy_test_matrix} generates a matrix of p-values for dummy t-chi-test. The null hypothesis(H0) is that the missing mechanism is MCAR.
#' The position [i,j] of this matrix shows the p-value of the test that the missigness in Yi does not depend on the value of Yj.
#'
#' We note Yj_1 as the part of Yj where Yi is missing, and Yj_0 as the part of Yj where Yi is observed.
#' Mj_1 and Mj_0 correspond to the mask of missingness where Yi is missing or observed. Mi is the mask of missingness for Yi.
#' For example, if Yi[3] is missing and Yj[3] is observed, then Mj_1[3]=0, Mi[3]=1.
#'
#' There are four situations:
#' \itemize{
#'  \item Yj is completely missing. In this case, no test will be done.
#'  \item Yj is partially observed, but Yj_1 (or Yj_0) is completely missing.
#'  In this case, a t-test is performed to test if the mean of Mj_0 (or Mj_1) is 1.
#'  \item Yj is numerical, Yj_1 and Yj_0 are both partially observed.
#'  In this case, a paired t-test is performed to test if Yj_1 and Yj_0 have the same mean.
#'  \item Yj is categorical, Yj_1 and Yj_0 are both partially observed.
#'  In this case, a chi-squared test is performed to test if Yj and Mi are independent.
#' }
#' @param df An incomplete dataframe.
#' @param col_cat The categorical columns index.
#' @return \code{p.matrix} A matrix of p-value, where the position [i,j] shows the p-value of the test that the missigness in Yi does not depend on the value of Yj.
#' @return \code{dof} Degree of freedom for the chi-squared statistics in Fisher's method.
#' @return \code{chi2stat} Chi-squared statistics by Fisher's method.
#' @return \code{p.value} Combined p-value for the MCAR test.
#' @export
#' @references
#' Missing value analysis & Data imputation, G. David Garson, 2015
#' @examples
#' n <- 10000
#' mu.X <- c(1, 2, 3)
#' Sigma.X <- matrix(c(9, 3, 2, 3, 4, 0, 2, 0, 1), nrow = 3)
#' X.complete.cont <- MASS::mvrnorm(n, mu.X, Sigma.X)
#' rs <- generate_miss(X.complete.cont, 0.5, mechanism = "MNAR2")
#' dummy_test(rs$X.incomp, c())
#'
#' # dummy_test(airquality, col_cat = c(5:6))
dummy_test <- function(df, col_cat = c()) {
  p_matrix <- dummy_test_matrix(df, col_cat)
  p_vector <- as.vector(p_matrix)
  p_vector <- p_vector[!is.na(p_vector)]
  dof <- 2 * length(p_vector)
  res <- -2 * sum(log(p_vector))
  p_dum <- stats::pchisq(res, df = dof, lower.tail = FALSE)
  return(list(
    "p.matrix" = p_matrix,
    "dof" = dof,
    "chi2stat" = res,
    "p.value" = p_dum
  ))
}
