#' result_mice: Multiple imputation with bayesian method
#'
#' @description \code{result_mice} is a function that retrieve the final imputed
#' dataset after running function \code{mice} in the 'mice' package.
#' As for categorical columns, both factor form and onehot probability vector
#' form are returned
#'
#' More details about the MICE implementation could be found in the
#' documentation of function \code{mice} from the 'mice' package.
#' @param res Result returned by \code{mice} function
#' @param col_cat Categorical columns index
#' @param impnum Number of multiple imputations
#' @export
#' @return \code{ximp} Final imputed data matrix, which is obtained with
#' Rubin's Rule.
#' @return \code{ximp.disj} Final imputed data matrix of same type as 'ximp' for
#' the numeric columns. For the categorical columns, the prediction of
#' probability for each category is shown in form of onehot probability vector.
#' @references Van Buuren, S., Groothuis-Oudshoorn, K. (2011). \code{mice}:
#' Multivariate Imputation by Chained Equations in \code{R}. \emph{Journal of
#' Statistical Software}, \bold{45}(3), 1-67.
#' \url{https://www.jstatsoft.org/v45/i03/}
#'
#' Van Buuren, S. (2018).
#' \href{https://stefvanbuuren.name/fimd/sec-FCS.html#sec:MICE}{\emph{Flexible
#' Imputation of Missing Data. Second Edition.}}
#' Chapman & Hall/CRC. Boca Raton, FL.
#'
#' Van Buuren, S., Brand, J.P.L., Groothuis-Oudshoorn C.G.M., Rubin, D.B. (2006)
#' Fully conditional specification in multivariate imputation.  \emph{Journal of
#' Statistical Computation and Simulation}, \bold{76}, 12, 1049--1064.
#'
#' Van Buuren, S. (2007) Multiple imputation of discrete and continuous data by
#' fully conditional specification.  \emph{Statistical Methods in Medical
#' Research}, \bold{16}, 3, 219--242.
#'
#' Van Buuren, S., Boshuizen, H.C., Knook, D.L. (1999) Multiple imputation of
#' missing blood pressure covariates in survival analysis.  \emph{Statistics in
#' Medicine}, \bold{18}, 681--694.
#'
#' Brand, J.P.L. (1999) \emph{Development, implementation and evaluation of
#' multiple imputation strategies for the statistical analysis of incomplete
#' data sets.} Dissertation. Rotterdam: Erasmus University.
#'
#' Statistical Analysis with Missing Data, by Little and Rubin, 2002
result_mice <- function(res, impnum, col_cat = c()) {
  is_abind_package_installed()
  exist_cat <- !all(c(0, col_cat) == c(0))
  if (exist_cat) {
    dict_cat <- dict_onehot(res$data, col_cat)
  }
  # Extract the imputations with categorical columns in onehot form
  res.disj <- list()
  idx <- c(1:impnum)
  for (j in idx) {
    data <- mice::complete(res, j)
    dummy <- dummyVars(" ~ .", data = data, sep = "_")
    data.disj <- data.frame(predict(dummy, newdata = data))
    res.disj[[j]] <- data.disj
  }
  # Final result
  df_new_merge <- abind::abind(res.disj, along = 3)
  ximp.all <- data.frame(apply(df_new_merge, c(1, 2), mean))
  ximp.disj <- ximp.all
  if (exist_cat) {
    names_cat <- names(dict_cat)

    for (name in names_cat) {
      ximp.all[[name]] <- apply(
        ximp.all[dict_cat[[name]]], 1,
        which_max_cat, name, dict_cat
      )
      ximp.all[[name]] <- unlist(ximp.all[[name]])
      ximp.all[[name]] <- factor(ximp.all[[name]])
    }
    ximp <- ximp.all[colnames(res$data)]
  } else {
    ximp <- ximp.disj
  }
  return(list(ximp = ximp, ximp.disj = ximp.disj))
}
