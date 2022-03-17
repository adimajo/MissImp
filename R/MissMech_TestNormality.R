#' TestMCARNormality: Testing Homoscedasticity, Multivariate Normality, and Missing Completely at Random
#'
#' @description
#' This is a function from \code{MissMech} package. The description of this function in the original package is as following:
#' The main purpose of this function is to test whether the missing data mechanism, for an incompletely observed data set,
#' is one of missing completely at random (MCAR). As a by product, however, this package has the capabilities of imputing incomplete data,
#' performing a test to determine whether data have a multivariate normal distribution, performing a test of equality of covariances for groups,
#' and obtaining normal-theory maximum likelihood estimates for mean and covariance when data are incomplete.
#' The test of MCAR follows the methodology proposed by Jamshidian and Jalal (2010).
#' It is based on testing equality of covariances between groups having identical missing data patterns.
#' The data are imputed, using two options of normality and distribution free, and the test of equality of covariances between
#' groups with identical missing data patterns is performed also with options of assuming normality (Hawkins test) or non-parametrically.
#' Users can optionally use their own method of data imputation as well. Multiple imputation is an additional feature of the program that can be used
#' as a diagnostic tool to help identify cases or variables that contribute to rejection of MCAR, when the MCAR test is rejecetd
#' (See Jamshidian and Jalal, 2010 for details). As explained in Jamshidian, Jalal, and Jansen (2014), this package can also be used for imputing missing data,
#'  test of multivariate normality, and test of equality of covariances between several groups when data are completly observed.
#'
#' @param data A matrix or data frame consisting of at least two columns. Values must be numerical with missing data indicated by NA.
#' @param del.lesscases Missing data patterns consisting of del.lesscases number of cases or less will be removed from the data set.
#' @param imputation.number Number of imputations to be used, if data are to be multiply imputed.
#' @param method method is an option that allows the user to select one of the methods of Hawkins or nonparametric for the test.
#' If the user is certain that data have multivariate normal distribution, the method="Hawkins" should be selected.
#' On the other hand if data are not normally distributed, then method="Nonparametric" should be used. If the user is unsure,
#' then the default value of method="Auto" will be used, in which case both the Hawkins and the nonparametric tests will be run,
#' and the default output follows the recommendation by Jamshidian and Jalal (2010) outlined in their flowchart given in Figure 7 of their paper.
#' @param imputation.method "Dist.Free": Missing data are imputed nonparametrically using the method of Sirvastava and Dolatabadi (2009);
#' also see Jamshidian and Jalal (2010).
#'
#' "normal": Missing data are imputed assuming that the data come from a multivariate normal distribution.
#' The maximum likelihood estimate of the mean and covariance obtained from Mls is used for generating imputed values.
#' The imputed values are based on the conditional distribution of the missing variables given the observed variables;
#' see Jamshidian and Jalal (2010) for more details.
#' @param nrep Number of replications used to simulate the Neyman distribution to determine the cut off value for the Neyman test in the program SimNey.
#'  Larger values increase the accuracy of the Neyman test.
#' @param n.min The minimum number of cases in a group that triggers the use of asymptotic Chi distribution in place of the emprical distribution in the Neyman test of uniformity.
#' @param seed An initial random number generator seed. The default is 110 that can be reset to a user selected number. If the value is set to NA, a system selected seed is used.
#' @param alpha The significance level at which tests are performed.
#' @param imputed.data The user can optionally provide an imputed data set. In this case the program will not impute the data and will use the imputed data set for the tests performed.
#' Note that the order of cases in the imputed data set should be the same as that of the incomplete data set.
#'
#' @return \code{analyzed.data} The data that were used in the analysis. If del.lesscases=0, this is the same as the orginal data inputted. If del.lesscases > 0, then this is the data with cases removed.
#' @return \code{imputed.data} The analyzed.data after imputation. If imputation.number > 1, the first imputed data set is returned.
#' @return \code{ordered.data} The analyzed.data ordered according to missing data pattern, usin the function OrderMissing.
#' @return \code{caseorder} A mapping of case number indices from ordered.data to the original data.
#' More specifically, the j-th row of the ordered.data is the caseorder[j]-th (the j-th element of caseorder) row of the original data.
#' @return \code{pnormality} p-value for the nonparametric test: When imputation.number > 1, this is a vector with each element corresponding to each of the imputed data sets.
#' @return \code{adistar} A matrix consisting of the Anderson-Darling test statistic for each group (columns) and each imputation (rows).
#' @return \code{adstar} Sum of adistar: When imputation.number >1, this is a vector with each element corresponding to each of the imputed data sets.
#' @return \code{pvalcomb} p-value for the Hawkins test: When imputation.number >1, this is a vector with each element corresponding to each of the imputed data sets.
#' @return \code{pvalsn} A matrix consisting of Hawkins test statistics for each group (columns) and each imputation (rows).
#' @return \code{g} Number of patterns used in the analysis.
#' @return \code{combp} Hawkins test statistic: When imputation.number > 1, this is a vector with each element corresponding to each of the imputed data sets.
#' @return \code{alpha} The significance level at which the hypothesis tests are performed.
#' @return \code{patcnt} A vector consisting the number of cases corresponding to each pattern in patused.
#' @return \code{patused} A matrix indicating the missing data patterns in the data set, using 1 and NA's.
#' @return \code{imputation.number} A value greater than or equal to 1. If a value larger than 1 is used, data will be imputed imputation.number times.
#' @return \code{mu} The normal-theory maximum likelihood estimate of the variables means.
#' @return \code{sigma} The normal-theory maximum likelihood estimate of the variables covariance matrix.
#'
#' @export
#'
#' @references
#' Jamshidian, M. and Bentler, P. M. (1999). ``ML estimation of mean and covariance structures with missing data using complete data routines.'' Journal of Educational and Behavioral Statistics, 24, 21-41.
#'
#' Jamshidian, M. and Jalal, S. (2010). ``Tests of homoscedasticity, normality, and missing at random for incomplete multivariate data,'' Psychometrika, 75, 649-674.
#'
#' Jamshidian, M. Jalal, S., and Jansen, C. (2014). `` MissMech: An R Package for Testing Homoscedasticity, Multivariate Normality, and Missing Completely at Random (MCAR),'' Journal of Statistical Software, 56(6), 1-31.
TestMCARNormality <- function(data, del.lesscases = 6, imputation.number = 1, method = "Auto",
                              imputation.method = "Dist.Free", nrep = 10000, n.min = 30,
                              seed = 110, alpha = 0.05, imputed.data = NA) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if (!is.na(imputed.data[1]) && imputation.number != 1) {
    stop("Warning: No multiple imputation allowed when imputed data is provided.\n")
  }
  if (!is.matrix(data)) {
    stop("Warning: Data is not a matrix or data frame.\n")
  }
  if (length(data) == 0) {
    stop("Warning: Data is empty.\n")
  }
  if (ncol(data) < 2) {
    stop("Warning: More than 1 variable is required.\n")
  }
  allempty <- which(apply(!is.na(data), 1, sum) == 0)
  if (length(allempty) != 0) {
    data <- data[apply(!is.na(data), 1, sum) != 0, ]
    cat("Warning:", length(allempty), "Cases with all variables missing have been removed \n
        from the data.\n")
  }
  newdata <- OrderMissing(data, del.lesscases)
  if (length(newdata$data) == 0) {
    stop("Warning: There are no data sets after deleting insufficient cases.\n")
  }

  if (newdata$g == 1) {
    stop("Warning: More than one missing data pattern should be present.\n")
  }
  if (sum(newdata$patcnt == 1) > 0) {
    stop("Warning: At least 2 cases needed in each missing data patterns.\n")
  }

  y <- newdata$data
  patused <- newdata$patused
  patcnt <- newdata$patcnt
  spatcnt <- newdata$spatcnt
  caseorder <- newdata$caseorder
  removedcases <- newdata$removedcases
  n <- nrow(y)
  p <- ncol(y)
  g <- newdata$g
  spatcntz <- c(0, spatcnt)
  pvalsn <- matrix(0, imputation.number, g)
  adistar <- matrix(0, imputation.number, g)
  pnormality <- c()
  x <- vector("list", g)
  n4sim <- vector("list", g)
  #------------------------------imputation-----------------------
  mu <- matrix(0, p, 1)
  sig <- diag(1, p)

  emest <- Mls(newdata, mu, sig, 1e-6)
  mu <- emest$mu
  sig <- emest$sig
  if (is.na(imputed.data[1])) {
    yimp <- y
    if (imputation.method == "Dist.Free") {
      iscomp <- apply(patused, 1, sum, na.rm = TRUE) == p

      cind <- which(iscomp)
      ncomp <- patcnt[cind]
      if (length(ncomp) == 0) ncomp <- 0
      use.normal <- FALSE
      if (ncomp >= 10 && ncomp >= 2 * p) {
        compy <- y[seq(spatcntz[cind] + 1, spatcntz[cind + 1]), ]
        ybar <- matrix(apply(compy, 2, mean))
        sbar <- stats::cov(compy)
        resid <- (ncomp / (ncomp - 1))^.5 *
          (compy - matrix(ybar, ncomp, p, byrow = TRUE))
      } else {
        cat("Warning: There is not sufficient number of complete cases.\n  Dist.Free imputation requires a least 10 complete cases\n  or 2*number of variables, whichever is bigger.\n  imputation.method = normal will be used instead.\n")
        use.normal <- TRUE
      }
    }
    for (k in 1:imputation.number)
    {
      #-----------------normal imputation--------------------------------
      if (imputation.method == "Normal" || use.normal) {
        yimp <- Impute(data = y, mu, sig, imputation.method = "Normal")
        yimp <- yimp$yimpOrdered
      }
      #-----------------distribution free imputation---------------------------------
      if (imputation.method == "Dist.Free" && !use.normal) {
        yimp <- Impute(data = y, ybar, sbar, imputation.method = "Dist.Free", resid)
        yimp <- yimp$yimpOrdered
      }
      if (k == 1) yimptemp <- yimp
      #--------------Hawkin's test on the completed data------------------
      templist <- Hawkins(yimp, spatcnt)
      fij <- templist$fij
      tail <- templist$a
      ni <- templist$ni
      if (method == "Auto" || method == "Hawkins") {
        # Neyman test of uniformity for each group
        for (i in 1:g)
        {
          if (ni[i] < n.min) {
            if (k == 1) {
              n4sim[[i]] <- SimNey(ni[i], nrep)
            }
          }
          templist <- TestUNey(tail[[i]], nrep, sim = n4sim[[i]], n.min)
          pn <- templist$pn
          n4 <- templist$n4
          pn <- pn + (pn == 0) / nrep
          pvalsn[k, i] <- pn
        }
      }
      #--------------Anderson darling test for equality of distribution
      if (method == "Auto" || method == "Nonparametric") {
        if (length(ni) < 2) {
          stop("Warning: Not enough groups for AndersonDarling test.")
        }
        templist <- AndersonDarling(fij, ni)
        p.ad <- templist$pn
        adistar[k, ] <- templist$adk.all
        pnormality <- c(pnormality, p.ad)
      }
    }
  } else {
    yimp <- imputed.data[caseorder, ]
    yimptemp <- yimp
    templist <- Hawkins(yimp, spatcnt)
    fij <- templist$fij
    tail <- templist$a
    ni <- templist$ni
    if (method == "Auto" || method == "Hawkins") {
      # Neyman test of uniformity for each group
      for (i in 1:g)
      {
        if (ni[i] < n.min) {
          n4sim[[i]] <- SimNey(ni[i], nrep)
        }
        templist <- TestUNey(tail[[i]], nrep, sim = n4sim[[i]], n.min)
        pn <- templist$pn
        n4 <- templist$n4
        pn <- pn + (pn == 0) / nrep
        pvalsn[1, i] <- pn
      }
    }
    if (method == "Auto" || method == "Nonparametric") {
      #--------------Anderson darling test for equality of distribution
      templist <- AndersonDarling(fij, ni)
      p.ad <- templist$pn
      adistar[1, ] <- templist$adk.all
      pnormality <- c(pnormality, p.ad)
    }
  }
  adstar <- apply(adistar, 1, sum)
  # combine p-values of test of uniformity
  combp <- -2 * apply(log(pvalsn), 1, sum)
  pvalcomb <- stats::pchisq(combp, 2 * g, lower.tail = FALSE)
  if (method == "Hawkins") {
    pnormality <- NULL
    adstar <- NULL
    adistar <- NULL
  }
  if (method == "Nonparametric") {
    pvalcomb <- NULL
    combp <- NULL
    pvalsn <- NULL
  }
  yimptemp <- yimptemp[order(caseorder), ]
  if (length(removedcases) == 0) {
    dataused <- data
  } else {
    dataused <- data[-1 * removedcases, ]
  }
  homoscedastic <- list(
    analyzed.data = dataused, imputed.data = yimptemp,
    ordered.data = y, caseorder = caseorder,
    pnormality = pnormality, adstar = adstar, adistar = adistar,
    pvalcomb = pvalcomb, combp = combp, pvalsn = pvalsn, g = g, alpha = alpha,
    patused = patused, patcnt = patcnt, imputation.number = imputation.number, mu = mu, sigma = sig
  )
  homoscedastic$call <- match.call()
  class(homoscedastic) <- "testhomosc"
  homoscedastic
}
#---------------------------------------------------------------------
# testmcar <- function(x, ...) UseImputationMethod("testmcar")
# testmcar.default <- function(data, ncases = 6, imputation.number = 10,
#                             imputation.method = "Normal", nrep = 10000)
# {
# test <- TestMCARNormality(data, ncases = 6, imputation.number = 10,
#                         imputation.method = "Normal", nrep = 10000)
# test$call <- match.call()
# class(test) <- "testmcar"
# test
# }
#---------------------------------------------------------------------
# printing format for the class "testhomosc"
print.testhomosc <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  # cat("\nNumber of imputation:\n")
  # print(x$imputation.number)
  ni <- x$patcnt
  cat("\nNumber of Patterns: ", x$g, "\n\nTotal number of cases used in the analysis: ", sum(ni), "\n")
  cat("\n Pattern(s) used:\n")
  alpha <- x$alpha
  disp.patt <- cbind(x$patused, ni)
  colnames(disp.patt)[ncol(disp.patt)] <- "Number of cases"
  rownames(disp.patt) <- rownames(disp.patt, do.NULL = FALSE, prefix = "group.")
  print(disp.patt, print.gap = 3)
  method <- "Auto"
  if (is.null(x$pnormality)) method <- "Hawkins"
  if (is.null(x$pvalcomb)) method <- "Nonparametric"
  cat("\n\n    Test of normality and Homoscedasticity:\n  -------------------------------------------\n")
  if (method == "Auto") {
    cat("\nHawkins Test:\n")
    cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", x$pvalcomb[1], "\n")
    if (x$pvalcomb[1] > alpha) {
      cat("\n    There is not sufficient evidence to reject normality
          or MCAR at", alpha, "significance level\n")
    } else {
      cat("\n    Either the test of multivariate normality or homoscedasticity (or both) is rejected.\n    Provided that normality can be assumed, the hypothesis of MCAR is
          rejected at", alpha, "significance level. \n")
      cat("\nNon-Parametric Test:\n")
      cat("\n    P-value for the non-parametric test of homoscedasticity: ", x$pnormality[1], "\n")
      if (x$pnormality[1] > alpha) {
        cat("\n    Reject Normality at", alpha, "significance level.
            There is not sufficient evidence to reject MCAR at", alpha, "significance level.\n")
      } else {
        cat("\n    Hypothesis of MCAR is rejected at ", alpha, "significance level.
            The multivariate normality test is inconclusive. \n")
      }
    }
  }
  if (method == "Hawkins") {
    cat("\nHawkins Test:\n")
    cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", x$pvalcomb[1], "\n")
  }
  if (method == "Nonparametric") {
    cat("\nNon-Parametric Test:\n")
    cat("\n    P-value for the non-parametric test of homoscedasticity: ", x$pnormality[1], "\n")
  }
}
#----------------------------------------------------------------------------
summary.testhomosc <- function(object, ...) {
  ni <- object$patcnt
  cat("\nNumber of imputation: ", object$imputation.number, "\n")
  cat("\nNumber of Patterns: ", object$g, "\n\nTotal number of cases used in the analysis: ", sum(ni), "\n")
  cat("\n Pattern(s) used:\n")
  alpha <- object$alpha
  disp.patt <- cbind(object$patused, ni)
  colnames(disp.patt)[ncol(disp.patt)] <- "Number of cases"
  rownames(disp.patt) <- rownames(disp.patt, do.NULL = FALSE, prefix = "group.")
  print(disp.patt, print.gap = 3)
  method <- "Auto"
  if (is.null(object$pnormality)) method <- "Hawkins"
  if (is.null(object$pvalcomb)) method <- "Nonparametric"
  cat("\n\n    Test of normality and Homoscedasticity:\n  -------------------------------------------\n")
  if (method == "Auto") {
    cat("\nHawkins Test:\n")
    cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", object$pvalcomb[1], "\n")
    cat("\nNon-Parametric Test:\n")
    cat("\n    P-value for the non-parametric test of homoscedasticity: ", object$pnormality[1], "\n")
  }
  if (method == "Hawkins") {
    cat("\nHawkins Test:\n")
    cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", object$pvalcomb[1], "\n")
  }
  if (method == "Nonparametric") {
    cat("\nNon-Parametric Test:\n")
    cat("\n    P-value for the non-parametric test of homoscedasticity: ", object$pnormality[1], "\n")
  }
}
#-----------------------------------------------------------------------------
# Plot "testhomosc"
boxplot.testhomosc <- function(x, ...) {
  if (is.null(x$pnormality)) {
    graphics::par(bg = "cornsilk")
    graphics::boxplot(x$pvalsn, col = "lightcyan", border = "blue", medlwd = .5, medcol = "red")
    graphics::title(
      main = "Boxplots of p-values corresponding to each set of the missing data patterns\n for the Neyman test of Uniformity",
      xlab = "Missing data pattern group", ylab = "P-value", font.main = 4,
      col.main = "blue4", cex.main = 1, font.lab = 4, cex.lab = 0.8,
      col.lab = "blue4"
    )
    graphics::abline(h = x$alpha / x$g, col = "red", lty = 2)
  }
  if (is.null(x$pvalcomb)) {
    graphics::par(bg = "cornsilk")
    graphics::boxplot(x$adistar, col = "lightcyan", border = "blue", medlwd = .5, medcol = "red")
    graphics::title(
      main = "Boxplots of the T-value test statistics corresponding to each set of missing\n data patterns for the non-parametric test",
      xlab = "Missing data pattern group", ylab = expression(T[i]),
      font.main = 4, col.main = "blue4", cex.main = 1, font.lab = 4,
      cex.lab = 0.8, col.lab = "blue4"
    )
  }
  if (!is.null(x$pvalcomb) && !is.null(x$pnormality)) {
    graphics::par(mfrow = c(2, 1), bg = "cornsilk")
    graphics::boxplot(x$pvalsn, col = "lightcyan", border = "blue", medlwd = .5, medcol = "red")
    graphics::title(
      main = "Boxplots of p-values corresponding to each set of the missing data patterns\n for the Neyman test of Uniformity",
      xlab = "Missing data pattern group", ylab = "P-value", font.main = 4,
      col.main = "blue4", cex.main = 1, font.lab = 4, cex.lab = 0.8,
      col.lab = "blue4"
    )
    graphics::abline(h = x$alpha / x$g, col = "red", lty = 2)
    graphics::boxplot(x$adistar, col = "lightcyan", border = "blue", medlwd = .5, medcol = "red")
    graphics::title(
      main = "Boxplots of the T-value test statistics corresponding to each set of missing\n data patterns for the non-parametric test",
      xlab = "Missing data pattern group", ylab = expression(T[i]),
      font.main = 4, col.main = "blue4", cex.main = 1, font.lab = 4,
      cex.lab = 0.8, col.lab = "blue4"
    )
  }
}

TestUNey <- function(x, nrep = 10000, sim = NA, n.min = 30) {
  # This routine tests whether the values in each row of x are unif(0,1). It
  # uses the Neyman's smooth test (see e.g., Ledwina 1994, TAS)
  # x is a vector
  # P-values are computed based on a
  # resampling method from unif(0,1).
  # All values of $x$ are between 0 and 1
  n <- length(x)
  pi <- LegNorm(x)
  n4 <- (apply(pi$p1, 2, sum)^2 + apply(pi$p2, 2, sum)^2 +
    apply(pi$p3, 2, sum)^2 + apply(pi$p4, 2, sum)^2) / n
  if (n < n.min) {
    if (is.na(sim[1])) {
      sim <- SimNey(n, nrep)
    }
    pn <- length(which(sim > n4)) / nrep
  } else {
    pn <- stats::pchisq(n4, 4, lower.tail = FALSE)
  }
  list(pn = pn, n4 = n4)
}
SimNey <- function(n, nrep) {
  x <- matrix(stats::runif(nrep * n), ncol = nrep)
  pi <- LegNorm(x)
  n4sim <- (apply(pi$p1, 2, sum)^2 + apply(pi$p2, 2, sum)^2 +
    apply(pi$p3, 2, sum)^2 + apply(pi$p4, 2, sum)^2) / n
  n4sim <- sort(n4sim)
}

AndersonDarling <- function(data, number.cases) {
  # Not adjusted for ties
  # x is data vector
  # ni is number of cases in each group
  x <- data
  ni <- number.cases
  if (length(ni) < 2) {
    cat("Warning: Not enough groups for AndersonDarling test.")
    stop("")
  }

  k <- length(ni)
  ni.z <- c(0, cumsum(ni))
  n <- length(x)
  x.sort <- sort(x)
  x.sort <- x.sort[1:(n - 1)]
  ind <- which(duplicated(x.sort) == 0)
  counts <- c(ind, length(x.sort) + 1) - c(0, ind)
  hj <- counts[2:(length(ind) + 1)]
  hn <- cumsum(hj)
  zj <- x.sort[ind]
  adk <- 0
  adk.all <- matrix(0, k, 1) # to keep contribution of kth group
  for (i in 1:k)
  {
    ind <- (ni.z[i] + 1):ni.z[i + 1]
    templist <- expand.grid(zj, x[ind])
    b <- templist[, 1] == templist[, 2]
    fij <- apply(matrix(b, length(zj)), 1, sum)
    mij <- cumsum(fij)
    num <- (n * mij - ni[i] * hn)^2
    dem <- hn * (n - hn)
    adk.all[i] <- (1 / ni[i] * sum(hj * (num / dem)))
    adk <- adk + adk.all[i]
  }
  adk <- (1 / n) * adk
  adk.all <- adk.all / n
  # Exact sample variance of the k-sample Anderson-Darling
  # Finding Variance of the statistics
  j <- sum(1 / ni)
  i <- seq(1:(n - 1))
  h <- sum(1 / i)
  g <- 0
  for (i in 1:(n - 2)) {
    g <- g + (1 / (n - i)) * sum(1 / seq((i + 1), (n - 1)))
  }
  a <- (4 * g - 6) * (k - 1) + (10 - 6 * g) * j
  b <- (2 * g - 4) * k^2 + 8 * h * k +
    (2 * g - 14 * h - 4) * j - 8 * h + 4 * g - 6
  c <- (6 * h + 2 * g - 2) * k^2 + (4 * h - 4 * g + 6) * k +
    (2 * h - 6) * j + 4 * h
  d <- (2 * h + 6) * k^2 - 4 * h * k
  var.adk <- ((a * n^3) + (b * n^2) + (c * n) + d) /
    ((n - 1) * (n - 2) * (n - 3))
  if (var.adk < 0) var.adk <- 0
  adk.s <- (adk - (k - 1)) / sqrt(var.adk)
  # k-sample Anderson-Darling P-value calculation by an extrapolate-interpolate
  # procedure
  a0 <- c(0.25, 0.10, 0.05, 0.025, 0.01)
  b0 <- c(0.675, 1.281, 1.645, 1.96, 2.326)
  b1 <- c(-0.245, 0.25, 0.678, 1.149, 1.822)
  b2 <- c(-0.105, -0.305, -0.362, -0.391, -0.396)
  c0 <- log((1 - a0) / a0)
  qnt <- b0 + b1 / sqrt(k - 1) + b2 / (k - 1)
  if (adk.s <= qnt[3]) {
    ind <- seq(1:4)
  } else {
    ind <- seq(2:5)
  }
  yy <- stats::spline(qnt[ind], c0[ind], xout = adk.s)$y
  p <- 1 / (1 + exp(yy))
  list(pn = p, adk.all = adk.all, adk = adk, var.sdk = var.adk)
} # end function

OrderMissing <- function(data, del.lesscases = 0) {
  # case order has the order of the original data
  y <- data
  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }
  if (!is.matrix(y)) {
    stop("Warning data is not a matrix or data frame")
  }
  if (length(y) == 0) {
    stop("Warning: data is empty")
  }
  names <- colnames(y)
  n <- nrow(y)
  pp <- ncol(y)
  yfinal <- c()
  patused <- c()
  patcnt <- c()
  caseorder <- c()
  removedcases <- c()
  ordertemp <- c(1:n)
  ntemp <- n
  ptemp <- pp
  done <- FALSE
  yatone <- FALSE
  while (!done) {
    pattemp <- is.na(y[1, ])
    indin <- c()
    indout <- c()
    done <- TRUE
    for (i in 1:ntemp)
    {
      if (all(is.na(y[i, ]) == pattemp)) {
        indout <- c(indout, i)
      } else {
        indin <- c(indin, i)
        done <- FALSE
      }
    }
    if (length(indin) == 1) yatone <- TRUE
    yfinal <- rbind(yfinal, y[indout, ])
    y <- y[indin, ]
    caseorder <- c(caseorder, ordertemp[indout])
    ordertemp <- ordertemp[indin]
    patcnt <- c(patcnt, length(indout))
    patused <- rbind(patused, pattemp)
    if (yatone) {
      pattemp <- is.na(y)
      yfinal <- rbind(yfinal, matrix(y, ncol = pp))
      y <- c()
      indin <- c()
      indout <- c(1)
      caseorder <- c(caseorder, ordertemp[indout])
      ordertemp <- ordertemp[indin]
      patcnt <- c(patcnt, length(indout))
      patused <- rbind(patused, pattemp)
      done <- TRUE
    }

    if (!done) ntemp <- nrow(y)
  }
  # yfinal <- rbind(yfinal, y)
  caseorder <- c(caseorder, ordertemp)
  patused <- ifelse(patused, NA, 1)
  rownames(patused) <- NULL
  colnames(patused) <- names
  spatcnt <- cumsum(patcnt)
  dataorder <- list(
    data = yfinal, patused = patused, patcnt = patcnt,
    spatcnt = spatcnt, g = length(patcnt),
    caseorder = caseorder, removedcases = removedcases
  )
  dataorder$call <- match.call()
  class(dataorder) <- "orderpattern"
  if (del.lesscases > 0) {
    dataorder <- DelLessData(dataorder, del.lesscases)
  }
  dataorder$patused <- matrix(dataorder$patused, ncol = pp)
  colnames(dataorder$patused) <- names
  dataorder
}
# Order <- function(x, ...) UseMethod("Order")
# Order.default <- function(x, ...) {
# temp <- OrderMissing(x)
# temp$call <- match.call()
# class(temp) <- "ordered"
# temp
# }
print.orderpattern <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nNumber of Patterns: ", x$g, "\n")
  cat("\nPattern used:\n")
  ni <- x$patcnt
  disp.patt <- cbind(x$patused, ni)
  colnames(disp.patt)[ncol(disp.patt)] <- "Number of cases"
  rownames(disp.patt) <- rownames(disp.patt, do.NULL = FALSE, prefix = "group.")
  print(disp.patt, print.gap = 3)
}

summary.orderpattern <- function(object, ...) {
  summary(object$data)
}

DelLessData <- function(data, ncases = 0) {
  # This function deletes cases of a missing pattern with less than or equal to ncases
  if (length(data) == 0) {
    stop("Warning: data is empty")
  }
  if (is.matrix(data)) {
    data <- OrderMissing(data)
  }
  n <- nrow(data$data)
  p <- ncol(data$data)
  ind <- which(data$patcnt <= ncases)
  spatcntz <- c(0, data$spatcnt)
  rm <- c()
  removedcases <- c()
  if (length(ind) != 0) {
    # cat("cases with insufficient number of observations were removed")
    for (i in 1:length(ind))
    {
      rm <- c(rm, seq(spatcntz[ind[i]] + 1, spatcntz[ind[i] + 1]))
    }
    y <- data$data[-1 * rm, ]
    removedcases <- data$caseorder[rm]
    patused <- data$patused[-1 * ind, ]
    patcnt <- data$patcnt[-1 * ind]
    caseorder <- data$caseorder[-1 * rm]
    spatcnt <- cumsum(patcnt)
  } else {
    patused <- data$patused
    patcnt <- data$patcnt
    spatcnt <- data$spatcnt
    caseorder <- data$caseorder
    y <- data$data
  }

  newdata <- list(
    data = y, patused = patused, patcnt = patcnt,
    spatcnt = spatcnt, g = length(patcnt), caseorder = caseorder,
    removedcases = removedcases
  )
  #  colnames(newdata)<-colnames(data)
  class(newdata) <- "orderpattern"
  newdata
}

Mls <- function(data, mu = NA, sig = NA, tol = 1e-6, Hessian = FALSE) {
  # mu is estimate of the mean
  # sig is estimate of the covariance
  if (!is.matrix(data) && class(data) != "orderpattern") {
    stop("Warning: data must have the classes of matrix or orderpattern.\n")
  }
  if (is.matrix(data)) {
    allempty <- which(apply(!is.na(data), 1, sum) == 0)
    if (length(allempty) != 0) {
      data <- data[apply(!is.na(data), 1, sum) != 0, ]
      cat("Warning:", length(allempty), "Cases with all variables missing have been removed
          from the data.\n")
    }
    data <- OrderMissing(data)
  }
  if (class(data) == "orderpattern") {
    allempty <- which(apply(!is.na(data$data), 1, sum) == 0)
    if (length(allempty) != 0) {
      data <- data$data
      data <- data[apply(!is.na(data), 1, sum) != 0, ]
      cat("Warning:", length(allempty), "Cases with all variables missing have been removed
          from the data.\n")
      data <- OrderMissing(data)
    }
  }
  if (length(data$data) == 0) {
    stop("Warning: Data is empty")
  }
  if (ncol(data$data) < 2) {
    stop("Warning: More than 1 variable is required.\n")
  }
  y <- data$data
  patused <- data$patused
  spatcnt <- data$spatcnt
  if (is.na(mu[1])) {
    mu <- matrix(0, ncol(y), 1)
    sig <- diag(1, ncol(y))
  }
  itcnt <- 0
  em <- 0
  repeat {
    emtemp <- Sexpect(y, mu, sig, patused, spatcnt)
    ysbar <- emtemp$ysbar
    sstar <- emtemp$sstar
    em <- max(abs(sstar - mu %*% t(mu) - sig), abs(mu - ysbar))
    mu <- ysbar
    sig <- sstar - mu %*% t(mu)
    itcnt <- itcnt + 1
    if (!(em > tol || itcnt < 2)) break()
  } # end repeat
  rownames(mu) <- colnames(y)
  colnames(sig) <- colnames(y)
  if (Hessian) {
    templist <- Ddf(y, mu, sig)
    hessian <- templist$dd
    stderror <- templist$se
    return(list(mu = mu, sig = sig, hessian = hessian, stderror = stderror, iteration = itcnt))
  }
  return(list(mu = mu, sig = sig, iteration = itcnt))
}
#------------------------------------------------------------------
Sexpect <- function(y, mu, sig, patused, spatcnt) {
  n <- nrow(y)
  pp <- ncol(y)
  sstar <- matrix(0, pp, pp)
  a <- nrow(mu)
  b <- ncol(mu)
  ysbar <- matrix(0, a, b)
  first <- 1
  for (i in 1:length(spatcnt)) {
    ni <- spatcnt[i] - first + 1
    stemp <- matrix(0, pp, pp)
    indm <- which(is.na(patused[i, ]))
    indo <- which(!is.na(patused[i, ]))
    yo <- matrix(y[first:spatcnt[i], indo], ni, length(indo))
    first <- spatcnt[i] + 1
    muo <- mu[indo]
    mum <- mu[indm]
    sigoo <- sig[indo, indo]
    sigooi <- solve(sigoo)
    soo <- t(yo) %*% yo
    stemp[indo, indo] <- soo
    ystemp <- matrix(0, ni, pp)
    ystemp[, indo] <- yo
    if (length(indm) != 0) {
      sigmo <- matrix(sig[indm, indo], length(indm), length(indo))
      sigmm <- sig[indm, indm]
      temp1 <- matrix(mum, ni, length(indm), byrow = TRUE)
      temp2 <- yo - matrix(muo, ni, length(indo), byrow = TRUE)
      ym <- temp1 + temp2 %*% sigooi %*% t(sigmo)
      som <- t(yo) %*% ym
      smm <- ni * (sigmm - sigmo %*% sigooi %*% t(sigmo)) + t(ym) %*% ym
      stemp[indo, indm] <- som
      stemp[indm, indo] <- t(som)
      stemp[indm, indm] <- smm
      ystemp[, indm] <- ym
    } # end if
    sstar <- sstar + stemp
    if (ni == 1) {
      ysbar <- t(ystemp) + ysbar
    } else {
      ysbar <- apply(ystemp, 2, sum) + ysbar
    }
  } # end for
  ysbar <- (1 / n) * ysbar
  sstar <- (1 / n) * sstar
  sstar <- (sstar + t(sstar)) / 2

  return(list(ysbar = ysbar, sstar = sstar))
}

Impute <- function(data, mu = NA, sig = NA, imputation.method = "Normal", resid = NA) { # Check if data is not ordered change it to ordered form
  if (!is.matrix(data) && class(data) != "orderpattern") {
    stop("Warning: data must have the classes of matrix or orderpattern.\n")
  }
  if (is.matrix(data)) {
    allempty <- which(apply(!is.na(data), 1, sum) == 0)
    if (length(allempty) != 0) {
      data <- data[apply(!is.na(data), 1, sum) != 0, ]
      cat("Warning:", length(allempty), "Cases with all variables missing have been removed
          from the data.\n")
    }
    data <- OrderMissing(data)
  }
  if (class(data) == "orderpattern") {
    allempty <- which(apply(!is.na(data$data), 1, sum) == 0)
    if (length(allempty) != 0) {
      data <- data$data
      data <- data[apply(!is.na(data), 1, sum) != 0, ]
      cat("Warning:", length(allempty), "Cases with all variables missing have been removed
          from the data.\n")
      data <- OrderMissing(data)
    }
  }
  if (length(data$data) == 0) {
    stop("Warning: data is empty")
  }
  if (ncol(data$data) < 2) {
    stop("More than 1 variable is required.\n")
  }
  y <- data$data
  patused <- data$patused
  spatcnt <- data$spatcnt
  patcnt <- data$patcnt
  g <- data$g
  caseorder <- data$caseorder
  spatcntz <- c(0, spatcnt)
  p <- ncol(y)
  n <- nrow(y)
  yimp <- y
  use.normal <- TRUE

  #---------impute the missing data with Servestava method(simple Imputation)--------
  if (imputation.method == "Dist.Free") {
    if (is.na(mu[1])) {
      ybar <- matrix(0, p, 1)
      sbar <- diag(1, p)
      iscomp <- apply(patused, 1, sum, na.rm = TRUE) == p

      cind <- which(iscomp)
      ncomp <- patcnt[cind]
      use.normal <- FALSE
      if (ncomp >= 10 && ncomp >= 2 * p) {
        compy <- y[seq(spatcntz[cind] + 1, spatcntz[cind + 1]), ]
        ybar <- matrix(apply(compy, 2, mean))
        sbar <- stats::cov(compy)
        if (is.na(resid[1])) {
          resid <- (ncomp / (ncomp - 1))^.5 *
            (compy - matrix(ybar, ncomp, p, byrow = TRUE))
        }
      } else {
        stop("Warning: There is not sufficient number of complete cases.\n  Dist.free imputation requires a least 10 complete cases\n  or 2*number of variables, whichever is bigger.\n")
      }
    }
    if (!is.na(mu[1])) {
      ybar <- mu
      sbar <- sig
      iscomp <- apply(patused, 1, sum, na.rm = TRUE) == p
      cind <- which(iscomp)
      ncomp <- patcnt[cind]
      use.normal <- FALSE
      compy <- y[seq(spatcntz[cind] + 1, spatcntz[cind + 1]), ]
      if (is.na(resid[1])) {
        resid <- (ncomp / (ncomp - 1))^.5 *
          (compy - matrix(ybar, ncomp, p, byrow = TRUE))
      }
    }
    indsample <- sample(ncomp, n - ncomp, replace = TRUE)
    resstar <- resid[indsample, ]
    indres1 <- 1
    for (i in 1:g) {
      if (sum(patused[i, ], na.rm = TRUE) != p) # choose a pattern not completely obsered
        {
          test <- y[(spatcntz[i] + 1):spatcntz[i + 1], ]
          indres2 <- indres1 + patcnt[i] - 1
          test <- MimputeS(
            matrix(test, ncol = p), patused[i, ], ybar, sbar,
            matrix(resstar[indres1:indres2, ], ncol = p)
          )
          indres1 <- indres2 + 1

          yimp[(spatcntz[i] + 1):spatcntz[i + 1], ] <- test
        }
    } # end loop
  }
  if (imputation.method == "Normal" | use.normal) {
    #-----------impute the missing data with normal assumption------------
    if (is.na(mu[1])) {
      emest <- Mls(data, tol = 1e-6)
      mu <- emest$mu
      sig <- emest$sig
    }
    for (i in 1:g) {
      if (sum(patused[i, ], na.rm = TRUE) != p) # choose a pattern not completely obsered
        {
          test <- y[(spatcntz[i] + 1):spatcntz[i + 1], ]
          test <- Mimpute(matrix(test, ncol = p), patused[i, ], mu, sig)
          yimp[(spatcntz[i] + 1):spatcntz[i + 1], ] <- test
        }
    } # end loop
  }
  yimpord <- yimp
  yimp <- yimp[order(caseorder), ]
  imputed <- list(
    yimp = yimp, yimpOrdered = yimpord, caseorder = caseorder,
    patused = patused, patcnt = patcnt
  )

  imputed
}

#--------------------------------------------------------------------------
Mimpute <- function(data, patused, mu, sig) {
  # This function imputes the missing data base on multivariate-normal, given the
  # observed and mu and sig for a single pattern it goeas through each patterns and uses the
  # conditional distribution of missing given observed and mu and sig to
  # impute from the appropriate posterior distribution
  ni <- nrow(data)
  pp <- ncol(data)
  indm <- which(is.na(patused))
  indo <- which(!is.na(patused))
  pm <- length(indm)
  po <- length(indo)
  muo <- mu[indo]
  mum <- mu[indm]
  sigooi <- solve(sig[indo, indo])
  sigmo <- matrix(sig[indm, indo], pm, po)
  sigmm <- matrix(sig[indm, indm], pm, pm)
  ss1 <- sigmo %*% sigooi
  varymiss <- sigmm - ss1 %*% t(sigmo)
  expymiss <- matrix(mum, ni, pm, byrow = TRUE) +
    (data[, indo] - matrix(muo, ni, po, byrow = TRUE)) %*% t(ss1)
  if (pm == 1) {
    a <- sqrt(varymiss)
  } else {
    svdvar <- svd(varymiss)
    a <- diag(sqrt(svdvar$d)) %*% t(svdvar$u)
  }
  data[, indm] <- matrix(stats::rnorm(ni * pm), ni, pm) %*% a + expymiss
  data
}
#---------------------------------------------------------------------------
MimputeS <- function(data, patused, y1, s1, e) {
  # This function imputes the missing data by Srivastava method,
  # given the observed and ybar and s from complete data set
  # for a single pattern it goeas through each patterns and uses the
  # linear regresion to predict missing given obsrved data and add the
  # residual to impute missing data

  ni <- nrow(data)
  pp <- ncol(data)
  indm <- which(is.na(patused))
  indo <- which(!is.na(patused))
  pm <- length(indm)
  po <- length(indo)
  a <- matrix(s1[indm, indo], pm, po) %*% solve(s1[indo, indo])
  dif <- data[, indo] - matrix(y1[indo], ni, po, byrow = TRUE)
  z <- matrix(y1[indm], ni, pm, byrow = TRUE) + dif %*% t(a)
  etta <- matrix(e[, indm], ni, pm) - matrix(e[, indo], ni, po) %*% t(a)
  zij <- z + etta
  data[, indm] <- zij
  data
}

Hawkins <- function(data, spatcnt) {
  # This function performs the Hawkin's method for testing equality of
  # covariances among groups. It is assumed that the data in y is ordered so
  # that the cases from the same group are adjacent.
  # gind is a vector indicating the end of the cases for each group. For
  # example [8 23 30] means that there are three groups with cases 1-8 in one
  # group, 9-23 and 24-30 in another.

  # Also in the output, it gives the statistic computed for each group in a
  # cell array, called A. This statistic should be tested for uniformity.
  # also ni(i) gives the number of components of each [a[i]]
  y <- data
  n <- nrow(y)
  p <- ncol(y)
  g <- length(spatcnt)
  spool <- matrix(0, p, p)
  gind <- c(0, spatcnt)
  ygc <- matrix(0, n, p)
  ni <- matrix(0, g, 1)
  for (i in 1:g)
  {
    yg <- y[seq(gind[i] + 1, gind[i + 1]), ]
    ni[i] <- nrow(yg)
    spool <- spool + (ni[i] - 1) * stats::cov(yg)
    ygmean <- apply(yg, 2, mean)
    ygc[seq(gind[i] + 1, gind[i + 1]), ] <-
      yg - matrix(ygmean, ni[i], p, byrow = TRUE)
  }
  spool <- spool / (n - g)
  spool <- solve(spool)
  f <- matrix(0, n, 1)
  nu <- n - g - 1
  a <- vector("list", g)
  for (i in 1:g)
  {
    vij <- ygc[seq(gind[i] + 1, gind[i + 1]), ]
    vij <- apply(vij %*% spool * vij, 1, sum)
    vij <- vij * ni[i]
    f[seq(gind[i] + 1, gind[i + 1])] <- ((n - g - p) * vij) /
      (p * ((ni[i] - 1) * (n - g) - vij))
    a[[i]] <- 1 - stats::pf(f[seq(gind[i] + 1, gind[i + 1])], p, (nu - p + 1))
  }
  list(fij = f, a = a, ni = ni)
}

LegNorm <- function(x) {
  # This function evaluates Legendre polynomials on [0,1] (note not [-1,1] at
  # a value x, and the polynomial are such that they have norm 1. It only
  # calculates p1,p2,p3, and p4 [Note p0=1, so it is not reurned]
  # x can be a vector or a matrix.
  # in return, each element of P1, p2, p3, p4 the values of the poly at each
  # component of x

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  } # end if
  x <- 2 * x - 1 # Transforming the legenre's to 0,1, and we will divide by the norm in each case
  p0 <- matrix(1, nrow(x), ncol(x))
  p1 <- x
  p2 <- (3 * x * p1 - p0) / 2
  p3 <- (5 * x * p2 - 2 * p1) / 3
  p4 <- (7 * x * p3 - 3 * p2) / 4
  p1 <- sqrt(3) * p1
  p2 <- sqrt(5) * p2
  p3 <- sqrt(7) * p3
  p4 <- 3 * p4
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

Ddf <- function(data, mu, sig) {
  y <- data
  n <- nrow(y)
  p <- ncol(y)
  ns <- p * (p + 1) / 2
  nparam <- ns + p
  ddss <- matrix(0, ns, ns)
  ddmm <- matrix(0, p, p)
  ddsm <- matrix(0, ns, p)
  for (i in 1:n) {
    obs <- which(!is.na(y[i, ]))
    lo <- length(obs)
    tmp <- cbind(rep(obs, 1, each = lo), rep(obs, lo))
    tmp <- matrix(tmp[tmp[, 1] >= tmp[, 2], ], ncol = 2)
    lolo <- lo * (lo + 1) / 2
    subsig <- sig[obs, obs]
    submu <- mu[obs, ]
    temp <- matrix(y[i, obs] - submu, nrow = 1)
    a <- solve(subsig)
    b <- a %*% (2 * t(temp) %*% temp - subsig) * a
    d <- temp %*% a
    ddimm <- 2 * a
    # ==================== DD(mu, mu)
    ddmm[obs, obs] <- ddmm[obs, obs] + ddimm
    # ==================== DD(sig, mu)
    rcnt <- 0
    ddism <- matrix(0, lolo, lo)
    for (k in 1:lo) {
      for (l in 1:k) {
        rcnt <- rcnt + 1
        ccnt <- 0
        for (kk in 1:lo) {
          ccnt <- ccnt + 1
          ddism[rcnt, ccnt] <- 2 * (1 - 0.5 * (k == l)) *
            (a[kk, l] %*% d[k] + a[kk, k] %*% d[l])
        }
      }
    }
    for (k in 1:lolo) {
      par1 <- tmp[k, 1] * (tmp[k, 1] - 1) / 2 + tmp[k, 2]
      for (j in 1:lo) {
        ddsm[par1, obs[j]] <- ddsm[par1, obs[j]] + ddism[k, j]
      }
    }
    #------------------------ Test part
    ssi <- matrix(0, lolo, lolo)
    for (m in 1:lolo) {
      u <- which(obs == tmp[m, 1])
      v <- which(obs == tmp[m, 2])
      for (q in 1:m) {
        k <- which(obs == tmp[q, 1])
        l <- which(obs == tmp[q, 2])
        ssi[m, q] <- (b[v, k] * a[l, u] + b[v, l] * a[k, u] + b[u, k] * a[l, v] +
          b[u, l] * a[k, v]) * (1 - 0.5 * (u == v)) * (1 - 0.5 * (k == l))
      }
    }
    for (k in 1:lolo) {
      par1 <- tmp[k, 1] * (tmp[k, 1] - 1) / 2 + tmp[k, 2]
      for (l in 1:k) {
        par2 <- tmp[l, 1] * (tmp[l, 1] - 1) / 2 + tmp[l, 2]
        ddss[par1, par2] <- ddss[par1, par2] + ssi[k, l]
        ddss[par2, par1] <- ddss[par1, par2]
      }
    }
  }
  dd <- -1 * rbind(cbind(ddmm, t(ddsm)), cbind(ddsm, ddss)) / 2
  se <- -solve(dd)
  list(dd = dd, se = se)
}
