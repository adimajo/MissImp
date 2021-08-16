#' Fast Imputation of Missing Values by Chained Random Forests
#'
#' @description \code{missRanger_mod} is a modified imputation function of \code{missRanger} in 'missRanger' package.
#' The only difference is that the disjunctive imputed dataset is also returned (with categorical columns in form of onehot probability vector).
#'
#' Uses the "ranger" package (Wright & Ziegler) to do fast missing value imputation by chained random forests, see Stekhoven & Buehlmann and Van Buuren & Groothuis-Oudshoorn.
#' Between the iterative model fitting, it offers the option of predictive mean matching. This firstly avoids imputation with values not present in the original data (like a value 0.3334 in a 0-1 coded variable). Secondly, predictive mean matching tries to raise the variance in the resulting conditional distributions to a realistic level. This allows to do multiple imputation when repeating the call to missRanger().
#' The iterative chaining stops as soon as \code{maxiter} is reached or if the average out-of-bag estimate of performance stops improving. In the latter case, except for the first iteration, the second last (i.e. best) imputed data is returned.
#'
#' A note on `mtry`: Be careful when passing a non-default `mtry` to `ranger()` because the number of available covariables might be growing during the first iteration, depending on the missing pattern. Values \code{NULL} (default) and 1 are safe choices. Additionally, recent versions of `ranger()` allow `mtry` to be a single-argument function of the number of available covariables, e.g. `mtry = function(m) max(1, m %/% 3)`.
#'
#' @importFrom stats var reformulate terms.formula predict setNames
#'
#' @param data A \code{data.frame} or \code{tibble} with missing values to impute.
#' @param formula A two-sided formula specifying variables to be imputed (left hand side) and variables used to impute (right hand side). Defaults to . ~ ., i.e. use all variables to impute all variables.
#' If e.g. all variables (with missings) should be imputed by all variables except variable "ID", use . ~ . - ID. Note that a "." is evaluated separately for each side of the formula. Further note that variables
#' with missings must appear in the left hand side if they should be used on the right hand side.
#' @param pmm.k Number of candidate non-missing values to sample from in the predictive mean matching steps. 0 to avoid this step.
#' @param maxiter Maximum number of chaining iterations.
#' @param seed Integer seed to initialize the random generator.
#' @param verbose Controls how much info is printed to screen. 0 to print nothing. 1 (default) to print a "." per iteration and variable, 2 to print the OOB prediction error per iteration and variable (1 minus R-squared for regression).
#' Furthermore, if \code{verbose} is positive, the variables used for imputation are listed as well as the variables to be imputed (in the imputation order). This will be useful to detect if some variables are unexpectedly skipped.
#' @param returnOOB Logical flag. If TRUE, the final average out-of-bag prediction error is added to the output as attribute "oob". This does not work in the special case when the variables are imputed univariately.
#' @param case.weights Vector with non-negative case weights.
#' @param col_cat Column index of categorical variables.
#' @param ... Arguments passed to \code{ranger()}. If the data set is large, better use less trees (e.g. \code{num.trees = 20}) and/or a low value of \code{sample.fraction}.
#' The following arguments are e.g. incompatible with \code{ranger}: \code{write.forest}, \code{probability}, \code{split.select.weights}, \code{dependent.variable.name}, and \code{classification}.
#'
#' @return \code{ximp} An imputed \code{data.frame}.
#' @return \code{ximp.disj} A disjunctive imputed \code{data.frame}.
#' For example if Y7 has levels a and b, then in \code{ximp.disj} the column Y7_1 corresponds to the probability of Y7=a.
#' @references
#' \enumerate{
#'   \item Wright, M. N. & Ziegler, A. (2016). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. Journal of Statistical Software, in press. <arxiv.org/abs/1508.04409>.
#'   \item Stekhoven, D.J. and Buehlmann, P. (2012). 'MissForest - nonparametric missing value imputation for mixed-type data', Bioinformatics, 28(1) 2012, 112-118. https://doi.org/10.1093/bioinformatics/btr597.
#'   \item Van Buuren, S., Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67. http://www.jstatsoft.org/v45/i03/
#' }
#' @export
#'

missRanger_mod <- function(data, formula = . ~ ., pmm.k = 0L, maxiter = 10L, seed = NULL,
                           verbose = 1, returnOOB = FALSE, case.weights = NULL, col_cat = c(), ...) {
  is_ranger_package_installed()
  if (verbose) {
    cat("\nMissing value imputation by random forests\n")
  }

  ## Add:
  exist_cat <- !all(c(0, col_cat) == c(0))
  if (exist_cat) {
    name_cat <- colnames(data)[col_cat]
    # Deal with the problem that nlevels(df[[col]]) > length(unique(df[[col]]))
    for (col in name_cat) {
      data[[col]] <- factor(as.character(data[[col]]))
    }
    # remember the levels for each categorical column
    dict_lev <- dict_level(data, col_cat)
    # preserve colnames for ximp.disj
    dummy <- dummyVars(" ~ .", data = data, sep = "_")
    col_names.disj <- colnames(data.frame(predict(dummy, newdata = data)))
    # represent the factor columns with their ordinal levels
    data <- factor_ordinal_encode(data, col_cat)
    # Create dict_cat with categroical columns
    dict_cat <- dict_onehot(data, col_cat)
  }
  ## Add: Last iteration will be used to predict the onehot probability for the categorical columns
  maxiter <- maxiter - 1

  # 1) INITIAL CHECKS

  stopifnot(
    is.data.frame(data), dim(data) >= 1L,
    inherits(formula, "formula"),
    length(formula <- as.character(formula)) == 3L,
    is.numeric(pmm.k), length(pmm.k) == 1L, pmm.k >= 0L,
    is.numeric(maxiter), length(maxiter) == 1L, maxiter >= 1L,
    !(c(
      "write.forest", "probability", "split.select.weights",
      "dependent.variable.name", "classification"
    ) %in% names(list(...)))
  )

  if (!is.null(case.weights)) {
    stopifnot(length(case.weights) == nrow(data), !anyNA(case.weights))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # 2) SELECT AND CONVERT VARIABLES TO IMPUTE

  # Extract lhs and rhs from formula
  relevantVars <- lapply(formula[2:3], function(z) {
    attr(terms.formula(
      reformulate(z),
      data = data[1, ]
    ), "term.labels")
  })

  # Pick variables from lhs with some but not all missings
  toImpute <- relevantVars[[1]][vapply(data[, relevantVars[[1]], drop = FALSE],
    FUN.VALUE = TRUE, function(z) anyNA(z) && !all(is.na(z))
  )]

  # Try to convert special variables to numeric/factor in order to be safely predicted by ranger
  converted <- convert(data[, toImpute, drop = FALSE], check = TRUE)
  data[, toImpute] <- converted$X

  # Remove variables that cannot be safely converted
  visitSeq <- setdiff(toImpute, converted$bad)

  if (verbose) {
    cat("\n  Variables to impute:\t\t")
    cat(visitSeq, sep = ", ")
  }

  if (!length(visitSeq)) {
    if (verbose) {
      cat("\n")
    }
    return(data)
  }

  # Get missing indicators and order variables by number of missings
  dataNA <- is.na(data[, visitSeq, drop = FALSE])
  visitSeq <- names(sort(colSums(dataNA)))

  # 3) SELECT VARIABLES USED TO IMPUTE

  # Variables on the rhs should either appear in "visitSeq" or do not contain any missings
  imputeBy <- relevantVars[[2]][relevantVars[[2]] %in% visitSeq |
    !vapply(data[, relevantVars[[2]], drop = FALSE], anyNA, TRUE)]
  completed <- setdiff(imputeBy, visitSeq)

  if (verbose) {
    cat("\n  Variables used to impute:\t")
    cat(imputeBy, sep = ", ")
  }

  # 4) IMPUTATION

  # Initialization
  j <- 1L # iterator
  crit <- TRUE # criterion on OOB prediction error to keep iterating
  verboseDigits <- 4L # formatting of OOB prediction errors (if verbose = 2)
  predError <- setNames(rep(1, length(visitSeq)), visitSeq)

  if (verbose >= 2) {
    cat("\n", abbreviate(visitSeq, minlength = verboseDigits + 2L), sep = "\t")
  }

  dataLast <- data
  # Looping over iterations and variables to impute
  while (crit && j <= maxiter) {
    if (verbose) {
      cat("\niter ", j, ":\t", sep = "")
    }
    dataLast2 <- dataLast
    dataLast <- data
    predErrorLast <- predError

    for (v in visitSeq) {
      v.na <- dataNA[, v]

      if (length(completed) == 0L) {
        data[[v]] <- imputeUnivariate(data[[v]])
      } else {
        fit <- ranger::ranger(
          formula = reformulate(completed, response = v),
          data = data[!v.na, union(v, completed), drop = FALSE],
          case.weights = case.weights[!v.na]
        )
        pred <- predict(fit, data[v.na, completed, drop = FALSE])$predictions
        data[v.na, v] <- if (pmm.k) {
          pmm(
            xtrain = fit$predictions,
            xtest = pred,
            ytrain = data[[v]][!v.na],
            k = pmm.k
          )
        } else {
          pred
        }
        predError[[v]] <- fit$prediction.error / (if (fit$treetype == "Regression") var(data[[v]][!v.na]) else 1)

        if (is.nan(predError[[v]])) {
          predError[[v]] <- 0
        }
      }

      if (j == 1L && (v %in% imputeBy)) {
        completed <- union(completed, v)
      }

      if (verbose == 1) {
        cat(".")
      } else if (verbose >= 2) {
        cat(format(round(predError[[v]], verboseDigits), nsmall = verboseDigits), "\t")
      }
    }

    j <- j + 1L
    crit <- mean(predError) < mean(predErrorLast)
  }

  if (verbose) {
    cat("\n")
  }


  ##### Add: Get the onehot probability result for categorical columns
  if (exist_cat) {
    dummy <- dummyVars(" ~ .", data = data, sep = "_")
    data.disj <- data.frame(predict(dummy, newdata = data))
  }

  if (verbose) {
    cat("Last iter:\t", sep = "")
  }

  if (j == maxiter + 1) {
    # Last iteration
    dataLast <- data
    predErrorLast <- predError
  } else {
    # if crit, use dataLast2 to redo the iteration
    data <- dataLast2
  }

  for (v in visitSeq) {
    v.na <- dataNA[, v]

    if (length(completed) == 0L) {
      data[[v]] <- imputeUnivariate(data[[v]])
    } else {
      ## Add: prediction one onehot form for categorical columns
      fit <- ranger::ranger(
        formula = reformulate(completed, response = v),
        data = data[!v.na, union(v, completed), drop = FALSE],
        case.weights = case.weights[!v.na]
      )
      pred <- predict(fit, data[v.na, completed, drop = FALSE])$predictions
      if (exist_cat) {
        if (v %in% name_cat) {
          fit.disj <- ranger::ranger(
            formula = reformulate(completed, response = v),
            data = data[!v.na, union(v, completed), drop = FALSE],
            case.weights = case.weights[!v.na], probability = TRUE
          )
          pred.disj <- predict(fit.disj, data[v.na, completed, drop = FALSE])$predictions
          if (pmm.k) {
            data.disj[v.na, dict_cat[[v]]] <- pmm(
              xtrain = fit.disj$predictions,
              xtest = pred.disj,
              ytrain = data.disj[[v]][!v.na],
              k = pmm.k
            )
          } else if (ncol(data.disj[v.na, dict_cat[[v]]]) == ncol(pred.disj)) {
            data.disj[v.na, dict_cat[[v]]] <- pred.disj
          } else { # When there are dropped unused levels
            colnames(pred.disj) <- paste0(v, "_", colnames(pred.disj))
            data.disj[v.na, colnames(pred.disj)] <- pred.disj
          }
        } else {
          data.disj[v.na, v] <- if (pmm.k) {
            pmm(
              xtrain = fit$predictions,
              xtest = pred,
              ytrain = data[[v]][!v.na],
              k = pmm.k
            )
          } else {
            pred
          }
        }
      }
      data[v.na, v] <- if (pmm.k) {
        pmm(
          xtrain = fit$predictions,
          xtest = pred,
          ytrain = data[[v]][!v.na],
          k = pmm.k
        )
      } else {
        pred
      }
      predError[[v]] <- fit$prediction.error / (if (fit$treetype == "Regression") var(data[[v]][!v.na]) else 1)


      if (is.nan(predError[[v]])) {
        predError[[v]] <- 0
      }
    }

    if (j == 1L && (v %in% imputeBy)) {
      completed <- union(completed, v)
    }

    if (verbose == 1) {
      cat(".")
    } else if (verbose >= 2) {
      cat(format(round(predError[[v]], verboseDigits), nsmall = verboseDigits), "\t")
    }
  }

  j <- j + 1L
  crit <- mean(predError) < mean(predErrorLast)
  ######################################

  if (verbose) {
    cat("\n")
  }

  if (j == 2L || (j == maxiter && crit)) {
    dataLast <- data
    predErrorLast <- predError
  }

  if (returnOOB) {
    attr(dataLast, "oob") <- predErrorLast
  }

  # Revert the conversions
  if (!exist_cat) {
    data.disj <- data
  }
  ximp <- revert(converted, X = data)

  # change back to original levels
  if (exist_cat) {
    for (col in name_cat) {
      levels(ximp[[col]]) <- dict_lev[[col]]
    }
    colnames(data.disj) <- col_names.disj
  }
  return(list(ximp = ximp, ximp.disj = data.disj))
}







#' A version of \code{typeof} internally used by \code{missRanger}.
#'
#' @description Returns either "numeric" (double or integer), "factor", "character", "logical", "special" (mode numeric, but neither double nor integer) or "" (otherwise).
#' \code{missRanger} requires this information to deal with response types not natively supported by \code{ranger}.
#'
#' @author Michael Mayer
#'
#' @param object Any object.
#'
#' @return A string.
typeof2 <- function(object) {
  if (is.numeric(object)) {
    "numeric"
  } else
  if (is.factor(object)) {
    "factor"
  } else
  if (is.character(object)) {
    "character"
  } else
  if (is.logical(object)) {
    "logical"
  } else
  if (mode(object) == "numeric") "special" else ""
}

#' Conversion of non-factor/non-numeric variables.
#'
#' @description Converts non-factor/non-numeric variables in a data frame to factor/numeric. Stores information to revert back.
#'
#' @author Michael Mayer
#'
#' @param X A data frame.
#' @param check If \code{TRUE}, the function checks if the converted columns can be reverted without changes.
#'
#' @return A list with the following elements: \code{X} is the converted data frame, \code{vars}, \code{types}, \code{classes} are the names, types and classes of the converted variables. Finally, \code{bad} names variables in \code{X} that should have been converted but could not.
convert <- function(X, check = FALSE) {
  stopifnot(is.data.frame(X))

  if (!ncol(X)) {
    return(list(
      X = X, bad = character(0), vars = character(0),
      types = character(0), classes = character(0)
    ))
  }

  types <- vapply(X, typeof2, FUN.VALUE = "")
  bad <- types == "" | if (check) {
    mapply(function(a, b) {
      isFALSE(all.equal(a, b))
    }, X, revert(convert(X)))
  } else {
    FALSE
  }
  types <- types[!(types %in% c("numeric", "factor") | bad)]
  vars <- names(types)
  classes <- lapply(X[, vars, drop = FALSE], class)

  X[, vars] <- lapply(X[, vars, drop = FALSE], function(v) {
    if (is.character(v) || is.logical(v)) as.factor(v) else as.numeric(v)
  })

  list(X = X, bad = names(X)[bad], vars = vars, types = types, classes = classes)
}

#' Revert conversion.
#'
#' @description Reverts conversions done by \code{convert}.
#'
#' @author Michael Mayer
#'
#' @param con A list returned by \code{convert}.
#' @param X A data frame with some columns to be converted back according to the information stored in \code{converted}.
#'
#' @return A data frame.
revert <- function(con, X = con$X) {
  stopifnot(c("vars", "types", "classes") %in% names(con), is.data.frame(X))

  if (!length(con$vars)) {
    return(X)
  }

  f <- function(v, ty, cl) {
    switch(ty,
      logical = as.logical(v),
      character = as.character(v),
      special = {
        class(v) <- cl
        v
      },
      v
    )
  }
  X[, con$vars] <- Map(f, X[, con$vars, drop = FALSE], con$types, con$classes)
  X
}

#' Univariate Imputation
#'
#' Fills missing values of a vector, matrix or data frame by sampling with replacement from the non-missing values. For data frames, this sampling is done within column.
#'
#' @param x A vector, matrix or data frame.
#' @param v A character vector of column names to impute (only relevant if \code{x} is a data frame). The default \code{NULL} imputes all columns.
#' @param seed An integer seed.
#'
#' @return \code{x} with imputed values.
#' @export
#'

imputeUnivariate <- function(x, v = NULL, seed = NULL) {
  stopifnot(is.atomic(x) || is.data.frame(x))

  if (!is.null(seed)) {
    set.seed(seed)
  }

  imputeVec <- function(z) {
    na <- is.na(z)
    if ((s <- sum(na))) {
      if (s == length(z)) {
        stop("No non-missing elements to sample from.")
      }
      z[na] <- sample(z[!na], s, replace = TRUE)
    }
    z
  }

  # vector or matrix
  if (is.atomic(x)) {
    return(imputeVec(x))
  }

  # data frame
  v <- if (is.null(v)) names(x) else intersect(v, names(x))
  x[, v] <- lapply(x[, v, drop = FALSE], imputeVec)

  x
}


#' Predictive Mean Matching
#'
#' For each value in the prediction vector \code{xtest}, one of the closest \code{k} values in the prediction vector \code{xtrain} is randomly chosen and its observed value in \code{ytrain} is returned.
#'
#' @importFrom stats rmultinom
#'
#' @param xtrain Vector with predicted values in the training data. Can be of type logical, numeric, character, or factor.
#' @param xtest Vector as \code{xtrain} with predicted values in the test data. Missing values are not allowed.
#' @param ytrain Vector of the observed values in the training data. Must be of same length as \code{xtrain}. Missing values in either of \code{xtrain} or \code{ytrain} will be dropped in a pairwise manner.
#' @param k Number of nearest neighbours to sample from.
#' @param seed Integer random seed.
#'
#' @return Vector of the same length as \code{xtest} with values from \code{xtrain}.
#' @export
#'
#' @examples
#' pmm(xtrain = c(0.2, 0.2, 0.8), xtest = 0.3, ytrain = c(0, 0, 1)) # 0
#' pmm(xtrain = c(TRUE, FALSE, TRUE), xtest = FALSE, ytrain = c(2, 0, 1)) # 0
#' pmm(xtrain = c(0.2, 0.8), xtest = 0.3, ytrain = c("A", "B"), k = 2) # "A" or "B"
#' pmm(xtrain = c("A", "A", "B"), xtest = "A", ytrain = c(2, 2, 4), k = 2) # 2
#' pmm(xtrain = factor(c("A", "B")), xtest = factor("C"), ytrain = 1:2) # 2
pmm <- function(xtrain, xtest, ytrain, k = 1L, seed = NULL) {
  stopifnot(
    length(xtrain) == length(ytrain),
    sum(ok <- !is.na(xtrain) & !is.na(ytrain)) >= 1L,
    (nt <- length(xtest)) >= 1L, !anyNA(xtest),
    mode(xtrain) %in% c("logical", "numeric", "character"),
    mode(xtrain) == mode(xtest),
    k >= 1L
  )
  is_FNN_package_installed()
  xtrain <- xtrain[ok]
  ytrain <- ytrain[ok]

  # Handle trivial case
  if (length(u <- unique(ytrain)) == 1L) {
    return(rep(u, nt))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # STEP 1: Turn xtrain and xtest into numbers.
  # Handles the case of inconsistent factor levels of xtrain and xtest.
  if (is.factor(xtrain) && (nlevels(xtrain) != nlevels(xtest) ||
    !all(levels(xtrain) == levels(xtest)))) {
    xtrain <- as.character(xtrain)
    xtest <- as.character(xtest)
  }

  # Turns character vectors into factors.
  if (is.character(xtrain)) {
    lvl <- unique(c(xtrain, xtest))
    xtrain <- factor(xtrain, levels = lvl)
    xtest <- factor(xtest, levels = lvl)
  }

  # Turns everything into numbers.
  if (!is.numeric(xtrain) && mode(xtrain) %in% c("logical", "numeric")) {
    xtrain <- as.numeric(xtrain)
    xtest <- as.numeric(xtest)
  }

  # STEP 2: PMM based on k-nearest neightbour.
  k <- min(k, length(xtrain))
  nn <- FNN::knnx.index(xtrain, xtest, k)
  take <- t(rmultinom(nt, 1L, rep(1L, k)))
  ytrain[rowSums(nn * take)]
}
