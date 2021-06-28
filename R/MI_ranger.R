missRanger_mod_draw <- function(data, formula = . ~ ., pmm.k = 0L, maxiter = 10L, seed = NULL,
                                verbose = 1, returnOOB = FALSE, case.weights = NULL, col_cat = c(), ...) {
  if (verbose) {
    cat("\nMissing value imputation by random forests\n")
  }

  ## Add: Create dict_cat with categroical columns
  exist_cat <- !all(c(0, col_cat) == c(0))
  if (exist_cat) {
    dict_cat <- dict_onehot(data, col_cat)
    name_cat <- colnames(data)[col_cat]
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
        # pred <- predict(fit, data[v.na, completed, drop = FALSE])$predictions
        ### MI###
        pred1 <- predict(fit, data[v.na, completed, drop = FALSE], predict.all = TRUE)$predictions
        pred_draw <- apply(pred1, 1, sample, 1)
        #######
        data[v.na, v] <- if (pmm.k) {
          pmm(
            xtrain = fit$predictions,
            xtest = pred_draw, # pred,
            ytrain = data[[v]][!v.na],
            k = pmm.k
          )
        } else {
          # pred
          pred_draw
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
    ## TODO: Change crit####
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
      ### MI###
      pred1 <- predict(fit, data[v.na, completed, drop = FALSE], predict.all = TRUE)$predictions
      pred_draw <- apply(pred1, 1, sample, 1)
      #######
      if (exist_cat) {
        if (v %in% name_cat) {
          fit.disj <- ranger::ranger(
            formula = reformulate(completed, response = v),
            data = data[!v.na, union(v, completed), drop = FALSE],
            case.weights = case.weights[!v.na], probability = TRUE
          )
          # pred.disj <- predict(fit.disj, data[v.na, completed, drop = FALSE])$predictions
          ### MI###
          pred1.disj <- predict(fit.disj, data[v.na, completed, drop = FALSE], predict.all = TRUE)$predictions
          tree_chosen <- sample(c(1:fit.disj$num.trees), dim(pred1.disj)[1], replace = TRUE) # choose tree result for each ligne
          # Draw from tress results
          df_tmp <- data.frame(tree_chosen)
          df_tmp$tree_chosen <- factor(df_tmp$tree_chosen)
          levels(df_tmp$tree_chosen) <- c(1:fit.disj$num.trees)
          dummy <- dummyVars(" ~ .", data = df_tmp, sep = "_")
          df.disj <- data.frame(predict(dummy, newdata = df_tmp))
          mask.disj <- array(as.matrix(df.disj), c(dim(df.disj), dim(pred1.disj)[2]))
          mask.disj <- aperm(mask.disj, c(1, 3, 2))
          pred_draw.disj <- apply(mask.disj * pred1.disj, c(1, 2), sum)
          #######
          data.disj[v.na, dict_cat[[v]]] <- if (pmm.k) {
            pmm(
              xtrain = fit.disj$predictions,
              xtest = pred.disj,
              ytrain = data.disj[[v]][!v.na],
              k = pmm.k
            )
          } else {
            pred_draw.disj
          }
        }
        else {
          data.disj[v.na, v] <- if (pmm.k) {
            pmm(
              xtrain = fit$predictions,
              xtest = pred_draw,
              ytrain = data[[v]][!v.na],
              k = pmm.k
            )
          } else {
            pred_draw
          }
        }
      }
      data[v.na, v] <- if (pmm.k) {
        pmm(
          xtrain = fit$predictions,
          xtest = pred_draw,
          ytrain = data[[v]][!v.na],
          k = pmm.k
        )
      } else {
        pred_draw
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
  return(list(ximp = revert(converted, X = data), ximp.disj = data.disj))
}


MI_missRanger <- function(data, formula = . ~ ., pmm.k = 0L, maxiter = 10L, seed = NULL,
                          verbose = 1, returnOOB = FALSE, case.weights = NULL, col_cat = c(), num_mi = 5, ...) {
  imputations <- list()
  imputations.disj <- list()
  for (i in seq(num_mi)) {
    res <- missRanger_mod_draw(data, formula, pmm.k, maxiter, seed, verbose, returnOOB, case.weights, col_cat)
    imputations[[i]] <- res$ximp
    imputations.disj[[i]] <- res$ximp.disj
  }
  exist_cat <- !all(c(0, col_cat) == c(0))
  df_new_merge <- Reduce(function(dtf1, dtf2) {
    abind::abind(dtf1, dtf2, along = 3)
  }, imputations.disj)
  final.imp.disj <- data.frame(apply(df_new_merge, c(1, 2), mean))
  final.imp <- final.imp.disj
  if (exist_cat) {
    dict_cat <- dict_onehot(data, col_cat)
    names_cat <- names(dict_cat)
    for (name in names_cat) {
      final.imp[[name]] <- apply(final.imp.disj[dict_cat[[name]]], 1, which_max_cat, name, dict_cat)
      final.imp[[name]] <- unlist(final.imp[[name]])
      final.imp[[name]] <- factor(final.imp[[name]])
      final.imp <- final.imp[, !names(final.imp) %in% dict_cat[[name]]]
    }
  } else {
    final.imp <- final.imp.disj
  }
  return(list(ls_imputations = imputations, ls_imputations.disj = imputations.disj, ximp = final.imp, ximp.disj = final.imp.disj))
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
#' @examples
#' imputeUnivariate(c(NA, 0, 1, 0, 1))
#' imputeUnivariate(c("A", "A", NA))
#' imputeUnivariate(as.factor(c("A", "A", NA)))
#' head(imputeUnivariate(generateNA(iris)))
#' head(imputeUnivariate(generateNA(iris), v = "Species"))
#' head(imputeUnivariate(generateNA(iris), v = c("Species", "Petal.Length")))
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
