#' missRanger_mod_draw
#'
#' @description \code{missRanger_mod_draw} create one imputation result for multiple imputation with \code{missRanger} method.
#' Please find the detailed explanation of \code{missRanger} single imputation method in  the documentation of \code{missRanger} in 'missRanger' package.
#' In this document, only the differences will be explained.
#'
#' \code{missRanger} is an imputation method based on random forest. In \code{missRanger_mod_draw}, during the last iteration, for a certain prediction,
#'  instead of taking average of the prediction result from each tree of the random forest, we draw one result from the empirical
#'  distribution constructed by predictions of trees. The other steps of the imputation are identical as those of \code{missRanger}.
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
#' @param ... Arguments passed to \code{ranger()}. If the data set is large, better use less trees (e.g. \code{num.trees = 20}) and/or a low value of \code{sample.fraction}.
#' The following arguments are e.g. incompatible with \code{ranger}: \code{write.forest}, \code{probability}, \code{split.select.weights}, \code{dependent.variable.name}, and \code{classification}.
#' @param col_cat Indices of categorical columns.
#' @param num_mi Number of multiple imputations.
#' @export
#' @return \code{ls_ximp} List of imputed datasets for multiple imputation in \code{MI_missRanger}.
#' @return \code{ls_ximp.disj} List of disjunctive imputed datasets for multiple imputation in \code{MI_missRanger}.
missRanger_mod_draw <- function(data, formula = . ~ ., pmm.k = 0L, maxiter = 10L, seed = NULL,
                                verbose = 1, returnOOB = FALSE, case.weights = NULL, col_cat = c(), num_mi = 5, ...) {
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
    cat("Last iter and multiple imputation:\t", sep = "")
  }

  if (j == maxiter + 1) {
    # Last iteration
    dataLast <- data
    predErrorLast <- predError
  } else {
    # if crit, use dataLast2 to redo the iteration
    data <- dataLast2
  }

  #### Multiple Imputation ####
  data_keep <- data
  ls_ximp <- list()
  ls_ximp.disj <- list()
  for (i in seq(num_mi)) {
    data <- data_keep
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
        # pred <- predict(fit, data[v.na, completed, drop = FALSE])$predictions
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
            if (pmm.k) {
              data.disj[v.na, dict_cat[[v]]] <- pmm(
                xtrain = fit.disj$predictions,
                xtest = pred_draw.disj,
                ytrain = data.disj[[v]][!v.na],
                k = pmm.k
              )
            } else if (ncol(data.disj[v.na, dict_cat[[v]]]) == ncol(pred_draw.disj)) {
              data.disj[v.na, dict_cat[[v]]] <- pred_draw.disj
            } else { # When there are dropped unused levels
              colnames(pred_draw.disj) <- paste0(v, "_", colnames(pred_draw.disj))
              data.disj[v.na, colnames(pred_draw.disj)] <- pred_draw.disj
            }
          } else {
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
    ls_ximp[[i]] <- revert(converted, X = data)
    ls_ximp.disj[[i]] <- data.disj
  }

  return(list(ls_ximp = ls_ximp, ls_ximp.disj = ls_ximp.disj))
}


#' MI_missRanger
#'
#' @description \code{MI_missRanger} is a function of multiple imputation with \code{missRanger} method.
#'
#'  In \code{missRanger_mod_draw}, for a certain prediction,
#'  instead of taking average of the prediction result from each tree of the random forest, during the last iteration, we draw one result from the empirical
#'  distribution constructed by predictions of trees. The other steps of the imputation are identical as those of \code{missRanger}
#'  from 'missRanger' package.
#'
#'  \code{MI_missRanger} takes all the imputation results from \code{missRanger_mod_draw} and combine them with Rubin's Rule
#'  to generate the final imputed data set.
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
#' @param ... Arguments passed to \code{ranger()}. If the data set is large, better use less trees (e.g. \code{num.trees = 20}) and/or a low value of \code{sample.fraction}.
#' The following arguments are e.g. incompatible with \code{ranger}: \code{write.forest}, \code{probability}, \code{split.select.weights}, \code{dependent.variable.name}, and \code{classification}.
#' @param col_cat Indices of categorical columns
#' @param num_mi Number of multiple imputation
#'
#' @export
#' @return \code{ximp} Final imputed dataset.
#' @return \code{ximp.disj} Final disjunctive imputed dataset.
#' @return \code{ls_imputations} List of imputed dataset from multiple imputation.
#' @return \code{ls_imputations.disj} List of disjunctive imputed dataset from multiple imputation.
MI_missRanger <- function(data, formula = . ~ ., pmm.k = 0L, maxiter = 10L, seed = NULL,
                          verbose = 1, returnOOB = FALSE, case.weights = NULL, col_cat = c(), num_mi = 5, ...) {
  is_ranger_package_installed()
  is_abind_package_installed()
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
    # Create the dictionary for disjunctive variable names
    dummy <- dummyVars(" ~ .", data = data, sep = "_")
    col_names.disj.new <- colnames(data.frame(predict(dummy, newdata = data)))
    dict_name_disj <- list()
    for (i in seq(length(col_names.disj.new))) {
      dict_name_disj[[col_names.disj.new[i]]] <- col_names.disj[i]
    }
    # Create dict_cat with categroical columns
    dict_cat <- dict_onehot(data, col_cat)
  }
  imputations <- list()
  imputations.disj <- list()
  res <- missRanger_mod_draw(data, formula, pmm.k, maxiter, seed, verbose, returnOOB, case.weights, col_cat, num_mi)
  imputations <- res$ls_ximp
  imputations.disj <- res$ls_ximp.disj
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
    for (col in name_cat) {
      levels(final.imp[[col]]) <- dict_lev[[col]]
    }
    colnames(final.imp.disj) <- dict_name_disj[colnames(final.imp.disj)]
    for (i in seq(length(imputations))) {
      for (col in name_cat) {
        levels(imputations[[i]][[col]]) <- dict_lev[[col]]
      }
      colnames(imputations.disj[[i]]) <- dict_name_disj[colnames(imputations.disj[[i]])]
    }
  } else {
    final.imp <- final.imp.disj
  }
  return(list(ls_imputations = imputations, ls_imputations.disj = imputations.disj, ximp = final.imp, ximp.disj = final.imp.disj))
}
