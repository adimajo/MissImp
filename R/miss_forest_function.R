##
## MissForest - nonparametric missing value imputation for mixed-type data
##
## This R script contains the actual missForest function.
##
## Author: Daniel Stekhoven, stekhoven@stat.math.ethz.ch
##
## Acknowledgement: Steve Weston for input regarding parallel execution (2012)
##############################################################################
#' missForest_mod: modified missForest with onehot probability
#'
#' @description \code{missForest_mod} is a modified version of the function \code{missForest} by Daniel Stekhoven.
#' Please find the detailed documentation of \code{missForest} in the missForest package. Only the modifications are explained on this page.
#' The original \code{missForest} function returns the final imputation result after convergence or \code{maxiter} iterations.
#' The results of categorical columns are returned in form of vector. In \code{missForest_mod} function, during the last iteration,
#' not only the final result, but also the onehot probability for each category is returned.
#' @param xmis data matrix with missing values.
#' @param maxiter stop after how many iterations (default = 10).
#' @param ntree how many trees are grown in the forest (default = 100).
#' @param col_cat index of categorical columns.
#' @param variablewise (boolean) return OOB errors for each variable separately.
#' @param decreasing (boolean) if TRUE the columns are sorted with decreasing amount of missing values.
#' @param verbose (boolean) if TRUE then missForest returns error estimates, runtime and if available true error during iterations.
#' @param mtry how many variables should be tried randomly at each node.
#' @param replace (boolean) if TRUE bootstrap sampling (with replacements) is performed, else subsampling (without replacements).
#' @param classwt list of priors of the classes in the categorical variables.
#' @param cutoff list of class cutoffs for each categorical variable.
#' @param strata list of (factor) variables used for stratified sampling.
#' @param sampsize list of size(s) of sample to draw
#' @param nodesize minimum size of terminal nodes, vector of length 2, with
#' number for continuous variables in the first entry and
#' number for categorical variables in the second entry.
#' @param maxnodes maximum number of terminal nodes for individual trees
#' @param xtrue complete data matrix
#' @param parallelize TODO
#' @export
#' @import foreach
#' @return \code{ximp} imputed data matrix of same type as 'xmis'.
#' @return \code{ximp.disj} imputed data matrix of same type as 'xmis' for the numeric columns.
#'  For the categorical columns, the prediction of probability for each category is shown in form of onehot vector.
#' @return \code{OOBerror} estimated OOB imputation error. For the set of continuous variables in 'xmis' the NRMSE and
#' for the set of categorical variables the proportion of falsely classified entries is returned.
#' See Details for the exact definition of these error measures.
#' If 'variablewise' is set to 'TRUE' then this will be a vector of length 'p' where 'p' is the number of variables and the entries will be the OOB error for each variable separately.
#' @return \code{error} true imputation error. This is only available if 'xtrue' was supplied.
#' The error measures are the same as for 'OOBerror'.
#'
missForest_mod <- function(xmis, maxiter = 10, ntree = 100, variablewise = FALSE,
                           decreasing = FALSE, verbose = FALSE,
                           mtry = floor(sqrt(ncol(xmis))), replace = TRUE,
                           classwt = NULL, cutoff = NULL, strata = NULL,
                           sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                           xtrue = NA, parallelize = c("no", "variables", "forests"),
                           col_cat = c()) { ## ----------------------------------------------------------------------
  ## Arguments:
  ## xmis         = data matrix with missing values
  ## maxiter      = stop after how many iterations (default = 10)
  ## ntree        = how many trees are grown in the forest (default = 100)
  ## variablewise = (boolean) return OOB errors for each variable separately
  ## decreasing   = (boolean) if TRUE the columns are sorted with decreasing
  ##                amount of missing values
  ## verbose      = (boolean) if TRUE then missForest returns error estimates,
  ##                runtime and if available true error during iterations
  ## mtry         = how many variables should be tried randomly at each node
  ## replace      = (boolean) if TRUE bootstrap sampling (with replacements)
  ##                is performed, else subsampling (without replacements)
  ## classwt      = list of priors of the classes in the categorical variables
  ## cutoff       = list of class cutoffs for each categorical variable
  ## strata       = list of (factor) variables used for stratified sampling
  ## sampsize     = list of size(s) of sample to draw
  ## nodesize     = minimum size of terminal nodes, vector of length 2, with
  ##                number for continuous variables in the first entry and
  ##                number for categorical variables in the second entry
  ## maxnodes     = maximum number of terminal nodes for individual trees
  ## xtrue        = complete data matrix
  ## [add]col_cat      = index of categorical columns
  ## ----------------------------------------------------------------------
  ##

  ## Add: Create dict_cat with categroical columns
  exist_cat <- !all(c(0, col_cat) == c(0))
  if (exist_cat) {
    name_cat <- colnames(xmis)[col_cat]
    # Deal with the problem that nlevels(df[[col]]) > length(unique(df[[col]]))
    for (col in name_cat) {
      xmis[[col]] <- factor(as.character(xmis[[col]]))
    }
    # remember the levels for each categorical column
    dict_lev <- dict_level(xmis, col_cat)
    # preserve colnames for ximp.disj
    dummy <- dummyVars(" ~ .", data = xmis, sep = "_")
    col_names.disj <- colnames(data.frame(predict(dummy, newdata = xmis)))
    # represent the factor columns with their ordinal levels
    xmis <- factor_ordinal_encode(xmis, col_cat)
    # Create dict_cat with categroical columns
    dict_cat <- dict_onehot(xmis, col_cat)
  }

  ## Add: Last iteration will be used to predict the onehot probability for the categorical columns
  maxiter <- maxiter - 1

  ## stop in case of wrong inputs passed to randomForest

  n <- nrow(xmis)
  p <- ncol(xmis)
  ls_colname <- colnames(xmis)
  if (!is.null(classwt)) {
    stopifnot(length(classwt) == p, typeof(classwt) == "list")
  }
  if (!is.null(cutoff)) {
    stopifnot(length(cutoff) == p, typeof(cutoff) == "list")
  }
  if (!is.null(strata)) {
    stopifnot(length(strata) == p, typeof(strata) == "list")
  }
  if (!is.null(nodesize)) {
    stopifnot(length(nodesize) == 2)
  }

  ## remove completely missing variables
  if (any(apply(is.na(xmis), 2, sum) == n)) {
    indCmis <- which(apply(is.na(xmis), 2, sum) == n)
    xmis <- xmis[, -indCmis]
    p <- ncol(xmis)
    cat(
      "  removed variable(s)", indCmis,
      "due to the missingness of all entries\n"
    )
  }

  ## return feedback on parallelization setup
  parallelize <- match.arg(parallelize)
  if (parallelize %in% c("variables", "forests")) {
    if (getDoParWorkers() == 1) {
      stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
    } else if (verbose) {
      if (parallelize == "variables") {
        cat("  parallelizing over the variables of the input data matrix 'xmis'\n")
      } else {
        cat("  parallelizing computation of the random forest model objects\n")
      }
    }
    if (getDoParWorkers() > p) {
      stop("The number of parallel cores should not exceed the number of variables (p=", p, ")")
    }
  }

  ## perform initial S.W.A.G. on xmis (mean imputation)
  ximp <- xmis
  xAttrib <- lapply(xmis, attributes)
  varType <- character(p)
  for (t.co in 1:p) {
    if (is.null(xAttrib[[t.co]])) {
      varType[t.co] <- "numeric"
      ximp[is.na(xmis[, t.co]), t.co] <- mean(xmis[, t.co], na.rm = TRUE)
    } else {
      varType[t.co] <- "factor"
      ## take the level which is more 'likely' (majority vote)
      max.level <- max(table(ximp[, t.co]))
      ## if there are several classes which are major, sample one at random
      class.assign <- sample(names(which(max.level == summary(ximp[, t.co]))), 1)
      ## it shouldn't be the NA class
      if (class.assign != "NA's") {
        ximp[is.na(xmis[, t.co]), t.co] <- class.assign
      } else {
        while (class.assign == "NA's") {
          class.assign <- sample(names(which(max.level ==
            summary(ximp[, t.co]))), 1)
        }
        ximp[is.na(xmis[, t.co]), t.co] <- class.assign
      }
    }
  }

  ## extract missingness pattern
  NAloc <- is.na(xmis) # where are missings
  noNAvar <- apply(NAloc, 2, sum) # how many are missing in the vars
  sort.j <- order(noNAvar) # indices of increasing amount of NA in vars
  if (decreasing) {
    sort.j <- rev(sort.j)
  }
  sort.noNAvar <- noNAvar[sort.j]

  ## compute a list of column indices for variable parallelization
  nzsort.j <- sort.j[sort.noNAvar > 0]
  if (parallelize == "variables") {
    "%cols%" <- get("%dopar%")
    idxList <- as.list(isplitVector(nzsort.j, chunkSize = getDoParWorkers()))
  }
  #   else {
  #     ## force column loop to be sequential
  #     '%cols%' <- get('%do%')
  #     idxList <- nzsort.j
  #   }

  ## output
  Ximp <- vector("list", maxiter)

  ## initialize parameters of interest
  iter <- 0
  k <- length(unique(varType))
  convNew <- rep(0, k)
  convOld <- rep(Inf, k)
  OOBerror <- numeric(p)
  names(OOBerror) <- varType

  ## setup convergence variables w.r.t. variable types
  if (k == 1) {
    if (unique(varType) == "numeric") {
      names(convNew) <- c("numeric")
    } else {
      names(convNew) <- c("factor")
    }
    convergence <- c()
    OOBerr <- numeric(1)
  } else {
    names(convNew) <- c("numeric", "factor")
    convergence <- matrix(NA, ncol = 2)
    OOBerr <- numeric(2)
  }

  ## function to yield the stopping criterion in the following 'while' loop
  stopCriterion <- function(varType, convNew, convOld, iter, maxiter) {
    k <- length(unique(varType))
    if (k == 1) {
      (convNew < convOld) & (iter < maxiter)
    } else {
      ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)
    }
  }

  ## iterate missForest
  while (stopCriterion(varType, convNew, convOld, iter, maxiter)) {
    if (iter != 0) {
      convOld <- convNew
      OOBerrOld <- OOBerr
    }
    cat("  missForest iteration", iter + 1, "in progress...")
    t.start <- proc.time()
    ximp.old <- ximp

    if (parallelize == "variables") {
      for (idx in idxList) {
        results <- foreach(varInd = idx, .packages = "randomForest") %cols% {
          obsi <- !NAloc[, varInd] # which i's are observed
          misi <- NAloc[, varInd] # which i's are missing
          obsY <- ximp[obsi, varInd] # training response
          obsX <- ximp[obsi, seq(1, p)[-varInd]] # training variables
          misX <- ximp[misi, seq(1, p)[-varInd]] # prediction variables
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            RF <- randomForest::randomForest(
              x = obsX,
              y = obsY,
              ntree = ntree,
              mtry = mtry,
              replace = replace,
              sampsize = if (!is.null(sampsize)) {
                sampsize[[varInd]]
              } else
              if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
              nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
              maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
            )
            ## record out-of-bag error
            oerr <- RF$mse[ntree]
            #           }
            ## predict missing values in column varInd
            misY <- stats::predict(RF, misX)
          } else { # if Y is categorical
            obsY <- factor(obsY) ## remove empty classes
            summarY <- summary(obsY)
            if (length(summarY) == 1) { ## if there is only one level left
              oerr <- 0
              misY <- factor(rep(names(summarY), length(misi)))
            } else {
              RF <- randomForest::randomForest(
                x = obsX,
                y = obsY,
                ntree = ntree,
                mtry = mtry,
                replace = replace,
                classwt = if (!is.null(classwt)) {
                  classwt[[varInd]]
                } else {
                  rep(1, nlevels(obsY))
                },
                cutoff = if (!is.null(cutoff)) {
                  cutoff[[varInd]]
                } else {
                  rep(1 / nlevels(obsY), nlevels(obsY))
                },
                strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                sampsize = if (!is.null(sampsize)) {
                  sampsize[[varInd]]
                } else
                if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
              )
              ## record out-of-bag error
              oerr <- RF$err.rate[[ntree, 1]]
              #             }
              ## predict missing values in column varInd
              misY <- stats::predict(RF, misX)
            }
          }
          list(varInd = varInd, misY = misY, oerr = oerr)
        }
        ## update the master copy of the data
        for (res in results) {
          misi <- NAloc[, res$varInd]
          ximp[misi, res$varInd] <- res$misY
          OOBerror[res$varInd] <- res$oerr
        }
      }
    } else { # if parallelize != "variables"
      for (s in 1:p) {
        varInd <- sort.j[s]
        if (noNAvar[[varInd]] != 0) {
          obsi <- !NAloc[, varInd]
          misi <- NAloc[, varInd]
          obsY <- ximp[obsi, varInd]
          obsX <- ximp[obsi, seq(1, p)[-varInd]]
          misX <- ximp[misi, seq(1, p)[-varInd]]
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            if (parallelize == "forests") {
              xntree <- NULL
              RF <- foreach(
                xntree = idiv(ntree, chunks = getDoParWorkers()),
                .combine = "combine", .multicombine = TRUE,
                .packages = "randomForest"
              ) %dopar% {
                randomForest::randomForest(
                  x = obsX,
                  y = obsY,
                  ntree = xntree,
                  mtry = mtry,
                  replace = replace,
                  sampsize = if (!is.null(sampsize)) {
                    sampsize[[varInd]]
                  } else
                  if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                  nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                  maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
                )
              }
              ## record out-of-bag error
              OOBerror[varInd] <- mean((predict(RF) - RF$y)^2, na.rm = TRUE)
              #               OOBerror[varInd] <- RF$mse[ntree]
            } else {
              RF <- randomForest::randomForest(
                x = obsX,
                y = obsY,
                ntree = ntree,
                mtry = mtry,
                replace = replace,
                sampsize = if (!is.null(sampsize)) {
                  sampsize[[varInd]]
                } else
                if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
              )
              ## record out-of-bag error
              OOBerror[varInd] <- RF$mse[ntree]
            }
            misY <- stats::predict(RF, misX)
          } else {
            obsY <- factor(obsY)
            summarY <- summary(obsY)
            if (length(summarY) == 1) {
              misY <- factor(rep(names(summarY), sum(misi)))
            } else {
              if (parallelize == "forests") {
                RF <- foreach(
                  xntree = idiv(ntree, chunks = getDoParWorkers()),
                  .combine = "combine", .multicombine = TRUE,
                  .packages = "randomForest"
                ) %dopar% {
                  randomForest::randomForest(
                    x = obsX,
                    y = obsY,
                    ntree = xntree,
                    mtry = mtry,
                    replace = replace,
                    classwt = if (!is.null(classwt)) {
                      classwt[[varInd]]
                    } else {
                      rep(1, nlevels(obsY))
                    },
                    cutoff = if (!is.null(cutoff)) {
                      cutoff[[varInd]]
                    } else {
                      rep(1 / nlevels(obsY), nlevels(obsY))
                    },
                    strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                    sampsize = if (!is.null(sampsize)) {
                      sampsize[[varInd]]
                    } else
                    if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                    nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                    maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
                  )
                }
                ## record out-of-bag error
                ne <- as.integer(predict(RF)) != as.integer(RF$y)
                ne <- ne[!is.na(ne)]
                OOBerror[varInd] <- sum(ne) / length(ne)
              } else {
                RF <- randomForest::randomForest(
                  x = obsX,
                  y = obsY,
                  ntree = ntree,
                  mtry = mtry,
                  replace = replace,
                  classwt = if (!is.null(classwt)) {
                    classwt[[varInd]]
                  } else {
                    rep(1, nlevels(obsY))
                  },
                  cutoff = if (!is.null(cutoff)) {
                    cutoff[[varInd]]
                  } else {
                    rep(1 / nlevels(obsY), nlevels(obsY))
                  },
                  strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                  sampsize = if (!is.null(sampsize)) {
                    sampsize[[varInd]]
                  } else
                  if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                  nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                  maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
                )
                ## record out-of-bag error
                OOBerror[varInd] <- RF$err.rate[[ntree, 1]]
              }
              ## predict missing parts of Y
              misY <- stats::predict(RF, misX)
            }
          }
          ximp[misi, varInd] <- misY
        }
      }
    }
    cat("done!\n")

    iter <- iter + 1
    Ximp[[iter]] <- ximp

    t.co2 <- 1
    ## check the difference between iteration steps
    for (t.type in names(convNew)) {
      t.ind <- which(varType == t.type)
      if (t.type == "numeric") {
        convNew[t.co2] <- sum((ximp[, t.ind] - ximp.old[, t.ind])^2) / sum(ximp[, t.ind]^2)
      } else {
        dist <- sum(as.character(as.matrix(ximp[, t.ind])) != as.character(as.matrix(ximp.old[, t.ind])))
        convNew[t.co2] <- dist / (n * sum(varType == "factor"))
      }
      t.co2 <- t.co2 + 1
    }

    ## compute estimated imputation error
    if (!variablewise) {
      NRMSE <- sqrt(mean(OOBerror[varType == "numeric"]) /
        var(as.vector(as.matrix(xmis[, varType == "numeric"])),
          na.rm = TRUE
        ))
      PFC <- mean(OOBerror[varType == "factor"])
      if (k == 1) {
        if (unique(varType) == "numeric") {
          OOBerr <- NRMSE
          names(OOBerr) <- "NRMSE"
        } else {
          OOBerr <- PFC
          names(OOBerr) <- "PFC"
        }
      } else {
        OOBerr <- c(NRMSE, PFC)
        names(OOBerr) <- c("NRMSE", "PFC")
      }
    } else {
      OOBerr <- OOBerror
      names(OOBerr)[varType == "numeric"] <- "MSE"
      names(OOBerr)[varType == "factor"] <- "PFC"
    }

    if (any(!is.na(xtrue))) {
      err <- suppressWarnings(mixError(ximp, xmis, xtrue))
    }

    ## return status output, if desired
    if (verbose) {
      delta.start <- proc.time() - t.start
      if (any(!is.na(xtrue))) {
        cat("    error(s):", err, "\n")
      }
      cat("    estimated error(s):", OOBerr, "\n")
      cat("    difference(s):", convNew, "\n")
      cat("    time:", delta.start[3], "seconds\n\n")
    }
  } # end while((convNew<convOld)&(iter<maxiter)){


  # Add: Calculate the one-hot probability from randomforest
  if (iter == maxiter) {
    ximp <- Ximp[[iter]]
  } else {
    ximp <- Ximp[[iter - 1]]
  }

  dummy <- dummyVars(" ~ .", data = ximp, sep = "_")
  ximp.disj <- data.frame(predict(dummy, newdata = ximp))

  if (iter != 0) {
    convOld <- convNew
    OOBerrOld <- OOBerr
  }
  cat("  missForest last iteration in progress...")
  t.start <- proc.time()
  ximp.old <- ximp

  if (parallelize == "variables") {
    for (idx in idxList) {
      results <- foreach(varInd = idx, .packages = "randomForest") %cols% {
        obsi <- !NAloc[, varInd] # which i's are observed
        misi <- NAloc[, varInd] # which i's are missing
        obsY <- ximp[obsi, varInd] # training response
        obsX <- ximp[obsi, seq(1, p)[-varInd]] # training variables
        misX <- ximp[misi, seq(1, p)[-varInd]] # prediction variables
        typeY <- varType[varInd]
        if (typeY == "numeric") {
          RF <- randomForest::randomForest(
            x = obsX,
            y = obsY,
            ntree = ntree,
            mtry = mtry,
            replace = replace,
            sampsize = if (!is.null(sampsize)) {
              sampsize[[varInd]]
            } else
            if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
            nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
            maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
          )
          ## record out-of-bag error
          oerr <- RF$mse[ntree]
          #           }
          ## predict missing values in column varInd
          misY <- stats::predict(RF, misX)
        } else { # if Y is categorical
          obsY <- factor(obsY) ## remove empty classes
          summarY <- summary(obsY)
          if (length(summarY) == 1) { ## if there is only one level left
            oerr <- 0
            misY <- factor(rep(names(summarY), length(misi)))
          } else {
            RF <- randomForest::randomForest(
              x = obsX,
              y = obsY,
              ntree = ntree,
              mtry = mtry,
              replace = replace,
              classwt = if (!is.null(classwt)) {
                classwt[[varInd]]
              } else {
                rep(1, nlevels(obsY))
              },
              cutoff = if (!is.null(cutoff)) {
                cutoff[[varInd]]
              } else {
                rep(1 / nlevels(obsY), nlevels(obsY))
              },
              strata = if (!is.null(strata)) strata[[varInd]] else obsY,
              sampsize = if (!is.null(sampsize)) {
                sampsize[[varInd]]
              } else
              if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
              nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
              maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
            )
            ## record out-of-bag error
            oerr <- RF$err.rate[[ntree, 1]]
            #             }
            ## predict missing values in column varInd
            misY <- stats::predict(RF, misX)
            misY.disj <- predict.randomForest(RF, misX, type = "prob")
          }
        }
        list(varInd = varInd, misY = misY, oerr = oerr, misY.disj = misY.disj)
      }
      ## update the master copy of the data
      for (res in results) {
        misi <- NAloc[, res$varInd]
        ximp[misi, res$varInd] <- res$misY
        OOBerror[res$varInd] <- res$oerr
        if (exist_cat) {
          ximp.disj[misi, dict_cat[[ls_colname[res$varInd]]]] <- res$misY.disj
        }
      }
    }
  } else { # if parallelize != "variables"
    for (s in 1:p) {
      varInd <- sort.j[s]
      if (noNAvar[[varInd]] != 0) {
        obsi <- !NAloc[, varInd]
        misi <- NAloc[, varInd]
        obsY <- ximp[obsi, varInd]
        obsX <- ximp[obsi, seq(1, p)[-varInd]]
        misX <- ximp[misi, seq(1, p)[-varInd]]
        typeY <- varType[varInd]
        if (typeY == "numeric") {
          if (parallelize == "forests") {
            xntree <- NULL
            RF <- foreach(
              xntree = idiv(ntree, chunks = getDoParWorkers()),
              .combine = "combine", .multicombine = TRUE,
              .packages = "randomForest"
            ) %dopar% {
              randomForest::randomForest(
                x = obsX,
                y = obsY,
                ntree = xntree,
                mtry = mtry,
                replace = replace,
                sampsize = if (!is.null(sampsize)) {
                  sampsize[[varInd]]
                } else
                if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
              )
            }
            ## record out-of-bag error
            OOBerror[varInd] <- mean((predict(RF) - RF$y)^2, na.rm = TRUE)
            #               OOBerror[varInd] <- RF$mse[ntree]
          } else {
            RF <- randomForest::randomForest(
              x = obsX,
              y = obsY,
              ntree = ntree,
              mtry = mtry,
              replace = replace,
              sampsize = if (!is.null(sampsize)) {
                sampsize[[varInd]]
              } else
              if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
              nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
              maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
            )
            ## record out-of-bag error
            OOBerror[varInd] <- RF$mse[ntree]
          }
          misY <- predict(RF, misX)
        } else {
          obsY <- factor(obsY)
          summarY <- summary(obsY)
          if (length(summarY) == 1) {
            misY <- factor(rep(names(summarY), sum(misi)))
          } else {
            if (parallelize == "forests") {
              RF <- foreach(
                xntree = idiv(ntree, chunks = getDoParWorkers()),
                .combine = "combine", .multicombine = TRUE,
                .packages = "randomForest"
              ) %dopar% {
                randomForest::randomForest(
                  x = obsX,
                  y = obsY,
                  ntree = xntree,
                  mtry = mtry,
                  replace = replace,
                  classwt = if (!is.null(classwt)) {
                    classwt[[varInd]]
                  } else {
                    rep(1, nlevels(obsY))
                  },
                  cutoff = if (!is.null(cutoff)) {
                    cutoff[[varInd]]
                  } else {
                    rep(1 / nlevels(obsY), nlevels(obsY))
                  },
                  strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                  sampsize = if (!is.null(sampsize)) {
                    sampsize[[varInd]]
                  } else
                  if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                  nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                  maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
                )
              }
              ## record out-of-bag error
              ne <- as.integer(predict(RF)) != as.integer(RF$y)
              ne <- ne[!is.na(ne)]
              OOBerror[varInd] <- sum(ne) / length(ne)
            } else {
              RF <- randomForest::randomForest(
                x = obsX,
                y = obsY,
                ntree = ntree,
                mtry = mtry,
                replace = replace,
                classwt = if (!is.null(classwt)) {
                  classwt[[varInd]]
                } else {
                  rep(1, nlevels(obsY))
                },
                cutoff = if (!is.null(cutoff)) {
                  cutoff[[varInd]]
                } else {
                  rep(1 / nlevels(obsY), nlevels(obsY))
                },
                strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                sampsize = if (!is.null(sampsize)) {
                  sampsize[[varInd]]
                } else
                if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                maxnodes = if (!is.null(maxnodes)) maxnodes else NULL
              )
              ## record out-of-bag error
              OOBerror[varInd] <- RF$err.rate[[ntree, 1]]
            }
            ## predict missing parts of Y
            misY <- stats::predict(RF, misX)
            if (exist_cat) {
              misY.disj <- predict.randomForest(RF, misX, type = "prob")
              if (ncol(ximp.disj[misi, dict_cat[[ls_colname[varInd]]]]) == ncol(misY.disj)) {
                ximp.disj[misi, dict_cat[[ls_colname[varInd]]]] <- misY.disj
              } else { # When there are dropped unused levels
                colnames(misY.disj) <- paste0(v, "_", colnames(misY.disj))
                ximp.disj[misi, colnames(misY.disj)] <- misY.disj
              }
            }
          }
        }
        ximp[misi, varInd] <- misY
      }
    }
  }
  cat("done!\n")

  iter <- iter + 1
  Ximp_final <- ximp

  t.co2 <- 1

  ## compute estimated imputation error
  if (!variablewise) {
    NRMSE <- sqrt(mean(OOBerror[varType == "numeric"]) /
      var(as.vector(as.matrix(xmis[, varType == "numeric"])),
        na.rm = TRUE
      ))
    PFC <- mean(OOBerror[varType == "factor"])
    if (k == 1) {
      if (unique(varType) == "numeric") {
        OOBerr <- NRMSE
        names(OOBerr) <- "NRMSE"
      } else {
        OOBerr <- PFC
        names(OOBerr) <- "PFC"
      }
    } else {
      OOBerr <- c(NRMSE, PFC)
      names(OOBerr) <- c("NRMSE", "PFC")
    }
  } else {
    OOBerr <- OOBerror
    names(OOBerr)[varType == "numeric"] <- "MSE"
    names(OOBerr)[varType == "factor"] <- "PFC"
  }

  if (any(!is.na(xtrue))) {
    err <- suppressWarnings(mixError(ximp, xmis, xtrue))
  }

  ## return status output, if desired
  if (verbose) {
    delta.start <- proc.time() - t.start
    if (any(!is.na(xtrue))) {
      cat("    error(s):", err, "\n")
    }
    cat("    estimated error(s):", OOBerr, "\n")
    cat("    difference(s):", convNew, "\n")
    cat("    time:", delta.start[3], "seconds\n\n")
  }



  ## Add: return the original levels
  if (exist_cat) {
    for (col in name_cat) {
      levels(Ximp_final[[col]]) <- dict_lev[[col]]
    }
    colnames(ximp.disj) <- col_names.disj
  }

  ## Add: if there is no categorical columns, ximp.disj = ximp
  if (!exist_cat) {
    ximp.disj <- Ximp_final
  }

  ## produce output w.r.t. stopping rule
  if (iter == maxiter + 1) {
    if (any(is.na(xtrue))) {
      out <- list(ximp = Ximp_final, OOBerror = OOBerr, ximp.disj = ximp.disj)
    } else {
      out <- list(ximp = Ximp_final, OOBerror = OOBerr, ximp.disj = ximp.disj, error = err)
    }
  } else {
    if (any(is.na(xtrue))) {
      out <- list(ximp = Ximp_final, ximp.disj = ximp.disj, OOBerror = OOBerrOld)
    } else {
      out <- list(
        ximp = Ximp_final, ximp.disj = ximp.disj, OOBerror = OOBerrOld,
        error = suppressWarnings(mixError(Ximp_final, xmis, xtrue))
      )
    }
  }
  class(out) <- "missForest"
  return(out)
}


predict.randomForest <-
  function(object, newdata, type = "response", norm.votes = TRUE,
           predict.all = FALSE, proximity = FALSE, nodes = FALSE, cutoff, ...) {
    if (!inherits(object, "randomForest")) {
      stop("object not of class randomForest")
    }
    if (is.null(object$forest)) stop("No forest component in the object")
    out.type <- charmatch(
      tolower(type),
      c("response", "prob", "vote", "class")
    )
    if (is.na(out.type)) {
      stop("type must be one of 'response', 'prob', 'vote'")
    }
    if (out.type == 4) out.type <- 1
    if (out.type != 1 && object$type == "regression") {
      stop("'prob' or 'vote' not meaningful for regression")
    }
    if (out.type == 2) {
      norm.votes <- TRUE
    }
    if (missing(newdata)) {
      p <- if (!is.null(object$na.action)) {
        napredict(object$na.action, object$predicted)
      } else {
        object$predicted
      }
      if (object$type == "regression") {
        return(p)
      }
      if (proximity & is.null(object$proximity)) {
        warning("cannot return proximity without new data if random forest object does not already have proximity")
      }
      if (out.type == 1) {
        if (proximity) {
          return(list(
            pred = p,
            proximity = object$proximity
          ))
        } else {
          return(p)
        }
      }
      v <- object$votes
      if (!is.null(object$na.action)) v <- napredict(object$na.action, v)
      if (norm.votes) {
        t1 <- t(apply(v, 1, function(x) {
          x / sum(x)
        }))
        class(t1) <- c(class(t1), "votes")
        if (proximity) {
          return(list(pred = t1, proximity = object$proximity))
        } else {
          return(t1)
        }
      } else {
        if (proximity) {
          return(list(pred = v, proximity = object$proximity))
        } else {
          return(v)
        }
      }
    }
    if (missing(cutoff)) {
      cutoff <- object$forest$cutoff
    } else {
      if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 0) ||
        length(cutoff) != length(object$classes)) {
        stop("Incorrect cutoff specified.")
      }
      if (!is.null(names(cutoff))) {
        if (!all(names(cutoff) %in% object$classes)) {
          stop("Wrong name(s) for cutoff")
        }
        cutoff <- cutoff[object$classes]
      }
    }

    if (object$type == "unsupervised") {
      stop("Can't predict unsupervised forest.")
    }

    if (inherits(object, "randomForest.formula")) {
      newdata <- as.data.frame(newdata)
      rn <- row.names(newdata)
      Terms <- delete.response(object$terms)
      x <- model.frame(Terms, newdata, na.action = stats::na.omit)
      keep <- match(row.names(x), rn)
    } else {
      if (is.null(dim(newdata))) {
        dim(newdata) <- c(1, length(newdata))
      }
      x <- newdata
      if (nrow(x) == 0) {
        stop("newdata has 0 rows")
      }
      if (any(is.na(x))) {
        stop("missing values in newdata")
      }
      keep <- 1:nrow(x)
      rn <- rownames(x)
      if (is.null(rn)) rn <- keep
    }
    vname <- if (is.null(dim(object$importance))) {
      names(object$importance)
    } else {
      rownames(object$importance)
    }
    if (is.null(colnames(x))) {
      if (ncol(x) != length(vname)) {
        stop("number of variables in newdata does not match that in the training data")
      }
    } else {
      if (any(!vname %in% colnames(x))) {
        stop("variables in the training data missing in newdata")
      }
      x <- x[, vname, drop = FALSE]
    }
    if (is.data.frame(x)) {
      isFactor <- function(x) is.factor(x) & !is.ordered(x)
      xfactor <- which(sapply(x, isFactor))
      if (length(xfactor) > 0 && "xlevels" %in% names(object$forest)) {
        for (i in xfactor) {
          if (any(!levels(x[[i]]) %in% object$forest$xlevels[[i]])) {
            stop("New factor levels not present in the training data")
          }
          x[[i]] <-
            factor(x[[i]],
              levels = levels(x[[i]])[match(levels(x[[i]]), object$forest$xlevels[[i]])]
            )
        }
      }
      cat.new <- sapply(x, function(x) {
        if (is.factor(x) && !is.ordered(x)) {
          length(levels(x))
        } else {
          1
        }
      })
      if (!all(object$forest$ncat == cat.new)) {
        stop("Type of predictors in new data do not match that of the training data.")
      }
    }
    mdim <- ncol(x)
    ntest <- nrow(x)
    ntree <- object$forest$ntree
    maxcat <- max(object$forest$ncat)
    nclass <- object$forest$nclass
    nrnodes <- object$forest$nrnodes
    ## get rid of warning:
    op <- options(warn = -1)
    on.exit(options(op))
    x <- t(data.matrix(x))

    if (predict.all) {
      treepred <- if (object$type == "regression") {
        matrix(double(ntest * ntree), ncol = ntree)
      } else {
        matrix(integer(ntest * ntree), ncol = ntree)
      }
    } else {
      treepred <- numeric(ntest)
    }
    proxmatrix <- if (proximity) matrix(0, ntest, ntest) else numeric(1)
    nodexts <- if (nodes) integer(ntest * ntree) else integer(ntest)

    if (object$type == "regression") {
      if (!is.null(object$forest$treemap)) {
        object$forest$leftDaughter <-
          object$forest$treemap[, 1, , drop = FALSE]
        object$forest$rightDaughter <-
          object$forest$treemap[, 2, , drop = FALSE]
        object$forest$treemap <- NULL
      }

      keepIndex <- "ypred"
      if (predict.all) keepIndex <- c(keepIndex, "treepred")
      if (proximity) keepIndex <- c(keepIndex, "proximity")
      if (nodes) keepIndex <- c(keepIndex, "nodexts")
      ## Ensure storage mode is what is expected in C.
      if (!is.integer(object$forest$leftDaughter)) {
        storage.mode(object$forest$leftDaughter) <- "integer"
      }
      if (!is.integer(object$forest$rightDaughter)) {
        storage.mode(object$forest$rightDaughter) <- "integer"
      }
      if (!is.integer(object$forest$nodestatus)) {
        storage.mode(object$forest$nodestatus) <- "integer"
      }
      if (!is.double(object$forest$xbestsplit)) {
        storage.mode(object$forest$xbestsplit) <- "double"
      }
      if (!is.double(object$forest$nodepred)) {
        storage.mode(object$forest$nodepred) <- "double"
      }
      if (!is.integer(object$forest$bestvar)) {
        storage.mode(object$forest$bestvar) <- "integer"
      }
      if (!is.integer(object$forest$ndbigtree)) {
        storage.mode(object$forest$ndbigtree) <- "integer"
      }
      if (!is.integer(object$forest$ncat)) {
        storage.mode(object$forest$ncat) <- "integer"
      }

      ans <- .C("regForest",
        as.double(x),
        ypred = double(ntest),
        as.integer(mdim),
        as.integer(ntest),
        as.integer(ntree),
        object$forest$leftDaughter,
        object$forest$rightDaughter,
        object$forest$nodestatus,
        nrnodes,
        object$forest$xbestsplit,
        object$forest$nodepred,
        object$forest$bestvar,
        object$forest$ndbigtree,
        object$forest$ncat,
        as.integer(maxcat),
        as.integer(predict.all),
        treepred = as.double(treepred),
        as.integer(proximity),
        proximity = as.double(proxmatrix),
        nodes = as.integer(nodes),
        nodexts = as.integer(nodexts),
        # DUP=FALSE,
        PACKAGE = "randomForest"
      )[keepIndex]
      ## Apply bias correction if needed.
      yhat <- rep(NA, length(rn))
      names(yhat) <- rn
      if (!is.null(object$coefs)) {
        yhat[keep] <- object$coefs[1] + object$coefs[2] * ans$ypred
      } else {
        yhat[keep] <- ans$ypred
      }
      if (predict.all) {
        treepred <- matrix(NA, length(rn), ntree,
          dimnames = list(rn, NULL)
        )
        treepred[keep, ] <- ans$treepred
      }
      if (!proximity) {
        res <- if (predict.all) {
          list(aggregate = yhat, individual = treepred)
        } else {
          yhat
        }
      } else {
        res <- list(
          predicted = yhat,
          proximity = structure(ans$proximity,
            dim = c(ntest, ntest), dimnames = list(rn, rn)
          )
        )
      }
      if (nodes) {
        attr(res, "nodes") <- matrix(ans$nodexts, ntest, ntree,
          dimnames = list(rn[keep], 1:ntree)
        )
      }
    } else {
      countts <- matrix(0, ntest, nclass)
      t1 <- .C("classForest",
        mdim = as.integer(mdim),
        ntest = as.integer(ntest),
        nclass = as.integer(object$forest$nclass),
        maxcat = as.integer(maxcat),
        nrnodes = as.integer(nrnodes),
        jbt = as.integer(ntree),
        xts = as.double(x),
        xbestsplit = as.double(object$forest$xbestsplit),
        pid = object$forest$pid,
        cutoff = as.double(cutoff),
        countts = as.double(countts),
        treemap = as.integer(aperm(
          object$forest$treemap,
          c(2, 1, 3)
        )),
        nodestatus = as.integer(object$forest$nodestatus),
        cat = as.integer(object$forest$ncat),
        nodepred = as.integer(object$forest$nodepred),
        treepred = as.integer(treepred),
        jet = as.integer(numeric(ntest)),
        bestvar = as.integer(object$forest$bestvar),
        nodexts = as.integer(nodexts),
        ndbigtree = as.integer(object$forest$ndbigtree),
        predict.all = as.integer(predict.all),
        prox = as.integer(proximity),
        proxmatrix = as.double(proxmatrix),
        nodes = as.integer(nodes),
        # DUP=FALSE,
        PACKAGE = "randomForest"
      )
      if (out.type > 1) {
        out.class.votes <- t(matrix(t1$countts, nrow = nclass, ncol = ntest))
        if (norm.votes) {
          out.class.votes <-
            sweep(out.class.votes, 1, rowSums(out.class.votes), "/")
        }
        z <- matrix(NA, length(rn), nclass,
          dimnames = list(rn, object$classes)
        )
        z[keep, ] <- out.class.votes
        class(z) <- c(class(z), "votes")
        res <- z
      } else {
        out.class <- factor(rep(NA, length(rn)),
          levels = 1:length(object$classes),
          labels = object$classes
        )
        out.class[keep] <- object$classes[t1$jet]
        names(out.class)[keep] <- rn[keep]
        res <- out.class
      }
      if (predict.all) {
        treepred <- matrix(object$classes[t1$treepred],
          nrow = length(keep), dimnames = list(rn[keep], NULL)
        )
        res <- list(aggregate = res, individual = treepred)
      }
      if (proximity) {
        res <- list(predicted = res, proximity = structure(t1$proxmatrix,
          dim = c(ntest, ntest),
          dimnames = list(rn[keep], rn[keep])
        ))
      }
      if (nodes) {
        attr(res, "nodes") <- matrix(t1$nodexts, ntest, ntree,
          dimnames = list(rn[keep], 1:ntree)
        )
      }
    }
    res
  }
