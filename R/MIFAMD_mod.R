#' MIFAMD_mod: modified multiple imputation with FAMD
#'
#' @description \code{MIFAMD_mod} is a modified multiple imputation function with FAMD (Factorial
#' Analysis of Mixed Data)
#' that returns categorical columns results both in factor and in onehot probability vector form.
#' Please find the detailed documentation of \code{MIFAMD} in the 'missMDA' package.
#' Only the modifications are explained on this page.
#'
#' With \code{MIFAMD_mod}, not only the multiple imputation results are returned,
#' but the disjunctive multiple imputation results are also returned (The categorical columns
#'  are in form of onehot probability vector). Besides, instead of returning the final imputed dataset by
#'  performing one time FAMD imputation, \code{MIFAMD_mod} returns the final imputed dataset by combining
#'   the multiple imputation results with Rubin's Rule.
#' @param X Data frame with missing values.
#' @param ncp  Number of components used to reconstruct data with the FAMD reconstruction formular.
#' @param method "Regularized" by default or "EM"
#' @param coeff.ridge 1 by default to perform the regularized imputeFAMD algorithm.
#' Other regularization terms can be implemented by setting the value to less than 1
#' in order to regularized less (to get closer to the results of an EM method) or
#' more than 1 to regularized more (to get closer to the results of the proportion imputation).
#' @param threshold Threshold for the criterion convergence.
#' @param seed integer, by default seed = NULL implies that missing values are initially imputed by
#' the mean of each variable for the continuous variables and by the proportion of the category
#' for the categorical variables coded with indicator matrices of dummy variables.
#' Other values leads to a random initialization.
#' @param maxiter Maximum number of iterations for the algorithm.
#' @param nboot Number of multiple imputations.
#' @param verbose verbose=TRUE for screen printing of iteration numbers.
#' @export
#' @return \code{res.MI} A list of imputed dataset after mutiple imputation.
#' @return \code{res.MI.disj} A list of disjunctive imputed dataset after mutiple imputation.
#' @return \code{ximp} Final imputed dataset by combining \code{res.MI.disj} with Rubin's Rule.
#' @return \code{ximp.disj} Disjunctive imputed data matrix of same type as 'ximp' for the numeric columns.
#'  For the categorical columns, the prediction of probability for each category is shown in form of onehot probability vector.
#' @return \code{res.imputeFAMD} Output obtained with the function imputeFAMD (single imputation).
#' @return \code{call} The matched call.
#' @references
#' Audigier, V., Husson, F. & Josse, J. (2015). A principal components method to impute mixed data. Advances in Data Analysis and Classification, 10(1), 5-26. <doi:10.1007/s11634-014-0195-1>
#'
#' Audigier, V., Husson, F., Josse, J. (2017). MIMCA: Multiple imputation for categorical variables with multiple correspondence analysis. <doi:10.1007/s11222-016-9635-4>
#'
#' Little R.J.A., Rubin D.B. (2002) Statistical Analysis with Missing Data. Wiley series in probability and statistics, New-York

MIFAMD_mod <-
  function(X,
           ncp = 2,
           method = c("Regularized", "EM"),
           coeff.ridge = 1,
           threshold = 1e-06,
           seed = NULL,
           maxiter = 1000,
           nboot = 20,
           verbose = T) {

    # intern functions

    estim.sigma2 <- function(Xquanti, Xquali, M, Zhat, ncp, WW, D) {
      tab.disjonctif.NA <- function(tab) {
        if (ncol(tab) == 0) {
          return(NULL)
        }
        tab <- as.data.frame(tab)
        modalite.disjonctif <- function(i) {
          moda <- tab[, i]
          ######### Change#######
          if (is.numeric(moda)) {
            return(moda)
          }
          ######################
          nom <- names(tab)[i]
          n <- length(moda)
          moda <- as.factor(moda)
          x <- matrix(0, n, length(levels(moda)))
          ind <- (1:n) + n * (unclass(moda) - 1)
          indNA <- which(is.na(ind))
          x[(1:n) + n * (unclass(moda) - 1)] <- 1
          x[indNA, ] <- NA
          if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(
            1:nlevels(moda),
            "n", "N", "y", "Y"
          ))) {
            dimnames(x) <- list(row.names(tab), paste(nom,
              levels(moda),
              sep = "_"
            ))
          } else {
            dimnames(x) <- list(row.names(tab), levels(moda))
          }
          return(x)
        }
        if (ncol(tab) == 1) {
          res <- modalite.disjonctif(1)
        } else {
          res <- lapply(1:ncol(tab), modalite.disjonctif)
          res <- as.matrix(data.frame(res, check.names = FALSE))
        }
        return(res)
      }

      reconst.FAMD <- function(xxquanti, xxquali, M = NULL, D = NULL, ncp, coeff.ridge = 1, method = "em") {
        zz <- cbind.data.frame(xxquanti, xxquali)
        if (is.null(M)) {
          M <- c(1 / apply(xxquanti, 2, var), colMeans(zz)[-c(1:ncol(xxquanti))])
        }
        if (is.null(D)) {
          D <- rep(1 / nrow(zz), nrow(zz))
        }
        moy <- colMeans(zz)
        zzimp <- sweep(zz, MARGIN = 2, FUN = "-", STATS = moy)
        res.svd <- FactoMineR::svd.triplet(zzimp, col.w = M, row.w = D, ncp = ncp)
        tmp <- seq(ncol(zz) - ncol(xxquali))
        if (nrow(zz) > length(tmp)) {
          moyeig <- mean(res.svd$vs[tmp[-seq(ncp)]]^2)
        } else {
          moyeig <- mean(res.svd$vs[-c(1:ncp)]^2)
        }
        moyeig <- min(moyeig * coeff.ridge, res.svd$vs[ncp + 1]^2)
        moyeigret <- moyeig
        if (method == "em") {
          moyeig <- 0
        }
        if (ncp > 1) {
          eig.shrunk <- (res.svd$vs[1:ncp]^2 - moyeig) / res.svd$vs[1:ncp]
        } else if (ncp == 1) {
          eig.shrunk <- matrix((res.svd$vs[1:ncp]^2 - moyeig) / res.svd$vs[1:ncp], 1, 1)
        }
        zzhat <- tcrossprod(res.svd$U %*% diag(eig.shrunk), res.svd$V[which(apply(is.finite(res.svd$V), 1, any)), , drop = FALSE])
        zzhat <- sweep(zzhat, MARGIN = 2, FUN = "+", STATS = moy)
        return(list(zzhat = zzhat, moyeig = moyeigret, res.svd = res.svd, M = M))
      }


      nb.obs <- sum(WW[, -cumsum(sapply(Xquali, nlevels))])
      nb.obs.quanti <- sum(WW[, seq(ncol(Xquanti))])
      Zhat2 <- reconst.FAMD(Zhat[, seq(ncol(Xquanti))], Zhat[, -seq(ncol(Xquanti))], ncp = ncp, D = D, M = M)
      Residu <- as.matrix(Zhat - Zhat2$zzhat) %*% diag(M)^{
        1 / 2
      }
      Residu[is.na(cbind.data.frame(Xquanti, tab.disjonctif.NA(Xquali)))] <- 0
      sigma2 <- sum(WW[, seq(ncol(Xquanti))] * ((Residu[, seq(ncol(Xquanti))])^2)) / (nb.obs.quanti - (ncol(Xquanti) + ncp * (sum(D) - 1) + ncol(Xquanti) - ncp))
      return(sigma2)
    }


    imputeFAMD.stoch <- function(don,
                                 ncp = 4,
                                 method = c("Regularized", "EM"),
                                 row.w = NULL,
                                 coeff.ridge = 1,
                                 threshold = 1e-06,
                                 seed = NULL,
                                 maxiter = 1000,
                                 nboot,
                                 verbose = FALSE) {
      normtdc <- function(tab.disj, data.na) {
        tdc <- tab.disj
        tdc[tdc < 0] <- 0
        tdc[tdc > 1] <- 1
        col.suppr <- cumsum(sapply(data.na, function(x) {
          nlevels(x)
        }))
        tdc <- t(apply(tdc, 1, FUN = function(x, col.suppr) {
          if (sum(x[1:col.suppr[1]]) != 1) {
            x[1:col.suppr[1]] <- x[1:col.suppr[1]] / sum(x[1:col.suppr[1]])
          }
          if (length(col.suppr) > 1) {
            for (i in 2:length(col.suppr)) {
              x[(col.suppr[i - 1] + 1):(col.suppr[i])] <- x[(col.suppr[i -
                1] + 1):(col.suppr[i])] / sum(x[(col.suppr[i -
                1] + 1):col.suppr[i]])
            }
          }
          return(x)
        }, col.suppr = col.suppr))
        return(tdc)
      }
      draw <- function(tabdisj, Don) {
        nbdummy <- rep(1, ncol(Don))
        is.quali <- which(!sapply(Don, is.numeric))
        # change
        nbdummy[is.quali] <- sapply(Don[, is.quali, drop = FALSE], nlevels)
        ##
        vec <- c(0, cumsum(nbdummy))
        Donres <- Don
        for (i in is.quali) {
          Donres[, i] <- as.factor(levels(Don[, i])[apply(tabdisj[
            ,
            (vec[i] + 1):vec[i + 1]
          ], 1, function(x) {
            sample(1:length(x), size = 1, prob = x)
          })])
          Donres[, i] <- factor(Donres[, i], levels(Don[, is.quali, drop = FALSE][
            ,
            i
          ]))
        }
        return(don.imp = Donres)
      }

      quanti <- which(sapply(don, is.numeric))
      Ncol <- length(quanti) + sum(sapply(don[, -quanti], nlevels))

      W <- matrix(0, nrow(don), Ncol)
      W[!is.na(don)] <- 1

      if (is.null(row.w)) {
        Row.w <- as.integer(rep(1, nrow(don)))
      } else {
        Row.w <- row.w
      }

      if (is.integer(Row.w)) {
        D <- Row.w / sum(Row.w)
        D[which(Row.w == 0)] <- 1 / (1000 * nrow(don)) # D without 0
        WW <- diag(Row.w) %*% W
      } else {
        stop(paste0("row.w of class ", class(Row.w), ", it needs to be an integer"))
      }


      res.imp <- imputeFAMD_mod(
        X = don, ncp = ncp, method = method, row.w = D,
        coeff.ridge = coeff.ridge, threshold = threshold, seed = seed, maxiter = maxiter
      )

      var_homo <- estim.sigma2(
        Xquanti = don[, quanti, drop = FALSE],
        Xquali = don[, -quanti, drop = FALSE],
        M = 1 / apply(res.imp$tab.disj, 2, var),
        Zhat = res.imp$tab.disj,
        ncp = ncp,
        D = D, WW = WW
      )
      sigma2 <- var_homo / (1 / apply(res.imp$tab.disj, 2, var)[quanti])
      res.imp$fittedX <- res.imp$tab.disj
      res.imp$quanti.act <- which(sapply(don, is.numeric))



      classvar <- unlist(lapply(lapply(don, class), "[", 1)) # quand le type est "ordered" il y a 2 classes pour la variable
      if ("integer" %in% classvar) {
        classvar[classvar == "integer"] <- "numeric"
      }
      if ("ordered" %in% classvar) {
        classvar[classvar == "ordered"] <- "factor"
      }

      donimp <- don
      donimp[, which(classvar == "numeric")] <- res.imp$tab.disj[, seq(length(res.imp$quanti.act))] # les quanti active sont les premi?res variables du tdc
      missing.quanti <- is.na(don[, res.imp$quanti.act])
      res.MI <- vector("list", length = nboot)
      names(res.MI) <- paste("nboot=", 1:nboot, sep = "")
      res.MI.disj <- vector("list", length = nboot)
      names(res.MI.disj) <- paste("nboot=", 1:nboot, sep = "")
      for (i in seq(nboot)) {
        if (verbose) {
          cat(paste(i, "...", sep = ""))
        }
        donimp.tmp <- donimp
        if (any("factor" %in% classvar)) {
          tdc.imp <- res.imp$tab.disj[, (length(c(res.imp$quanti.act, res.imp$quanti.sup)) + 1):(ncol(res.imp$tab.disj) - length(res.imp$quali.sup)), drop = FALSE]
          tdc.norm <- normtdc(tab.disj = tdc.imp, data.na = don[, which(classvar == "factor"), drop = FALSE])
          donimp.quali <- draw(tdc.norm, don[, which(classvar == "factor"), drop = FALSE])
          donimp.tmp[, which(classvar == "factor")] <- donimp.quali[, names(which(classvar == "factor"))]
        }
        donimp.tmp[, res.imp$quanti.act][missing.quanti] <- res.imp$fittedX[, res.imp$quanti.act][missing.quanti] + mvtnorm::rmvnorm(nrow(don), sigma = diag(sigma2))[missing.quanti]
        res.MI[[paste("nboot=", i, sep = "")]] <- donimp.tmp
        ## change##
        donimp.disj <- donimp.tmp[, res.imp$quanti.act]
        donimp.disj[colnames(tdc.norm)] <- tdc.norm
        res.MI.disj[[paste("nboot=", i, sep = "")]] <- donimp.disj
        #########
      }
      if (verbose) {
        cat("\ndone!\n")
      }
      res.MI <- list(res.MI = res.MI, sigma2 = sigma2, res.MI.disj = res.MI.disj)
      class(res.MI) <- "MIFAMD"
      return(res.MI)
    }

    imputeFAMD.stoch.print <- function(don, ncp, method = c(
                                         "Regularized",
                                         "EM"
                                       ), row.w = NULL, coeff.ridge = 1, threshold = 1e-06,
                                       seed = NULL, maxiter = 1000, verbose, printm) {
      if (verbose) {
        cat(paste(printm, "...", sep = ""))
      }
      res <- imputeFAMD.stoch(
        don = don, ncp = ncp, method = method,
        row.w = row.w, coeff.ridge = coeff.ridge, threshold = threshold,
        seed = seed, maxiter = maxiter, nboot = 1
      )
      rs <- list(res = res$res.MI[[1]], res.disj = res$res.MI.disj[[1]])
      return(rs)
    }


    if (sum(sapply(X, is.numeric)) == ncol(X)) {
      rs <- missMDA::MIPCA(
        X = X,
        ncp = ncp,
        method = method,
        threshold = threshold,
        nboot = nboot,
        method.mi = "Boot",
        verbose = verbose
      )
      df_new_merge <- abind::abind(rs$res.MI, along = 3)
      ximp <- data.frame(apply(df_new_merge, c(1, 2), mean))
      res <- list(
        res.MI = rs$res.MI,
        res.MI.disj = rs$res.MI,
        ximp.disj = ximp,
        ximp = ximp,
        res.imputeFAMD = rs$res.imputePCA,
        call = list(X = X, nboot = nboot, ncp = ncp, coeff.ridge = coeff.ridge, threshold = threshold, seed = seed, maxiter = maxiter)
      )
      return(res)
    }
    # }else if(sum(sapply(X,is.numeric))==0){
    #   rs <- missMDA::imputeMCA(don = X,
    #             ncp = ncp,
    #             method = method,
    #             threshold = threshold,
    #             coeff.ridge= coeff.ridge,
    #             seed = seed,
    #             maxiter = maxiter)
    #   res <- list(
    #     res.MI = rs$res.MI,
    #     res.MI.disj = rs$res.MI,
    #     ximp.disj = rs$res.imputePCA,
    #     ximp = rs$res.imputePCA,
    #     res.imputeFAMD = imputePCA(X, ncp = ncp, coeff.ridge = coeff.ridge, method = method, threshold = threshold, maxiter = maxiter, seed = seed),
    #     call = list(X = X, nboot = nboot, ncp = ncp, coeff.ridge = coeff.ridge, threshold = threshold, seed = seed, maxiter = maxiter)
    #   )}

    ## Add:
    col_cat <- which(!sapply(X, is.numeric))
    exist_cat <- !all(c(0, col_cat) == c(0))
    if (exist_cat) {
      name_cat <- colnames(X)[col_cat]
      # Deal with the problem that nlevels(df[[col]]) > length(unique(df[[col]]))
      for (col in name_cat) {
        X[[col]] <- factor(as.character(X[[col]]))
      }
      # remember the levels for each categorical column
      dict_lev <- dict_level(X, col_cat)
      # preserve colnames for ximp.disj
      dummy <- dummyVars(" ~ .", data = X, sep = "_")
      col_names.disj <- colnames(data.frame(predict(dummy, newdata = X)))
      # represent the factor columns with their ordinal levels
      X <- factor_ordinal_encode(X, col_cat)
      # Create dict_cat with categroical columns
      dict_cat <- dict_onehot(X, col_cat)
      dummy <- dummyVars(" ~ .", data = X, sep = "_")
      col_names.disj.after <- colnames(data.frame(predict(dummy, newdata = X)))
      # Build a dictionary between Y1_1 and Y1_a
      dict_disj <- list()
      for (i in seq(length(col_names.disj.after))) {
        dict_disj[[col_names.disj.after[i]]] <- col_names.disj[i]
      }
    }



    # variables are ordered
    don <- X[, c(which(sapply(X, is.numeric)), which(!sapply(X, is.numeric)))]
    # print
    temp <- if (coeff.ridge == 1) {
      "regularized"
    } else if ((coeff.ridge == 0) | (method == "EM")) {
      "EM"
    } else {
      paste("coeff.ridge=", coeff.ridge)
    }

    if (verbose) {
      cat(
        "Multiple Imputation using", temp, "FAMD using", nboot,
        "imputed arrays", "\n"
      )
    }

    # multiple imputation
    n <- nrow(don)
    Boot <- matrix(sample(1:n, size = nboot * n, replace = T), n, nboot)
    Boot <- lapply(as.data.frame(Boot), FUN = function(xx) {
      yy <- as.factor(xx)
      levels(yy) <- c(levels(yy), xx[which(!xx %in% as.numeric(as.character(levels(yy))))])
      return(yy)
    })
    Weight <- as.data.frame(matrix(0, n, nboot, dimnames = list(1:n, paste("nboot=", 1:nboot, sep = ""))))
    Boot.table <- lapply(Boot, table)
    for (i in 1:nboot) {
      Weight[names(Boot.table[[i]]), i] <- Boot.table[[i]]
    }
    Weight <- do.call(cbind.data.frame, lapply(Weight, as.integer))
    res.imp <- mapply(Weight,
      FUN = imputeFAMD.stoch.print,
      MoreArgs = list(
        don = don, ncp = ncp,
        coeff.ridge = coeff.ridge,
        method = method,
        threshold = threshold,
        maxiter = maxiter,
        verbose = verbose,
        seed = NULL
      ), printm = as.character(1:nboot),
      SIMPLIFY = FALSE
    )
    ## Change##
    res <- list(res.MI = res.imp)
    res.MI <- lapply(res$res.MI, function(xx, nom) {
      return(xx$res[, nom])
    }, nom = colnames(don))
    res.MI.disj <- lapply(res$res.MI, function(xx) {
      return(xx$res.disj)
    })
    df_new_merge <- abind::abind(res.MI.disj, along = 3)
    ximp.all <- data.frame(apply(df_new_merge, c(1, 2), mean))
    ximp.disj <- ximp.all

    if (any(!sapply(don, is.numeric))) {
      names_cat <- names(dict_cat)
      for (name in names_cat) {
        ximp.all[[name]] <- apply(ximp.all[dict_cat[[name]]], 1, which_max_cat, name, dict_cat)
        ximp.all[[name]] <- unlist(ximp.all[[name]])
        ximp.all[[name]] <- factor(ximp.all[[name]])
      }
      ximp <- ximp.all[colnames(X)]
    } else {
      ximp <- ximp.disj
    }

    # Add: change back to original levels
    if (exist_cat) {
      for (col in name_cat) {
        levels(ximp[[col]]) <- dict_lev[[col]]
      }
      colnames(ximp.disj) <- dict_disj[colnames(ximp.disj)]
      for (i in seq(length(res.MI))) {
        for (col in name_cat) {
          levels(res.MI[[i]][[col]]) <- dict_lev[[col]]
        }
        colnames(res.MI.disj[[i]]) <- dict_disj[colnames(res.MI.disj[[i]])]
      }
    }
    ###########
    res <- list(
      res.MI = res.MI,
      res.MI.disj = res.MI.disj,
      ximp.disj = ximp.disj,
      ximp = ximp,
      res.imputeFAMD = imputeFAMD_mod(X, ncp = ncp, coeff.ridge = coeff.ridge, method = method, threshold = threshold, maxiter = maxiter, seed = seed),
      call = list(X = X, nboot = nboot, ncp = ncp, coeff.ridge = coeff.ridge, threshold = threshold, seed = seed, maxiter = maxiter)
    )
    class(res) <- c("MIFAMD", "list")
    if (verbose) {
      cat("\ndone!\n")
    }
    return(res)
  }
