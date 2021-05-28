
#################### Original function with print added#############################

estim_ncpFAMD_1 <- function(don, ncp.min = 0, ncp.max = 5, method = c("Regularized", "EM"),
                            method.cv = c("Kfold", "loo"), nbsim = 100, pNA = 0.05, ind.sup = NULL, sup.var = NULL,
                            threshold = 1e-04, verbose = TRUE) {
  tab.disjonctif.NA <- function(tab) {
    print("In tab.disjonctif.NA")
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
      moda <- tab[, i]
      nom <- names(tab)[i]
      print(paste("nom:", nom))
      ## Add to stop the algorithm faster####
      if (nom == "Y5") {
        stop("Problem in tab.disjonctif.NA")
      }
      ######################################
      n <- length(moda)
      print(paste("n:", n))
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind <- (1:n) + n * (unclass(moda) - 1)
      indNA <- which(is.na(ind))
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA, ] <- NA
      if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), "n", "N", "y", "Y"))) {
        dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), sep = "."))
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
    print("Done tab.disjonctif.NA")
    return(res)
  }

  # from missForest package
  prodna <- function(x, noNA) {
    print("In prodna")
    n <- nrow(x)
    p <- ncol(x)
    NAloc <- rep(FALSE, n * p)
    NAloc[sample(n * p, floor(n * p * noNA))] <- TRUE
    x[matrix(NAloc, nrow = n, ncol = p)] <- NA
    print("Done prodna")
    return(x)
  }

  print("In estim_ncpFAMD")
  #  if(!("numeric"%in%(lapply(don,class)) & "factor"%in%(lapply(don,class)))){stop("Your data set must contain mixed data.")}
  don <- as.data.frame(don)
  if (!is.null(ind.sup)) don <- don[-ind.sup, ]
  if (!is.null(sup.var)) don <- don[, -sup.var]
  if ((sum(sapply(don, is.numeric)) == 0) || (sum(!sapply(don, is.numeric)) == 0)) {
    stop("Your data set must contain mixed data.")
  }
  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"), several.ok = T)[1]
  method.cv <- match.arg(method.cv, c("loo", "Kfold", "kfold", "LOO"), several.ok = T)[1]
  method <- tolower(method)
  method.cv <- tolower(method.cv)
  # reagencement des variables
  jeu <- don[, c(which(sapply(don, is.numeric)), which(sapply(don, is.factor))), drop = F]
  nbquanti <- sum(sapply(don, is.numeric))
  jeu[, 1:nbquanti] <- lapply(jeu[, 1:nbquanti, drop = FALSE], as.double)

  # suppression niveaux non pris
  jeu <- droplevels(jeu)

  vrai.tab <- cbind(jeu[, 1:nbquanti, drop = F], tab.disjonctif.NA(jeu[, (nbquanti + 1):ncol(jeu), drop = F]))
  if (method.cv == "kfold") {
    print("In kfold")
    res <- matrix(NA, ncp.max - ncp.min + 1, nbsim)
    if (verbose) pb <- txtProgressBar(min = 1 / nbsim * 100, max = 100, style = 3)
    for (sim in 1:nbsim) {
      print(paste("sim:", sim))
      continue <- TRUE
      while (continue) {
        jeuNA <- prodna(jeu, pNA)
        continue <- continue <- (sum(unlist(sapply(as.data.frame(droplevels(jeuNA[, -c(1:nbquanti), drop = F])), nlevels))) != sum(unlist(sapply(jeu, nlevels))))
      }


      for (nbaxes in ncp.min:ncp.max) {
        print(paste("nbaxes:", nbaxes))
        tab.disj.comp <- imputeFAMD(as.data.frame(jeuNA), ncp = nbaxes, method = method, threshold = threshold, maxiter = maxiter)$tab.disj
        if (sum(is.na(jeuNA)) != sum(is.na(jeu))) {
          res[nbaxes - ncp.min + 1, sim] <- sum((tab.disj.comp - vrai.tab)^2, na.rm = TRUE) / (sum(is.na(tab.disjonctif.NA(jeuNA))) - sum(is.na(tab.disjonctif.NA(jeu))))
        }
      }
      if (verbose) setTxtProgressBar(pb, sim / nbsim * 100)
    }
    crit <- apply(res, 1, mean, na.rm = TRUE)
    names(crit) <- c(ncp.min:ncp.max)
    if (verbose) close(pb)
    result <- list(ncp = as.integer(which.min(crit) + ncp.min - 1), criterion = crit)
    print("Done estim_ncpFAMD")
    return(result)
  }

  if (method.cv == "loo") {
    if (verbose) pb <- txtProgressBar(min = 0, max = 100, style = 3)
    crit <- NULL
    tab.disj.hat <- vrai.tab
    col.in.indicator <- c(0, rep(1, nbquanti), sapply(jeu[, (nbquanti:ncol(jeu)), drop = F], nlevels))
    for (nbaxes in ncp.min:ncp.max) {
      for (i in 1:nrow(jeu)) {
        for (j in 1:ncol(jeu)) {
          if (!is.na(jeu[i, j])) {
            jeuNA <- as.matrix(jeu)
            jeuNA[i, j] <- NA
            if (!any(summary(jeuNA[, j] == 0))) {
              tab.disj.hat[i, (cumsum(col.in.indicator)[j] + 1):(cumsum(col.in.indicator)[j + 1])] <- imputeFAMD(as.data.frame(jeuNA),
                ncp = nbaxes, method = method, threshold = threshold
              )$tab.disj[
                i,
                (cumsum(col.in.indicator)[j] + 1):(cumsum(col.in.indicator)[j + 1])
              ]
            }
          }
        }
        if (verbose) setTxtProgressBar(pb, round((((1:length(ncp.min:ncp.max))[which(nbaxes == (ncp.min:ncp.max))] - 1) * nrow(jeu) + i) / (length(ncp.min:ncp.max) * nrow(jeu)) * 100))
      }
      crit <- c(crit, mean((tab.disj.hat - vrai.tab)^2, na.rm = TRUE))
    }
    if (verbose) close(pb)
    names(crit) <- c(ncp.min:ncp.max)
    print("Done estim_ncpFAMD")
    return(list(ncp = as.integer(which.min(crit) + ncp.min -
      1), criterion = crit))
  }
}

#################### Original function with MSE part changed and maxiter added#############################

estim_ncpFAMD_2 <- function(don, ncp.min = 0, ncp.max = 5, method = c("Regularized", "EM"),
                            method.cv = c("Kfold", "loo"), nbsim = 100, pNA = 0.05, ind.sup = NULL, sup.var = NULL,
                            threshold = 1e-04, verbose = TRUE, maxiter = 1000) {
  tab.disjonctif.NA <- function(tab) {
    print("In tab.disjonctif.NA")
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
      moda <- tab[, i]
      nom <- names(tab)[i]
      print(paste("nom:", nom))
      n <- length(moda)
      print(paste("n:", n))
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind <- (1:n) + n * (unclass(moda) - 1)
      indNA <- which(is.na(ind))
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA, ] <- NA
      if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), "n", "N", "y", "Y"))) {
        dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), sep = "."))
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
    print("Done tab.disjonctif.NA")
    return(res)
  }

  # from missForest package
  prodna <- function(x, noNA) {
    print("In prodna")
    n <- nrow(x)
    p <- ncol(x)
    NAloc <- rep(FALSE, n * p)
    NAloc[sample(n * p, floor(n * p * noNA))] <- TRUE
    x[matrix(NAloc, nrow = n, ncol = p)] <- NA
    print("Done prodna")
    return(x)
  }

  print("In estim_ncpFAMD")
  #  if(!("numeric"%in%(lapply(don,class)) & "factor"%in%(lapply(don,class)))){stop("Your data set must contain mixed data.")}
  don <- as.data.frame(don)
  if (!is.null(ind.sup)) don <- don[-ind.sup, ]
  if (!is.null(sup.var)) don <- don[, -sup.var]
  if ((sum(sapply(don, is.numeric)) == 0) || (sum(!sapply(don, is.numeric)) == 0)) {
    stop("Your data set must contain mixed data.")
  }
  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"), several.ok = T)[1]
  method.cv <- match.arg(method.cv, c("loo", "Kfold", "kfold", "LOO"), several.ok = T)[1]
  method <- tolower(method)
  method.cv <- tolower(method.cv)
  # reagencement des variables
  jeu <- don[, c(which(sapply(don, is.numeric)), which(sapply(don, is.factor))), drop = F]
  nbquanti <- sum(sapply(don, is.numeric))
  jeu[, 1:nbquanti] <- lapply(jeu[, 1:nbquanti, drop = FALSE], as.double)

  # suppression niveaux non pris
  jeu <- droplevels(jeu)

  vrai.tab <- cbind(jeu[, 1:nbquanti, drop = F], tab.disjonctif.NA(jeu[, (nbquanti + 1):ncol(jeu), drop = F]))
  if (method.cv == "kfold") {
    print("In kfold")
    res <- matrix(NA, ncp.max - ncp.min + 1, nbsim)
    if (verbose) pb <- txtProgressBar(min = 1 / nbsim * 100, max = 100, style = 3)
    for (sim in 1:nbsim) {
      print(paste("sim:", sim))
      continue <- TRUE
      while (continue) {
        jeuNA <- prodna(jeu, pNA)
        continue <- continue <- (sum(unlist(sapply(as.data.frame(droplevels(jeuNA[, -c(1:nbquanti), drop = F])), nlevels))) != sum(unlist(sapply(jeu, nlevels))))
      }
      ######### Changed part 1##############
      jeuNA.tab <- cbind(jeuNA[, 1:nbquanti, drop = F], tab.disjonctif.NA(jeuNA[, (nbquanti + 1):ncol(jeuNA), drop = F]))
      #####################################
      for (nbaxes in ncp.min:ncp.max) {
        # print(paste("nbaxes:", nbaxes))
        # Add the option maxiter in the impute function
        tab.disj.comp <- imputeFAMD(as.data.frame(jeuNA), ncp = nbaxes, method = method, threshold = threshold, maxiter = maxiter)$tab.disj
        if (sum(is.na(jeuNA)) != sum(is.na(jeu))) {
          ######### Changed part 2##############
          # res[nbaxes - ncp.min + 1, sim] <- sum((tab.disj.comp - vrai.tab)^2, na.rm = TRUE)/(sum(is.na(tab.disjonctif.NA(jeuNA))) - sum(is.na(tab.disjonctif.NA(jeu))))
          res[nbaxes - ncp.min + 1, sim] <- sum((tab.disj.comp - vrai.tab)^2, na.rm = TRUE) / (sum(is.na(jeuNA.tab)) - sum(is.na(vrai.tab)))
          #####################################
        }
      }
      if (verbose) setTxtProgressBar(pb, sim / nbsim * 100)
    }
    crit <- apply(res, 1, mean, na.rm = TRUE)
    names(crit) <- c(ncp.min:ncp.max)
    if (verbose) close(pb)
    result <- list(ncp = as.integer(which.min(crit) + ncp.min - 1), criterion = crit)
    print("Done estim_ncpFAMD")
    return(result)
  }

  if (method.cv == "loo") {
    if (verbose) pb <- txtProgressBar(min = 0, max = 100, style = 3)
    crit <- NULL
    tab.disj.hat <- vrai.tab
    col.in.indicator <- c(0, rep(1, nbquanti), sapply(jeu[, (nbquanti:ncol(jeu)), drop = F], nlevels))
    for (nbaxes in ncp.min:ncp.max) {
      for (i in 1:nrow(jeu)) {
        for (j in 1:ncol(jeu)) {
          if (!is.na(jeu[i, j])) {
            jeuNA <- as.matrix(jeu)
            jeuNA[i, j] <- NA
            if (!any(summary(jeuNA[, j] == 0))) {
              tab.disj.hat[i, (cumsum(col.in.indicator)[j] + 1):(cumsum(col.in.indicator)[j + 1])] <- imputeFAMD(as.data.frame(jeuNA),
                ncp = nbaxes, method = method, threshold = threshold
              )$tab.disj[
                i,
                (cumsum(col.in.indicator)[j] + 1):(cumsum(col.in.indicator)[j + 1])
              ]
            }
          }
        }
        if (verbose) setTxtProgressBar(pb, round((((1:length(ncp.min:ncp.max))[which(nbaxes == (ncp.min:ncp.max))] - 1) * nrow(jeu) + i) / (length(ncp.min:ncp.max) * nrow(jeu)) * 100))
      }
      crit <- c(crit, mean((tab.disj.hat - vrai.tab)^2, na.rm = TRUE))
    }
    if (verbose) close(pb)
    names(crit) <- c(ncp.min:ncp.max)
    print("Done estim_ncpFAMD")
    return(list(ncp = as.integer(which.min(crit) + ncp.min -
      1), criterion = crit))
  }
}

################### estim_ncpFAMD_2 with some modifications on impute function##############################
#' estim_ncpFAMD_mod: modified estim_ncpFAMD from missMDA
#' 
#' @description \code{estim_ncpFAMD_mod} is a modified version of the function \code{estim_ncpFAMD} by Francois Husson and Julie Josse.
#' Please find the detailed documentation of \code{estim_ncpFAMD} in the missMDA package. 
#' The modified version of \code{estim_ncpFAMD} ameliorates the convergence problem by modifying the function \code{impute}, and resolve the memory error problem.
#' @export
estim_ncpFAMD_mod <- function(don, ncp.min = 0, ncp.max = 5, method = c("Regularized", "EM"),
                            method.cv = c("Kfold", "loo"), nbsim = 100, pNA = 0.05, ind.sup = NULL, sup.var = NULL,
                            threshold = 1e-04, verbose = TRUE, maxiter = 1000) {
  tab.disjonctif.NA <- function(tab) {
    # print("In tab.disjonctif.NA")
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
      moda <- tab[, i]
      nom <- names(tab)[i]
      # print(paste("nom:", nom))
      n <- length(moda)
      # print(paste("n:", n))
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind <- (1:n) + n * (unclass(moda) - 1)
      indNA <- which(is.na(ind))
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA, ] <- NA
      if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), "n", "N", "y", "Y"))) {
        dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), sep = "_"))
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
    # print("Done tab.disjonctif.NA")
    return(res)
  }

  # from missForest package
  prodna <- function(x, noNA) {
    # print("In prodna")
    n <- nrow(x)
    p <- ncol(x)
    NAloc <- rep(FALSE, n * p)
    NAloc[sample(n * p, floor(n * p * noNA))] <- TRUE
    x[matrix(NAloc, nrow = n, ncol = p)] <- NA
    # print("Done prodna")
    return(x)
  }

  # print("In estim_ncpFAMD")
  #  if(!("numeric"%in%(lapply(don,class)) & "factor"%in%(lapply(don,class)))){stop("Your data set must contain mixed data.")}
  don <- as.data.frame(don)
  if (!is.null(ind.sup)) don <- don[-ind.sup, ]
  if (!is.null(sup.var)) don <- don[, -sup.var]
  if ((sum(sapply(don, is.numeric)) == 0) || (sum(!sapply(don, is.numeric)) == 0)) {
    stop("Your data set must contain mixed data.")
  }
  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"), several.ok = T)[1]
  method.cv <- match.arg(method.cv, c("loo", "Kfold", "kfold", "LOO"), several.ok = T)[1]
  method <- tolower(method)
  method.cv <- tolower(method.cv)
  # reagencement des variables
  jeu <- don[, c(which(sapply(don, is.numeric)), which(sapply(don, is.factor))), drop = F]
  nbquanti <- sum(sapply(don, is.numeric))
  jeu[, 1:nbquanti] <- lapply(jeu[, 1:nbquanti, drop = FALSE], as.double)

  # suppression niveaux non pris
  jeu <- droplevels(jeu)

  vrai.tab <- cbind(jeu[, 1:nbquanti, drop = F], tab.disjonctif.NA(jeu[, (nbquanti + 1):ncol(jeu), drop = F]))
  if (method.cv == "kfold") {
    # print("In kfold")
    res <- matrix(NA, ncp.max - ncp.min + 1, nbsim)
    if (verbose) pb <- txtProgressBar(min = 1 / nbsim * 100, max = 100, style = 3)
    for (sim in 1:nbsim) {
      # print(paste("sim:", sim))
      continue <- TRUE
      while (continue) {
        jeuNA <- prodna(jeu, pNA)
        continue <- continue <- (sum(unlist(sapply(as.data.frame(droplevels(jeuNA[, -c(1:nbquanti), drop = F])), nlevels))) != sum(unlist(sapply(jeu, nlevels))))
      }
      ######### Changed part 1##############
      jeuNA.tab <- cbind(jeuNA[, 1:nbquanti, drop = F], tab.disjonctif.NA(jeuNA[, (nbquanti + 1):ncol(jeuNA), drop = F]))
      #####################################
      for (nbaxes in ncp.min:ncp.max) {
        # print(paste("nbaxes:", nbaxes))
        tab.disj.comp <- imputeFAMD_mod(as.data.frame(jeuNA), ncp = nbaxes, method = method, threshold = threshold, maxiter = maxiter)$tab.disj
        if (sum(is.na(jeuNA)) != sum(is.na(jeu))) {
          ######### Changed part 2##############
          # res[nbaxes - ncp.min + 1, sim] <- sum((tab.disj.comp - vrai.tab)^2, na.rm = TRUE)/(sum(is.na(tab.disjonctif.NA(jeuNA))) - sum(is.na(tab.disjonctif.NA(jeu))))
          res[nbaxes - ncp.min + 1, sim] <- sum((tab.disj.comp - vrai.tab)^2, na.rm = TRUE) / (sum(is.na(jeuNA.tab)) - sum(is.na(vrai.tab)))
          #####################################
        }
      }
      if (verbose) setTxtProgressBar(pb, sim / nbsim * 100)
    }
    crit <- apply(res, 1, mean, na.rm = TRUE)
    names(crit) <- c(ncp.min:ncp.max)
    if (verbose) close(pb)
    result <- list(ncp = as.integer(which.min(crit) + ncp.min - 1), criterion = crit)
    # print("Done estim_ncpFAMD")
    return(result)
  }

  if (method.cv == "loo") {
    if (verbose) pb <- txtProgressBar(min = 0, max = 100, style = 3)
    crit <- NULL
    tab.disj.hat <- vrai.tab
    col.in.indicator <- c(0, rep(1, nbquanti), sapply(jeu[, (nbquanti:ncol(jeu)), drop = F], nlevels))
    for (nbaxes in ncp.min:ncp.max) {
      for (i in 1:nrow(jeu)) {
        for (j in 1:ncol(jeu)) {
          if (!is.na(jeu[i, j])) {
            jeuNA <- as.matrix(jeu)
            jeuNA[i, j] <- NA
            if (!any(summary(jeuNA[, j] == 0))) {
              tab.disj.hat[i, (cumsum(col.in.indicator)[j] + 1):(cumsum(col.in.indicator)[j + 1])] <- imputeFAMD(as.data.frame(jeuNA),
                ncp = nbaxes, method = method, threshold = threshold, maxiter = maxiter
              )$tab.disj[
                i,
                (cumsum(col.in.indicator)[j] + 1):(cumsum(col.in.indicator)[j + 1])
              ]
            }
          }
        }
        if (verbose) setTxtProgressBar(pb, round((((1:length(ncp.min:ncp.max))[which(nbaxes == (ncp.min:ncp.max))] - 1) * nrow(jeu) + i) / (length(ncp.min:ncp.max) * nrow(jeu)) * 100))
      }
      crit <- c(crit, mean((tab.disj.hat - vrai.tab)^2, na.rm = TRUE))
    }
    if (verbose) close(pb)
    names(crit) <- c(ncp.min:ncp.max)
    # print("Done estim_ncpFAMD")
    return(list(ncp = as.integer(which.min(crit) + ncp.min -
      1), criterion = crit))
  }
}

#' imputeFAMD_mod: modified imputeFAMD from missMDA
#' 
#' @description \code{imputeFAMD_mod} is a modified version of the function \code{imputeFAMD} by Francois Husson and Julie Josse.
#' Please find the detailed documentation of \code{imputeFAMD} in the missMDA package. 
#' The modified version of \code{imputeFAMD} ameliorates the convergence problem by modifying the function \code{impute}.
#' @export
imputeFAMD_mod <- function(X, ncp = 2, method = c("Regularized", "EM"), row.w = NULL, coeff.ridge = 1, threshold = 1e-6,
                         ind.sup = NULL, sup.var = NULL, seed = NULL, maxiter = 1000, ...) {
  X <- as.data.frame(X)
  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"), several.ok = T)[1]
  method <- tolower(method)
  type <- rep("s", ncol(X))
  type[!sapply(X, is.numeric)] <- "n"
  res <- imputeMFA_mod(X = X, group = rep(1, ncol(X)), type = type, ncp = ncp, method = method, row.w = row.w, coeff.ridge = coeff.ridge, ind.sup = ind.sup, num.group.sup = sup.var, threshold = threshold, seed = seed, maxiter = maxiter)
  return(res)
}

#' impute_mod: modified impute from missMDA
#' 
#' @description \code{impute_mod} is a modified version of the function \code{impute} by Francois Husson and Julie Josse.
#' Please find the detailed documentation of \code{impute} in the missMDA package. 
#' The modified version of \code{impute} ameliorates the convergence problem by restricting the onehot probabilities after PCA being positive.
#' @export
impute_mod <- function(X, group, ncp = 2, type = rep("s", length(group)),
                     method = NULL, threshold = 1e-06, ind.sup = NULL, num.group.sup = NULL, seed = NULL, maxiter = 1000,
                     row.w = NULL, coeff.ridge = 1, ...) {
  # print("in impute")
  moy.p <- function(V, poids) {
    res <- sum(V * poids, na.rm = TRUE) / sum(poids[!is.na(V)])
  }
  ec <- function(V, poids) {
    res <- sqrt(sum(V^2 * poids, na.rm = TRUE) / sum(poids[!is.na(V)]))
  }
  tab.disjonctif.NA <- function(tab) {
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
      moda <- tab[, i]
      nom <- names(tab)[i]
      n <- length(moda)
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind <- (1:n) + n * (unclass(moda) - 1)
      indNA <- which(is.na(ind))
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA, ] <- NA
      if ((ncol(tab) != 1) & (levels(moda)[1] %in%
        c(1:nlevels(moda), "n", "N", "y", "Y"))) {
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
  find.category <- function(X, tabdisj) {
    nbdummy <- rep(1, ncol(X))
    is.quali <- which(!unlist(lapply(X, is.numeric)))
    nbdummy[is.quali] <- unlist(lapply(X[, is.quali,
      drop = FALSE
    ], nlevels))
    vec <- c(0, cumsum(nbdummy))
    Xres <- X
    for (i in 1:ncol(X)) {
      if (i %in% is.quali) {
        Xres[, i] <- as.factor(levels(X[, i])[apply(tabdisj[
          ,
          (vec[i] + 1):vec[i + 1]
        ], 1, which.max)])
      } else {
        Xres[, i] <- tabdisj[, vec[i] + 1]
      }
    }
    return(Xres)
  }

  method <- match.arg(method, c(
    "Regularized", "regularized",
    "EM", "em"
  ), several.ok = T)[1]
  method <- tolower(method)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  X <- as.data.frame(X)
  X <- droplevels(X)
  if ("n" %in% type) {
    niveau <- NULL
    for (j in 1:ncol(X)) {
      if (!is.numeric(X[, j])) {
        niveau <- c(niveau, levels(X[, j]))
      }
    }
    # Construct Y7_1,...
    for (j in 1:ncol(X)) {
      if (!is.numeric(X[, j])) {
        if (sum(niveau %in% levels(X[, j])) != nlevels(X[, j])) {
          levels(X[, j]) <- paste(colnames(X)[j], levels(X[, j]), sep = "_")
        }
      }
    }
  }
  group.mod <- group
  Xhat <- matrix(0, nrow(X), 0)
  Xhat2 <- matrix(0, nrow(X), 0)
  MM <- vector(mode = "list", length = length(group))
  ET <- vector(mode = "list", length = length(group))
  tab.disj.comp <- vector(mode = "list", length = length(group))
  ponderation <- rep(1, length(group))

  idx_obs <- list()

  # case ncp=0
  if (ncp == 0) {
    result <- list()

    if (sum(unlist(lapply(X, is.numeric))) > 0) {
      ind.quanti <- which(unlist(lapply(X, is.numeric)))
    } else {
      ind.quanti <- NULL
    }
    if (sum(lapply(X, class) == "factor") > 0) {
      ind.quali <- which((lapply(X, class)) == "factor")
    } else {
      ind.quali <- NULL
    }
    # quali
    result$completeObs <- X
    if (!is.null(ind.quali)) result$completeObs[, ind.quali] <- find.category(X[, ind.quali, drop = F], FactoMineR::tab.disjonctif.prop(X[, ind.quali, drop = F], row.w = row.w))
    # quanti
    if (!is.null(ind.quanti)) {
      tab.disj <- X[, ind.quanti, drop = F]
      Moy <- matrix(colMeans(tab.disj, na.rm = T), nrow = nrow(tab.disj), ncol = ncol(tab.disj), byrow = T)
      tab.disj[is.na(tab.disj)] <- Moy[is.na(tab.disj)]
      result$completeObs[, ind.quanti] <- tab.disj
    }

    # tdc
    nbdummy <- rep(1, ncol(X))
    is.quali <- which(!unlist(lapply(X, is.numeric)))
    if (length(is.quali) != 0) {
      nbdummy[is.quali] <- unlist(lapply(X[, is.quali, drop = FALSE], nlevels))
      tabdisj <- matrix(NA, nrow(X), ncol = sum(nbdummy))
      tabdisj[, cumsum(nbdummy)[which(nbdummy == 1)]] <- as.matrix(result$completeObs[, ind.quanti, drop = F])
      auxQuali <- FactoMineR::tab.disjonctif.prop(X[, ind.quali, drop = F], row.w = row.w)
      tabdisj[, -cumsum(nbdummy)[which(nbdummy == 1)]] <- auxQuali
      rownames(tabdisj) <- rownames(X)
      colnames(tabdisj) <- paste0("v", 1:ncol(tabdisj))
      colnames(tabdisj)[cumsum(nbdummy)[which(nbdummy == 1)]] <- colnames(result$completeObs[, ind.quanti, drop = F])
      colnames(tabdisj)[-cumsum(nbdummy)[which(nbdummy == 1)]] <- colnames(auxQuali)
      result$tab.disj <- tabdisj
    }

    # ind.var, group.mod

    for (g in 1:length(group)) {
      if (g == 1) {
        aux.base <- X[, 1:group[1], drop = FALSE]
      } else {
        aux.base <- X[, (cumsum(group)[g - 1] + 1):cumsum(group)[g],
          drop = FALSE
        ]
      }

      if (type[g] == "n") {
        tab.disj <- FactoMineR::tab.disjonctif.prop(aux.base, seed, row.w = row.w)
        tab.disj.comp[[g]] <- tab.disj
        group.mod[g] <- ncol(tab.disj)
      }
    }
    ind.var <- vector(mode = "list", length = length(group))
    ind.var[[1]] <- 1:group.mod[1]
    for (g in 2:length(group)) ind.var[[g]] <- (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]

    result$call$group.mod <- group.mod
    ind.var <- vector(mode = "list", length = length(group))
    ind.var[[1]] <- 1:result$call$group.mod[1]
    for (g in 2:length(group)) ind.var[[g]] <- (cumsum(result$call$group.mod)[g - 1] + 1):cumsum(result$call$group.mod)[g]
    result$call$ind.var <- ind.var
    return(result)
  }


  for (g in 1:length(group)) {
    if (g == 1) {
      aux.base <- X[, 1:group[1], drop = FALSE]
    } else {
      aux.base <- X[, (cumsum(group)[g - 1] + 1):cumsum(group)[g], drop = FALSE]
    }
    if (type[g] == "s") {
      Xhat2 <- cbind.data.frame(Xhat2, aux.base)
      MM[[g]] <- apply(as.data.frame(aux.base), 2, moy.p, row.w)
      aux.base <- t(t(as.matrix(aux.base)) - MM[[g]])
      ET[[g]] <- apply(as.data.frame(aux.base), 2, ec, row.w)
      aux.base <- t(t(as.matrix(aux.base)) / ET[[g]])
      missing <- which(is.na(as.matrix(aux.base)))
      if (any(is.na(aux.base))) aux.base[missing] <- 0
      ponderation[g] <- FactoMineR::svd.triplet(aux.base, ncp = 1, row.w = row.w)$vs[1]
      Xhat <- cbind.data.frame(Xhat, aux.base / ponderation[g])
      if (!is.null(seed) & (length(missing) != 0)) Xhat[missing, ] <- rnorm(length(missing))
    }
    if (type[g] == "c") {
      Xhat2 <- cbind.data.frame(Xhat2, aux.base)
      MM[[g]] <- apply(
        as.data.frame(aux.base), 2,
        moy.p, row.w
      )
      aux.base <- t(t(as.matrix(aux.base)) - MM[[g]])
      missing <- which(is.na(as.matrix(aux.base)))
      if (any(is.na(aux.base))) aux.base[missing] <- 0
      ponderation[g] <- FactoMineR::svd.triplet(aux.base, ncp = 1, row.w = row.w)$vs[1]
      Xhat <- cbind.data.frame(Xhat, aux.base / ponderation[g])
      if (!is.null(seed) & (length(missing) != 0)) Xhat[missing, ] <- rnorm(length(missing))
    }
    if (type[g] == "n") {
      tab.disj <- FactoMineR::tab.disjonctif.prop(data.frame(aux.base), seed, row.w = row.w)
      ####### Change: add the index of observation############
      obs_cat <- data.frame(!is.na(tab.disjonctif.NA(data.frame(aux.base))))
      idx_obs[[g]] <- which(apply(obs_cat, 1, any))
      #######################################################
      tab.disj.comp[[g]] <- tab.disj
      group.mod[g] <- ncol(tab.disj)
      MM[[g]] <- apply(tab.disj, 2, moy.p, row.w) / ncol(aux.base)


      ####### Change: equivalent code of the calculation############
      # Original code:
      # Z = t(t(tab.disj)/apply(tab.disj, 2, moy.p, row.w))
      # Z = t(t(Z)- apply(Z, 2, moy.p, row.w))
      # Zscale = t(t(Z) * sqrt(MM[[g]]))

      # The three lines above are equivalent to these two lines:
      Z <- t(t(tab.disj) / sqrt(apply(tab.disj, 2, moy.p, row.w)))
      Zscale <- t(t(Z) - sqrt(MM[[g]]))
      # Here we could use -MM[[g]] instead of -sqrt(MM[[g]])
      #######################################################

      ponderation[g] <- FactoMineR::svd.triplet(Zscale, row.w = row.w)$vs[1]
      Xhat <- cbind.data.frame(Xhat, Zscale / ponderation[g])
      Xhat2 <- cbind.data.frame(Xhat2, as.data.frame(tab.disjonctif.NA(aux.base)))
    }
  }
  ind.var <- vector(mode = "list", length = length(group))
  ind.var[[1]] <- 1:group.mod[1]
  for (g in 2:length(group)) ind.var[[g]] <- (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]
  fittedX <- Xhat <- as.matrix(Xhat)
  if (ncp >= min(nrow(Xhat) - 2, ncol(Xhat) - 1)) stop("ncp is too large")
  ncp <- min(ncp, ncol(X) - 1, nrow(X) - 2)
  missing <- which(is.na(as.matrix(Xhat2)))
  nb.iter <- 1
  old <- Inf

  while (nb.iter > 0) {
    # print(paste("nb.iter:",nb.iter))
    for (g in 1:length(group)) {
      if (g == 1) {
        aux.base <- Xhat[, 1:group.mod[1], drop = FALSE]
      } else {
        aux.base <- Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g], drop = FALSE]
      }
      aux.base <- aux.base * ponderation[g]
      if (type[g] == "s") {
        aux.base <- t((t(aux.base) * ET[[g]]) + MM[[g]])
        MM[[g]] <- apply(aux.base, 2, moy.p, row.w)
        aux.base <- t(t(aux.base) - MM[[g]])
        ET[[g]] <- apply(aux.base, 2, ec, row.w)
        aux.base <- t(t(aux.base) / ET[[g]])
        ponderation[g] <- FactoMineR::svd.triplet(aux.base, ncp = 1, row.w = row.w)$vs[1]
      }
      if (type[g] == "c") {
        aux.base <- t(t(aux.base) + MM[[g]])
        MM[[g]] <- apply(aux.base, 2, moy.p, row.w)
        aux.base <- t(t(aux.base) - MM[[g]])
        ponderation[g] <- FactoMineR::svd.triplet(aux.base, ncp = 1, row.w = row.w)$vs[1]
      }
      if (type[g] == "n") {


        ####### Change: equivalent code of the calculation############
        # Original code:
        # tab.disj = t(t(aux.base)/sqrt(MM[[g]])) + matrix(1, nrow(aux.base), ncol(aux.base))
        # tab.disj = t(t(tab.disj) * apply(tab.disj.comp[[g]], 2, moy.p, row.w))
        # The two lines above are equivalent to these two lines:
        tab.disj <- t(t(aux.base) + sqrt(MM[[g]]))
        tab.disj <- t(t(tab.disj) * sqrt(MM[[g]] * ncol(aux.base)))
        #######################################################


        ####### correct some computational error e-17######
        tab.disj[idx_obs[[g]], ] <- tab.disj.comp[[g]][idx_obs[[g]], ]
        ##################################################


        tab.disj.comp[[g]] <- tab.disj
        MM[[g]] <- apply(tab.disj, 2, moy.p, row.w) / ncol(aux.base)
        ############# Add threshold#######################
        minMM <- 1e-04 / ncol(aux.base)
        if (any(MM[[g]] < minMM)) {
          # Sometimes, one category could have a extremely small probability that is near 0
          # We could limit it to a threshold
          MM[[g]][which(MM[[g]] < minMM)] <- minMM
          # stop(paste("The algorithm fails to converge. Choose a number of components (ncp) less or equal than ",
          #            ncp - 1, " or a number of iterations (maxiter) less or equal than ",
          #            maxiter - 1, sep = ""))
        }
        ##################################################



        ####### Change: equivalent code of the calculation############
        # Original code:
        # Z = t(t(tab.disj)/apply(tab.disj, 2, moy.p, row.w))
        # Z = t(t(Z) - apply(Z, 2, moy.p, row.w))
        # aux.base = t(t(Z) * sqrt(MM[[g]]))
        # The three lines above are equivalent to these two lines:
        Z <- t(t(tab.disj) / sqrt(MM[[g]] * ncol(aux.base)))
        aux.base <- t(t(Z) - sqrt(MM[[g]]))
        ##################################################

        ponderation[g] <- FactoMineR::svd.triplet(aux.base, row.w = row.w, ncp = 1)$vs[1]
      }
      if (g == 1) {
        Xhat[, 1:group.mod[1]] <- aux.base / ponderation[g]
      } else {
        Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] <- aux.base / ponderation[g]
      }
    }
    if (!is.null(num.group.sup)) {
      for (g in num.group.sup) {
        if (g == 1) {
          Xhat[, 1:group.mod[1]] <- Xhat[, 1:group.mod[1]] * 1e-08
        } else {
          Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] <- Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] * 1e-08
        }
      }
    }

    svd.res <- FactoMineR::svd.triplet(Xhat, row.w = row.w, ncp = ncp)
    if (length(num.group.sup) > 0) {
      sigma2 <- mean(svd.res$vs[-c(1:ncp, (ncol(Xhat) - sum(group.mod[num.group.sup]) + 1):ncol(Xhat))]^2)
    } else {
      sigma2 <- mean(svd.res$vs[-c(1:ncp)]^2)
    }
    sigma2 <- min(sigma2 * coeff.ridge, svd.res$vs[ncp + 1]^2)
    if (method == "em") {
      sigma2 <- 0
    }
    lambda.shrinked <- (svd.res$vs[1:ncp]^2 - sigma2) / svd.res$vs[1:ncp]
    if (ncp == 1) {
      fittedX <- tcrossprod((svd.res$U[, 1, drop = FALSE] * row.w) * lambda.shrinked, svd.res$V[, 1, drop = FALSE])
    } else {
      fittedX <- tcrossprod(t(t(svd.res$U[, 1:ncp] * row.w) * lambda.shrinked), svd.res$V[, 1:ncp])
    }
    fittedX <- fittedX / row.w
    diff <- Xhat - fittedX
    diff[missing] <- 0
    objective <- sum(diff^2 * row.w)
    criterion <- abs(1 - objective / old)
    old <- objective
    nb.iter <- nb.iter + 1
    Xhat[missing] <- fittedX[missing]
    if (!is.null(num.group.sup)) {
      for (g in num.group.sup) {
        if (g == 1) {
          Xhat[, 1:group.mod[1]] <- Xhat[, 1:group.mod[1]] * 1e+08
        } else {
          Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] <- Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] * 1e+08
        }
      }
    }

    if (!is.nan(criterion)) {
      if ((criterion < threshold) && (nb.iter > 5)) nb.iter <- 0
      if ((objective < threshold) && (nb.iter > 5)) nb.iter <- 0
    }
    if (nb.iter > maxiter) {
      nb.iter <- 0
      warning(paste("Stopped after ", maxiter, " iterations\n"))
    }
  }



  # Xhat[missing] <- fittedX[missing]  ## A ajouter ?
  for (g in 1:length(group)) {
    if (g == 1) {
      aux.base <- Xhat[, 1:group.mod[1], drop = FALSE]
    } else {
      aux.base <- Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g], drop = FALSE]
    }
    aux.base <- aux.base * ponderation[g]
    if (type[g] == "s") {
      aux.base <- t(t(aux.base) * ET[[g]])
      aux.base <- t(t(aux.base) + MM[[g]])
    }
    if (type[g] == "c") aux.base <- sweep(aux.base, 2, MM[[g]], FUN = "+")
    if (type[g] == "n") {
      tab.disj <- t(t(aux.base) / sqrt(MM[[g]])) + matrix(1, nrow(aux.base), ncol(aux.base))
      aux.base <- t(t(tab.disj) * apply(tab.disj.comp[[g]], 2, moy.p, row.w))
    }
    if (g == 1) {
      Xhat[, 1:group.mod[1]] <- aux.base
    } else {
      Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] <- aux.base
    }
  }
  completeObs <- as.matrix(Xhat2)
  completeObs[missing] <- Xhat[missing]
  result <- list()
  result$tab.disj <- completeObs
  result$completeObs <- find.category(X, completeObs)
  result$call$group.mod <- group.mod
  result$call$ind.var <- ind.var
  return(result)
}


#' imputeMFA_mod: modified imputeMFA from missMDA
#' 
#' @description \code{imputeMFA_mod} is a modified version of the function \code{imputeMFA} by Francois Husson and Julie Josse.
#' Please find the detailed documentation of \code{imputeMFA} in the missMDA package. 
#' The modified version of \code{impute} ameliorates the convergence problem by modifying the function \code{impute}.
#' @export
imputeMFA_mod <- function(X, group, ncp = 2, type = rep("s", length(group)),
                        method = c("Regularized", "EM"), row.w = NULL, coeff.ridge = 1,
                        threshold = 1e-06, ind.sup = NULL, num.group.sup = NULL, seed = NULL, maxiter = 1000, ...) {
  # print("in imputeMFA")
  obj <- Inf
  method <- tolower(method)
  if (is.null(row.w)) row.w <- rep(1, nrow(X)) / nrow(X)
  if (length(ind.sup) > 0) row.w[ind.sup] <- row.w[ind.sup] * 1e-08
  if (!any(is.na(X))) stop("no missing values in X, this function is not useful. Perform MFA on X.")
  res.impute <- impute_mod(X,
    group = group, ncp = ncp, type = type,
    method = method, threshold = threshold, seed = seed,
    maxiter = maxiter, row.w = row.w, ind.sup = ind.sup,
    num.group.sup = num.group.sup, coeff.ridge = coeff.ridge
  )
  return(res.impute)
}
