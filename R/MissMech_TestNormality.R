TestMCARNormality <- function(data, del.lesscases = 6, imputation.number = 1, method = "Auto", 
                              imputation.method = "Dist.Free", nrep = 10000, n.min = 30, 
                              seed = 110, alpha = 0.05, imputed.data = NA)
{
  if (!is.na(seed))
    set.seed(seed) 
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  } 
  if(!is.na(imputed.data[1]) && imputation.number!=1)
  {
    cat("Warning: No multiple imputation allowed when imputed data is provided.\n")
  }
  if(!is.matrix(data))
  {
    cat("Warning: Data is not a matrix or data frame.\n")
    stop("")
  }
  if(length(data)==0)
  {
    cat("Warning: Data is empty.\n")
    stop("")
  }
  if(ncol(data)<2)
  {
    cat("Warning: More than 1 variable is required.\n")
    stop("")
  }
  allempty <- which(apply(!is.na(data),1,sum) == 0)
  if (length(allempty) != 0) {
    data <- data[apply(!is.na(data), 1, sum) != 0, ]
    cat("Warning:", length(allempty), "Cases with all variables missing have been removed \n
        from the data.\n")
  }
  newdata <- OrderMissing(data, del.lesscases)
  if(length(newdata$data)==0)
  {
    cat("Warning: There are no data sets after deleting insufficient cases.\n")
    stop("")
  }
  
  if(newdata$g == 1)
  {
    cat("Warning: More than one missing data pattern should be present.\n")
    stop("")
  }
  if(sum(newdata$patcnt==1) > 0)
  {
    cat("Warning: At least 2 cases needed in each missing data patterns.\n")
    stop("")
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
  n4sim <- vector("list",g)
  #------------------------------imputation-----------------------
  mu <- matrix(0, p, 1)
  sig <- diag(1, p)
  
  emest <- Mls(newdata, mu, sig, 1e-6)
  mu <- emest$mu
  sig <- emest$sig
  if(is.na(imputed.data[1]))
  {
    yimp <- y
    if (imputation.method == "Dist.Free") {
      iscomp <- apply(patused, 1, sum, na.rm = TRUE) == p
      
      cind <- which(iscomp)
      ncomp <- patcnt[cind]
      if (length(ncomp) == 0) ncomp <- 0
      use.normal <- FALSE
      if (ncomp >= 10 && ncomp>=2*p){
        compy <- y[seq(spatcntz[cind] + 1, spatcntz[cind + 1]), ]
        ybar <- matrix(apply(compy, 2, mean))
        sbar <- cov(compy)
        resid <- (ncomp / (ncomp - 1)) ^ .5 * 
          (compy - matrix(ybar, ncomp, p, byrow = TRUE))
      } else {
        cat("Warning: There is not sufficient number of complete cases.\n  Dist.Free imputation requires a least 10 complete cases\n  or 2*number of variables, whichever is bigger.\n  imputation.method = normal will be used instead.\n")
        use.normal <- TRUE
      }
    }
    for(k in 1:imputation.number)
    {
      #-----------------normal imputation--------------------------------
      if (imputation.method == "Normal" || use.normal){
        yimp <- Impute(data = y, mu, sig, imputation.method = "Normal")
        yimp <- yimp$yimpOrdered
      }
      #-----------------distribution free imputation---------------------------------
      if (imputation.method == "Dist.Free" && !use.normal){
        yimp <- Impute(data = y, ybar, sbar, imputation.method = "Dist.Free", resid)
        yimp <- yimp$yimpOrdered
      }
      if (k == 1) yimptemp <- yimp
      #--------------Hawkin's test on the completed data------------------
      templist <- Hawkins(yimp,spatcnt)
      fij <- templist$fij
      tail <- templist$a
      ni <- templist$ni
      if (method == "Auto" || method == "Hawkins") {
        #Neyman test of uniformity for each group
        for(i in 1:g)
        {
          if (ni[i] < n.min){
            if (k == 1) {
              n4sim[[i]] <- SimNey(ni[i], nrep)
            }
          }
          templist <- TestUNey(tail[[i]], nrep, sim = n4sim[[i]], n.min)
          pn <- templist$pn
          n4 <- templist$n4
          pn <- pn + (pn == 0) / nrep
          pvalsn[k,i] <- pn
        }
      }
      #--------------Anderson darling test for equality of distribution
      if (method == "Auto" || method == "Nonparametric") {
        if(length(ni)<2)
        {
          cat("Warning: Not enough groups for AndersonDarling test.")
          stop("")
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
    templist <- Hawkins(yimp,spatcnt)
    fij <- templist$fij
    tail <- templist$a
    ni <- templist$ni
    if (method == "Auto" || method == "Hawkins") {
      #Neyman test of uniformity for each group
      for(i in 1:g)
      {
        if (ni[i] < n.min){
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
  adstar <- apply(adistar,1,sum)
  #combine p-values of test of uniformity
  combp <- -2 * apply(log(pvalsn), 1, sum)
  pvalcomb <- pchisq(combp, 2*g, lower.tail = FALSE)
  if (method == "Hawkins") {
    pnormality <- NULL
    adstar <- NULL
    adistar <- NULL
  }
  if (method == "Nonparametric") {
    pvalcomb = NULL
    combp = NULL
    pvalsn = NULL
  }
  yimptemp <- yimptemp[order(caseorder), ]
  if (length(removedcases) == 0) {
    dataused <- data
  }else {dataused <- data[-1 * removedcases, ]}
  homoscedastic <- list(analyzed.data = dataused, imputed.data = yimptemp,
                        ordered.data =  y, caseorder = caseorder,
                        pnormality = pnormality, adstar = adstar, adistar = adistar,  
                        pvalcomb = pvalcomb, combp = combp, pvalsn = pvalsn, g = g, alpha = alpha,
                        patused = patused, patcnt = patcnt, imputation.number = imputation.number, mu = mu, sigma = sig)
  homoscedastic$call <- match.call()
  class(homoscedastic) <- "testhomosc"
  homoscedastic
  }
#---------------------------------------------------------------------
#testmcar <- function(x, ...) UseImputationMethod("testmcar")
#testmcar.default <- function(data, ncases = 6, imputation.number = 10,
#                             imputation.method = "Normal", nrep = 10000)
#{
#test <- TestMCARNormality(data, ncases = 6, imputation.number = 10,
#                         imputation.method = "Normal", nrep = 10000)
#test$call <- match.call()
#class(test) <- "testmcar"
#test
#}
#---------------------------------------------------------------------
# printing format for the class "testhomosc"
print.testhomosc <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  #cat("\nNumber of imputation:\n")
  #print(x$imputation.number)
  ni <- x$patcnt
  cat("\nNumber of Patterns: ", x$g,"\n\nTotal number of cases used in the analysis: ", sum(ni),"\n")
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
    cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", x$pvalcomb[1],"\n")
    if (x$pvalcomb[1] > alpha){
      cat("\n    There is not sufficient evidence to reject normality
          or MCAR at", alpha,"significance level\n")
    }else {
      cat("\n    Either the test of multivariate normality or homoscedasticity (or both) is rejected.\n    Provided that normality can be assumed, the hypothesis of MCAR is 
          rejected at",alpha,"significance level. \n")
      cat("\nNon-Parametric Test:\n")
      cat("\n    P-value for the non-parametric test of homoscedasticity: ", x$pnormality[1],"\n")
      if (x$pnormality[1] > alpha){
        cat("\n    Reject Normality at",alpha,"significance level.
            There is not sufficient evidence to reject MCAR at",alpha,"significance level.\n")
      }else {
        cat("\n    Hypothesis of MCAR is rejected at ",alpha,"significance level.
            The multivariate normality test is inconclusive. \n")
      }
      }
    }
  if (method == "Hawkins"){
    cat("\nHawkins Test:\n")
    cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", x$pvalcomb[1],"\n")
  }
  if (method == "Nonparametric"){
    cat("\nNon-Parametric Test:\n")
    cat("\n    P-value for the non-parametric test of homoscedasticity: ", x$pnormality[1],"\n")
  }
    }
#----------------------------------------------------------------------------
summary.testhomosc <- function(object, ...) {
  ni <- object$patcnt
  cat("\nNumber of imputation: ", object$imputation.number,"\n")
  cat("\nNumber of Patterns: ", object$g,"\n\nTotal number of cases used in the analysis: ", sum(ni),"\n")
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
    cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", object$pvalcomb[1],"\n")
    cat("\nNon-Parametric Test:\n")
    cat("\n    P-value for the non-parametric test of homoscedasticity: ", object$pnormality[1],"\n")
  }
  if (method == "Hawkins"){
    cat("\nHawkins Test:\n")
    cat("\n    P-value for the Hawkins test of normality and homoscedasticity: ", object$pvalcomb[1],"\n")
  }
  if (method == "Nonparametric"){
    cat("\nNon-Parametric Test:\n")
    cat("\n    P-value for the non-parametric test of homoscedasticity: ", object$pnormality[1],"\n")
  }
}
#-----------------------------------------------------------------------------
# Plot "testhomosc"
boxplot.testhomosc <- function(x, ...) {
  if (is.null(x$pnormality)) {
    par(bg = "cornsilk")
    boxplot(x$pvalsn, col="lightcyan", border = "blue", medlwd = .5, medcol = "red")
    title(main = "Boxplots of p-values corresponding to each set of the missing data patterns\n for the Neyman test of Uniformity",
          xlab = "Missing data pattern group", ylab = "P-value", font.main = 4, 
          col.main = "blue4", cex.main = 1, font.lab = 4, cex.lab = 0.8, 
          col.lab = "blue4")
    abline(h = x$alpha / x$g, col = "red", lty = 2)
  }
  if (is.null(x$pvalcomb)) {
    par(bg = "cornsilk")
    boxplot(x$adistar, col="lightcyan", border = "blue", medlwd = .5, medcol = "red")
    title(main = "Boxplots of the T-value test statistics corresponding to each set of missing\n data patterns for the non-parametric test",
          xlab = "Missing data pattern group", ylab = expression(T[i]), 
          font.main = 4, col.main = "blue4", cex.main = 1, font.lab = 4, 
          cex.lab = 0.8, col.lab = "blue4")
  }
  if (!is.null(x$pvalcomb) && !is.null(x$pnormality)) {
    par(mfrow=c(2,1), bg = "cornsilk")
    boxplot(x$pvalsn, col="lightcyan", border = "blue", medlwd = .5, medcol = "red")
    title(main = "Boxplots of p-values corresponding to each set of the missing data patterns\n for the Neyman test of Uniformity",
          xlab = "Missing data pattern group", ylab = "P-value", font.main = 4, 
          col.main = "blue4", cex.main = 1, font.lab = 4, cex.lab = 0.8, 
          col.lab = "blue4")
    abline(h = x$alpha / x$g, col = "red", lty = 2)
    boxplot(x$adistar, col="lightcyan", border = "blue", medlwd = .5, medcol = "red")
    title(main = "Boxplots of the T-value test statistics corresponding to each set of missing\n data patterns for the non-parametric test",
          xlab = "Missing data pattern group", ylab = expression(T[i]), 
          font.main = 4, col.main = "blue4", cex.main = 1, font.lab = 4, 
          cex.lab = 0.8, col.lab = "blue4")
  }
}

TestUNey <- function(x, nrep = 10000, sim = NA, n.min = 30)
{
  # This routine tests whether the values in each row of x are unif(0,1). It
  # uses the Neyman's smooth test (see e.g., Ledwina 1994, TAS)
  # x is a vector
  # P-values are computed based on a
  # resampling method from unif(0,1).
  # All values of $x$ are between 0 and 1
  
  n <- length(x)
  pi <- LegNorm(x)
  n4 <- (apply(pi$p1, 2, sum) ^ 2 + apply(pi$p2, 2, sum) ^ 2 + 
           apply(pi$p3, 2, sum) ^ 2 + apply(pi$p4, 2, sum) ^ 2) / n 
  if (n < n.min){
    if(is.na(sim[1])) {
      sim <- SimNey(n, nrep)
    }
    pn <- length(which(sim > n4)) / nrep
  } else {
    pn <- pchisq(n4, 4, lower.tail = FALSE)
  }
  list(pn = pn, n4 = n4)
}
SimNey <- function(n, nrep)
{
  x <- matrix(runif(nrep * n), ncol = nrep)
  pi <- LegNorm(x)
  n4sim <- (apply(pi$p1, 2, sum) ^ 2 + apply(pi$p2, 2, sum) ^ 2 + 
              apply(pi$p3, 2, sum) ^ 2 + apply(pi$p4, 2, sum) ^ 2) / n 
  n4sim <- sort(n4sim)
}

AndersonDarling <- function(data, number.cases) { 
  # Not adjusted for ties
  # x is data vector
  # ni is number of cases in each group
  x <- data
  ni <- number.cases
  if(length(ni)<2)
  {
    cat("Warning: Not enough groups for AndersonDarling test.")
    stop("")
  }
  
  k <- length(ni)
  ni.z <- c(0, cumsum(ni))
  n <- length(x)
  x.sort <- sort(x)
  x.sort <- x.sort[1:(n-1)]
  ind <- which(duplicated(x.sort) == 0)
  counts <- c(ind, length(x.sort) + 1) - c(0, ind)
  hj <- counts[2 : (length(ind) + 1)]
  hn <- cumsum(hj)
  zj <- x.sort[ind]
  adk <- 0
  adk.all <- matrix(0, k, 1) # to keep contribution of kth group
  for(i in 1:k)
  {
    ind <- (ni.z[i] + 1) : ni.z[i + 1]
    templist <- expand.grid(zj, x[ind])
    b <- templist[, 1] == templist[, 2]
    fij <- apply(matrix(b, length(zj)), 1, sum)
    mij <- cumsum(fij)
    num <- (n * mij - ni[i] * hn) ^ 2
    dem <- hn*(n - hn)
    adk.all[i] <- (1 / ni[i] * sum(hj * (num / dem)))
    adk <- adk + adk.all[i]
  }
  adk <- (1 / n) * adk
  adk.all <- adk.all / n
  #Exact sample variance of the k-sample Anderson-Darling
  #Finding Variance of the statistics
  j <- sum(1 / ni)
  i <- seq(1:(n - 1))
  h <- sum(1 / i)
  g <- 0
  for (i in 1:(n - 2)) {
    g <- g + (1 / (n - i)) * sum(1 / seq((i + 1), (n - 1)))
  }
  a <- (4 * g - 6) * (k - 1) + (10 - 6 * g) * j
  b <- (2 * g - 4) * k ^ 2 + 8 * h * k + 
    (2 * g - 14 * h - 4) * j - 8 * h + 4 * g - 6
  c <- (6 * h + 2 * g - 2) * k ^ 2 + (4 * h - 4 * g + 6) * k +
    (2 * h - 6) * j + 4 * h
  d <- (2 * h + 6) * k ^ 2 - 4 * h * k
  var.adk <- ((a * n ^ 3) + (b * n ^ 2) + (c * n) + d) /
    ((n - 1) * (n - 2) * (n - 3))
  if(var.adk<0) var.adk=0
  adk.s <- (adk - (k - 1)) / sqrt(var.adk) 
  #k-sample Anderson-Darling P-value calculation by an extrapolate-interpolate
  #procedure
  a0 <- c(0.25, 0.10, 0.05, 0.025, 0.01)
  b0 <- c(0.675, 1.281, 1.645, 1.96, 2.326)
  b1 <- c(-0.245, 0.25, 0.678 ,1.149, 1.822)
  b2 <- c(-0.105, -0.305, -0.362, -0.391, -0.396)
  c0 <- log((1 - a0) / a0)
  qnt <- b0 + b1 / sqrt(k - 1) + b2 / (k - 1)
  if (adk.s <= qnt[3]) {
    ind <- seq(1:4)
  } else {
    ind <- seq(2:5)
  }
  yy <- spline(qnt[ind], c0[ind], xout = adk.s)$y
  p <- 1 / (1 + exp(yy))
  list(pn = p, adk.all = adk.all, adk = adk, var.sdk = var.adk)
}#end function

OrderMissing <- function(data, del.lesscases = 0)
{
  # case order has the order of the original data
  y <- data
  if (is.data.frame(y)) {
    y <- as.matrix(y)
  } 
  if(!is.matrix(y))
  {
    cat("Warning data is not a matrix or data frame")
    stop("")
  }
  if(length(y)==0)
  {
    cat("Warning: data is empty")
    return
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
  while(!done)
  {
    pattemp <- is.na(y[1, ])
    indin <- c()
    indout <- c()
    done <- TRUE
    for(i in 1:ntemp)
    {
      if(all(is.na(y[i, ]) == pattemp))
      {
        indout <- c(indout, i)
      }  else {
        indin <- c(indin, i)
        done <- FALSE
      }
    }
    if(length(indin) == 1) yatone = TRUE
    yfinal <- rbind(yfinal, y[indout, ])
    y <- y[indin, ]
    caseorder <- c(caseorder, ordertemp[indout])
    ordertemp <- ordertemp[indin]
    patcnt <- c(patcnt, length(indout))
    patused <- rbind(patused, pattemp)
    if(yatone)
    {
      pattemp <- is.na(y)
      yfinal <- rbind(yfinal, matrix(y,ncol=pp))
      y <- c()
      indin <- c()
      indout <- c(1)
      caseorder <- c(caseorder, ordertemp[indout])
      ordertemp <- ordertemp[indin]
      patcnt <- c(patcnt, length(indout))
      patused <- rbind(patused, pattemp)
      done <- TRUE
    }
    
    if(!done) ntemp <- nrow(y)
  }
  #yfinal <- rbind(yfinal, y)
  caseorder <- c(caseorder, ordertemp)
  patused <- ifelse(patused, NA, 1)
  rownames(patused) <- NULL
  colnames(patused) <- names
  spatcnt <- cumsum(patcnt)
  dataorder <- list(data = yfinal, patused = patused, patcnt = patcnt, 
                    spatcnt = spatcnt, g = length(patcnt), 
                    caseorder = caseorder, removedcases = removedcases)
  dataorder$call <- match.call()
  class(dataorder) <- "orderpattern"
  if(del.lesscases > 0)
  {
    dataorder <- DelLessData(dataorder, del.lesscases)
  }   
  dataorder$patused <- matrix(dataorder$patused, ncol = pp)
  colnames(dataorder$patused) <- names
  dataorder
}
#Order <- function(x, ...) UseMethod("Order")
#Order.default <- function(x, ...) {
# temp <- OrderMissing(x)
# temp$call <- match.call()
# class(temp) <- "ordered"
# temp
#}
print.orderpattern <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nNumber of Ptterns: ", x$g, "\n")
  cat("\nPttern used:\n")
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
  #This function deletes cases of a missing pattern with less than or equal to ncases 
  if(length(data)==0)
  {
    cat("Warning: data is empty")
    return
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
  if(length(ind) != 0){
    #cat("cases with insufficient number of observations were removed")
    for(i in 1:length(ind))
    {
      rm <- c(rm, seq(spatcntz[ind[i]] + 1, spatcntz[ind[i] + 1]));
    }
    y <- data$data[-1 * rm, ]
    removedcases <- data$caseorder[rm]
    patused <- data$patused[-1 * ind, ]
    patcnt <- data$patcnt[-1 * ind]
    caseorder <- data$caseorder[-1 * rm]
    spatcnt <- cumsum(patcnt)
  }else {
    patused <- data$patused
    patcnt <- data$patcnt
    spatcnt <- data$spatcnt
    caseorder <- data$caseorder
    y <- data$data
  }
  
  newdata <- list(data = y, patused = patused, patcnt = patcnt,
                  spatcnt = spatcnt, g = length(patcnt), caseorder = caseorder, 
                  removedcases = removedcases)
  #  colnames(newdata)<-colnames(data)
  class(newdata) <- "orderpattern"
  newdata
}

Mls  <- function(data, mu = NA, sig = NA, tol = 1e-6, Hessian = FALSE)
{
  # mu is estimate of the mean
  # sig is estimate of the covariance
  if (!is.matrix(data) && class(data) != "orderpattern") {
    cat("Warning: data must have the classes of matrix or orderpattern.\n")
    stop("")
  }
  if (is.matrix(data)) {
    allempty <- which(apply(!is.na(data),1,sum) == 0)
    if (length(allempty) != 0) {
      data <- data[apply(!is.na(data), 1, sum) != 0, ]
      cat("Warning:", length(allempty), "Cases with all variables missing have been removed
          from the data.\n")
    }
    data <- OrderMissing(data)
    }
  if (class(data) == "orderpattern") {
    allempty <- which(apply(!is.na(data$data),1,sum) == 0)
    if (length(allempty) != 0) {
      data <- data$data
      data <- data[apply(!is.na(data), 1, sum) != 0, ]
      cat("Warning:", length(allempty), "Cases with all variables missing have been removed
          from the data.\n")
      data <- OrderMissing(data)
    }
  }
  if (length(data$data)==0)
  {
    cat("Warning: Data is empty")
    stop("")
  }
  if(ncol(data$data)<2)
  {
    cat("Warning: More than 1 variable is required.\n")
    stop("")
  }
  y <- data$data
  patused <- data$patused
  spatcnt <- data$spatcnt
  if (is.na(mu[1])){
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
    if(!(em > tol || itcnt < 2)) break()
  }#end repeat
  rownames(mu) <- colnames(y)
  colnames(sig) <- colnames(y)
  if(Hessian)
  {
    templist <- Ddf(y,mu,sig)
    hessian <- templist$dd
    stderror <- templist$se
    return (list(mu = mu, sig = sig, hessian = hessian, stderror = stderror, iteration = itcnt))
  }
  return(list(mu = mu, sig = sig, iteration = itcnt))
  }
#------------------------------------------------------------------
Sexpect <- function(y, mu, sig, patused, spatcnt)
{
  
  n <-  nrow(y)
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
    if (length(indm)!= 0) {
      sigmo <- matrix(sig[indm, indo], length(indm), length(indo))
      sigmm <- sig[indm, indm]
      temp1 <- matrix(mum, ni, length(indm), byrow = TRUE)
      temp2 <- yo - matrix(muo, ni, length(indo), byrow = TRUE)
      ym <- temp1 + temp2 %*% sigooi %*% t(sigmo)
      som <- t(yo) %*% ym
      smm <- ni * (sigmm - sigmo %*% sigooi %*% t(sigmo))+ t(ym)%*%ym
      stemp[indo, indm] <- som
      stemp[indm, indo] <- t(som)
      stemp[indm, indm] <- smm
      ystemp[, indm] <- ym
    }# end if 
    sstar <- sstar + stemp;
    if (ni == 1){
      ysbar <- t(ystemp) + ysbar
    }else { 
      ysbar <- apply(ystemp, 2, sum) + ysbar
    }
  }#end for
  ysbar <- (1 / n) * ysbar
  sstar <- (1 / n) * sstar
  sstar <- (sstar + t(sstar))/2
  
  return(list(ysbar = ysbar, sstar = sstar))
}

Impute <- function(data, mu = NA, sig = NA, imputation.method = "Normal", resid = NA) 
{ # Check if data is not ordered change it to ordered form
  if (!is.matrix(data) && class(data) != "orderpattern") {
    cat("Warning: data must have the classes of matrix or orderpattern.\n")
    stop("")
  }
  if (is.matrix(data)) {
    allempty <- which(apply(!is.na(data),1,sum) == 0)
    if (length(allempty) != 0) {
      data <- data[apply(!is.na(data), 1, sum) != 0, ]
      cat("Warning:", length(allempty), "Cases with all variables missing have been removed
          from the data.\n")
    }
    data <- OrderMissing(data)
    }
  if (class(data) == "orderpattern") {
    allempty <- which(apply(!is.na(data$data),1,sum) == 0)
    if (length(allempty) != 0) {
      data <- data$data
      data <- data[apply(!is.na(data), 1, sum) != 0, ]
      cat("Warning:", length(allempty), "Cases with all variables missing have been removed
          from the data.\n")
      data <- OrderMissing(data)
    }
  }
  if(length(data$data)==0)
  {
    cat("Warning: data is empty")
    return
  }
  if(ncol(data$data)<2)
  {
    cat("More than 1 variable is required.\n")
    stop("")
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
      if (ncomp >= 10 && ncomp>=2*p){
        compy <- y[seq(spatcntz[cind] + 1, spatcntz[cind + 1]), ]
        ybar <- matrix(apply(compy, 2, mean))
        sbar <- cov(compy)
        if (is.na(resid[1])){ 
          resid <- (ncomp / (ncomp - 1)) ^ .5 * 
            (compy - matrix(ybar, ncomp, p, byrow = TRUE))
        }
      } else {
        cat("Warning: There is not sufficient number of complete cases.\n  Dist.free imputation requires a least 10 complete cases\n  or 2*number of variables, whichever is bigger.\n")
        return
      }
    }
    if (!is.na(mu[1])) {
      ybar <- mu
      sbar <- sig
      iscomp <-  apply(patused, 1, sum, na.rm = TRUE) == p
      cind <- which(iscomp)
      ncomp <- patcnt[cind]
      use.normal <- FALSE
      compy <- y[seq(spatcntz[cind] + 1, spatcntz[cind + 1]), ]
      if (is.na(resid[1])){ 
        resid <- (ncomp / (ncomp - 1)) ^ .5 * 
          (compy - matrix(ybar, ncomp, p, byrow = TRUE))
      }
    }
    indsample <- sample(ncomp, n - ncomp, replace = TRUE)
    resstar <- resid[indsample, ]
    indres1 <- 1
    for (i in 1:g) {
      if (sum(patused[i, ], na.rm = TRUE) != p) # choose a pattern not completely obsered 
      {  test <- y[(spatcntz[i] + 1) : spatcntz[i + 1], ] 
      indres2 <- indres1 + patcnt[i] - 1
      test <- MimputeS(matrix(test,ncol=p), patused[i, ], ybar, sbar,
                       matrix(resstar[indres1:indres2, ],ncol=p))
      indres1 <- indres2 + 1
      
      yimp[(spatcntz[i] + 1) : spatcntz[i + 1], ] <- test
      }
    }#end loop
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
      {  test <- y[(spatcntz[i] + 1) : spatcntz[i + 1], ] 
      test <- Mimpute(matrix(test,ncol=p), patused[i, ], mu, sig)
      yimp[(spatcntz[i] + 1) : spatcntz[i + 1], ] <- test
      }
    }#end loop
  }
  yimpord <- yimp
  yimp <- yimp[order(caseorder), ]
  imputed <- list(yimp = yimp, yimpOrdered = yimpord, caseorder = caseorder, 
                  patused = patused, patcnt = patcnt)
  
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
  data[, indm]= matrix(rnorm(ni * pm), ni, pm) %*% a + expymiss
  data
}
#---------------------------------------------------------------------------
MimputeS <- function(data, patused, y1, s1, e)
{
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
  dif <-  data[, indo] - matrix(y1[indo], ni, po, byrow = TRUE)
  z <- matrix(y1[indm], ni, pm, byrow = TRUE) + dif %*% t(a)
  etta <- matrix(e[, indm], ni, pm) - matrix(e[, indo], ni, po) %*% t(a)
  zij <- z + etta
  data[, indm] <- zij 
  data
}


Hawkins <- function(data, spatcnt)
{
  # This function performs the Hawkin's method for testing equality of
  # covariances among groups. It is assumed that the data in y is ordered so
  # that the cases from the same group are adjacent.
  #gind is a vector indicating the end of the cases for each group. For
  #example [8 23 30] means that there are three groups with cases 1-8 in one
  #group, 9-23 and 24-30 in another.
  
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
  for(i in 1:g)
  {
    yg <- y[seq(gind[i] + 1, gind[i + 1]), ]
    ni[i] <- nrow(yg)
    spool <- spool + (ni[i] - 1) * cov(yg)
    ygmean <- apply(yg, 2, mean)
    ygc[seq(gind[i]+1, gind[i + 1]), ] <-
      yg - matrix(ygmean, ni[i], p, byrow = TRUE)
  }
  spool <- spool / (n - g)
  spool <- solve(spool)
  f <- matrix(0, n, 1)
  nu <- n - g - 1
  a <- vector("list",g)
  for(i in 1:g)
  {
    vij <- ygc [seq(gind[i] + 1, gind[i + 1]), ]
    vij <- apply(vij %*% spool * vij, 1, sum)
    vij <-  vij*ni[i]
    f[seq(gind[i] + 1, gind[i + 1])] <- ((n - g - p) * vij)/
      (p * ((ni[i] -1 ) * (n - g) - vij))
    a[[i]] <- 1 - pf(f[seq(gind[i] + 1, gind[i + 1])], p, (nu - p + 1))
  }
  list(fij = f, a = a, ni = ni)
} 


LegNorm <- function(x)
{
  # This function evaluates Legendre polynomials on [0,1] (note not [-1,1] at
  # a value x, and the polynomial are such that they have norm 1. It only
  # calculates p1,p2,p3, and p4 [Note p0=1, so it is not reurned]
  # x can be a vector or a matrix. 
  # in return, each element of P1, p2, p3, p4 the values of the poly at each
  # component of x
  
  if (!is.matrix(x))
  {
    x <- as.matrix(x)
  }#end if
  x <- 2 * x - 1 # Transforming the legenre's to 0,1, and we will divide by the norm in each case
  p0 <- matrix(1,nrow(x), ncol(x))
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
