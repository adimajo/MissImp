single.complete <- function(data, where, imp, ell) {
  if (ell == 0L) {
    return(data)
  }
  if (is.null(where)) {
    where <- is.na(data)
  }
  idx <- seq_len(ncol(data))[apply(where, 2, any)]
  for (j in idx) {
    if (is.null(imp[[j]])) {
      data[where[, j], j] <- NA
    } else {
      data[where[, j], j] <- imp[[j]][, ell]
    }
  }
  data
}

result_mice <- function(res, impnum, col_cat = c()) {
  exist_cat <- !all(c(0, col_cat) == c(0))
  if (exist_cat) {
    dict_cat <- dict_onehot(res$data, col_cat)
  }
  # Extract the imputations with categorical columns in onehot form
  res.disj <- list()
  idx <- c(1:imp_num)
  for (j in idx) {
    data <- single.complete(res$data, res$where, res$imp, j)
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
      ximp.all[[name]] <- apply(ximp.all[dict_cat[[name]]], 1, which_max_cat, name, dict_cat)
      ximp.all[[name]] <- unlist(ximp.all[[name]])
      ximp.all[[name]] <- factor(ximp.all[[name]])
    }
    ximp <- ximp.all[colnames(res$data)]
  } else {
    ximp <- ximp.disj
  }
  return(list(ximp = ximp, ximp.disj = ximp.disj))
}
