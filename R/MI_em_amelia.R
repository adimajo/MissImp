MI_em_amelia <- function(df_with_mv, col_num, col_cat=NULL, num_imp=5){
  exist_cat <- !all(c(0, col_cat) == c(0))
  if(exist_cat){
    dict_name_cat <- dict_onehot(df_with_mv, col_cat)
    # imputation
    imp_amelia <- Amelia::amelia(df_with_mv, m = num_imp, p2s = 0, noms = col_cat, boot.type='none')
    imp_amelia_disj <- list()
    i <- 1
    for(imp in imp_amelia$imputations){
      imp$index <- as.numeric(row.names(imp))
      imp_num <- imp[,col_num]
      dummy <- dummyVars(" ~ .", data = imp[,col_cat], sep = "_")
      imp_cat <- data.frame(predict(dummy, newdata = imp[,col_cat]))
      imp_num$index <- as.numeric(row.names(imp_num))
      imp_cat$index <- as.numeric(row.names(imp_cat))
      imp_amelia_disj[[i]] <- merge(x = imp_num, y = imp_cat, by = "index", all = TRUE)
      i <- i+1
    }
    imp_merge_disj <- Reduce(function(dtf1, dtf2) {
      rbind(dtf1, dtf2)
    }, imp_amelia_disj)
    
    ximp.disj <- stats::aggregate(. ~ index, data = imp_merge_disj, mean)
    ximp.all <- ximp.disj
    
    names_cat <- names(dict_name_cat)
    which_max_cat <- function(x, name) {
      return(dict_name_cat[[name]][which.max(x)])
    }
    
    for (name in names_cat) {
      ximp.all[[name]] <- apply(ximp.all[dict_name_cat[[name]]], 1, which_max_cat, name)
      ximp.all[[name]] <- unlist(ximp.all[[name]])
      ximp.all[[name]] <- factor(ximp.all[[name]])
    }
    ximp <- ximp.all[colnames(df_with_mv)]
    ximp.disj <- ximp.disj[-c(1)]
  }
  else{
    imp_amelia <- Amelia::amelia(df_with_mv, m = num_imp, p2s = 0, boot.type='none')
    imp_amelia_disj <- list()
    i <- 1
    for(imp in imp_amelia$imputations){
      imp$index <- as.numeric(row.names(imp))
      imp_amelia_disj[[i]] <- imp
      i <- i+1
    }
    
    imp_merge <- Reduce(function(dtf1, dtf2) {
      rbind(dtf1, dtf2)
    }, imp_amelia_disj)
    ximp <- stats::aggregate(. ~ index, data = imp_merge, mean)
    ximp.disj <- ximp
  }
  
  
  return(list(ximp=ximp,ximp.disj=ximp.disj, ls_ximp=imp_amelia$imputations, lx_ximp.disj = imp_amelia_disj))
}