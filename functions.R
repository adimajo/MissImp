
#Bootstrap
#Input: 1 incomplete dataset
#Output: num_sample incomplete dataset
#Suggest that num_sample > Log(n), where n is the number of rows
bootsample = function(df,num_sample){
  num_row = nrow(df)
  if(num_sample<log(num_row)){
    warning("Warning: We suggest that the bootstrap number > log(number of rows), if not, all the rows may not be covered in the bootstrap samples")
  }
  ls_df_new = list()
  i = 1
  while(i<=num_sample){
    chosen_idx = resample(c(1:num_row), num_row, replace = TRUE)
    df_new = df[chosen_idx,]
    #df_new[["old_idx"]]=chosen_idx
    ls_df_new[[i]] = df_new
    i = i + 1
  }
  return(ls_df_new)
}

#Jackknife
#Input: 1 incomplete dataset
#Output: num_sample incomplete dataset
jacksample = function(df, num_sample){
  num_row = nrow(df)
  row_tranch = num_row %/% num_sample
  rest = num_row %% num_sample
  ls_df_new = list()
  i = 1
  start = 1
  while(i <= num_sample-rest){ # num_sample-rest samples with row_tranch rows
    end = start + row_tranch -1
    if(start==1){
      df_new = df[(end+1):num_row,]
    }
    else if(end==num_row){
      df_new = df[1:(start-1),]
    }
    else{
      df_new = df[c(1:(start-1),(end+1):num_row),]
    }
    #df_new[["old_idx"]] = c(start:end)
    ls_df_new[[i]] = df_new
    i = i + 1
    start = end + 1
  }
  
  while(i <= num_sample){ # rest samples with row_tranch+1 rows
    end = start + (row_tranch + 1) -1
    df_new = df[c(1:(start-1),(end+1):num_row),]
    #df_new[["old_idx"]] = c(start:end)
    ls_df_new[[i]] = df_new
    i = i + 1
    start = end + 1
  }
  return(ls_df_new)
}


# The most frequent result for one categorical variable
# Used in combine_boot with method = 'factor'
Mode_cat = function(x) {
  ux = unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# Calculate the Wilcox's VarNC for a vector of categorical values
VA_fact=function(x){
  lev_miss = length(levels(x))-length(unique(x))
  freq = table(x)/length(x)
  return(VA(c(freq, rep(0,lev_miss))))
}

# Combine imputed bootstrap datasets
# Coded for onehot categorical variables and the factor encoded categorical variables
# Input: a list of imputed bootstrapped dataset
# Output: the final combined imputed dataset and the variance for each imputation
combine_boot = function(ls_df, col_con, col_dis, col_cat, num_row_origin, method='onehot',dict_cat=NULL, var_cat='unalike'){
  method = match.arg(method, c("onehot","factor"))
  var_cat = match.arg(var_cat, c("unalike","wilcox_va"))
  is_unalike = (var_cat=='unalike')
  is_onehot = (method=='onehot')
  if(is_onehot & is.null(dict_cat)){
    stop("dict_cat is needed when the method is onehot.")
  }
  ls_df_new = list()
  
  ls_col_name = colnames(ls_df[[1]])
  col_num = c(col_con,col_dis)
  exist_dis = !all(c(0,col_dis)==c(0))
  if(exist_dis){
    col_name_dis = ls_col_name[col_dis]
  }
  exist_cat =  !all(c(0,col_cat)==c(0))
  if(exist_cat){
    col_name_cat = ls_col_name[col_cat]
    if(is_onehot){
      col_name_cat = ls_col_name[col_cat]
      #We treat everything as numeric if the categorical variables are encoded by onehot probability
      col_num = c(col_num,col_cat)
    }
  }
  col_name_num = ls_col_name[col_num]
  
  # Deal with doublets in each imputed dataset
  i = 1
  for(df in ls_df){
    df$index <- as.numeric(row.names(df))
    df$index = floor(df$index)
    df_num = aggregate(.~index , data =df[c("index",col_name_num)], mean)
    if(exist_cat && !is_onehot){
      df_cat = aggregate(.~index , data =df[c("index",col_name_cat)], Mode_cat)
      ls_df_new[[i]] = merge(df_num,df_cat,by="index")
    }
    else{
      ls_df_new[[i]] = df_num
    }
    
    # the expected number of rows is n(1-(n-1)^k/n^k) where k = n
    # expectation of unique rows for sampling with replacement
    # We need a varaince matrix here to pass ? No, see https://stats.stackexchange.com/questions/399382/bootstrap-rubins-rules-and-uncertainty-of-sub-estimates
    
    i = i + 1
  }
  
  # Combine all imputed datasets together
  df_new_merge = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),ls_df_new)
  #Add the combined result for categorical variables deducted from the onehot result
  if(exist_cat && is_onehot){
    names_cat = names(dict_cat)
    which_max_cat = function(x, name){
      return(dict_cat[[name]][which.max(x)])
    }
    for(name in names_cat){
      df_new_merge[[name]]=apply(df_new_merge[dict_cat[[name]]],1,which_max_cat,name)
      df_new_merge[[name]]=unlist(df_new_merge[[name]])
      df_new_merge[[name]]=factor(df_new_merge[[name]])
    }
  }
  
  #By Bootstrap combining rules
  df_new_mean_num = aggregate(.~index , data =df_new_merge[c("index",col_name_num)], mean)
  df_new_num_var = aggregate(.~index , data =df_new_merge[c("index",col_name_num)], var)
  if(exist_dis){
    df_new_mean_num[col_name_dis] = round(df_new_mean_num[col_name_dis]) # for the discret variables
  }
  
  if(exist_cat && !is_onehot){
    df_new_merge=factor_encode(df_new_merge,col_cat+1)
    df_new_mode_cat = aggregate(.~index , data =df_new_merge[c("index",col_name_cat)], Mode_cat)
    df_new = merge(df_new_mean_num,df_new_mode_cat,by="index")
    df_new = factor_encode(df_new, col_cat+1) # +1 is for the index column
    if(is_unalike){
      df_new_cat_var = aggregate(.~index , data =df_new_merge[c("index",col_name_cat)], unalike)
    }
    else{
      df_new_cat_var = df_new_merge[c("index",col_name_cat)] %>% group_by(index) %>% summarise(across(col_name_cat, VA_fact))
      #df_new_cat_var = aggregate(.~index , data =df_new_merge[c("index",col_name_cat)], VA_fact)
    }
    df_new_var = merge(df_new_num_var,df_new_cat_var,by="index")
    # We need to use unalikeability to mesure the uncertainty of the categorical variables
  }
  else{
    df_new = df_new_mean_num
    df_new_var = df_new_num_var
  }
  #Add the combined result for categorical variables deducted from the onehot result to df_new
  if(is_onehot && exist_cat){
    for(name in names_cat){
      df_new[[name]]=apply(df_new[dict_cat[[name]]],1,which_max_cat,name)
      df_new[[name]]=unlist(df_new[[name]])
    }
    if(is_unalike){
      df_new_cat_var = aggregate(.~index , data =df_new_merge[c("index",names_cat)], unalike)
    }
    else{
      df_new_cat_var = df_new_merge[c("index",names_cat)] %>% group_by(index) %>% summarise(across(names_cat, VA_fact))
      #df_new_cat_var = aggregate(.~index , simplify=FALSE, data =df_new_merge[c("index",names_cat)], VA_fact)
    }
    df_new_var = merge(df_new_num_var,df_new_cat_var,by="index")
  }

  # Add the uncovered rows, fill in with NA
  num_row_new = nrow(df_new)
  if(num_row_new<num_row_origin){
    print(paste0("Covered number of rows: ", num_row_new))
    print(paste0("Original number of rows: ", num_row_origin))
    warning("Warning: All the rows in the original dataset is not selected in the bootstrap samples. An increase on the number of bootstrap samples is suggested.\n")
  }
  df_idx = data.frame("index"=c(1:num_row_origin))
  df_result_with_idx = merge(x = df_idx, y = df_new, by = "index", all = TRUE)
  df_result_var_with_idx = merge(x = df_idx, y = df_new_var, by = "index", all = TRUE)
  
  #remove index column and the combined categorical columns
  if(is_onehot && exist_cat){
    df_result_disj=df_result_with_idx[,!(names(df_result_with_idx) %in% c(names_cat,'index'))]
    df_result_var_disj=df_result_var_with_idx[,!(names(df_result_with_idx) %in% c(names_cat,'index'))]
  }
  else{
    df_result_disj=df_result_with_idx[-c(1)]
    df_result_var_disj=df_result_var_with_idx[-c(1)]
  }
  
 
  
  #Final result for the categorical variables
  if(is_onehot){
    df_result = df_result_with_idx[-c(1,col_cat+1)] #remove index and onehot columns, +1 is for the adjustment
    df_result_var = df_result_var_with_idx[-c(1,col_cat+1)]
  }
  else{
    df_result = df_result_disj
    df_result_var = df_result_var_disj
  }
  
  
  
  return(list("imp.disj"=df_result_disj,"uncertainty.disj"=df_result_var_disj, "imp"=df_result, "uncertainty"=df_result_var))
  
}


# Combine imputed jackknife datasets
# Only coded for onehot categorical variables
# Input: a list of imputed jackknifed dataset and a imputed whole dataset
# Output: the final combined imputed dataset and the variance for each imputation
combine_jack = function(ls_df, df_full, col_con, col_dis, col_cat, method='onehot',dict_cat=NULL,var_cat='unalike'){
  ls_df_minus = list()
  ls_col_name = colnames(ls_df[[1]])
  exist_dis = !all(c(0,col_dis)==c(0))
  if(exist_dis){
    col_name_dis = ls_col_name[col_dis]
  }
  exist_cat = !all(c(0,col_cat)==c(0))
  if(exist_cat){
    col_name_cat = ls_col_name[col_cat]
    var_cat = match.arg(var_cat, c("unalike","wilcox_va"))
    is_unalike = (var_cat=='unalike')
    if(is.null(dict_cat)){
      stop("dict_cat is needed when there are onehot categorical columns")
    }
  }
  i = 1
  n_sample = length(ls_df)
  df_full$index <- as.numeric(row.names(df_full))
  while(i<=n_sample){
    ls_df[[i]]$index <- as.numeric(row.names(ls_df[[i]]))
    #only for numeric columns
    ls_df_minus[[i]]= n_sample * df_full[which(df_full$index %in% ls_df[[i]]$index),] - (n_sample-1)*ls_df[[i]]
    ls_df_minus[[i]]$index = ls_df[[i]]$index
    i = i + 1
  }
  # Put all imputed datasets together
  df_new_merge = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),ls_df_minus)
  #Add the combined result for categorical variables deducted from the onehot result
  if(exist_cat){
    names_cat = names(dict_cat)
    which_max_cat = function(x, name){
      return(dict_cat[[name]][which.max(x)])
    }
    for(name in names_cat){
      df_new_merge[[name]]=apply(df_new_merge[dict_cat[[name]]],1,which_max_cat,name)
      df_new_merge[[name]]=unlist(df_new_merge[[name]])
      df_new_merge[[name]]=factor(df_new_merge[[name]])
    }
  }
 

  
  #By Jackknife combining rules for disjunctive part
  df_new = aggregate(.~index , data =df_new_merge[c("index",ls_col_name)], mean)# Here for each indice, there are n_sample-1 values
  if(exist_dis){
    df_new[col_name_dis] = round(df_new[col_name_dis]) # round for the discret variables
  }
  df_new_var_disj = aggregate(.~index , data =df_new_merge[c("index",ls_col_name)], var)
  df_new_var_disj[c(ls_col_name)] = df_new_var_disj[c(ls_col_name)]/(n_sample-1)
  
  
  #Final result according to the mean of onehot probability
  if(exist_cat){
    for(name in names_cat){
      df_new[[name]]=apply(df_new[dict_cat[[name]]],1,which_max_cat,name)
      df_new[[name]]=unlist(df_new[[name]])
      df_new[[name]]=factor(df_new[[name]])
    }
    df_result = df_new[-c(1,col_cat+1)] #remove the index and the onehot columns
    if(is_unalike){
      df_new_cat_var = aggregate(.~index , data =df_new_merge[c("index",names_cat)], unalike)
    }
    else{
      df_new_cat_var = df_new_merge[c("index",names_cat)] %>% group_by(index) %>% summarise(across(names_cat, VA_fact))
    }
    
    df_new_var = merge(df_new_var_disj,df_new_cat_var,by="index")
    df_result_var = df_new_var[-c(1,col_cat+1)] #remove the index and the onehot columns
    
    
    #remove index column and the combined categorical columns
    df_result_disj=df_new[,!(names(df_new) %in% c(names_cat,'index'))]
    df_result_var_disj=df_new_var[,!(names(df_new) %in% c(names_cat,'index'))]
    return(list("imp.disj"=df_result_disj,"uncertainty.disj"=df_result_var_disj,"imp"=df_result, "uncertainty"=df_result_var))
    
  }
  else{
    return(list("imp"=df_new[-c(1)], "uncertainty"=df_new_var_disj[-c(1)]))
    
  }
  
  
 }






# Calculate the missing proportion in MAR3, used in 'generate_miss'
# Solve (1-x)^p + (1-m)*p*x -1 = 0 in (0, 1)
# where m = miss_perc, p = num_co
monot_quantil = function(miss_perc,num_col){
  m = miss_perc
  p = num_col
  tmp_result = c()
  tempt = linspace(0.01,1,n=1000)
  i=1
  for(x in tempt){
    tmp_result[i]= abs((1-x)^p + (1-m)*p*x -1)
    i=i+1
  }
  return(tempt[which.min(tmp_result)])
  
}

# Generate missing values with different mechanisms
# Input: A complete dataset, the mechanism and the proportion of missingness
# Output: A dataset with missing values according to the mechanism and the proportion of missingness
generate_miss = function( df,
                          miss_perc,
                          mechanism="MCAR", # c("MCAR", "MAR1", "MAR2","MAR3","MNAR1","MNAR2")
                          # We could add here the parameters for produce_NA function in MAR1 and MNAR1
                          mar2.col.ctrl = 1){
  
  mechanism <- match.arg(mechanism, c("MCAR", "MAR1", "MAR2","MAR3","MNAR1","MNAR2"))
  ls_col_name = colnames(df)
  num_col = length(ls_col_name)
  if(mechanism=="MCAR"){
    # Bernouilli
    mcar = produce_NA(df, mechanism="MCAR", perc.missing = miss_perc) 
    X.mcar = mcar$data.incomp
    R.mcar = data.frame(mcar$idx_newNA)
    real_miss_perc = sum(R.mcar*1)/prod(dim(R.mcar*1))
    return(list("X.incomp"=X.mcar, "R.mask"=R.mcar, "real_miss_perc" = real_miss_perc))
  }
  else if(mechanism=="MAR1"){
    # Logistic regression to determinate the missingness
    # The options in this produce_NA function could be added to the main function
    mar1 = produce_NA(df, mechanism="MAR", perc.missing = miss_perc, by.patterns = F, logit.model = 'MID')
    X.mar1 = mar1$data.incomp
    R.mar1 = data.frame(mar1$idx_newNA)
    real_miss_perc = sum(R.mar1*1)/prod(dim(R.mar1*1))
    return(list("X.incomp"=X.mar1, "R.mask"=R.mar1, "real_miss_perc" = real_miss_perc))
  }
  else if(mechanism=="MAR2"){
    # Censoring algorithm. Everything depends on the quantile of one specified complete column.
    # For example, the Y2[i] will be removed if Y1[i]<q(30%) of Y1)
    # For the categorical variable cat, the quantile is taken on the levels(cat)
    idx_ctrl = mar2.col.ctrl
    X.mar2 = df
    for(coll in ls_col_name[1:num_col]){
      if(coll == ls_col_name[idx_ctrl]){
        next
      }
      X.mar2[,coll] = delete_MAR_censoring(X.mar2, miss_perc, 
                                           coll, cols_ctrl =ls_col_name[idx_ctrl])[,coll]
    }
    R.mar2 = data.frame(is.na(X.mar2))
    real_miss_perc = sum(R.mar2*1)/prod(dim(R.mar2*1))
    return(list("X.incomp"=X.mar2, "R.mask"=R.mar2, "real_miss_perc" = real_miss_perc))
  }
  else if(mechanism=="MAR3"){
    # Monotone with censoring mechanism, each column depends on the quantile of the observed data of the column before
    X.mar3 = df
    if(miss_perc *num_col >= (num_col-1)){
      stop("Error: MAR3 mechanism cannot work with this miss_perc")
    }
    perc = monot_quantil(miss_perc = miss_perc, num_col = num_col)
    i=1
    while (i < num_col) {
      ls_row = which(!is.na(X.mar3[,ls_col_name[i]]))
      # if(i != num_col-1){
        for(coll in ls_col_name[(i+1):num_col]){
          X.mar3[ls_row,coll] = delete_MAR_censoring(X.mar3[ls_row,], perc, coll, 
                                                     cols_ctrl = ls_col_name[i])[ls_row,coll]
        }
      # }
      # else{ # the last column adjust the missing percentage to approach the target missing percentage
      #   R = data.frame(is.na(X.mar3))
      #   p_adjust = (prod(dim(R*1)) * miss_perc - sum(R*1))/length(ls_row)
      #   X.mar3[ls_row,ls_col_name[num_col]] = delete_MAR_censoring(X.mar3[ls_row,], p_adjust, ls_col_name[num_col], 
      #                                              cols_ctrl = ls_col_name[i])[ls_row,ls_col_name[num_col]]
      # }
      
      i = i+1
    }
    R.mar3 = data.frame(is.na(X.mar3))
    real_miss_perc = sum(R.mar3*1)/prod(dim(R.mar3*1))
    return(list("X.incomp"=X.mar3, "R.mask"=R.mar3, "real_miss_perc" = real_miss_perc))
  }
  else if(mechanism=="MNAR1"){
    #logistic regression to determinate the missingness, with num_patt_mnar random patterns
    mnar1 = produce_NA(df, mechanism="MNAR", perc.missing = miss_perc, 
                       by.patterns= F,logit.model = 'LEFT')
    X.mnar1 = mnar1$data.incomp
    R.mnar1 = data.frame(mnar1$idx_newNA)
    real_miss_perc = sum(R.mnar1*1)/prod(dim(R.mnar1*1))
    return(list("X.incomp"=X.mnar1, "R.mask"=R.mnar1, "real_miss_perc" = real_miss_perc))
  }
  else if(mechanism=="MNAR2"){
    # A missing value in "X", if the x-value is below the miss_perc % quantile of "the first column "X"
    X.mnar2 = df
    for(coll in ls_col_name){
      X.mnar2[,coll] = delete_MNAR_censoring(X.mnar2, miss_perc, coll)[,coll]
    }
    R.mnar2 = data.frame(is.na(X.mnar2))
    real_miss_perc = sum(R.mnar2*1)/prod(dim(R.mnar2*1))
    return(list("X.incomp"=X.mnar2, "R.mask"=R.mnar2, "real_miss_perc" = real_miss_perc))
    
  }
}



# Generate a list of incomplete data with all the mechanisms that are in generate_miss function
generate_miss_ls = function(df, miss_perc){
  return(list("mcar" = generate_miss(df, miss_perc, mechanism = "MCAR"),
              "mar1" = generate_miss(df, miss_perc, mechanism = "MAR1"),
              "mar2" = generate_miss(df, miss_perc, mechanism = "MAR2"),
              "mar3" = generate_miss(df, miss_perc, mechanism = "MAR3"),
              "mnar1" = generate_miss(df, miss_perc, mechanism = "MNAR1"),
              "mnar2" = generate_miss(df, miss_perc, mechanism = "MNAR2")))
}


# Encoding functions
ordinal_encode = function(df, idx_col_cat){
  for(j in idx_col_cat){
    df[,j] = as.numeric(factor(df[,j], levels=levels(df[,j])))
  }
  return(df)
}

factor_encode = function(df, idx_col_cat){
  for(j in idx_col_cat){
    df[,j] = factor(df[,j])
  }
  return(df)
}

normalize_num = function(df, idx_col_num){
  df[,idx_col_num] = data.frame(apply(df[,idx_col_num],2, scale))
  return(df)
}





# Create the matrix of p-value for dummy t-chi-test
dummy_test_matrix = function(df, col_cat){
  test_result_dummy = data.frame()
  df = factor_encode(df,col_cat)  # in case the categorical column is not a factor
  R_df = data.frame(is.na(df))
  ls_col_name = colnames(df)
  ls_row_name = c()
  num_col = length(ls_col_name)
  i = 1
  while(i<=num_col){
    col_ctr = ls_col_name[i] # control column
    df_1 = df[R_df[[col_ctr]]==1,]
    df_0 = df[R_df[[col_ctr]]==0,]
    R_1 = R_df[R_df[[col_ctr]]==1,] # this matrix is needed for the t-test in a special case
    R_0 = R_df[R_df[[col_ctr]]==0,] 
    row1 = length(row.names(df_1))
    row0 = length(row.names(df_1))
    ls_row_name[i] = col_ctr
    if(row1==0 | row0==0){ # if one column contains only NA or doesn't contain any NA
      i = i + 1
      next
    }
    j = 1
    while(j <= num_col){
      if(i != j){
        col_test = ls_col_name[j]
        # if the test column is all NA, the test is not done on this column
        if(all(is.na(df[[col_test]]))){ 
          test_result_dummy[i,col_test] =NA
          j = j + 1
          next
        }
        
        # if when the control column is NA, the test column is all NA
        # or when the control colmun is observed , the test column is all NA
        # the paired t-test is done on the indicator matrix R_1 an R_0
        critc_1 = all(is.na(df_1[[col_test]]))
        critc_0 = all(is.na(df_0[[col_test]]))
        if(critc_1 || critc_0){ #the situation of critc_1 && critic_0 is discussed before
          if(critc_1 && !critc_0){
            test_result_dummy[i,col_test] = t.test(R_0[[col_test]]*1,mu=1)$p.value
          }
          else{
            test_result_dummy[i,col_test] = t.test(R_1[[col_test]]*1,mu=1)$p.value
          }
          
          j = j + 1
          next
        }
        
        # if the test column is categorical, we use chi 2 test
        if(j %in% col_cat){
          test_table = data.frame(col_test=df[[col_test]], R=R_df[[col_ctr]])
          test_result_dummy[i,col_test] = suppressWarnings(with(test_table, chisq.test(col_test,R))$p.value)
          # warnings about the approximation of normal distribution if the sample is small
          j = j + 1
          next
        }
        
        test_result_dummy[i,col_test] = t.test(df_1[[col_test]],df_0[[col_test]])$p.value
      }
      j = j + 1
    }
    i = i+1
  }
  row.names(test_result_dummy) = ls_row_name
  return(test_result_dummy)
}

# The combined p-value for dummy t-chi-test
# (Assume that the tests are independent: a correct assumption under MCAR?)
# Input: an incomplete dataset and the index of the categorical columns
# Output: a matrix of p-values of the dummy tests and a final p-value combined by Fisher's method
dummy_test = function(df, col_cat){
  p_matrix = dummy_test_matrix(df, col_cat)
  p_vector = as.vector(p_matrix)
  p_vector = p_vector[!is.na(p_vector)]
  dof = 2*length(p_vector)
  res = -2*sum(log(p_vector))
  p_dum = pchisq(res, df=dof, lower.tail=FALSE)
  return(list("p.matrix"=p_matrix, "dof"=dof, "chi2stat"=res, "p.value"=p_dum))
}



# There are density comparing functions for some of the imputation method, but here we try to define the function that could be used for every method
# Input: two datasets with same column names
# Output: the comparison of the density distribution for each column
dens_comp = function(df_comp,df_imp){
  ls_col_name = colnames(df_comp)
  ls_p = list()
  for( ycol in ls_col_name){
    dat = data.frame(y = c(df_comp[[ycol]],df_imp[[ycol]])
                     , lines = rep(c("complete", "imputed"), each = 100))
    p = ggplot(dat, aes(x = y, fill = lines)) + geom_density(alpha = 0.3) +labs(x = ycol) 
    show(p)
  }
}

# MSE ( average over missing number)
# As we did a scale operation on the complete dataset, we may need to rescale the output imputed data set then calculate the MSE
# Input: Original dataset, list of imputed dataset, mask of missingness, column index for numerical columns
# Output: A list of MSE values for each imputation result, a average MSE value, variance of MSE
ls_MSE = function(df_comp, ls_df_imp, mask, col_num) {
  ls_mse_result = c()
  mask_num = mask[,col_num]*1
  df_comp_num = df_comp[,col_num]
  i = 1
  for(df_imp in ls_df_imp){
    df_imp_num = df_imp[,col_num]
    mse_result = (sqrt(sum((as.matrix(df_comp_num) * mask_num - as.matrix(df_imp_num) * mask_num) ^ 2) / sum(mask_num)))
    ls_mse_result[i] = mse_result
    i = i + 1
  }
  return(list("list_MSE"=ls_mse_result,
              "Mean_MSE"=mean(ls_mse_result),
              "Variance_MSE"=var(ls_mse_result)))
}


# F1
# Input: Original dataset, list of imputed dataset, mask of missingness, column index for categorical columns
# Output: A list of F1-scores for each imputation result, a average F1-score, variance of F1-score
ls_F1 =function(df_comp, ls_df_imp, mask, col_cat){
  ls_f1_result = c()
  i = 1
  df_comp_concat = as.vector(as.matrix(df_comp[,col_cat]))
  mask_concat = as.vector(as.matrix(mask[,col_cat]*1))
  df_comp_true = df_comp_concat[mask_concat==1]
  for(df_imp in ls_df_imp){
    df_imp_concat =  as.vector(as.matrix(df_imp[,col_cat]))
    df_imp_predict = df_imp_concat[mask_concat==1]
    #micro f1 score : https://sebastianraschka.com/faq/docs/multiclass-metric.html 
    f1_result=F1_Score_micro(factor(df_comp_true),factor(df_imp_predict))
    ls_f1_result[i] = mean(f1_result)
    i = i + 1
  }
  return(list("list_F1"=ls_f1_result,
              "Mean_F1"=mean(ls_f1_result),
              "Variance_F1"=var(ls_f1_result)))
}














# generate_miss_ls = function(miss_perc, df){ 
# # all the columns in df requires to be quantitative. For the categorical variable, 
# # an encoder is needed before passing df into this function
#   num_patt_mar=5
#   num_patt_mnar=3
#   ls_col_name = colnames(df)
#   num_col = length(ls_col_name)
#   
#   
#   #MCAR (Bernouilli)
#   mcar = produce_NA(df, mechanism="MCAR", perc.missing = miss_perc) 
#   X.mcar = mcar$data.incomp
#   R.mcar = data.frame(mcar$idx_newNA)
#   
#   
#   #MAR1 (logistic regression to determinate the missingness, with num_patt patterns)
#   # n=0
#   # patt=c()
#   # while (n<num_patt_mar) { #generate a set of patterns
#   #   pa = rbinom(num_col, 1, 0.7)
#   #   if(!all(pa == 0)){ # pa is not 0,0,0,...,0
#   #     n = n + 1
#   #     patt=c(patt, pa)
#   #   }
#   # }
#   # freq_patt = rexp(num_patt_mar)
#   # we chose the exponential distribution so that the frequence for each pattern won't be so similar
#   # freq_patt = freq_patt/sum(freq_patt) #frequency of each pattern
#   # mar1 = produce_NA(df, mechanism="MAR", perc.missing = miss_perc, by.patterns = T, 
#   #                       patterns  = matrix(patt, ncol = num_col, byrow=T),
#   #                       freq.patterns = freq_patt,logit.model = 'MID')
#   mar1 = produce_NA(df, mechanism="MAR", perc.missing = miss_perc, by.patterns = F, logit.model = 'MID')
#   X.mar1 = mar1$data.incomp
#   R.mar1 = data.frame(mar1$idx_newNA)
#   
#   
#   #MAR2 (Censoring algorithm. Everything depends on the quantile of the first column, for example, the Y2[i] will be removed if Y1[i]<q(30%) of Y1)
#   #For the categorical variable cat, the quantile is taken on the levels(cat)
#   X.mar2 = df
#   
#   for(coll in ls_col_name[2:num_col]){
#     X.mar2[,coll] = delete_MAR_censoring(X.mar2, miss_perc, coll, cols_ctrl =ls_col_name[1])[,coll]
#   }
#   R.mar2 = data.frame(is.na(X.mar2))
#   
#   
#   
#   #MAR3 (monotone with censoring mechanisme, each column depends on the quantile of the observed data of the column before)
#   X.mar3 = df
#   if(miss_perc *num_col >= (num_col-1)){
#     stop("This mechanism cannot work with this miss_perc")
#   }
#   perc = monot_quantil(miss_perc = miss_perc, num_col = num_col)
#   i=1
#   while (i < num_col) {
#     ls_row = which(!is.na(X.mar3[,ls_col_name[i]]))
#     for(coll in ls_col_name[(i+1):num_col]){
#       X.mar3[ls_row,coll] = delete_MAR_censoring(X.mar3[ls_row,], perc, coll, 
#                                                cols_ctrl = ls_col_name[i])[ls_row,coll]
#     }
#     i = i+1
#   }
#   R.mar3 = data.frame(is.na(X.mar3))
#   
#   
#   
#   #MNAR1(logistic regression to determinate the missingness, with 3 patterns)
#   n=0
#   patt=c()
#   while (n<num_patt_mnar) {
#     pa = rbinom(num_col, 1, 0.5)
#     if(!all(pa == 0)){
#       n = n + 1
#       patt=c(patt, pa)
#     }
#   }
#   freq_patt = rexp(num_patt_mnar)
#   freq_patt = freq_patt/sum(freq_patt)
#   mnar1 = produce_NA(df, mechanism="MNAR", perc.missing = miss_perc, by.patterns= T,logit.model = 'LEFT',
#                     patterns  = matrix(patt, ncol = num_col, byrow=T),
#                    freq.patterns = freq_patt)
#   X.mnar1 = mnar1$data.incomp
#   R.mnar1 = data.frame(mnar1$idx_newNA)
#   
#   #MNAR2(A missing value in "X", if the x-value is below the miss_perc % quantile of "the first column "X")
#   X.mnar2 = df
#   for(coll in ls_col_name){
#     X.mnar2[,coll] = delete_MNAR_censoring(X.mnar2, miss_perc, coll)[,coll]
#   }
#   R.mnar2 = data.frame(is.na(X.mnar2))
#   
#   return(list("X.mcar"=X.mcar, "R.mcar"=R.mcar, "X.mar1"=X.mar1, "R.mar1"=R.mar1,
#               "X.mar2"=X.mar2, "R.mar2"=R.mar2, "X.mnar1"=X.mnar1, "R.mnar1"=R.mnar1,
#               "X.mnar2"=X.mnar2, "R.mnar2"=R.mnar2, "X.mar3"=X.mar3, "R.mar3"=R.mar3))
# }
