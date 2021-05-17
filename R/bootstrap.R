#Bootstrap
#Input: 1 incomplete dataset
#Output: num_sample incomplete dataset
#Suggest that num_sample > Log(n), where n is the number of rows
bootsample = function(df, num_sample) {
  num_row = nrow(df)
  if (num_sample < log(num_row)) {
    warning(
      "Warning: We suggest that the bootstrap number > log(number of rows), if not, all the rows may not be covered in the bootstrap samples"
    )
  }
  ls_df_new = list()
  i = 1
  while (i <= num_sample) {
    chosen_idx = resample(c(1:num_row), num_row, replace = TRUE)
    df_new = df[chosen_idx, ]
    #df_new[["old_idx"]]=chosen_idx
    ls_df_new[[i]] = df_new
    i = i + 1
  }
  return(ls_df_new)
}


# Combine imputed bootstrap datasets
# Coded for onehot categorical variables and the factor encoded categorical variables
# Input: a list of imputed bootstrapped dataset
# Output: the final combined imputed dataset and the variance for each imputation
combine_boot = function(ls_df,
                        col_con,
                        col_dis,
                        col_cat,
                        num_row_origin,
                        method = 'onehot',
                        dict_cat = NULL,
                        var_cat = 'unalike') {
  method = match.arg(method, c("onehot", "factor"))
  var_cat = match.arg(var_cat, c("unalike", "wilcox_va"))
  is_unalike = (var_cat == 'unalike')
  is_onehot = (method == 'onehot')
  if (is_onehot & is.null(dict_cat)) {
    stop("dict_cat is needed when the method is onehot.")
  }
  ls_df_new = list()
  
  ls_col_name = colnames(ls_df[[1]])
  col_num = c(col_con, col_dis)
  exist_dis = !all(c(0, col_dis) == c(0))
  if (exist_dis) {
    col_name_dis = ls_col_name[col_dis]
  }
  exist_cat =  !all(c(0, col_cat) == c(0))
  if (exist_cat) {
    col_name_cat = ls_col_name[col_cat]
    if (is_onehot) {
      col_name_cat = ls_col_name[col_cat]
      #We treat everything as numeric if the categorical variables are encoded by onehot probability
      col_num = c(col_num, col_cat)
    }
  }
  col_name_num = ls_col_name[col_num]
  
  # Deal with doublets in each imputed dataset
  i = 1
  for (df in ls_df) {
    df$index <- as.numeric(row.names(df))
    df$index = floor(df$index)
    df_num = aggregate(. ~ index , data = df[c("index", col_name_num)], mean)
    if (exist_cat && !is_onehot) {
      df_cat = aggregate(. ~ index , data = df[c("index", col_name_cat)], Mode_cat)
      ls_df_new[[i]] = merge(df_num, df_cat, by = "index")
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
  df_new_merge = Reduce(function(dtf1, dtf2)
    rbind(dtf1, dtf2), ls_df_new)
  #Add the combined result for categorical variables deducted from the onehot result
  if (exist_cat && is_onehot) {
    names_cat = names(dict_cat)
    which_max_cat = function(x, name) {
      return(dict_cat[[name]][which.max(x)])
    }
    for (name in names_cat) {
      df_new_merge[[name]] = apply(df_new_merge[dict_cat[[name]]], 1, which_max_cat, name)
      df_new_merge[[name]] = unlist(df_new_merge[[name]])
      df_new_merge[[name]] = factor(df_new_merge[[name]])
    }
  }
  
  #By Bootstrap combining rules
  df_new_mean_num = aggregate(. ~ index , data = df_new_merge[c("index", col_name_num)], mean)
  df_new_num_var = aggregate(. ~ index , data = df_new_merge[c("index", col_name_num)], var)
  if (exist_dis) {
    df_new_mean_num[col_name_dis] = round(df_new_mean_num[col_name_dis]) # for the discret variables
  }
  
  if (exist_cat && !is_onehot) {
    df_new_merge = factor_encode(df_new_merge, col_cat + 1)
    df_new_mode_cat = aggregate(. ~ index , data = df_new_merge[c("index", col_name_cat)], Mode_cat)
    df_new = merge(df_new_mean_num, df_new_mode_cat, by = "index")
    df_new = factor_encode(df_new, col_cat + 1) # +1 is for the index column
    if (is_unalike) {
      df_new_cat_var = aggregate(. ~ index , data = df_new_merge[c("index", col_name_cat)], unalike)
    }
    else{
      df_new_cat_var = df_new_merge[c("index", col_name_cat)] %>% group_by(index) %>% summarise(across(col_name_cat, VA_fact))
      #df_new_cat_var = aggregate(.~index , data =df_new_merge[c("index",col_name_cat)], VA_fact)
    }
    df_new_var = merge(df_new_num_var, df_new_cat_var, by = "index")
    # We need to use unalikeability to mesure the uncertainty of the categorical variables
  }
  else{
    df_new = df_new_mean_num
    df_new_var = df_new_num_var
  }
  #Add the combined result for categorical variables deducted from the onehot result to df_new
  if (is_onehot && exist_cat) {
    for (name in names_cat) {
      df_new[[name]] = apply(df_new[dict_cat[[name]]], 1, which_max_cat, name)
      df_new[[name]] = unlist(df_new[[name]])
    }
    if (is_unalike) {
      df_new_cat_var = aggregate(. ~ index , data = df_new_merge[c("index", names_cat)], unalike)
    }
    else{
      df_new_cat_var = df_new_merge[c("index", names_cat)] %>% group_by(index) %>% summarise(across(names_cat, VA_fact))
      #df_new_cat_var = aggregate(.~index , simplify=FALSE, data =df_new_merge[c("index",names_cat)], VA_fact)
    }
    df_new_var = merge(df_new_num_var, df_new_cat_var, by = "index")
  }
  
  # Add the uncovered rows, fill in with NA
  num_row_new = nrow(df_new)
  if (num_row_new < num_row_origin) {
    print(paste0("Covered number of rows: ", num_row_new))
    print(paste0("Original number of rows: ", num_row_origin))
    warning(
      "Warning: All the rows in the original dataset is not selected in the bootstrap samples. An increase on the number of bootstrap samples is suggested.\n"
    )
  }
  df_idx = data.frame("index" = c(1:num_row_origin))
  df_result_with_idx = merge(x = df_idx,
                             y = df_new,
                             by = "index",
                             all = TRUE)
  df_result_var_with_idx = merge(x = df_idx,
                                 y = df_new_var,
                                 by = "index",
                                 all = TRUE)
  
  #remove index column and the combined categorical columns
  if (is_onehot && exist_cat) {
    df_result_disj = df_result_with_idx[, !(names(df_result_with_idx) %in% c(names_cat, 'index'))]
    df_result_var_disj = df_result_var_with_idx[, !(names(df_result_with_idx) %in% c(names_cat, 'index'))]
  }
  else{
    df_result_disj = df_result_with_idx[-c(1)]
    df_result_var_disj = df_result_var_with_idx[-c(1)]
  }
  
  
  
  #Final result for the categorical variables
  if (is_onehot) {
    df_result = df_result_with_idx[-c(1, col_cat + 1)] #remove index and onehot columns, +1 is for the adjustment
    df_result_var = df_result_var_with_idx[-c(1, col_cat + 1)]
  }
  else{
    df_result = df_result_disj
    df_result_var = df_result_var_disj
  }
  
  
  
  return(
    list(
      "imp.disj" = df_result_disj,
      "uncertainty.disj" = df_result_var_disj,
      "imp" = df_result,
      "uncertainty" = df_result_var
    )
  )
  
}
