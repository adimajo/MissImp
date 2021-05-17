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

# MSE (average over missing number)
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
