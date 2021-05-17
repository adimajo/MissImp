#' Create the matrix of p-value for dummy t-chi-test
#' 
#' @description
#' \code{dummy_test_matrix} TODO
#' @param df TODO
#' @param col_cat TODO
#' @return TODO
#' @export
#' @references
#' toto et al.
#' @examples
#' print(toto)
dummy_test_matrix = function(df, col_cat) {
  test_result_dummy = data.frame()
  df = factor_encode(df, col_cat)  # in case the categorical column is not a factor
  R_df = data.frame(is.na(df))
  ls_col_name = colnames(df)
  ls_row_name = c()
  num_col = length(ls_col_name)
  i = 1
  while (i <= num_col) {
    col_ctr = ls_col_name[i] # control column
    df_1 = df[R_df[[col_ctr]] == 1, ]
    df_0 = df[R_df[[col_ctr]] == 0, ]
    R_1 = R_df[R_df[[col_ctr]] == 1, ] # this matrix is needed for the t-test in a special case
    R_0 = R_df[R_df[[col_ctr]] == 0, ]
    row1 = length(row.names(df_1))
    row0 = length(row.names(df_1))
    ls_row_name[i] = col_ctr
    if (row1 == 0 |
        row0 == 0) {
      # if one column contains only NA or doesn't contain any NA
      i = i + 1
      next
    }
    j = 1
    while (j <= num_col) {
      if (i != j) {
        col_test = ls_col_name[j]
        # if the test column is all NA, the test is not done on this column
        if (all(is.na(df[[col_test]]))) {
          test_result_dummy[i, col_test] = NA
          j = j + 1
          next
        }
        
        # if when the control column is NA, the test column is all NA
        # or when the control column is observed , the test column is all NA
        # the paired t-test is done on the indicator matrix R_1 an R_0
        critc_1 = all(is.na(df_1[[col_test]]))
        critc_0 = all(is.na(df_0[[col_test]]))
        if (critc_1 ||
            critc_0) {
          #the situation of critic_1 && critic_0 is discussed before
          if (critc_1 && !critc_0) {
            test_result_dummy[i, col_test] = stats::t.test(R_0[[col_test]] * 1, mu = 1)$p.value
          }
          else{
            test_result_dummy[i, col_test] = t.test(R_1[[col_test]] * 1, mu = 1)$p.value
          }
          
          j = j + 1
          next
        }
        
        # if the test column is categorical, we use chi 2 test
        if (j %in% col_cat) {
          test_table = data.frame(col_test = df[[col_test]], R = R_df[[col_ctr]])
          test_result_dummy[i, col_test] = suppressWarnings(with(test_table, chisq.test(col_test, R))$p.value)
          # warnings about the approximation of normal distribution if the sample is small
          j = j + 1
          next
        }
        
        test_result_dummy[i, col_test] = t.test(df_1[[col_test]], df_0[[col_test]])$p.value
      }
      j = j + 1
    }
    i = i + 1
  }
  row.names(test_result_dummy) = ls_row_name
  return(test_result_dummy)
}

# The combined p-value for dummy t-chi-test
# (Assume that the tests are independent: a correct assumption under MCAR?)
# Input: an incomplete dataset and the index of the categorical columns
# Output: a matrix of p-values of the dummy tests and a final p-value combined by Fisher's method
dummy_test = function(df, col_cat) {
  p_matrix = dummy_test_matrix(df, col_cat)
  p_vector = as.vector(p_matrix)
  p_vector = p_vector[!is.na(p_vector)]
  dof = 2 * length(p_vector)
  res = -2 * sum(log(p_vector))
  p_dum = pchisq(res, df = dof, lower.tail = FALSE)
  return(list(
    "p.matrix" = p_matrix,
    "dof" = dof,
    "chi2stat" = res,
    "p.value" = p_dum
  ))
}
