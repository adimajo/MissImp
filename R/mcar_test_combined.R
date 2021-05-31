#' MCAR Test
#' @description \code{mcar_test_combined} shows if a missing data mechanism is MCAR by performing three classical statistic tests:
#'  the dummy variable test \code{dummy_test}, Little's MCAR test \code{mcar_test} from package \code{naniar} and 
#'  the normality tests (Hawkin's test and non-parametric test included) \code{TestMCARNormality} from package \code{MissMech}. 
#'  Our null hypothesis H0 is that the missing data mechanism is MCAR. 
#'  These three tests focus on testing different aspects of MCAR mechanism. 
#'  The dummy variable test is to verify if the value of one variable is correlated with the missingnness of another variable. 
#'  Little's MCAR tet focus on testing the if the mean vector for each missig pattern is part of the total mean vector. And the normality tests 
#'  focus on testing the normality and homoscedasticity according to the variance matrix. 
#'  The conclusion of these three tests as well as their statistics will be returned.
#' @param df Incomplete dataframe.
#' @param col_cat Categorical columns index.
#' @param p_val Significance level.
#' @return \code{test_results} The final results from the three tests mentioned in 'Description'. 
#' \code{True} means that we don't have enough evidence to reject our H0: the missing data mechanism is MCAR. 
#' And \code{False} means that H0 is rejected, the missing data mechanism is not MCAR.
#' @return \code{p_values} The p-values of our three tests.
#' @return \code{result_dummy_test} The detailed result from dummy variable test.
#' @return \code{result_little_test} The detailed result from Little's MCAR test.
#' @return \code{result_missmech} The detailed result from the normality tests (Hawkin's test and non-parametric test).
#' @export
#' @examples
#'  
#' mcar_test_combined(airquality, col_cat = c(5:6))$test_result
#' mcar_test_combined(oceanbuoys)$test_result
#' 
mcar_test_combined <- function(df, col_cat=c(), p_val=0.1){
  
  out_dummy <- dummy_test(df,col_cat)
  out_little <- naniar::mcar_test(df)
  exist_cat <- !all(c(0, col_cat) == c(0))
  if(exist_cat){
    warning("The normality test is written for numerical dataframes. Here the test is performed on only the numerical part of the input dataframe.\n")
    df = subset(df,select = -col_cat)
  }
  out_missmech <- TestMCARNormality(df)
  p_values <- data.frame("Dummy variable test" = out_dummy$p.value,
                               "Little's MCAR test" = out_little$p.value,
                               "Hawkin's test" = out_missmech$pvalcomb,
                               "Non parametric test" = out_missmech$pnormality)
  test_results <- data.frame("Dummy variable test" = (out_dummy$p.value>p_val),
                            "Little's MCAR test" = (out_little$p.value>p_val),
                            "MissMech test" = ((out_missmech$pvalcomb<p_val)*(out_missmech$pnormality<p_val)
                                               == 0))
  return(list(test_results=test_results ,p_values=p_values,
              result_dummy_test=out_dummy, result_little_test=out_little, result_missmech=out_missmech ))
}

