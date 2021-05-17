

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
