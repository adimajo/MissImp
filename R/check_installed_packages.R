is_ranger_package_installed <- function() {
  installed_ranger <- requireNamespace("ranger", quietly = TRUE)
  if (!installed_ranger) {
    stop("Package ranger not installed.")
  }
  return(installed_ranger)
}

is_foreach_package_installed <- function() {
  installed_foreach <- requireNamespace("foreach", quietly = TRUE)
  if (!installed_foreach) {
    stop("Package foreach not installed.")
  }
  return(installed_foreach)
}

is_glmnet_package_installed <- function() {
  installed_glmnet <- requireNamespace("glmnet", quietly = TRUE)
  if (!installed_glmnet) {
    stop("Package glmnet not installed.")
  }
  return(installed_glmnet)
}

is_LiblineaR_package_installed <- function() {
  installed_LiblineaR <- requireNamespace("LiblineaR", quietly = TRUE)
  if (!installed_LiblineaR) {
    stop("Package LiblineaR not installed.")
  }
  return(installed_LiblineaR)
}

is_norm_package_installed <- function() {
  installed_norm <- requireNamespace("norm", quietly = TRUE)
  if (!installed_norm) {
    stop("Package norm not installed.")
  }
  return(installed_norm)
}

is_missMethods_package_installed <- function() {
  installed_missMethods <-
    requireNamespace("missMethods", quietly = TRUE)
  if (!installed_missMethods) {
    stop("Package missMethods not installed.")
  }
  return(installed_missMethods)
}

is_pracma_package_installed <- function() {
  installed_pracma <- requireNamespace("pracma", quietly = TRUE)
  if (!installed_pracma) {
    stop("Package pracma not installed.")
  }
  return(installed_pracma)
}

is_mltools_package_installed <- function() {
  installed_mltools <- requireNamespace("mltools", quietly = TRUE)
  if (!installed_mltools) {
    stop("Package mltools not installed.")
  }
  return(installed_mltools)
}

is_uwo4419_package_installed <- function() {
  installed_uwo4419 <- requireNamespace("uwo4419", quietly = TRUE)
  if (!installed_uwo4419) {
    stop("Package uwo4419 not installed.")
  }
  return(installed_uwo4419)
}

is_randomForest_package_installed <- function() {
  installed_randomForest <-
    requireNamespace("randomForest", quietly = TRUE)
  if (!installed_randomForest) {
    stop("Package randomForest not installed.")
  }
  return(installed_randomForest)
}

is_FactoMineR_package_installed <- function() {
  installed_FactoMineR <-
    requireNamespace("FactoMineR", quietly = TRUE)
  if (!installed_FactoMineR) {
    stop("Package FactoMineR not installed.")
  }
  return(installed_FactoMineR)
}

is_qualvar_package_installed <- function() {
  installed_qualvar <- requireNamespace("qualvar", quietly = TRUE)
  if (!installed_qualvar) {
    stop("Package qualvar not installed.")
  }
  return(installed_qualvar)
}

is_VIM_package_installed <- function() {
  installed_VIM <- requireNamespace("VIM", quietly = TRUE)
  if (!installed_VIM) {
    stop("Package VIM not installed.")
  }
  return(installed_VIM)
}

is_data.table_package_installed <- function() {
  installed_data.table <-
    requireNamespace("data.table", quietly = TRUE)
  if (!installed_data.table) {
    stop("Package data.table not installed.")
  }
  return(installed_data.table)
}

is_mix_package_installed <- function() {
  installed_mix <- requireNamespace("mix", quietly = TRUE)
  if (!installed_mix) {
    stop("Package mix not installed.")
  }
  return(installed_mix)
}

is_Amelia_package_installed <- function() {
  installed_Amelia <- requireNamespace("Amelia", quietly = TRUE)
  if (!installed_Amelia) {
    stop("Package Amelia not installed.")
  }
  return(installed_Amelia)
}

is_missMDA_package_installed <- function() {
  installed_missMDA <- requireNamespace("missMDA", quietly = TRUE)
  if (!installed_missMDA) {
    stop("Package missMDA not installed.")
  }
  return(installed_missMDA)
}

is_abind_package_installed <- function() {
  installed_abind <- requireNamespace("abind", quietly = TRUE)
  if (!installed_abind) {
    stop("Package abind not installed.")
  }
  return(installed_abind)
}

is_mvtnorm_package_installed <- function() {
  installed_mvtnorm <- requireNamespace("mvtnorm", quietly = TRUE)
  if (!installed_mvtnorm) {
    stop("Package mvtnorm not installed.")
  }
  return(installed_mvtnorm)
}

is_gdata_package_installed <- function() {
  installed_gdata <- requireNamespace("gdata", quietly = TRUE)
  if (!installed_gdata) {
    stop("Package gdata not installed.")
  }
  return(installed_gdata)
}

is_naniar_package_installed <- function() {
  installed_naniar <- requireNamespace("naniar", quietly = TRUE)
  if (!installed_naniar) {
    stop("Package naniar not installed.")
  }
  return(installed_naniar)
}

is_FNN_package_installed <- function() {
  installed_FNN <- requireNamespace("FNN", quietly = TRUE)
  if (!installed_FNN) {
    stop("Package FNN not installed.")
  }
  return(installed_FNN)
}

is_missForest_package_installed <- function() {
  installed_missForest <-
    requireNamespace("missForest", quietly = TRUE)
  if (!installed_missForest) {
    stop("Package missForest not installed.")
  }
  return(installed_missForest)
}

is_iterators_package_installed <- function() {
  installed_iterators <- requireNamespace("iterators", quietly = TRUE)
  if (!installed_iterators) {
    stop("Package iterators not installed.")
  }
  return(installed_iterators)
}

is_itertools_package_installed <- function() {
  installed_itertools <- requireNamespace("itertools", quietly = TRUE)
  if (!installed_itertools) {
    stop("Package itertools not installed.")
  }
  return(installed_itertools)
}

is_rlang_package_installed <- function() {
  installed_rlang <- requireNamespace("rlang", quietly = TRUE)
  if (!installed_rlang) {
    stop("Package rlang not installed.")
  }
  return(installed_rlang)
}
