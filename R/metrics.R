#' @title Confusion Matrix
#'
#' @description
#' Compute confusion matrix to evaluate the accuracy of a classification.
#'
#' @param y_pred Predicted labels vector, as returned by a classifier
#' @param y_true Ground truth (correct) 0-1 labels vector
#' @return a table of Confusion Matrix
#' @examples
#' data(cars)
#' logreg <- glm(
#'   formula = vs ~ hp + wt,
#'   family = binomial(link = "logit"), data = mtcars
#' )
#' pred <- ifelse(logreg$fitted.values < 0.5, 0, 1)
#' ConfusionMatrix(y_pred = pred, y_true = mtcars$vs)
#' @export
ConfusionMatrix <- function(y_pred, y_true) {
  Confusion_Mat <- table(y_true, y_pred)
  return(Confusion_Mat)
}

#' @title Confusion Matrix (Data Frame Format)
#'
#' @description
#' Compute data frame format confusion matrix for internal usage.
#'
#' @param y_pred Predicted labels vector, as returned by a classifier
#' @param y_true Ground truth (correct) 0-1 labels vector
#' @return a data.frame of Confusion Matrix
#' @examples
#' data(cars)
#' logreg <- glm(
#'   formula = vs ~ hp + wt,
#'   family = binomial(link = "logit"), data = mtcars
#' )
#' pred <- ifelse(logreg$fitted.values < 0.5, 0, 1)
#' ConfusionDF(y_pred = pred, y_true = mtcars$vs)
#' @keywords internal
#' @export
ConfusionDF <- function(y_pred, y_true) {
  Confusion_DF <- transform(as.data.frame(ConfusionMatrix(y_pred, y_true)),
    y_true = as.character(y_true),
    y_pred = as.character(y_pred),
    Freq = as.integer(Freq)
  )
  return(Confusion_DF)
}
utils::globalVariables("Freq")

#' @title Precision (micro averaged)
#'
#' @description
#' Compute the precision score of multi-class problem using the "micro" average.
#' details: https://sebastianraschka.com/faq/docs/multiclass-metric.html
#'
#' @param y_pred Predicted labels vector, as returned by a classifier
#' @param y_true Ground truth (correct) labels vector
#' @param labels An optional vector containing the list of the existent
#'   (unique) labels.
#' @return Precision (micro averaged)
#' @examples
#' labels <- c("Q1", "Q2", "Q3", "Q4")
#' truth <- sample(labels, 10, replace = TRUE)
#' pred <- sample(labels, 10, replace = TRUE)
#' Precision_micro(y_pred = pred, y_true = truth, labels)
#' @export
Precision_micro <- function(y_true, y_pred, labels = NULL) {
  Confusion_DF <- ConfusionDF(y_pred, y_true)

  if (is.null(labels) == TRUE) labels <- unique(c(y_true, y_pred))
  # this is not bulletproof since there might be labels missing (in strange cases)
  # in strange cases where they existed in training set but are missing from test ground truth and predictions.

  TP <- c()
  FP <- c()
  for (i in c(1:length(labels))) {
    positive <- labels[i]

    # it may happen that a label is never predicted (missing from y_pred) but exists in y_true
    # in this case ConfusionDF will not have these lines and thus the simplified code crashes
    # TP[i] <- as.integer(Confusion_DF[which(Confusion_DF$y_true==positive & Confusion_DF$y_pred==positive), "Freq"])
    # FP[i] <- as.integer(sum(Confusion_DF[which(Confusion_DF$y_true!=positive & Confusion_DF$y_pred==positive), "Freq"]))

    # workaround:
    # i don't want to change ConfusionDF since i don't know if the current behaviour is a feature or a bug.
    tmp <- Confusion_DF[which(Confusion_DF$y_true == positive & Confusion_DF$y_pred == positive), "Freq"]
    TP[i] <- if (length(tmp) == 0) 0 else as.integer(tmp)

    tmp <- Confusion_DF[which(Confusion_DF$y_true != positive & Confusion_DF$y_pred == positive), "Freq"]
    FP[i] <- if (length(tmp) == 0) 0 else as.integer(sum(tmp))
  }
  Precision_micro <- sum(TP) / (sum(TP) + sum(FP))
  return(Precision_micro)
}

#' @title Recall (micro averaged)
#'
#' @description
#' Compute the recall score of multi-class problem using the "micro" average.
#' details: https://sebastianraschka.com/faq/docs/multiclass-metric.html
#'
#' @param y_pred Predicted labels vector, as returned by a classifier
#' @param y_true Ground truth (correct) labels vector
#' @param labels An optional vector containing the list of the existent
#'   (unique) labels.
#' @return Recall (micro averaged)
#' @examples
#' labels <- c("Q1", "Q2", "Q3", "Q4")
#' truth <- sample(labels, 10, replace = TRUE)
#' pred <- sample(labels, 10, replace = TRUE)
#' Recall_micro(y_pred = pred, y_true = truth, labels)
#' @export
Recall_micro <- function(y_true, y_pred, labels = NULL) {
  Confusion_DF <- ConfusionDF(y_pred, y_true)

  if (is.null(labels) == TRUE) labels <- unique(c(y_true, y_pred))
  # this is not bulletproof since there might be labels missing (in strange cases)
  # in strange cases where they existed in training set but are missing from test ground truth and predictions.

  TP <- c()
  FN <- c()
  for (i in c(1:length(labels))) {
    positive <- labels[i]

    # short version, comment out due to bug or feature of Confusion_DF
    # TP[i] <- as.integer(Confusion_DF[which(Confusion_DF$y_true==positive & Confusion_DF$y_pred==positive), "Freq"])
    # FP[i] <- as.integer(sum(Confusion_DF[which(Confusion_DF$y_true==positive & Confusion_DF$y_pred!=positive), "Freq"]))

    # workaround:
    tmp <- Confusion_DF[which(Confusion_DF$y_true == positive & Confusion_DF$y_pred == positive), "Freq"]
    TP[i] <- if (length(tmp) == 0) 0 else as.integer(tmp)

    tmp <- Confusion_DF[which(Confusion_DF$y_true == positive & Confusion_DF$y_pred != positive), "Freq"]
    FN[i] <- if (length(tmp) == 0) 0 else as.integer(sum(tmp))
  }
  Recall_micro <- sum(TP) / (sum(TP) + sum(FN))
  return(Recall_micro)
}

#' @title F1 Score (micro averaged)
#'
#' @description
#' Compute the F1 Score of multi-class problem using the "micro" average.
#' details: https://sebastianraschka.com/faq/docs/multiclass-metric.html
#'
#' @param y_pred Predicted labels vector, as returned by a classifier
#' @param y_true Ground truth (correct) labels vector
#' @param labels An optional vector containing the list of the existent
#'   (unique) labels.
#' @return F1 Score (micro averaged)
#' @examples
#' labels <- c("Q1", "Q2", "Q3", "Q4")
#' truth <- sample(labels, 10, replace = TRUE)
#' pred <- sample(labels, 10, replace = TRUE)
#' F1_Score_micro(y_pred = pred, y_true = truth, labels)
#' @export
F1_Score_micro <- function(y_true, y_pred, labels = NULL) {
  if (is.null(labels) == TRUE) {
    labels <- unique(c(y_true, y_pred))
  } # possible problems if labels are missing from y_*
  Precision <- Precision_micro(y_true, y_pred, labels)
  Recall <- Recall_micro(y_true, y_pred, labels)
  F1_Score_micro <- 2 * (Precision * Recall) / (Precision + Recall)
  return(F1_Score_micro)
}
