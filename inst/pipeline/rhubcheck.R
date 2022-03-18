library(MissImp)
if (!exists("arguments")) {
  suppressPackageStartupMessages(library("argparse"))
  parser <- ArgumentParser()
  parser$add_argument("-p", "--platform", type="character",
                      help="platform")
  parser$add_argument("-t", "--token", type="character",
                      help="token")
  arguments <- parser$parse_args()
}

if (length(arguments) < 2L) {
  stop("Incorrect number of args, needs 2: platform (string), token (string)")
}

platform <- arguments$platform
token <- arguments$token
if (!is.element(platform, rhub::platforms()[[1L]])) {
  stop(paste(platform, "not in rhub::platforms()[[1L]]"))
}

rhub::validate_email(
  email = substr(utils::maintainer(pkg = "MissImp"),
                 regexec("<", utils::maintainer(pkg = "MissImp"))[[1]][1] + 1,
                 nchar(utils::maintainer(pkg = "MissImp")) - 1),
  token = token
)
cr <- rhub::check(platform = platform, show_status = TRUE)
statuses <- cr[[".__enclos_env__"]][["private"]][["status_"]]

res <- do.call(rbind, lapply(statuses, function(thisStatus) {
  data.frame(
    plaform = thisStatus[["platform"]][["name"]],
    errors = length(thisStatus[["result"]][["errors"]]),
    warnings = length(thisStatus[["result"]][["warnings"]]),
    notes = length(thisStatus[["result"]][["notes"]]),
    stringsAsFactors = FALSE
  )
}))
print(res)

if (any(colSums(res[2L:3L]) > 0)) {
  stop("Some checks with errors, or warnings.")
}
if (colSums(res[4L]) > 0) {
  print("Some checks with notes")
}
