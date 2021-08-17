# Lint
docOrg <- lapply(dir("R", pattern = "*.R", recursive = TRUE, full.names = TRUE), readLines)
styler::style_pkg()
docNew <- lapply(dir("R", pattern = "*.R", recursive = TRUE, full.names = TRUE), readLines)
if (!identical(docOrg, docNew)) {
  for (doc_new in seq(length(docOrg))) {
    if (!identical(docOrg[[doc_new]], docNew[[doc_new]])) {
      for (line_new in seq(1, length(docOrg[[doc_new]]))) {
        if (!identical(docOrg[[doc_new]][[line_new]], docNew[[doc_new]][[line_new]])) {
          cat(paste("docOrg:", docOrg[[doc_new]][[line_new]], "\n"))
          cat(paste("docNew:", docNew[[doc_new]][[line_new]], "\n\n"))
        }
      }
    }
  }
  stop("Lint was not performed, try running styler::style_pkg().")
}
