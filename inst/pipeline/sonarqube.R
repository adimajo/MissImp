hacked_sonarqube <- function(cov, filename = "sonarqube.xml") {
  loadNamespace("xml2")

  df <- covr::tally_coverage(cov, by = "line")

  d <- xml2::xml_new_document()

  top <- xml2::xml_add_child(d, "coverage", version = "1")

  files <- unique(df$filename)

  for (f in files) {
    file <- xml2::xml_add_child(top, "file",
                                path = paste(as.character(f), sep = ""))

    for (fun_name in unique(na.omit(df[df$filename == f, "functions"]))) {
      fun_lines <- which(df$functions == fun_name)
      for (i in fun_lines) {
        line <- df[i, ]
        xml2::xml_add_child(file, "lineToCover",
          lineNumber = as.character(line$line),
          covered = tolower(as.character(line$value > 0))
        )
      }
    }
  }

  xml2::write_xml(d, file = filename)

  invisible(d)
}

the_coverage <- covr::package_coverage(quiet = FALSE, clean = FALSE)

print(the_coverage)

hacked_sonarqube(the_coverage, filename = "coverage.xml")
