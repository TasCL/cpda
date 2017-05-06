bib <- function() {
  sub("\\.bib$", "", system.file("bib", "cpda.bib", package = "cpda"))
}
