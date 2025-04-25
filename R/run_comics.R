#' Run the COMICS Shiny Application
#'
#' This function launches the COMICS (Combined Outlier Method for Identifying
#' Candidate Signals) Shiny application. The application allows users to
#' analyze genomic data using ICS (Invariant Coordinate Selection) and
#' visualize results.
#'
#' @export
#' @importFrom shiny runApp
#' @examples
#' \dontrun{
#' run_comics()
#' }
run_comics <- function() {
  app_dir <- system.file("shiny", "comics", package = "COMICS")
  if (app_dir == "") {
    stop("Could not find app directory. Try reinstalling `COMICS`.", call. = FALSE)
  }
  
  shiny::runApp(app_dir, display.mode = "normal")
}
