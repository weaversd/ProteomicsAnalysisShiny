# launch.R
if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny", repos = "https://cloud.r-project.org", type = "binary")
}

shiny::runApp(".", launch.browser = TRUE)