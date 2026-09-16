# 1. Establish isolated application library
user_lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_lib)) {
  user_lib <- file.path(getwd(), "app_lib")
}
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}

.libPaths(c(user_lib, .libPaths()))

# Clean up stale locks if any crash occurred
lock_dirs <- list.dirs(user_lib, recursive = FALSE, full.names = TRUE)
lock_dirs <- lock_dirs[grepl("00LOCK", basename(lock_dirs))]
if (length(lock_dirs) > 0) unlink(lock_dirs, recursive = TRUE, force = TRUE)

# 2. Base package managers
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib, repos = "https://cloud.r-project.org", type = "binary")
}

# 3. Core dependencies to verify and install automatically on first run
cran_pkgs <- c("shiny", "bslib", "dtplyr", "dplyr", "tidyr", "stringr", 
               "ggplot2", "ggrepel", "plotly", "DT", "openxlsx", "jsonlite", "colourpicker")

bioc_pkgs <- c("QFeatures", "limma", "MsCoreUtils", "preprocessCore")

# Install missing CRAN packages
missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages into local library...")
  install.packages(missing_cran, lib = user_lib, repos = "https://cloud.r-project.org", type = "binary")
}

# Install missing Bioconductor packages
missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  message("Installing missing Bioconductor packages into local library...")
  BiocManager::install(missing_bioc, lib = user_lib, ask = FALSE, update = FALSE, type = "binary")
}

# 4. Launch Shiny App
shiny::runApp(".", launch.browser = TRUE)