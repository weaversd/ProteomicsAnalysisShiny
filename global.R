# global.R
# Define CRAN packages
cran_packages <- c(
  "shiny", "bslib", "dtplyr", "dplyr", "tidyr", 
  "stringr", "ggplot2", "ggrepel", "plotly", "DT", 
  "openxlsx", "jsonlite", "colourpicker", "glue"
)

# Define Bioconductor packages
bioc_packages <- c(
  "QFeatures", "limma", "MsCoreUtils"
)

# 1. Ensure BiocManager is installed (handles Bioconductor releases & dependencies)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# 2. Get vector of currently installed packages once (faster than repeated checks)
installed_pkgs <- installed.packages()[, "Package"]

# 3. Install missing CRAN packages
missing_cran <- setdiff(cran_packages, installed_pkgs)
if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# 4. Install missing Bioconductor packages
missing_bioc <- setdiff(bioc_packages, installed_pkgs)
if (length(missing_bioc) > 0) {
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# 5. Load all packages
all_packages <- c(cran_packages, bioc_packages)
invisible(lapply(all_packages, library, character.only = TRUE))

# Color Palette from baseline scripts
my_palette <- c(
  "#000000", "#FF0066", "#107F80", "#40007F", 
  "#AA66FF", "#66CCFE", "#FFAA00", "#00CC33", "#992600"
)

# Helper for QFeatures creation
readQFeatures2 <- function(table, ecol, fnames = "Protein.ID", name = "raw_proteins") {
  table_df <- as.data.frame(table)
  if (is.character(ecol)) ecol_idx <- match(ecol, colnames(table_df)) else ecol_idx <- ecol
  
  feature_names <- make.unique(as.character(table_df[[fnames]]))
  assay_mat <- as.matrix(table_df[, ecol_idx, drop = FALSE])
  rownames(assay_mat) <- feature_names
  
  row_data <- table_df[, -ecol_idx, drop = FALSE]
  rownames(row_data) <- feature_names
  
  se <- SummarizedExperiment(
    assays = setNames(list(assay_mat), name),
    rowData = row_data
  )
  QFeatures(setNames(list(se), name))
}

get_environment_audit <- function() {
  # Capture R Version
  r_ver <- R.version.string
  
  # Capture loaded package versions
  pkg_info <- installed.packages()[, "Version"]
  loaded_pkgs <- loadedNamespaces()
  
  # Filter for loaded packages and their active versions
  pkg_versions <- as.list(pkg_info[intersect(names(pkg_info), loaded_pkgs)])
  
  list(
    r_version = r_ver,
    package_versions = pkg_versions,
    platform = R.version$platform,
    timestamp = Sys.time()
  )
}

plot_normalization_qc <- function(object, 
                                  i = 1, 
                                  plot_type = c("density", "histogram"),
                                  methods = NULL) {
  
  plot_type <- match.arg(plot_type)
  
  # 1. Determine available normalization methods
  all_supported <- normalizeMethods()
  
  if (is.null(methods)) {
    methods <- c("none", "center.mean", "center.median", "div.mean", 
                 "div.median", "sum", "max", "quantiles", "quantiles.robust")
  }
  
  valid_methods <- intersect(methods, c("none", all_supported))
  
  if (length(valid_methods) == 0) {
    stop("None of the requested methods are supported by QFeatures::normalizeMethods().")
  }
  
  # 2. Extract and normalize data for each method
  plot_data_list <- list()
  
  for (method in valid_methods) {
    tryCatch({
      if (method == "none") {
        norm_matrix <- assay(object[[i]])
      } else {
        norm_matrix <- assay(normalize(object[[i]], method = method))
      }
      
      # Reshape matrix into long format for ggplot
      df_long <- norm_matrix %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column(var = "Feature") %>% 
        pivot_longer(
          cols = -Feature, 
          names_to = "Sample", 
          values_to = "Intensity"
        ) %>% 
        mutate(Method = method)
      
      plot_data_list[[method]] <- df_long
      
    }, error = function(e) {
      warning(paste("Method", method, "failed and was skipped:", e$message))
    })
  }
  
  # Combine results
  all_df <- bind_rows(plot_data_list)
  all_df$Method <- factor(all_df$Method, levels = valid_methods)
  
  # 3. Build ggplot
  p <- ggplot(all_df, aes(x = Intensity, color = Sample, group = Sample))
  
  if (plot_type == "density") {
    p <- p + geom_density(alpha = 0.3, linewidth = 0.6)
  } else {
    p <- p + geom_freqpoly(bins = 50, alpha = 0.7, linewidth = 0.6)
  }
  
  p <- p + 
    facet_wrap(~ Method, scales = "free") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92")
    ) +
    labs(
      title = paste0("Normalization Method Comparison (Assay: ", names(object)[i], ")"),
      subtitle = paste("Plot type:", plot_type),
      x = "Intensity Value",
      y = ifelse(plot_type == "density", "Density", "Count")
    )
  
  return(p)
}