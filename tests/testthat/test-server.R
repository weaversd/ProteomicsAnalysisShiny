# tests/testthat/test-server.R
library(testthat)
library(shiny)

source("../../global.R")
source("../../R/parsers.R")
source("../../server.R")

test_that("Full pipeline executes successfully inside testServer", {
  testServer(server, {
    # 1. Create a dummy input file
    mock_df <- data.frame(
      Protein.ID = paste0("PROT_", 1:10),
      Gene = paste0("GENE_", 1:10),
      Description = paste0("Desc_", 1:10),
      WT1.Intensity = rnorm(10, mean = 20, sd = 1),
      WT2.Intensity = rnorm(10, mean = 20, sd = 1),
      KO1.Intensity = rnorm(10, mean = 18, sd = 1),
      KO2.Intensity = rnorm(10, mean = 18, sd = 1)
    )
    
    tmp_path <- tempfile(fileext = ".tsv")
    write.table(mock_df, tmp_path, sep = "\t", row.names = FALSE)
    
    # 2. Simulate Upload & Data Source selection
    session$setInputs(
      data_source = "Generic Table (Long or Wide)",
      file_upload = list(datapath = tmp_path, name = "mock_protein.tsv")
    )
    
    expect_false(is.null(rv$raw_data))
    expect_equal(length(unique(rv$raw_data$Protein.ID)), 10)
    
    # 3. Simulate Sample Mapping Trigger
    session$setInputs(
      cond_WT1 = "WT", br_WT1 = "1",
      cond_WT2 = "WT", br_WT2 = "2",
      cond_KO1 = "KO", br_KO1 = "1",
      cond_KO2 = "KO", br_KO2 = "2",
      btn_process_import = 1
    )
    
    expect_equal(unique(rv$raw_data$condition), c("WT", "KO"))
    
    # 4. Simulate Normalization and Imputation Application
    session$setInputs(
      norm_method = "center.median",
      impute_method = "None",
      btn_apply_norm = 1
    )
    
    expect_false(is.null(rv$norm_data))
    expect_true("WT-KO" %in% names(rv$de_results) || "KO-WT" %in% names(rv$de_results))
    
    # 5. Verify Plotly & ggplot Evaluation
    comp_name <- names(rv$de_results)[1]
    session$setInputs(
      selected_comparison = comp_name,
      plot_type = "Volcano",
      adj_p_cutoff = 0.05,
      fc_cutoff = 1.0,
      show_labels = TRUE,
      label_size = 3,
      max_overlaps = 10,
      point_size = 2,
      text_size = 12,
      col_up = "#66CCFE",
      col_down = "#FF0066"
    )
    
    # Evaluate reactive ggplot output without throwing errors
    p <- final_ggplot_object()
    expect_s3_class(p, "ggplot")
    
    unlink(tmp_path)
  })
})