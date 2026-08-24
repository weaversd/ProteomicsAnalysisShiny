# server.R

server <- function(input, output, session) {
  
  # Reactive values store
  rv <- reactiveValues(
    raw_data = NULL,
    geneDict = NULL,
    sample_map = NULL,
    Qprot = NULL,
    norm_data = NULL,
    de_results = list(),
    audit_log = list(),
    impute_stats = NULL
  )
  
  # ----------------------------------------------------------------------------
  # 1. Parsing Input & Re-importing State
  # ----------------------------------------------------------------------------
  observeEvent(input$file_upload, {
    req(input$file_upload)
    
    # Wrap parsing inside tryCatch for graceful error handling
    res <- tryCatch({
      
      if (input$data_source == "Spectronaut") {
        parse_spectronaut(input$file_upload$datapath)
      } else if (input$data_source == "MSFragger") {
        parse_msfragger(input$file_upload$datapath)
      } else if (input$data_source == "Generic Table (Long or Wide)") {
        parse_generic(input$file_upload$datapath)
      } else if (input$data_source == "Re-import Exported RDS/State") {
        saved_state <- readRDS(input$file_upload$datapath)
        
        # Basic sanity check on imported state structure
        if (!all(c("raw_data", "geneDict") %in% names(saved_state))) {
          stop("Invalid state file structure.")
        }
        
        rv$raw_data   <- saved_state$raw_data
        rv$geneDict   <- saved_state$geneDict
        rv$norm_data  <- saved_state$norm_data
        rv$de_results <- saved_state$de_results
        rv$audit_log  <- saved_state$audit_log
        rv$sample_map <- saved_state$sample_map
        
        showNotification("Saved state successfully restored!", type = "message")
        return(NULL)
      }
      
    }, error = function(e) {
      showNotification(
        ui = paste0("File formatted incorrectly or unsupported format selected: ", e$message),
        type = "error",
        duration = 8
      )
      return(NULL)
    })
    
    # Exit early if an error occurred during parsing
    if (is.null(res)) return()
    
    # Successfully parsed: Store data into reactive state
    rv$raw_data <- res$data
    rv$geneDict <- res$geneDict
    
    # --- AUTO-INITIALIZE RAW QFEATURES OBJECT FOR QC PLOT ---
    wide_data_init <- rv$raw_data %>%
      select(Protein.ID, ID, LogInt) %>%
      pivot_wider(names_from = ID, values_from = LogInt)
    
    rv$Qprot <- readQFeatures2(
      wide_data_init, 
      ecol = 2:ncol(wide_data_init), 
      fnames = "Protein.ID", 
      name = "raw"
    )
    # --------------------------------------------------------
    # Initialize basic audit information
    rv$audit_log$import_time <- Sys.time()
    rv$audit_log$source_type <- input$data_source
    rv$audit_log$file_name   <- input$file_upload$name
    
    # Attach R and Package Versions
    if (exists("get_environment_audit")) {
      rv$audit_log$environment <- get_environment_audit()
    }
    
    showNotification("File loaded successfully! Please verify sample metadata below.", type = "message")
  })
  
  # Render Dynamic Mapping Controls for Each Detected Sample
  output$sample_mapping_ui <- renderUI({
    req(rv$raw_data)
    
    unique_ids <- unique(rv$raw_data$ID)
    
    mapping_rows <- lapply(unique_ids, function(id) {
      default_cond <- str_remove_all(str_extract(id, "x[0-9]|[A-Za-z]+"), "x")
      default_br   <- str_extract(id, "\\d+$")
      if (is.na(default_cond) || default_cond == "") default_cond <- "Cond"
      if (is.na(default_br)   || default_br == "")   default_br <- "1"
      
      fluidRow(
        column(4, tags$strong(id, style = "font-size: 11px; vertical-align: -20px;")),
        column(4, textInput(paste0("cond_", id), label = NULL, value = default_cond, placeholder = "Condition")),
        column(4, textInput(paste0("br_", id),   label = NULL, value = default_br,   placeholder = "Replicate"))
      )
    })
    
    tagList(
      fluidRow(
        column(4, tags$b("Detected ID")),
        column(4, tags$b("Condition")),
        column(4, tags$b("Replicate"))
      ),
      hr(style = "margin-top: 5px; margin-bottom: 10px;"),
      mapping_rows
    )
  })
  
  # Apply Sample Mapping Customizations
  observeEvent(input$btn_process_import, {
    req(rv$raw_data)
    
    unique_ids <- unique(rv$raw_data$ID)
    
    mapping_df <- do.call(rbind, lapply(unique_ids, function(id) {
      cond_val <- input[[paste0("cond_", id)]]
      br_val   <- input[[paste0("br_", id)]]
      
      if (is.null(cond_val) || cond_val == "") cond_val <- "Unspecified"
      if (is.null(br_val)   || br_val == "")   br_val   <- "1"
      
      data.frame(
        ID = id,
        user_condition = cond_val,
        user_BR = br_val,
        new_ID = paste0(cond_val, br_val),
        stringsAsFactors = FALSE
      )
    }))
    
    rv$raw_data <- rv$raw_data %>%
      select(-any_of(c("condition", "BR"))) %>%
      left_join(mapping_df, by = "ID") %>%
      mutate(
        condition = user_condition,
        BR = user_BR,
        ID = new_ID
      ) %>%
      select(Protein.ID, ID, Intensity, condition, BR, LogInt)
    
    rv$sample_map <- mapping_df
    rv$audit_log$sample_mapping <- mapping_df
    
    # --- UPDATE QFEATURES OBJECT WITH NEW SAMPLE NAMES ---
    wide_data_mapped <- rv$raw_data %>%
      select(Protein.ID, ID, LogInt) %>%
      pivot_wider(names_from = ID, values_from = LogInt)
    
    rv$Qprot <- readQFeatures2(
      wide_data_mapped, 
      ecol = 2:ncol(wide_data_mapped), 
      fnames = "Protein.ID", 
      name = "raw"
    )
    # -----------------------------------------------------
    
    showNotification("Sample metadata successfully mapped!", type = "message")
  })
  
  # Data Preview Table
  output$import_preview_table <- renderDT({
    req(rv$raw_data)
    datatable(head(rv$raw_data, 100), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # ----------------------------------------------------------------------------
  # 2. Normalization & Imputation Processing
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_apply_norm, {
    req(rv$raw_data)
    
    wide_data <- rv$raw_data %>%
      select(Protein.ID, ID, LogInt) %>%
      pivot_wider(names_from = ID, values_from = LogInt)
    
    qobj <- readQFeatures2(wide_data, ecol = 2:ncol(wide_data), fnames = "Protein.ID", name = "raw")
    
    if (input$norm_method != "none") {
      qobj <- addAssay(qobj, normalize(qobj[["raw"]], method = input$norm_method), name = "norm")
    } else {
      qobj <- addAssay(qobj, qobj[["raw"]], name = "norm")
    }
    
    norm_mat <- assay(qobj[["norm"]])
    total_nas <- sum(is.na(norm_mat))
    mar_cells_count <- 0
    mnar_cells_count <- 0
    
    if (input$impute_method == "Hybrid (MAR: KNN / MNAR: MinDet)") {
      global_mar <- MsCoreUtils::impute_matrix(norm_mat, method = "nbavg")
      cond_lookup <- setNames(rv$sample_map$user_condition, rv$sample_map$new_ID)
      conds <- unname(cond_lookup[colnames(norm_mat)])
      
      imputed_mat <- norm_mat
      
      for (cond in unique(conds)) {
        cols <- which(conds == cond)
        n_reps <- length(cols)
        
        sub_m <- norm_mat[, cols, drop = FALSE]
        present <- rowSums(!is.na(sub_m))
        mar_threshold <- ceiling(n_reps / 2)
        
        mar_proteins <- names(present[present >= mar_threshold & present < n_reps])
        if (length(mar_proteins) > 0) {
          # Count exact missing cells for MAR
          mar_cells_count <- mar_cells_count + sum(is.na(sub_m[mar_proteins, , drop = FALSE]))
          imputed_mat[mar_proteins, cols] <- global_mar[mar_proteins, cols]
        }
        
        mnar_proteins <- names(present[present < mar_threshold])
        if (length(mnar_proteins) > 0) {
          # Count exact missing cells for MNAR
          mnar_cells_count <- mnar_cells_count + sum(is.na(sub_m[mnar_proteins, , drop = FALSE]))
          mnar_sub_mat <- sub_m[mnar_proteins, , drop = FALSE]
          mnar_imputed_sub <- MsCoreUtils::impute_matrix(mnar_sub_mat, method = "MinDet")
          imputed_mat[mnar_proteins, cols] <- mnar_imputed_sub
        }
      }
      
      norm_mat <- imputed_mat
      rv$impute_stats <- list(
        method = input$impute_method,
        total_na = total_nas,
        mar = mar_cells_count,
        mnar = mnar_cells_count
      )
    } else if (input$impute_method != "None") {
      norm_mat <- MsCoreUtils::impute_matrix(norm_mat, method = input$impute_method)
      rv$impute_stats <- list(
        method = input$impute_method,
        total_na = total_nas,
        imputed_total = total_nas
      )
    } else {
      rv$impute_stats <- list(
        method = "None",
        total_na = total_nas
      )
    }
    
    rv$norm_data <- as.data.frame(norm_mat)
    rv$Qprot <- qobj
    
    run_limma_analysis()
    
    rv$audit_log$normalization <- input$norm_method
    rv$audit_log$imputation <- input$impute_method
    rv$audit_log$imputation_stats <- rv$impute_stats
  })
  
  output$imputation_stats_text <- renderText({
    if (is.null(rv$impute_stats)) {
      return("Imputation has not been applied yet. Select options and click 'Apply Transformation'.")
    }
    
    stats <- rv$impute_stats
    if (stats$method == "Hybrid (MAR: KNN / MNAR: MinDet)") {
      paste0(
        "Method: Hybrid Imputation\n",
        "Total Missing Values: ", stats$total_na, "\n",
        " - Missing at Random (MAR, >=50% present in condition -> KNN/nbavg): ", stats$mar, " values\n",
        " - Missing Not at Random (MNAR, <50% present in condition -> MinDet): ", stats$mnar, " values"
      )
    } else if (stats$method == "None") {
      paste0("Method: None\nRemaining Missing Values: ", stats$total_na)
    } else {
      paste0(
        "Method: ", stats$method, "\n",
        "Total Values Imputed: ", stats$total_na
      )
    }
  })
  
  # Render QC Plot for Normalization
  output$norm_qc_plot <- renderPlot({
    req(rv$Qprot)
    plot_normalization_qc(rv$Qprot, i = 1)
  })
  
  # Render Transformed Data Output Table
  output$transformed_data_table <- renderDT({
    req(rv$norm_data)
    datatable(rv$norm_data, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # ----------------------------------------------------------------------------
  # 3. Differential Expression (limma Engine)
  # ----------------------------------------------------------------------------
  run_limma_analysis <- function() {
    req(rv$norm_data, rv$sample_map)
    data <- rv$norm_data
    
    cond_lookup <- setNames(rv$sample_map$user_condition, rv$sample_map$new_ID)
    groups <- unname(cond_lookup[colnames(data)])
    
    design <- model.matrix(~0 + factor(groups))
    colnames(design) <- make.names(unique(groups))
    
    fit1 <- lmFit(data, design)
    
    unique_groups <- colnames(design)
    if (length(unique_groups) < 2) return()
    
    combos <- combn(unique_groups, 2, simplify = FALSE)
    contrast_strings <- sapply(combos, function(x) paste0(x[1], "-", x[2]))
    
    cm <- makeContrasts(contrasts = contrast_strings, levels = design)
    fit2 <- eBayes(contrasts.fit(fit1, cm))
    
    long_data <- data %>%
      mutate(Accession = rownames(.)) %>%
      pivot_longer(cols = -Accession, names_to = "Sample", values_to = "Quant") %>%
      mutate(Condition = make.names(unname(cond_lookup[Sample])))
    
    results_list <- list()
    for (comp in colnames(cm)) {
      dt <- topTable(fit2, coef = comp, number = Inf, adjust.method = "BH", confint = TRUE)
      dt$Accession <- rownames(dt)
      
      exp_cond <- str_trim(str_split(comp, "-")[[1]][1])
      ref_cond <- str_trim(str_split(comp, "-")[[1]][2])
      
      scatter_means <- long_data %>%
        filter(Condition %in% c(exp_cond, ref_cond)) %>%
        group_by(Accession) %>%
        summarise(
          ExpQuant = mean(Quant[Condition == exp_cond], na.rm = TRUE),
          RefQuant = mean(Quant[Condition == ref_cond], na.rm = TRUE),
          .groups = "drop"
        )
      
      dt <- dt %>% 
        left_join(scatter_means, by = "Accession") %>%
        left_join(rv$geneDict, by = "Accession")
      
      attr(dt, "exp_cond") <- exp_cond
      attr(dt, "ref_cond") <- ref_cond
      
      results_list[[comp]] <- dt
    }
    rv$de_results <- results_list
  }
  
  output$comparison_selector <- renderUI({
    req(rv$de_results)
    selectInput("selected_comparison", "Select Contrast Comparison:", choices = names(rv$de_results))
  })
  
  # ----------------------------------------------------------------------------
  # 4. Interactive Plot Engine (ggplot & Plotly)
  # ----------------------------------------------------------------------------
  base_ggplot <- reactive({
    req(rv$de_results, input$selected_comparison)
    df <- rv$de_results[[input$selected_comparison]]
    
    exp_name <- attr(df, "exp_cond") %||% "Experimental"
    ref_name <- attr(df, "ref_cond") %||% "Reference"
    comp_title <- input$selected_comparison
    
    # Flip sign and condition labels if enabled
    if (isTRUE(input$flip_direction)) {
      df$logFC <- -df$logFC
      if ("t" %in% names(df)) df$t <- -df$t
      
      tmp <- exp_name
      exp_name <- ref_name
      ref_name <- tmp
      comp_title <- paste0(exp_name, "-", ref_name)
    }
    
    raw_p_cutoff <- input$adj_p_cutoff
    log10_p_line <- -log10(raw_p_cutoff)
    
    df <- df %>%
      mutate(
        significant = abs(logFC) > input$fc_cutoff & adj.P.Val < raw_p_cutoff,
        Regulation = case_when(
          significant & logFC > 0  ~ "Upregulated",
          significant & logFC <= 0 ~ "Downregulated",
          TRUE                     ~ "Not Significant"
        ),
        Regulation = factor(Regulation, levels = c("Upregulated", "Downregulated", "Not Significant"))
      )
    
    if (input$plot_type == "Volcano") {
      p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), text = paste("Gene:", gene, "<br>Accession:", Accession))) +
        geom_point(aes(color = Regulation), size = input$point_size) +
        geom_hline(yintercept = log10_p_line, linetype = 2, color = "grey50") +
        geom_vline(xintercept = c(-input$fc_cutoff, input$fc_cutoff), linetype = 2, color = "grey50") +
        labs(
          title = paste("Volcano Plot:", comp_title), 
          x = paste0("Log2 Fold Change (", exp_name, " / ", ref_name, ")"), 
          y = "-Log10 Adjusted p-value"
        )
      
    } else if (input$plot_type == "MA") {
      p <- ggplot(df, aes(x = AveExpr, y = logFC, text = paste("Gene:", gene, "<br>Accession:", Accession))) +
        geom_point(aes(color = Regulation), size = input$point_size) +
        geom_hline(yintercept = 0, color = "grey30") +
        geom_hline(yintercept = c(-input$fc_cutoff, input$fc_cutoff), linetype = 2, color = "grey50") +
        labs(
          title = paste("MA Plot:", comp_title), 
          x = "Log2 Average Expression", 
          y = paste0("Log2 Fold Change (", exp_name, " / ", ref_name, ")")
        )
      
    } else if (input$plot_type == "Scatter") {
      # Swap axes if flipped
      x_val <- if (isTRUE(input$flip_direction)) "ExpQuant" else "RefQuant"
      y_val <- if (isTRUE(input$flip_direction)) "RefQuant" else "ExpQuant"
      
      p <- ggplot(df, aes(x = .data[[x_val]], y = .data[[y_val]], text = paste("Gene:", gene, "<br>Accession:", Accession))) +
        geom_point(aes(color = Regulation), size = input$point_size) +
        geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey50") +
        labs(
          title = paste("Scatter Plot:", comp_title), 
          x = paste("Log2", ref_name, "Average Abundance"), 
          y = paste("Log2", exp_name, "Average Abundance")
        )
    }
    
    p <- p + 
      scale_color_manual(
        name = "Expression Status",
        values = c(
          "Upregulated"     = input$col_up, 
          "Downregulated"   = input$col_down, 
          "Not Significant" = "grey70"
        ),
        drop = FALSE
      ) +
      theme_bw(base_size = input$text_size) +
      theme(
        panel.grid = element_blank(),
        legend.position = "right"
      )
    
    p
  })
  
  # Render Plotly Output
  output$plotly_view <- renderPlotly({
    req(base_ggplot())
    ggplotly(base_ggplot(), tooltip = "text") %>%
      layout(
        autosize = TRUE,
        legend = list(title = list(text = "Expression Status"))
      )
  })
  
  # Single Unified Reactive Expression for Custom ggplot with Labels
  final_ggplot_object <- reactive({
    show_lbls   <- input$show_labels
    lbl_size    <- input$label_size
    max_ovrlaps <- input$max_overlaps
    raw_p_cut   <- input$adj_p_cutoff
    fc_cut      <- input$fc_cutoff
    
    p <- base_ggplot()
    req(rv$de_results, input$selected_comparison)
    
    if (isTRUE(show_lbls)) {
      df <- rv$de_results[[input$selected_comparison]]
      
      if (isTRUE(input$flip_direction)) {
        df$logFC <- -df$logFC
      }
      
      df$significant <- abs(df$logFC) > fc_cut & df$adj.P.Val < raw_p_cut
      
      sig_df <- df %>% filter(significant & !is.na(gene) & gene != "")
      
      if (nrow(sig_df) > 0) {
        p <- p + geom_label_repel(
          data = sig_df,
          aes(label = gene, group = paste0(gene, "_", max_ovrlaps)), 
          size = lbl_size,
          max.overlaps = max_ovrlaps,
          show.legend = FALSE,
          inherit.aes = TRUE
        )
      }
    }
    
    return(p)
  })
  
  # Render ggplot Output
  output$ggplot_view <- renderPlot({
    final_ggplot_object()
  })
  
  # ----------------------------------------------------------------------------
  # 5. Audit Log & Export Handlers
  # ----------------------------------------------------------------------------
  output$audit_preview <- renderText({
    toJSON(rv$audit_log, pretty = TRUE)
  })
  
  output$download_excel <- downloadHandler(
    filename = function() { paste0(input$export_prefix, "_results.xlsx") },
    content = function(file) {
      wb <- createWorkbook()
      for (comp in names(rv$de_results)) {
        addWorksheet(wb, comp)
        writeData(wb, comp, rv$de_results[[comp]])
      }
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$download_audit <- downloadHandler(
    filename = function() { 
      paste0(input$export_prefix, "_audit_trail.json") 
    },
    content = function(file) {
      if (exists("get_environment_audit")) {
        rv$audit_log$environment <- get_environment_audit()
      }
      writeLines(toJSON(rv$audit_log, pretty = TRUE, auto_unbox = TRUE), file)
    }
  )
  
  output$download_state <- downloadHandler(
    filename = function() { paste0(input$export_prefix, "_state.rds") },
    content = function(file) {
      saveRDS(list(
        raw_data   = rv$raw_data,
        geneDict   = rv$geneDict,
        norm_data  = rv$norm_data,
        de_results = rv$de_results,
        audit_log  = rv$audit_log,
        sample_map = rv$sample_map
      ), file)
    }
  )
  
  output$download_plot_png <- downloadHandler(
    filename = function() {
      paste0(input$selected_comparison, "_", input$plot_type, ".png")
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = final_ggplot_object(),
        device = "png",
        width = input$plot_width,
        height = input$plot_height,
        dpi = 300
      )
    }
  )
  
  output$download_plot_pdf <- downloadHandler(
    filename = function() {
      paste0(input$selected_comparison, "_", input$plot_type, ".pdf")
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = final_ggplot_object(),
        device = "pdf",
        width = input$plot_width,
        height = input$plot_height
      )
    }
  )
  
  output$experiment_summary_table <- renderDT({
    req(rv$raw_data)
    
    # Calculate per condition metrics using unlogged Intensity
    summary_df <- rv$raw_data %>%
      group_by(condition) %>%
      summarise(
        Total_Samples = n_distinct(ID),
        Proteins_Identified_Any = n_distinct(Protein.ID[!is.na(Intensity) & Intensity > 0]),
        .groups = "drop"
      )
    
    # Calculate %CV per protein per condition on raw intensity scale
    cv_per_cond <- rv$raw_data %>%
      filter(!is.na(Intensity) & Intensity > 0) %>%
      group_by(condition, Protein.ID) %>%
      summarise(
        n_obs = n(),
        mean_int = mean(Intensity, na.rm = TRUE),
        sd_int = sd(Intensity, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_obs > 1) %>% # Need at least 2 replicates to compute CV
      mutate(CV = (sd_int / mean_int) * 100) %>%
      group_by(condition) %>%
      summarise(
        Median_CV_Percent = round(median(CV, na.rm = TRUE), 2),
        .groups = "drop"
      )
    
    final_summary <- summary_df %>%
      left_join(cv_per_cond, by = "condition") %>%
      mutate(Median_CV_Percent = ifelse(is.na(Median_CV_Percent), "N/A (<2 reps)", paste0(Median_CV_Percent, "%")))
    
    datatable(
      final_summary,
      colnames = c("Condition", "Number of Replicates", "Proteins Identified", "Median %CV"),
      options = list(dom = 't', paging = FALSE),
      rownames = FALSE
    )
  })
  
  output$download_script <- downloadHandler(
    filename = function() {
      paste0(input$export_prefix, "_analysis_pipeline.R")
    },
    content = function(file) {
      req(rv$raw_data, rv$sample_map)
      
      mode_choice       <- input$script_source_type
      data_source_type  <- input$data_source
      raw_file_name     <- input$file_upload$name %||% "raw_proteomics_data.tsv"
      norm_method_val   <- input$norm_method
      impute_method_val <- input$impute_method
      flip_val          <- isTRUE(input$flip_direction)
      adj_p_cut         <- input$adj_p_cutoff
      fc_cut            <- input$fc_cutoff
      pt_size           <- input$point_size
      txt_size          <- input$text_size
      lbl_size          <- input$label_size
      max_overlaps_val  <- input$max_overlaps
      show_lbls         <- isTRUE(input$show_labels)
      col_up_val        <- input$col_up
      col_down_val      <- input$col_down
      export_pfx        <- input$export_prefix
      
      # Serialize the sample mapping metadata as an inline R expression
      mapping_dput <- paste(capture.output(dput(rv$sample_map)), collapse = "\n")
      
      # Header & Libraries
      header_code <- glue::glue('
# ==============================================================================
# Automated Reproducible Proteomics Pipeline
# Generated from Proteomics Explorer Dashboard
# Timestamp: {Sys.time()}
# Execution Mode: {ifelse(mode_choice == "raw", "Raw Ingestion Pipeline", "State Ingestion Pipeline")}
# ==============================================================================

# 1. Load Required Packages
required_pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "ggrepel", "limma", 
  "QFeatures", "MsCoreUtils", "stringr", "readr"
)
for (p in required_pkgs) {{
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}}

readQFeatures2 <- function(table, ecol, fnames = "Protein.ID", name = "raw_proteins") {{
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
}}
')
  
  # Ingestion Block Selection
  if (mode_choice == "raw") {
    ingestion_code <- glue::glue('
# 2. Raw Ingestion & Parsers
raw_file_path <- "{raw_file_name}" # Adjust file path if located in a different working directory

parse_spectronaut <- function(path) {{
  df <- read.delim(path, sep = "\\t", check.names = FALSE)
  acc_col <- intersect(c("PG.ProteinAccessions", "Protein.ID", "ProteinAccessions"), names(df))[1]
  if (is.na(acc_col)) stop("Could not find a valid Protein Accession column in Spectronaut file.")
  
  gene_col <- intersect(c("PG.Genes", "Gene", "Gene.Name", "PG.GeneNames"), names(df))[1]
  desc_col <- intersect(c("PG.ProteinDescriptions", "Description", "PG.ProteinDescriptions"), names(df))[1]
  
  geneDict <- df %>% 
    mutate(
      Accession   = .data[[acc_col]],
      gene        = if (!is.na(gene_col) && gene_col %in% names(df)) .data[[gene_col]] else .data[[acc_col]],
      description = if (!is.na(desc_col) && desc_col %in% names(df)) .data[[desc_col]] else .data[[acc_col]]
    ) %>% 
    select(Accession, gene, description) %>% 
    mutate(
      gene = ifelse(is.na(gene) | gene == "", Accession, gene),
      description = ifelse(is.na(description) | description == "", Accession, description)
    ) %>%
    distinct(Accession, .keep_all = TRUE)
  
  df_proc <- df %>%
    mutate(
      Protein.ID = .data[[acc_col]],
      Intensity  = ifelse(PG.Quantity == 0, NA, PG.Quantity),
      condition  = as.character(R.Condition),
      BR         = as.character(R.Replicate),
      ID         = paste0(R.Condition, R.Replicate),
      LogInt     = log2(Intensity)
    ) %>%
    select(Protein.ID, ID, Intensity, condition, BR, LogInt)
  
  list(data = df_proc, geneDict = geneDict)
}}

parse_msfragger <- function(path) {{
  df <- read.delim(path, sep = "\\t", check.names = FALSE)
  acc_col <- intersect(c("Protein.ID", "Protein", "Protein ID", "Accession"), names(df))[1]
  if (is.na(acc_col)) stop("Could not find a valid Protein ID column in MSFragger file.")
  
  gene_col <- intersect(c("Gene", "PG.Genes", "Gene Name"), names(df))[1]
  desc_col <- intersect(c("Description", "PG.ProteinDescriptions"), names(df))[1]
  
  geneDict <- df %>% 
    mutate(
      Accession   = .data[[acc_col]],
      gene        = if (!is.na(gene_col) && gene_col %in% names(df)) .data[[gene_col]] else .data[[acc_col]],
      description = if (!is.na(desc_col) && desc_col %in% names(df)) .data[[desc_col]] else .data[[acc_col]]
    ) %>% 
    select(Accession, gene, description) %>% 
    mutate(
      gene = ifelse(is.na(gene) | gene == "", Accession, gene),
      description = ifelse(is.na(description) | description == "", Accession, description)
    ) %>%
    distinct(Accession, .keep_all = TRUE)
  
  int_cols <- names(df)[which(grepl(" Intensity$", names(df)) & !grepl("MaxLFQ", names(df)))]
  if (length(int_cols) == 0) {{
    int_cols <- names(df)[which(grepl(" MaxLFQ Intensity$", names(df)))]
  }}
  if (length(int_cols) == 0) stop("No valid intensity columns found in MSFragger file.")
  
  df_proc <- df %>%
    mutate(Protein.ID = .data[[acc_col]]) %>%
    select(Protein.ID, all_of(int_cols)) %>%
    pivot_longer(cols = -Protein.ID, names_to = "ID", values_to = "Intensity") %>%
    mutate(
      Intensity = ifelse(Intensity == 0, NA, Intensity),
      ID        = str_remove(ID, "\\\\.(MaxLFQ\\\\.)?Intensity$"),
      condition = str_remove_all(str_extract(ID, "x[0-9]+|[A-Za-z]+"), "x"),
      BR        = str_extract(ID, "\\\\d+$"),
      LogInt    = log2(Intensity)
    ) %>%
    select(Protein.ID, ID, Intensity, condition, BR, LogInt)
  
  list(data = df_proc, geneDict = geneDict)
}}

# Parse source file
parsed_res <- {ifelse(data_source_type == "Spectronaut", "parse_spectronaut(raw_file_path)", "parse_msfragger(raw_file_path)")}
raw_data <- parsed_res$data
geneDict <- parsed_res$geneDict

# Apply metadata mapping configured in Dashboard
sample_map <- {mapping_dput}

raw_data <- raw_data %>%
  select(-any_of(c("condition", "BR"))) %>%
  left_join(sample_map, by = "ID") %>%
  mutate(
    condition = user_condition,
    BR = user_BR,
    ID = new_ID
  ) %>%
  select(Protein.ID, ID, Intensity, condition, BR, LogInt)
')
  } else {
    ingestion_code <- glue::glue('
# 2. Load State Objects from RDS
saved_state <- readRDS("{export_pfx}_state.rds")
raw_data   <- saved_state$raw_data
geneDict   <- saved_state$geneDict
sample_map <- saved_state$sample_map
')
  }
  
  # Downstream Analysis & Plotting
  pipeline_code <- glue::glue('
# 3. Normalization & Imputation
wide_data <- raw_data %>%
  select(Protein.ID, ID, LogInt) %>%
  pivot_wider(names_from = ID, values_from = LogInt)

qobj <- readQFeatures2(wide_data, ecol = 2:ncol(wide_data), fnames = "Protein.ID", name = "raw")

# Normalization: {norm_method_val}
if ("{norm_method_val}" != "none") {{
  qobj <- addAssay(qobj, normalize(qobj[["raw"]], method = "{norm_method_val}"), name = "norm")
}} else {{
  qobj <- addAssay(qobj, qobj[["raw"]], name = "norm")
}}

norm_mat <- assay(qobj[["norm"]])

# Imputation: {impute_method_val}
if ("{impute_method_val}" == "Hybrid (MAR: KNN / MNAR: MinDet)") {{
  global_mar <- MsCoreUtils::impute_matrix(norm_mat, method = "nbavg")
  cond_lookup <- setNames(sample_map$user_condition, sample_map$new_ID)
  conds <- unname(cond_lookup[colnames(norm_mat)])
  
  imputed_mat <- norm_mat
  
  for (cond in unique(conds)) {{
    cols <- which(conds == cond)
    n_reps <- length(cols)
    sub_m <- norm_mat[, cols, drop = FALSE]
    present <- rowSums(!is.na(sub_m))
    mar_threshold <- ceiling(n_reps / 2)
    
    mar_proteins <- names(present[present >= mar_threshold & present < n_reps])
    if (length(mar_proteins) > 0) {{
      imputed_mat[mar_proteins, cols] <- global_mar[mar_proteins, cols]
    }}
    
    mnar_proteins <- names(present[present < mar_threshold])
    if (length(mnar_proteins) > 0) {{
      mnar_sub_mat <- sub_m[mnar_proteins, , drop = FALSE]
      mnar_imputed_sub <- MsCoreUtils::impute_matrix(mnar_sub_mat, method = "MinDet")
      imputed_mat[mnar_proteins, cols] <- mnar_imputed_sub
    }}
  }}
  norm_mat <- imputed_mat
}} else if ("{impute_method_val}" != "None") {{
  norm_mat <- MsCoreUtils::impute_matrix(norm_mat, method = "{impute_method_val}")
}}

norm_data <- as.data.frame(norm_mat)

# Export processed/normalized abundance matrix to TSV
norm_export_df <- norm_data %>%
  tibble::rownames_to_column(var = "Protein.ID") %>%
  left_join(geneDict, by = c("Protein.ID" = "Accession"))

readr::write_tsv(norm_export_df, file = "{export_pfx}_normalized_imputed_matrix.tsv")

# 4. Differential Expression (limma)
cond_lookup <- setNames(sample_map$user_condition, sample_map$new_ID)
groups <- unname(cond_lookup[colnames(norm_data)])

design <- model.matrix(~0 + factor(groups))
colnames(design) <- make.names(unique(groups))

fit1 <- lmFit(norm_data, design)
unique_groups <- colnames(design)
combos <- combn(unique_groups, 2, simplify = FALSE)
contrast_strings <- sapply(combos, function(x) paste0(x[1], "-", x[2]))

cm <- makeContrasts(contrasts = contrast_strings, levels = design)
fit2 <- eBayes(contrasts.fit(fit1, cm))

long_data <- norm_data %>%
  mutate(Accession = rownames(.)) %>%
  pivot_longer(cols = -Accession, names_to = "Sample", values_to = "Quant") %>%
  mutate(Condition = make.names(unname(cond_lookup[Sample])))

de_results <- list()
for (comp in colnames(cm)) {{
  dt <- topTable(fit2, coef = comp, number = Inf, adjust.method = "BH", confint = TRUE)
  dt$Accession <- rownames(dt)
  
  exp_cond <- str_trim(str_split(comp, "-")[[1]][1])
  ref_cond <- str_trim(str_split(comp, "-")[[1]][2])
  
  scatter_means <- long_data %>%
    filter(Condition %in% c(exp_cond, ref_cond)) %>%
    group_by(Accession) %>%
    summarise(
      ExpQuant = mean(Quant[Condition == exp_cond], na.rm = TRUE),
      RefQuant = mean(Quant[Condition == ref_cond], na.rm = TRUE),
      .groups = "drop"
    )
  
  dt <- dt %>% 
    left_join(scatter_means, by = "Accession") %>%
    left_join(geneDict, by = "Accession")
  
  attr(dt, "exp_cond") <- exp_cond
  attr(dt, "ref_cond") <- ref_cond
  
  de_results[[comp]] <- dt
}}

# 5. Generate and Save All TSV Results and Plots (Volcano, MA, Scatter)
flip_direction <- {flip_val}
adj_p_cutoff <- {adj_p_cut}
fc_cutoff <- {fc_cut}
show_labels <- {show_lbls}

for (comp_name in names(de_results)) {{
  df <- de_results[[comp_name]]
  exp_name <- attr(df, "exp_cond")
  ref_name <- attr(df, "ref_cond")
  curr_title <- comp_name
  
  if (flip_direction) {{
    df$logFC <- -df$logFC
    if ("t" %in% names(df)) df$t <- -df$t
    tmp <- exp_name; exp_name <- ref_name; ref_name <- tmp
    curr_title <- paste0(exp_name, "-", ref_name)
  }}
  
  log10_p_line <- -log10(adj_p_cutoff)
  
  df <- df %>%
    mutate(
      significant = abs(logFC) > fc_cutoff & adj.P.Val < adj_p_cutoff,
      Regulation = case_when(
        significant & logFC > 0  ~ "Upregulated",
        significant & logFC <= 0 ~ "Downregulated",
        TRUE                     ~ "Not Significant"
      ),
      Regulation = factor(Regulation, levels = c("Upregulated", "Downregulated", "Not Significant"))
    )
  
  # Export full statistics for this contrast to TSV
  readr::write_tsv(df, file = paste0("{export_pfx}_", comp_name, "_de_results.tsv"))
  
  sig_df <- df %>% filter(significant & !is.na(gene) & gene != "")
  
  # A. Volcano Plot
  p_volcano <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Regulation), size = {pt_size}) +
    geom_hline(yintercept = log10_p_line, linetype = 2, color = "grey50") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = 2, color = "grey50") +
    scale_color_manual(values = c("Upregulated" = "{col_up_val}", "Downregulated" = "{col_down_val}", "Not Significant" = "grey70"), drop = FALSE) +
    theme_bw(base_size = {txt_size}) +
    theme(panel.grid = element_blank()) +
    labs(title = paste("Volcano Plot:", curr_title), x = paste0("Log2 Fold Change (", exp_name, " / ", ref_name, ")"), y = "-Log10 Adjusted p-value")
  
  if (show_labels && nrow(sig_df) > 0) {{
    p_volcano <- p_volcano + geom_label_repel(data = sig_df, aes(label = gene), size = {lbl_size}, max.overlaps = {max_overlaps_val}, show.legend = FALSE)
  }}
  ggsave(paste0("{export_pfx}_", comp_name, "_volcano.png"), plot = p_volcano, width = 8, height = 6, dpi = 300)
  
  # B. MA Plot
  p_ma <- ggplot(df, aes(x = AveExpr, y = logFC)) +
    geom_point(aes(color = Regulation), size = {pt_size}) +
    geom_hline(yintercept = 0, color = "grey30") +
    geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = 2, color = "grey50") +
    scale_color_manual(values = c("Upregulated" = "{col_up_val}", "Downregulated" = "{col_down_val}", "Not Significant" = "grey70"), drop = FALSE) +
    theme_bw(base_size = {txt_size}) +
    theme(panel.grid = element_blank()) +
    labs(title = paste("MA Plot:", curr_title), x = "Log2 Average Expression", y = paste0("Log2 Fold Change (", exp_name, " / ", ref_name, ")"))
  
  if (show_labels && nrow(sig_df) > 0) {{
    p_ma <- p_ma + geom_label_repel(data = sig_df, aes(label = gene), size = {lbl_size}, max.overlaps = {max_overlaps_val}, show.legend = FALSE)
  }}
  ggsave(paste0("{export_pfx}_", comp_name, "_ma.png"), plot = p_ma, width = 8, height = 6, dpi = 300)
  
  # C. Scatter Plot
  x_col <- if (flip_direction) "ExpQuant" else "RefQuant"
  y_col <- if (flip_direction) "RefQuant" else "ExpQuant"
  
  p_scatter <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(color = Regulation), size = {pt_size}) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey50") +
    scale_color_manual(values = c("Upregulated" = "{col_up_val}", "Downregulated" = "{col_down_val}", "Not Significant" = "grey70"), drop = FALSE) +
    theme_bw(base_size = {txt_size}) +
    theme(panel.grid = element_blank()) +
    labs(title = paste("Scatter Plot:", curr_title), x = paste("Log2", ref_name, "Average Abundance"), y = paste("Log2", exp_name, "Average Abundance"))
  
  if (show_labels && nrow(sig_df) > 0) {{
    p_scatter <- p_scatter + geom_label_repel(data = sig_df, aes(label = gene), size = {lbl_size}, max.overlaps = {max_overlaps_val}, show.legend = FALSE)
  }}
  ggsave(paste0("{export_pfx}_", comp_name, "_scatter.png"), plot = p_scatter, width = 8, height = 6, dpi = 300)
}}

cat("Pipeline completed successfully! Normalized matrix, DE result tables (.tsv), and plots (.png) have been exported.\\n")
')
  
  writeLines(paste0(header_code, ingestion_code, pipeline_code), file)
    }
  )
}