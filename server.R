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
  
  # Helper to resolve flip sign input ID whether named 'flip_direction' or 'flip_sign'
  is_flipped <- reactive({
    isTRUE(input$flip_direction) || isTRUE(input$flip_sign)
  })
  
  # ----------------------------------------------------------------------------
  # 1. Parsing Input & Re-importing State
  # ----------------------------------------------------------------------------
  observeEvent(input$file_upload, {
    req(input$file_upload)
    
    res <- tryCatch({
      if (input$data_source == "Spectronaut") {
        parse_spectronaut(input$file_upload$datapath)
      } else if (input$data_source == "MSFragger") {
        parse_msfragger(input$file_upload$datapath)
      } else if (input$data_source == "DIA-NN") {
        parse_diann(input$file_upload$datapath)
      } else if (input$data_source == "PEAKS") {
        parse_peaks(input$file_upload$datapath)
      } else if (input$data_source == "Generic Table (Long or Wide)") {
        parse_generic(input$file_upload$datapath)
      } else if (input$data_source == "Re-import Exported RDS/State") {
        saved_state <- readRDS(input$file_upload$datapath)
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
    
    if (is.null(res)) return()
    
    rv$raw_data <- res$data
    rv$geneDict <- res$geneDict
    
    # Auto-initialize raw QFeatures container
    wide_data_init <- rv$raw_data %>%
      select(Protein.ID, ID, LogInt) %>%
      pivot_wider(names_from = ID, values_from = LogInt)
    
    rv$Qprot <- readQFeatures2(
      wide_data_init, 
      ecol = 2:ncol(wide_data_init), 
      fnames = "Protein.ID", 
      name = "raw"
    )
    
    rv$audit_log$app_version <- APP_VERSION
    rv$audit_log$import_time <- Sys.time()
    rv$audit_log$source_type <- input$data_source
    rv$audit_log$file_name   <- input$file_upload$name
    
    if (exists("get_environment_audit")) {
      rv$audit_log$environment <- get_environment_audit()
    }
    
    showNotification("File loaded successfully! Please verify sample metadata below.", type = "message")
  })
  
  # Dynamic Mapping Controls
  output$sample_mapping_ui <- renderUI({
    req(rv$raw_data)
    unique_ids <- unique(rv$raw_data$ID)
    show_tr <- isTRUE(input$has_tech_reps)
    
    mapping_rows <- lapply(unique_ids, function(id) {
      # 1. Guess Technical Replicate (e.g., _TR1, _inj2, or trailing number after delimiter)
      tr_match <- str_extract(id, "(?i)(?<=[-_](tr|inj|tech))\\d+")
      if (is.na(tr_match)) {
        # Check if double-numbered pattern like Cond1_BR1_1
        tr_match <- str_extract(id, "(?<=[_-])\\d+$")
      }
      default_tr <- if (!is.na(tr_match)) tr_match else "1"
      
      # 2. Clean ID string to infer Condition and BioRep
      id_clean <- str_remove(id, "(?i)[-_](tr|inj|tech)?\\d+$")
      
      default_cond <- str_remove_all(str_extract(id_clean, "x[0-9]|[A-Za-z]+"), "x")
      default_br   <- str_extract(id_clean, "\\d+$")
      
      if (is.na(default_cond) || default_cond == "") default_cond <- "Cond"
      if (is.na(default_br)   || default_br == "")   default_br   <- "1"
      
      if (show_tr) {
        fluidRow(
          column(3, tags$strong(id, style = "font-size: 11px; word-break: break-all;")),
          column(3, textInput(paste0("cond_", id), label = NULL, value = default_cond, placeholder = "Cond")),
          column(3, textInput(paste0("br_", id),   label = NULL, value = default_br,   placeholder = "BioRep")),
          column(3, textInput(paste0("tr_", id),   label = NULL, value = default_tr,   placeholder = "TechRep"))
        )
      } else {
        fluidRow(
          column(4, tags$strong(id, style = "font-size: 11px; word-break: break-all;")),
          column(4, textInput(paste0("cond_", id), label = NULL, value = default_cond, placeholder = "Condition")),
          column(4, textInput(paste0("br_", id),   label = NULL, value = default_br,   placeholder = "Replicate"))
        )
      }
    })
    
    header_row <- if (show_tr) {
      fluidRow(
        column(3, tags$b("Detected Run ID")),
        column(3, tags$b("Condition")),
        column(3, tags$b("Bio Rep")),
        column(3, tags$b("Tech Rep"))
      )
    } else {
      fluidRow(
        column(4, tags$b("Detected Run ID")),
        column(4, tags$b("Condition")),
        column(4, tags$b("Replicate"))
      )
    }
    
    tagList(
      header_row,
      hr(style = "margin-top: 5px; margin-bottom: 10px;"),
      mapping_rows
    )
  })
  
  # Apply Sample Mapping Customizations
  observeEvent(input$btn_process_import, {
    req(rv$raw_data)
    unique_ids <- unique(rv$raw_data$ID)
    show_tr    <- isTRUE(input$has_tech_reps)
    
    mapping_df <- do.call(rbind, lapply(unique_ids, function(id) {
      cond_val <- input[[paste0("cond_", id)]]
      br_val   <- input[[paste0("br_", id)]]
      tr_val   <- if (show_tr) input[[paste0("tr_", id)]] else "1"
      
      if (is.null(cond_val) || cond_val == "") cond_val <- "Unspecified"
      if (is.null(br_val)   || br_val == "")   br_val   <- "1"
      if (is.null(tr_val)   || tr_val == "")   tr_val   <- "1"
      
      data.frame(
        ID             = id,
        user_condition = cond_val,
        user_BR        = br_val,
        user_TR        = tr_val,
        # Unique biological sample ID
        bio_sample_ID  = paste0(cond_val, br_val),
        stringsAsFactors = FALSE
      )
    }))
    
    # Merge mappings with raw dataset
    joined_data <- rv$raw_data %>%
      select(-any_of(c("condition", "BR", "TR"))) %>%
      left_join(mapping_df, by = "ID")
    
    if (show_tr) {
      # Log transformation -> Average across technical replicates per biological unit
      processed_data <- joined_data %>%
        group_by(Protein.ID, user_condition, user_BR, bio_sample_ID) %>%
        summarise(
          # Average log2 intensities; if all TRs are NA, returns NaN -> cast to NA
          LogInt = {
            valid_vals <- LogInt[!is.na(LogInt) & is.finite(LogInt)]
            if (length(valid_vals) > 0) mean(valid_vals) else NA_real_
          },
          .groups = "drop"
        ) %>%
        mutate(
          Intensity = ifelse(is.na(LogInt), NA_real_, 2^LogInt),
          condition = user_condition,
          BR        = user_BR,
          ID        = bio_sample_ID
        ) %>%
        select(Protein.ID, ID, Intensity, condition, BR, LogInt)
      
      # Audit trail metadata
      rv$audit_log$technical_replicate_averaging <- list(
        enabled = TRUE,
        total_runs_ingested = length(unique_ids),
        consolidated_samples = length(unique(processed_data$ID))
      )
    } else {
      # Standard 1:1 run-to-sample mapping
      processed_data <- joined_data %>%
        mutate(
          condition = user_condition,
          BR        = user_BR,
          ID        = bio_sample_ID
        ) %>%
        select(Protein.ID, ID, Intensity, condition, BR, LogInt)
      
      rv$audit_log$technical_replicate_averaging <- list(enabled = FALSE)
    }
    
    rv$raw_data   <- processed_data
    rv$sample_map <- mapping_df %>%
      rename(new_ID = bio_sample_ID)
    rv$audit_log$sample_mapping <- rv$sample_map
    
    # Rebuild QFeatures container using consolidated biological sample columns
    wide_data_mapped <- rv$raw_data %>%
      select(Protein.ID, ID, LogInt) %>%
      pivot_wider(names_from = ID, values_from = LogInt)
    
    rv$Qprot <- readQFeatures2(
      wide_data_mapped, 
      ecol = 2:ncol(wide_data_mapped), 
      fnames = "Protein.ID", 
      name = "raw"
    )
    
    if (show_tr) {
      showNotification(
        paste0("Technical replicates averaged. Consolidated into ", length(unique(rv$raw_data$ID)), " biological samples."),
        type = "message",
        duration = 5
      )
    } else {
      showNotification("Sample metadata successfully mapped!", type = "message")
    }
  })
  
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
    
    # ----------------------------------------------------------------------------
    # 1. Reference Protein Normalization (Prior Step)
    # ----------------------------------------------------------------------------
    if (isTRUE(input$enable_ref_norm) && !is.null(rv$geneDict)) {
      query       <- trimws(input$ref_target_query)
      target_type <- input$ref_target_type
      
      if (target_type == "gene") {
        matched <- rv$geneDict %>% filter(toupper(gene) == toupper(query))
      } else {
        matched <- rv$geneDict %>% filter(toupper(Accession) == toupper(query))
      }
      
      if (nrow(matched) == 0) {
        showNotification(
          paste0("Reference protein '", query, "' not found. Skipping reference normalization."),
          type = "warning",
          duration = 6
        )
      } else {
        ref_protein_id <- matched$Accession[1]
        ref_row <- wide_data %>% filter(Protein.ID == ref_protein_id)
        
        if (nrow(ref_row) == 0) {
          showNotification("Reference protein has no quantification data.", type = "error")
        } else {
          sample_cols <- setdiff(names(wide_data), "Protein.ID")
          ref_vals <- as.numeric(ref_row[1, sample_cols])
          
          if (any(is.na(ref_vals))) {
            showNotification(
              "Warning: Reference protein contains missing values (NAs) in some samples. Offset may be incomplete.",
              type = "warning",
              duration = 6
            )
          }
          
          global_ref_mean <- mean(ref_vals, na.rm = TRUE)
          sample_offsets  <- ref_vals - global_ref_mean
          
          for (idx in seq_along(sample_cols)) {
            col_name <- sample_cols[idx]
            offset_val <- sample_offsets[idx]
            if (!is.na(offset_val)) {
              wide_data[[col_name]] <- wide_data[[col_name]] - offset_val
            }
          }
          
          rv$audit_log$reference_protein_norm <- list(
            enabled   = TRUE,
            target    = query,
            accession = ref_protein_id,
            gene      = matched$gene[1]
          )
        }
      }
    } else {
      rv$audit_log$reference_protein_norm <- list(enabled = FALSE)
    }
    
    # ----------------------------------------------------------------------------
    # 2. Global Normalization (QFeatures)
    # ----------------------------------------------------------------------------
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
    
    # ----------------------------------------------------------------------------
    # 3. Imputation
    # ----------------------------------------------------------------------------
    if (input$impute_method == "Hybrid (MAR: KNN / MNAR: MinDet)") {
      global_mar <- MsCoreUtils::impute_matrix(norm_mat, method = "nbavg")
      global_min <- min(norm_mat, na.rm = TRUE)
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
          mar_cells_count <- mar_cells_count + sum(is.na(sub_m[mar_proteins, , drop = FALSE]))
          imputed_mat[mar_proteins, cols] <- global_mar[mar_proteins, cols]
        }
        
        mnar_proteins <- names(present[present < mar_threshold])
        if (length(mnar_proteins) > 0) {
          mnar_cells_count <- mnar_cells_count + sum(is.na(sub_m[mnar_proteins, , drop = FALSE]))
          for (p in mnar_proteins) {
            row_vals <- sub_m[p, ]
            if (all(is.na(row_vals))) {
              imputed_mat[p, cols] <- global_min - 0.5
            } else {
              min_val <- min(row_vals, na.rm = TRUE)
              imputed_mat[p, cols] <- ifelse(is.na(row_vals), min_val - 0.5, row_vals)
            }
          }
        }
      }
      norm_mat <- imputed_mat
      rv$impute_stats <- list(method = input$impute_method, total_na = total_nas, mar = mar_cells_count, mnar = mnar_cells_count)
      
    } else if (!input$impute_method %in% c("None", "No Imputation (Show 1-Condition Dropouts on Margins)")) {
      norm_mat <- MsCoreUtils::impute_matrix(norm_mat, method = input$impute_method)
      rv$impute_stats <- list(method = input$impute_method, total_na = total_nas, imputed_total = total_nas)
    } else {
      rv$impute_stats <- list(method = "None", total_na = total_nas)
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
      paste0("Method: ", stats$method, "\nTotal Values Imputed: ", stats$total_na)
    }
  })
  
  output$norm_qc_plot <- renderPlot({
    req(rv$Qprot)
    plot_normalization_qc(rv$Qprot, i = 1)
  })
  
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
          ExpQuant  = mean(Quant[Condition == exp_cond], na.rm = TRUE),
          RefQuant  = mean(Quant[Condition == ref_cond], na.rm = TRUE),
          Exp_Valid = sum(!is.na(Quant[Condition == exp_cond])),
          Ref_Valid = sum(!is.na(Quant[Condition == ref_cond])),
          .groups   = "drop"
        )
      
      # 1. Join scatter means ONCE to the fitted limma table
      dt <- dt %>% left_join(scatter_means, by = "Accession")
      
      # 2. Append 1-condition dropouts if requested
      if (input$impute_method == "No Imputation (Show 1-Condition Dropouts on Margins)") {
        dropout_rows <- scatter_means %>%
          filter((Exp_Valid > 0 & Ref_Valid == 0) | (Exp_Valid == 0 & Ref_Valid > 0)) %>%
          filter(!Accession %in% dt$Accession) %>%
          mutate(
            logFC     = ifelse(Exp_Valid > 0, Inf, -Inf),
            AveExpr   = ifelse(Exp_Valid > 0, ExpQuant, RefQuant),
            t         = NA_real_,
            P.Value   = NA_real_,
            adj.P.Val = NA_real_,
            B         = NA_real_,
            CI.L      = NA_real_,
            CI.R      = NA_real_
          )
        
        dt <- bind_rows(dt, dropout_rows)
      }
      
      # 3. Cleanly attach gene dictionary metadata
      dt <- dt %>% left_join(rv$geneDict, by = "Accession")
      
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
  # server.R (Updated base_ggplot block)
  
  base_ggplot <- reactive({
    req(rv$de_results, input$selected_comparison)
    df <- rv$de_results[[input$selected_comparison]]
    
    exp_name <- attr(df, "exp_cond") %||% "Experimental"
    ref_name <- attr(df, "ref_cond") %||% "Reference"
    
    if (is_flipped()) {
      df <- df %>%
        mutate(
          logFC = -logFC,
          temp_exp = ExpQuant,
          ExpQuant = RefQuant,
          RefQuant = temp_exp
        ) %>%
        select(-temp_exp)
      
      temp_name <- exp_name
      exp_name  <- ref_name
      ref_name  <- temp_name
    }
    
    raw_p_cutoff <- input$adj_p_cutoff
    log10_p_line <- -log10(raw_p_cutoff)
    
    # 1. Base DE Status
    show_de <- isTRUE(input$show_de_colors)
    df <- df %>%
      mutate(
        significant = !is.na(adj.P.Val) & is.finite(logFC) & abs(logFC) > input$fc_cutoff & adj.P.Val < raw_p_cutoff,
        Regulation = case_when(
          !show_de                 ~ "Background",
          is.infinite(logFC)       ~ "Dropout (1-Condition)",
          significant & logFC > 0  ~ "Upregulated",
          significant & logFC <= 0 ~ "Downregulated",
          TRUE                     ~ "Not Significant"
        )
      )
    
    # 2. Assign Custom Protein Set Overrides
    df <- assign_custom_protein_groups(df, input, rv$num_custom_sets)
    
    # Define final plotting factor: Custom Set overrides DE Status if matched & color enabled
    df <- df %>%
      mutate(
        Final_Group = ifelse(!is.na(Custom_Group), Custom_Group, Regulation)
      )
    
    # 3. Assemble Dynamic Color Palette
    color_map <- c(
      "Upregulated"           = input$col_up,
      "Downregulated"         = input$col_down,
      "Dropout (1-Condition)" = "grey40",
      "Not Significant"       = "grey75",
      "Background"            = "grey75"
    )
    
    if (rv$num_custom_sets > 0) {
      for (i in seq_len(rv$num_custom_sets)) {
        set_name <- input[[paste0("custom_set_name_", i)]] %||% paste("Set", i)
        set_col  <- input[[paste0("custom_set_col_", i)]]
        if (!is.null(set_col)) color_map[set_name] <- set_col
      }
    }
    
    df$Final_Group <- factor(df$Final_Group, levels = c(names(color_map)))
    
    # ----------------------------------------------------------------------------
    # VOLCANO PLOT
    # ----------------------------------------------------------------------------
    if (input$plot_type == "Volcano") {
      plot_df <- df %>% filter(is.finite(logFC) & !is.na(adj.P.Val))
      
      p <- ggplot(plot_df, aes(x = logFC, y = -log10(adj.P.Val), text = paste("Gene:", gene, "<br>Accession:", Accession))) +
        geom_point(aes(color = Final_Group), size = input$point_size) +
        geom_hline(yintercept = log10_p_line, linetype = 2, color = "grey50") +
        geom_vline(xintercept = c(-input$fc_cutoff, input$fc_cutoff), linetype = 2, color = "grey50") +
        labs(
          title = paste("Volcano Plot:", input$selected_comparison), 
          x = paste0("Log2 Fold Change (", exp_name, " / ", ref_name, ")"), 
          y = "-Log10 Adjusted p-value"
        )
      
      # ----------------------------------------------------------------------------
      # MA PLOT
      # ----------------------------------------------------------------------------
    } else if (input$plot_type == "MA") {
      finite_fc <- df$logFC[!is.na(df$logFC) & is.finite(df$logFC)]
      y_cap <- if (length(finite_fc) > 0) max(abs(finite_fc), na.rm = TRUE) + 2 else 6
      
      plot_df <- df %>%
        mutate(
          has_exp = !is.na(ExpQuant) & is.finite(ExpQuant),
          has_ref = !is.na(RefQuant) & is.finite(RefQuant),
          AveExpr = case_when(
            has_exp & has_ref ~ ifelse(!is.na(AveExpr) & is.finite(AveExpr), AveExpr, (ExpQuant + RefQuant) / 2),
            has_exp & !has_ref ~ ExpQuant,
            !has_exp & has_ref ~ RefQuant,
            TRUE ~ NA_real_
          ),
          plot_logFC = case_when(
            has_exp & !has_ref ~ y_cap,
            !has_exp & has_ref ~ -y_cap,
            is.infinite(logFC) & logFC > 0 ~ y_cap,
            is.infinite(logFC) & logFC < 0 ~ -y_cap,
            TRUE ~ logFC
          )
        ) %>%
        filter(!is.na(AveExpr) & !is.na(plot_logFC))
      
      p <- ggplot(plot_df, aes(x = AveExpr, y = plot_logFC, text = paste("Gene:", gene, "<br>Accession:", Accession))) +
        geom_point(aes(color = Final_Group), size = input$point_size) +
        geom_hline(yintercept = 0, color = "grey30") +
        geom_hline(yintercept = c(-input$fc_cutoff, input$fc_cutoff), linetype = 2, color = "grey50") +
        geom_hline(yintercept = c(-y_cap, y_cap), linetype = 3, color = "grey70") +
        coord_cartesian(ylim = c(-y_cap - 0.5, y_cap + 0.5)) +
        labs(
          title = paste("MA Plot:", input$selected_comparison), 
          x = "Log2 Average Expression (Valid Group Mean)", 
          y = paste0("Log2 Fold Change (", exp_name, " / ", ref_name, ")")
        )
      
      # ----------------------------------------------------------------------------
      # SCATTER PLOT
      # ----------------------------------------------------------------------------
    } else if (input$plot_type == "Scatter") {
      all_quants <- c(df$RefQuant[is.finite(df$RefQuant)], df$ExpQuant[is.finite(df$ExpQuant)])
      axis_floor <- if (length(all_quants) > 0) min(all_quants, na.rm = TRUE) - 1 else 0
      
      plot_df <- df %>%
        mutate(
          plot_Ref = ifelse(is.na(RefQuant) | !is.finite(RefQuant), axis_floor, RefQuant),
          plot_Exp = ifelse(is.na(ExpQuant) | !is.finite(ExpQuant), axis_floor, ExpQuant)
        )
      
      p <- ggplot(plot_df, aes(x = plot_Ref, y = plot_Exp, text = paste("Gene:", gene, "<br>Accession:", Accession))) +
        geom_point(aes(color = Final_Group), size = input$point_size) +
        geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey50") +
        labs(
          title = paste("Scatter Plot:", input$selected_comparison), 
          x = paste("Log2", ref_name, "Average Abundance"), 
          y = paste("Log2", exp_name, "Average Abundance")
        )
    }
    
    p <- p + 
      scale_color_manual(name = "Group", values = color_map, drop = TRUE) +
      theme_bw(base_size = input$text_size) +
      theme(panel.grid = element_blank(), legend.position = "right")
    
    p
  })
  
  output$plotly_view <- renderPlotly({
    req(base_ggplot())
    ggplotly(base_ggplot(), tooltip = "text") %>%
      layout(
        autosize = TRUE,
        legend = list(title = list(text = "Expression Status"))
      )
  })
  
  # server.R (Updated final_ggplot_object with Dual Labeling Support)
  
  final_ggplot_object <- reactive({
    lbl_size    <- input$label_size
    max_ovrlaps <- input$max_overlaps
    raw_p_cut   <- input$adj_p_cutoff
    fc_cut      <- input$fc_cutoff
    show_de_lbl <- isTRUE(input$show_de_labels)
    
    p <- base_ggplot()
    req(rv$de_results, input$selected_comparison)
    
    df <- rv$de_results[[input$selected_comparison]]
    if (is_flipped()) {
      df <- df %>%
        mutate(
          logFC = -logFC,
          temp_exp = ExpQuant,
          ExpQuant = RefQuant,
          RefQuant = temp_exp
        ) %>%
        select(-temp_exp)
    }
    
    # Determine DE significant status
    df$de_sig <- !is.na(df$adj.P.Val) & is.finite(df$logFC) & abs(df$logFC) > fc_cut & df$adj.P.Val < raw_p_cut
    
    # Determine Custom set labels
    df <- assign_custom_protein_groups(df, input, rv$num_custom_sets)
    
    # Target proteins to label: (DE proteins IF toggled) OR (Custom proteins IF toggled)
    df <- df %>%
      mutate(to_label = (show_de_lbl & de_sig) | Custom_Label)
    
    sig_df <- df %>% filter(to_label & !is.na(gene) & gene != "")
    
    if (nrow(sig_df) > 0) {
      if (input$plot_type == "Volcano") {
        p <- p + geom_label_repel(
          data = sig_df %>% filter(is.finite(logFC) & !is.na(adj.P.Val)),
          aes(x = logFC, y = -log10(adj.P.Val), label = gene, group = paste0(gene, "_", max_ovrlaps)),
          size = lbl_size,
          max.overlaps = max_ovrlaps,
          show.legend = FALSE,
          inherit.aes = FALSE
        )
      } else if (input$plot_type == "MA") {
        p <- p + geom_label_repel(
          data = sig_df,
          aes(x = AveExpr, y = logFC, label = gene, group = paste0(gene, "_", max_ovrlaps)),
          size = lbl_size,
          max.overlaps = max_ovrlaps,
          show.legend = FALSE,
          inherit.aes = FALSE
        )
      } else if (input$plot_type == "Scatter") {
        p <- p + geom_label_repel(
          data = sig_df,
          aes(x = RefQuant, y = ExpQuant, label = gene, group = paste0(gene, "_", max_ovrlaps)),
          size = lbl_size,
          max.overlaps = max_ovrlaps,
          show.legend = FALSE,
          inherit.aes = FALSE
        )
      }
    }
    
    return(p)
  })
  
  output$ggplot_view <- renderPlot({
    final_ggplot_object()
  })
  
  # Contrasting default colors (Amber Orange, Vivid Green, Vivid Purple, Golden Yellow, Bright Teal)
  default_custom_colors <- c("#FF9900", "#33CC33", "#9933FF", "#FFCC00", "#00CCCC")
  
  # Initialize custom set counter inside rv
  rv$num_custom_sets <- 1
  
  # Add / Remove set event handlers
  observeEvent(input$btn_add_custom_set, {
    rv$num_custom_sets <- rv$num_custom_sets + 1
  })
  
  observeEvent(input$btn_remove_custom_set, {
    if (rv$num_custom_sets > 0) {
      rv$num_custom_sets <- rv$num_custom_sets - 1
    }
  })
  
  # Dynamic UI Renderer for Custom Sets
  # server.R (Updated output$custom_protein_sets_ui with input persistence)
  
  output$custom_protein_sets_ui <- renderUI({
    n <- rv$num_custom_sets
    if (n == 0) return(tags$em("No custom protein sets active. Click '+ Add Protein Set' to create one."))
    
    lapply(seq_len(n), function(i) {
      # 1. Grab previously entered values if they exist, otherwise use initial defaults
      existing_name <- isolate(input[[paste0("custom_set_name_", i)]])
      curr_name     <- if (!is.null(existing_name)) existing_name else paste("Set", i)
      
      existing_type <- isolate(input[[paste0("custom_set_type_", i)]])
      curr_type     <- if (!is.null(existing_type)) existing_type else "list"
      
      existing_input <- isolate(input[[paste0("custom_set_input_", i)]])
      curr_input     <- if (!is.null(existing_input)) existing_input else ""
      
      existing_col <- isolate(input[[paste0("custom_set_col_", i)]])
      def_col      <- default_custom_colors[((i - 1) %% length(default_custom_colors)) + 1]
      curr_col     <- if (!is.null(existing_col) && existing_col != "") existing_col else def_col
      
      existing_show_col <- isolate(input[[paste0("custom_set_show_col_", i)]])
      curr_show_col     <- if (!is.null(existing_show_col)) existing_show_col else TRUE
      
      existing_show_lbl <- isolate(input[[paste0("custom_set_show_lbl_", i)]])
      curr_show_lbl     <- if (!is.null(existing_show_lbl)) existing_show_lbl else TRUE
      
      # 2. Render UI Card with preserved values
      tags$div(
        style = "border: 1px solid #e3e3e3; border-radius: 5px; padding: 10px; margin-bottom: 10px; background-color: #fafafa;",
        fluidRow(
          column(8, tags$strong(paste("Set", i, "Name / Pattern:"))),
          column(4, textInput(paste0("custom_set_name_", i), label = NULL, value = curr_name, placeholder = "Set Name"))
        ),
        radioButtons(
          inputId  = paste0("custom_set_type_", i),
          label    = NULL,
          choices  = c("List (Genes/Accessions)" = "list", "Regex Match" = "regex"),
          inline   = TRUE,
          selected = curr_type
        ),
        textAreaInput(
          inputId     = paste0("custom_set_input_", i),
          label       = NULL,
          value       = curr_input,
          rows        = 2,
          placeholder = "e.g., ACTB, GAPDH, P04406 or ^HIST.*"
        ),
        fluidRow(
          column(6, colourpicker::colourInput(paste0("custom_set_col_", i), "Color:", value = curr_col)),
          column(6,
                 checkboxInput(paste0("custom_set_show_col_", i), "Highlight Color", value = curr_show_col),
                 checkboxInput(paste0("custom_set_show_lbl_", i), "Show Labels", value = curr_show_lbl)
          )
        )
      )
    })
  })
  
  output$ref_norm_match_status <- renderText({
    req(rv$raw_data, rv$geneDict, input$enable_ref_norm)
    query <- trimws(input$ref_target_query)
    if (query == "") return("Please enter a gene or accession.")
    
    target_type <- input$ref_target_type
    
    if (target_type == "gene") {
      matched <- rv$geneDict %>% filter(toupper(gene) == toupper(query))
    } else {
      matched <- rv$geneDict %>% filter(toupper(Accession) == toupper(query))
    }
    
    if (nrow(matched) == 0) {
      return(paste0("❌ No match found for: '", query, "'"))
    } else {
      return(paste0("✓ Matched: ", matched$gene[1], " (", matched$Accession[1], ")"))
    }
  })
  
  # ----------------------------------------------------------------------------
  # KEGG & GSEA Module Reactives & Observers
  # ----------------------------------------------------------------------------
  rv_enrich <- reactiveValues(
    raw_excel_sheets = list(),
    contrast_df = NULL,
    kegg_res = NULL,
    gsea_res = NULL
  )
  
  # 1. Read Excel if uploaded
  observeEvent(input$enrich_file_upload, {
    req(input$enrich_file_upload)
    sheets <- readxl::excel_sheets(input$enrich_file_upload$datapath)
    sheet_list <- lapply(sheets, function(s) {
      readxl::read_excel(input$enrich_file_upload$datapath, sheet = s)
    })
    names(sheet_list) <- sheets
    rv_enrich$raw_excel_sheets <- sheet_list
  })
  
  # 2. Dynamic Contrast / Sheet Selector
  output$enrich_contrast_selector <- renderUI({
    if (input$enrich_data_source == "app") {
      req(rv$de_results)
      selectInput("enrich_selected_contrast", "Select Contrast:", choices = names(rv$de_results))
    } else {
      req(length(rv_enrich$raw_excel_sheets) > 0)
      selectInput("enrich_selected_contrast", "Select Excel Sheet / Contrast:", choices = names(rv_enrich$raw_excel_sheets))
    }
  })
  
  # 3. Pull Contrast Data (Respects 'invert contrast direction' toggle)
  observe({
    req(input$enrich_selected_contrast)
    
    if (input$enrich_data_source == "app") {
      req(rv$de_results[[input$enrich_selected_contrast]])
      df <- rv$de_results[[input$enrich_selected_contrast]]
      
      # If inverted in the visualization tab, invert logFC and condition means
      if (is_flipped()) {
        df <- df %>%
          mutate(
            logFC = -logFC,
            temp_exp = ExpQuant,
            ExpQuant = RefQuant,
            RefQuant = temp_exp
          ) %>%
          select(-temp_exp)
      }
      
      rv_enrich$contrast_df <- df
    } else {
      req(rv_enrich$raw_excel_sheets[[input$enrich_selected_contrast]])
      rv_enrich$contrast_df <- rv_enrich$raw_excel_sheets[[input$enrich_selected_contrast]]
    }
  })
  
  # 4. Run Analysis Observer
  observeEvent(input$btn_run_enrichment, {
    req(rv_enrich$contrast_df)
    df <- as.data.frame(rv_enrich$contrast_df)
    
    # Identify key columns
    acc_col <- intersect(c("Accession", "Protein.ID", "UNIPROT"), names(df))[1]
    fc_col  <- intersect(c("logFC", "Log2FC"), names(df))[1]
    p_col   <- intersect(c("adj.P.Val", "FDR", "pvalue"), names(df))[1]
    
    if (is.na(acc_col) || is.na(fc_col)) {
      showNotification("Missing 'Accession' or 'logFC' columns in selected table.", type = "error")
      return()
    }
    
    # Select organism OrgDb
    org_code <- input$enrich_organism
    org_db <- if (org_code == "mmu") org.Mm.eg.db::org.Mm.eg.db else org.Hs.eg.db::org.Hs.eg.db
    
    # Clean UniProt IDs (remove isoform suffixes or pipes like sp|P12345|...)
    raw_acc <- as.character(df[[acc_col]])
    clean_acc <- gsub("^.*\\|([A-Za-z0-9]+)\\|.*$", "\\1", raw_acc)
    clean_acc <- gsub("-.*$", "", clean_acc)
    df$clean_accession <- clean_acc
    
    withProgress(message = "Running Pathway Enrichment...", value = 0.2, {
      
      # Map UniProt to Entrez ID
      incProgress(0.2, detail = "Mapping UniProt to Entrez IDs...")
      mapped_ids <- tryCatch({
        clusterProfiler::bitr(
          unique(df$clean_accession),
          fromType = "UNIPROT",
          toType   = "ENTREZID",
          OrgDb    = org_db
        )
      }, error = function(e) NULL)
      
      if (is.null(mapped_ids) || nrow(mapped_ids) == 0) {
        showNotification("Could not map Accessions to Entrez IDs for this organism.", type = "error")
        return()
      }
      
      df_mapped <- merge(df, mapped_ids, by.x = "clean_accession", by.y = "UNIPROT")
      
      # ------------------------------------------------------------------------
      # A. KEGG Over-Representation Analysis (ORA)
      # ------------------------------------------------------------------------
      if (input$enrich_method == "kegg") {
        incProgress(0.4, detail = "Calculating KEGG over-representation...")
        p_cut  <- input$enrich_p_cutoff
        fc_cut <- input$enrich_fc_cutoff
        
        # Direction filter
        if (input$enrich_kegg_dir == "up") {
          target_acc <- df$clean_accession[!is.na(df[[p_col]]) & df[[p_col]] < p_cut & df[[fc_col]] > fc_cut]
        } else if (input$enrich_kegg_dir == "down") {
          target_acc <- df$clean_accession[!is.na(df[[p_col]]) & df[[p_col]] < p_cut & df[[fc_col]] < -fc_cut]
        } else {
          target_acc <- df$clean_accession[!is.na(df[[p_col]]) & df[[p_col]] < p_cut & abs(df[[fc_col]]) > fc_cut]
        }
        
        target_entrez <- mapped_ids$ENTREZID[mapped_ids$UNIPROT %in% target_acc]
        universe_entrez <- unique(mapped_ids$ENTREZID)
        
        if (length(target_entrez) < 5) {
          showNotification("Fewer than 5 significant proteins mapped to Entrez IDs. Try relaxing cutoffs.", type = "warning")
        }
        
        kegg_out <- tryCatch({
          clusterProfiler::enrichKEGG(
            gene          = unique(target_entrez),
            universe      = universe_entrez,
            organism      = org_code,
            keyType       = "kegg",
            pAdjustMethod = "BH",
            pvalueCutoff  = input$enrich_kegg_pvalue_cut,
            qvalueCutoff  = 0.2
          )
        }, error = function(e) {
          showNotification(paste("KEGG enrichment error:", e$message), type = "error")
          NULL
        })
        
        rv_enrich$kegg_res <- kegg_out
        rv_enrich$gsea_res <- NULL
        
        # ------------------------------------------------------------------------
        # B. KEGG Gene Set Enrichment Analysis (GSEA)
        # ------------------------------------------------------------------------
      } else if (input$enrich_method == "gsea") {
        incProgress(0.4, detail = "Ranking genes and executing GSEA...")
        
        # Sort by absolute fold change to resolve duplicates, keeping top magnitude per Entrez
        df_ranked <- df_mapped[order(abs(df_mapped[[fc_col]]), decreasing = TRUE), ]
        df_ranked <- df_ranked[!duplicated(df_ranked$ENTREZID), ]
        
        gene_list <- df_ranked[[fc_col]]
        names(gene_list) <- df_ranked$ENTREZID
        gene_list <- sort(gene_list, decreasing = TRUE)
        
        gsea_out <- tryCatch({
          clusterProfiler::gseKEGG(
            geneList      = gene_list,
            organism      = org_code,
            pvalueCutoff  = input$enrich_gsea_pvalue_cut,
            pAdjustMethod = "BH",
            verbose       = FALSE
          )
        }, error = function(e) {
          showNotification(paste("GSEA execution error:", e$message), type = "error")
          NULL
        })
        
        rv_enrich$gsea_res <- gsea_out
        rv_enrich$kegg_res <- NULL
      }
      
      incProgress(0.2, detail = "Complete!")
    })
    
    showNotification("Analysis completed successfully!", type = "message")
  })
  
  # 5. Populate GSEA Pathway Dropdown
  output$gsea_pathway_selector <- renderUI({
    req(rv_enrich$gsea_res)
    res_df <- as.data.frame(rv_enrich$gsea_res)
    if (nrow(res_df) == 0) return(tags$em("No enriched pathways found at current cutoff."))
    
    choices_vec <- setNames(res_df$ID, paste0(res_df$Description, " (", res_df$ID, ")"))
    selectInput("gsea_selected_pathway", "Select Pathway to View:", choices = choices_vec)
  })
  
  # 6. Render Dot Plot (Reactively handles ORA vs GSEA)
  output$enrich_dotplot <- renderPlot({
    if (input$enrich_method == "kegg") {
      req(rv_enrich$kegg_res)
      if (nrow(as.data.frame(rv_enrich$kegg_res)) == 0) {
        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "No significant KEGG pathways found.")
        return()
      }
      enrichplot::dotplot(rv_enrich$kegg_res, showCategory = 15, title = "KEGG Pathway Over-Representation")
    } else {
      req(rv_enrich$gsea_res)
      if (nrow(as.data.frame(rv_enrich$gsea_res)) == 0) {
        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "No significant GSEA pathways found.")
        return()
      }
      enrichplot::dotplot(rv_enrich$gsea_res, showCategory = 10, split = ".sign") + 
        ggplot2::facet_grid(. ~ .sign) +
        ggplot2::labs(title = "KEGG GSEA Pathway Enrichment (Activated vs Suppressed)")
    }
  })
  
  # 7. Render Single Pathway GSEA Plot
  output$gsea_single_plot <- renderPlot({
    req(rv_enrich$gsea_res, input$gsea_selected_pathway)
    path_id <- input$gsea_selected_pathway
    res_df <- as.data.frame(rv_enrich$gsea_res)
    desc <- res_df$Description[res_df$ID == path_id][1]
    
    enrichplot::gseaplot2(rv_enrich$gsea_res, geneSetID = path_id, title = desc)
  })
  
  # 8. Render Results Table
  output$enrich_table <- renderDT({
    res_obj <- if (input$enrich_method == "kegg") rv_enrich$kegg_res else rv_enrich$gsea_res
    req(res_obj)
    df <- as.data.frame(res_obj)
    datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # 9. Download Handlers
  output$download_enrich_dotplot_png <- downloadHandler(
    filename = function() { paste0("Enrichment_dotplot_", input$enrich_method, ".png") },
    content = function(file) {
      p <- if (input$enrich_method == "kegg") {
        enrichplot::dotplot(rv_enrich$kegg_res, showCategory = 15, title = "KEGG Pathway Over-Representation")
      } else {
        enrichplot::dotplot(rv_enrich$gsea_res, showCategory = 10, split = ".sign") + ggplot2::facet_grid(. ~ .sign)
      }
      ggplot2::ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
    }
  )
  
  output$download_gsea_pathway_png <- downloadHandler(
    filename = function() { paste0("GSEA_Pathway_", input$gsea_selected_pathway, ".png") },
    content = function(file) {
      req(rv_enrich$gsea_res, input$gsea_selected_pathway)
      path_id <- input$gsea_selected_pathway
      res_df <- as.data.frame(rv_enrich$gsea_res)
      desc <- res_df$Description[res_df$ID == path_id][1]
      p <- enrichplot::gseaplot2(rv_enrich$gsea_res, geneSetID = path_id, title = desc)
      ggplot2::ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
    }
  )
  
  # ----------------------------------------------------------------------------
  # GSEA Pathway-Specific DE Plot (Volcano / MA / Scatter)
  # ----------------------------------------------------------------------------
  pathway_de_ggplot <- reactive({
    req(rv_enrich$gsea_res, input$gsea_selected_pathway, rv_enrich$contrast_df)
    
    gsea_obj <- rv_enrich$gsea_res
    path_id  <- input$gsea_selected_pathway
    df       <- as.data.frame(rv_enrich$contrast_df)
    
    # 1. Identify key columns in the active contrast dataframe
    acc_col <- intersect(c("Accession", "Protein.ID", "UNIPROT"), names(df))[1]
    fc_col  <- intersect(c("logFC", "Log2FC"), names(df))[1]
    p_col   <- intersect(c("adj.P.Val", "FDR", "pvalue"), names(df))[1]
    gene_col <- intersect(c("gene", "Gene", "Gene.Name", "PG.Genes"), names(df))[1]
    
    req(!is.na(acc_col), !is.na(fc_col))
    
    df$Accession_clean <- gsub("^.*\\|([A-Za-z0-9]+)\\|.*$", "\\1", as.character(df[[acc_col]]))
    df$Accession_clean <- gsub("-.*$", "", df$Accession_clean)
    
    df$gene_symbol <- if (!is.na(gene_col) && gene_col %in% names(df)) {
      as.character(df[[gene_col]])
    } else {
      df$Accession_clean
    }
    df$gene_symbol <- ifelse(is.na(df$gene_symbol) | df$gene_symbol == "", df$Accession_clean, df$gene_symbol)
    
    # 2. Extract Entrez IDs belonging to this specific pathway
    pathway_entrez <- gsea_obj@geneSets[[path_id]]
    req(length(pathway_entrez) > 0)
    
    # 3. Map pathway Entrez IDs back to UniProt using active OrgDb
    org_code <- input$enrich_organism
    org_db   <- if (org_code == "mmu") org.Mm.eg.db::org.Mm.eg.db else org.Hs.eg.db::org.Hs.eg.db
    
    mapped_back <- tryCatch({
      clusterProfiler::bitr(
        pathway_entrez,
        fromType = "ENTREZID",
        toType   = "UNIPROT",
        OrgDb    = org_db
      )
    }, error = function(e) NULL)
    
    pathway_uniprots <- if (!is.null(mapped_back)) unique(mapped_back$UNIPROT) else character(0)
    
    # 4. Mark pathway membership
    df <- df %>%
      mutate(
        is_pathway = Accession_clean %in% pathway_uniprots,
        plot_group = factor(
          ifelse(is_pathway, "Pathway Protein", "Other Proteins"),
          levels = c("Other Proteins", "Pathway Protein")
        )
      )
    
    # Extract path description for title
    res_df <- as.data.frame(gsea_obj)
    path_desc <- res_df$Description[res_df$ID == path_id][1]
    plot_title <- paste0(path_desc, " (", path_id, ")")
    
    # Inherit visual parameters from Tab 3 inputs (with safe fallbacks)
    pt_size     <- input$point_size %||% 2.5
    txt_size    <- input$text_size %||% 12
    lbl_size    <- input$label_size %||% 3.5
    max_ovrlaps <- input$max_overlaps %||% 15
    plot_layout <- input$plot_type %||% "Volcano"
    fc_cut      <- input$fc_cutoff %||% 1.0
    p_cut       <- input$adj_p_cutoff %||% 0.05
    
    # --------------------------------------------------------------------------
    # VOLCANO
    # --------------------------------------------------------------------------
    if (plot_layout == "Volcano") {
      plot_df <- df %>% filter(is.finite(.data[[fc_col]]) & !is.na(.data[[p_col]]))
      
      p <- ggplot(plot_df, aes(x = .data[[fc_col]], y = -log10(.data[[p_col]]))) +
        geom_point(aes(color = plot_group, size = plot_group, alpha = plot_group)) +
        geom_hline(yintercept = -log10(p_cut), linetype = 2, color = "grey60") +
        geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 2, color = "grey60") +
        labs(
          title = plot_title,
          x = "Log2 Fold Change",
          y = "-Log10 Adjusted p-value"
        )
      
      # --------------------------------------------------------------------------
      # MA PLOT
      # --------------------------------------------------------------------------
    } else if (plot_layout == "MA") {
      # Calculate AveExpr if missing
      exp_col <- intersect(names(df), c("ExpQuant", grep("_Mean$", names(df), value = TRUE)))[1]
      ref_col <- intersect(names(df), c("RefQuant", grep("_Mean$", names(df), value = TRUE)))[2]
      
      plot_df <- df %>%
        mutate(
          AveExpr = if ("AveExpr" %in% names(df) && !all(is.na(df$AveExpr))) {
            df$AveExpr
          } else if (!is.na(exp_col) && !is.na(ref_col)) {
            (as.numeric(.data[[exp_col]]) + as.numeric(.data[[ref_col]])) / 2
          } else {
            0
          }
        ) %>%
        filter(is.finite(.data[[fc_col]]) & is.finite(AveExpr))
      
      p <- ggplot(plot_df, aes(x = AveExpr, y = .data[[fc_col]])) +
        geom_point(aes(color = plot_group, size = plot_group, alpha = plot_group)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_hline(yintercept = c(-fc_cut, fc_cut), linetype = 2, color = "grey60") +
        labs(
          title = plot_title,
          x = "Average Log2 Expression",
          y = "Log2 Fold Change"
        )
      
      # --------------------------------------------------------------------------
      # SCATTER PLOT
      # --------------------------------------------------------------------------
    } else {
      exp_col <- intersect(names(df), c("ExpQuant", grep("_Mean$", names(df), value = TRUE)))[1]
      ref_col <- intersect(names(df), c("RefQuant", grep("_Mean$", names(df), value = TRUE)))[2]
      
      req(!is.na(exp_col), !is.na(ref_col))
      plot_df <- df %>% filter(is.finite(.data[[ref_col]]) & is.finite(.data[[exp_col]]))
      
      p <- ggplot(plot_df, aes(x = .data[[ref_col]], y = .data[[exp_col]])) +
        geom_point(aes(color = plot_group, size = plot_group, alpha = plot_group)) +
        geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey60") +
        labs(
          title = plot_title,
          x = "Reference Abundance",
          y = "Experimental Abundance"
        )
    }
    
    # 5. Styling: Distinct Color, Size, and Order
    p <- p +
      scale_color_manual(
        name   = "Status",
        values = c("Other Proteins" = "grey82", "Pathway Protein" = "#E69F00")
      ) +
      scale_size_manual(
        name   = "Status",
        values = c("Other Proteins" = pt_size * 0.8, "Pathway Protein" = pt_size * 1.5)
      ) +
      scale_alpha_manual(
        name   = "Status",
        values = c("Other Proteins" = 0.45, "Pathway Protein" = 1.0)
      ) +
      theme_bw(base_size = txt_size) +
      theme(
        panel.grid      = element_blank(),
        legend.position = "bottom",
        plot.title      = element_text(face = "bold", size = txt_size)
      )
    
    # 6. Add labels specifically for the pathway proteins
    pathway_pts <- plot_df %>% filter(is_pathway)
    if (nrow(pathway_pts) > 0) {
      if (plot_layout == "Volcano") {
        p <- p + ggrepel::geom_label_repel(
          data          = pathway_pts,
          aes(x = .data[[fc_col]], y = -log10(.data[[p_col]]), label = gene_symbol),
          size          = lbl_size,
          max.overlaps  = max_ovrlaps,
          box.padding   = 0.35,
          color         = "black",
          fill          = alpha("white", 0.85),
          show.legend   = FALSE
        )
      } else if (plot_layout == "MA") {
        p <- p + ggrepel::geom_label_repel(
          data          = pathway_pts,
          aes(x = AveExpr, y = .data[[fc_col]], label = gene_symbol),
          size          = lbl_size,
          max.overlaps  = max_ovrlaps,
          box.padding   = 0.35,
          color         = "black",
          fill          = alpha("white", 0.85),
          show.legend   = FALSE
        )
      } else {
        p <- p + ggrepel::geom_label_repel(
          data          = pathway_pts,
          aes(x = .data[[ref_col]], y = .data[[exp_col]], label = gene_symbol),
          size          = lbl_size,
          max.overlaps  = max_ovrlaps,
          box.padding   = 0.35,
          color         = "black",
          fill          = alpha("white", 0.85),
          show.legend   = FALSE
        )
      }
    }
    
    p
  })
  
  # Render plot output
  output$gsea_pathway_de_plot <- renderPlot({
    pathway_de_ggplot()
  })
  
  # Download Handler for the new plot
  output$download_gsea_de_plot_png <- downloadHandler(
    filename = function() { paste0("Pathway_DE_", input$gsea_selected_pathway, ".png") },
    content = function(file) {
      ggplot2::ggsave(file, plot = pathway_de_ggplot(), width = 8, height = 7, dpi = 300)
    }
  )
  
  # ----------------------------------------------------------------------------
  # 5. Audit Log & Export Handlers
  # ----------------------------------------------------------------------------
  output$audit_preview <- renderText({
    toJSON(rv$audit_log, pretty = TRUE)
  })
  
  # server.R (Updated download_excel handler)
  
  # Excel Workbook Download Handler
  output$download_excel <- downloadHandler(
    filename = function() { 
      paste0(input$export_prefix, "_results.xlsx") 
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    content = function(file) {
      wb <- openxlsx::createWorkbook()
      flip <- is_flipped()
      
      for (comp in names(rv$de_results)) {
        df <- rv$de_results[[comp]]
        exp_name <- attr(df, "exp_cond") %||% "Experimental"
        ref_name <- attr(df, "ref_cond") %||% "Reference"
        
        # 1. Handle Contrast Inversion (Sign & Columns)
        if (flip) {
          df <- df %>%
            mutate(
              logFC = -logFC,
              temp_exp = ExpQuant,
              ExpQuant = RefQuant,
              RefQuant = temp_exp,
              temp_v = Exp_Valid,
              Exp_Valid = Ref_Valid,
              Ref_Valid = temp_v
            ) %>%
            select(-any_of(c("temp_exp", "temp_v")))
          
          # Swap condition names and title
          temp_n   <- exp_name
          exp_name <- ref_name
          ref_name <- temp_n
          sheet_title <- paste0(exp_name, "-", ref_name)
        } else {
          sheet_title <- comp
        }
        
        # Excel sheet name character restrictions
        safe_sheet_title <- substr(gsub("[\\[\\]\\*\\?\\/\\\\]", "_", sheet_title), 1, 31)
        
        # 2. Base R renaming (bypasses the rlang `:=` operator requirement)
        rename_map <- c(
          "ExpQuant"  = paste0("Log2_", exp_name, "_Mean"),
          "RefQuant"  = paste0("Log2_", ref_name, "_Mean"),
          "Exp_Valid" = paste0(exp_name, "_Valid_Count"),
          "Ref_Valid" = paste0(ref_name, "_Valid_Count")
        )
        
        df_export <- df
        for (orig_col in names(rename_map)) {
          if (orig_col %in% names(df_export)) {
            names(df_export)[names(df_export) == orig_col] <- rename_map[[orig_col]]
          }
        }
        
        openxlsx::addWorksheet(wb, safe_sheet_title)
        openxlsx::writeData(wb, safe_sheet_title, df_export)
      }
      
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
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
    
    summary_df <- rv$raw_data %>%
      group_by(condition) %>%
      summarise(
        Total_Samples = n_distinct(ID),
        Proteins_Identified_Any = n_distinct(Protein.ID[!is.na(Intensity) & Intensity > 0]),
        .groups = "drop"
      )
    
    cv_per_cond <- rv$raw_data %>%
      filter(!is.na(Intensity) & Intensity > 0) %>%
      group_by(condition, Protein.ID) %>%
      summarise(
        n_obs = n(),
        mean_int = mean(Intensity, na.rm = TRUE),
        sd_int = sd(Intensity, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_obs > 1) %>%
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
      flip_val          <- is_flipped()
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
      
      mapping_dput <- paste(capture.output(dput(rv$sample_map)), collapse = "\n")
      
      header_code <- glue::glue('
# ==============================================================================
# Automated Reproducible Proteomics Pipeline
# Generated from Proteomics Explorer Dashboard v{APP_VERSION}
# Timestamp: {Sys.time()}
# Execution Mode: {ifelse(mode_choice == "raw", "Raw Ingestion Pipeline", "State Ingestion Pipeline")}
# ==============================================================================

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "limma", "QFeatures", "MsCoreUtils", "stringr", "readr")
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
  se <- SummarizedExperiment(assays = setNames(list(assay_mat), name), rowData = row_data)
  QFeatures(setNames(list(se), name))
}}
')
      
      if (mode_choice == "raw") {
        ingestion_code <- glue::glue('
raw_file_path <- "{raw_file_name}"

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

parsed_res <- {ifelse(data_source_type == "Spectronaut", "parse_spectronaut(raw_file_path)", "parse_msfragger(raw_file_path)")}
raw_data <- parsed_res$data
geneDict <- parsed_res$geneDict
sample_map <- {mapping_dput}

raw_data <- raw_data %>%
  select(-any_of(c("condition", "BR"))) %>%
  left_join(sample_map, by = "ID") %>%
  mutate(condition = user_condition, BR = user_BR, ID = new_ID) %>%
  select(Protein.ID, ID, Intensity, condition, BR, LogInt)
')
      } else {
        ingestion_code <- glue::glue('
saved_state <- readRDS("{export_pfx}_state.rds")
raw_data   <- saved_state$raw_data
geneDict   <- saved_state$geneDict
sample_map <- saved_state$sample_map
')
      }
      
      pipeline_code <- glue::glue('
wide_data <- raw_data %>%
  select(Protein.ID, ID, LogInt) %>%
  pivot_wider(names_from = ID, values_from = LogInt)

qobj <- readQFeatures2(wide_data, ecol = 2:ncol(wide_data), fnames = "Protein.ID", name = "raw")

if ("{norm_method_val}" != "none") {{
  qobj <- addAssay(qobj, normalize(qobj[["raw"]], method = "{norm_method_val}"), name = "norm")
}} else {{
  qobj <- addAssay(qobj, qobj[["raw"]], name = "norm")
}}

norm_mat <- assay(qobj[["norm"]])

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
      for (p in mnar_proteins) {{
        row_vals <- sub_m[p, ]
        if (all(is.na(row_vals))) {{
          imputed_mat[p, cols] <- min(norm_mat, na.rm = TRUE) - 0.5
        }} else {{
          min_val <- min(row_vals, na.rm = TRUE)
          imputed_mat[p, cols] <- ifelse(is.na(row_vals), min_val - 0.5, row_vals)
        }}
      }}
    }}
  }}
  norm_mat <- imputed_mat
}} else if (!("{impute_method_val}" %in% c("None", "No Imputation (Show 1-Condition Dropouts on Margins)"))) {{
  norm_mat <- MsCoreUtils::impute_matrix(norm_mat, method = "{impute_method_val}")
}}

norm_data <- as.data.frame(norm_mat)

norm_export_df <- norm_data %>%
  tibble::rownames_to_column(var = "Protein.ID") %>%
  left_join(geneDict, by = c("Protein.ID" = "Accession"))

readr::write_tsv(norm_export_df, file = "{export_pfx}_normalized_imputed_matrix.tsv")

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
      ExpQuant  = mean(Quant[Condition == exp_cond], na.rm = TRUE),
      RefQuant  = mean(Quant[Condition == ref_cond], na.rm = TRUE),
      Exp_Valid = sum(!is.na(Quant[Condition == exp_cond])),
      Ref_Valid = sum(!is.na(Quant[Condition == ref_cond])),
      .groups   = "drop"
    )
  
  dt <- dt %>% left_join(scatter_means, by = "Accession")
  
  if ("{impute_method_val}" == "No Imputation (Show 1-Condition Dropouts on Margins)") {{
    dropout_rows <- scatter_means %>%
      filter((Exp_Valid > 0 & Ref_Valid == 0) | (Exp_Valid == 0 & Ref_Valid > 0)) %>%
      filter(!Accession %in% dt$Accession) %>%
      mutate(
        logFC     = ifelse(Exp_Valid > 0, Inf, -Inf),
        AveExpr   = ifelse(Exp_Valid > 0, ExpQuant, RefQuant),
        t         = NA_real_,
        P.Value   = NA_real_,
        adj.P.Val = NA_real_,
        B         = NA_real_,
        CI.L      = NA_real_,
        CI.R      = NA_real_
      )
    dt <- bind_rows(dt, dropout_rows)
  }}
  
  dt <- dt %>% left_join(geneDict, by = "Accession")
  attr(dt, "exp_cond") <- exp_cond
  attr(dt, "ref_cond") <- ref_cond
  de_results[[comp]] <- dt
}}

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
    df <- df %>%
      mutate(
        logFC = -logFC,
        temp_exp = ExpQuant,
        ExpQuant = RefQuant,
        RefQuant = temp_exp,
        temp_v = Exp_Valid,
        Exp_Valid = Ref_Valid,
        Ref_Valid = temp_v
      ) %>%
      select(-any_of(c("temp_exp", "temp_v")))
    
    tmp <- exp_name
    exp_name <- ref_name
    ref_name <- tmp
    curr_title <- paste0(exp_name, "-", ref_name)
  }}
  
  log10_p_line <- -log10(adj_p_cutoff)
  df <- df %>%
    mutate(
      significant = !is.na(adj.P.Val) & is.finite(logFC) & abs(logFC) > fc_cutoff & adj.P.Val < adj_p_cutoff,
      Regulation = case_when(
        is.infinite(logFC)       ~ "Dropout (1-Condition)",
        significant & logFC > 0  ~ "Upregulated",
        significant & logFC <= 0 ~ "Downregulated",
        TRUE                     ~ "Not Significant"
      ),
      Regulation = factor(Regulation, levels = c("Upregulated", "Downregulated", "Dropout (1-Condition)", "Not Significant"))
    )
  
  # === INSERT STEP 2 HERE ===
  rename_map <- c(
    "ExpQuant"  = paste0("Log2_", exp_name, "_Mean"),
    "RefQuant"  = paste0("Log2_", ref_name, "_Mean"),
    "Exp_Valid" = paste0(exp_name, "_Valid_Count"),
    "Ref_Valid" = paste0(ref_name, "_Valid_Count")
  )
  
  export_df <- df
  for (orig_col in names(rename_map)) {{
    if (orig_col %in% names(export_df)) {{
      names(export_df)[names(export_df) == orig_col] <- rename_map[[orig_col]]
    }}
  }}
  
  readr::write_tsv(export_df, file = paste0("{export_pfx}_", curr_title, "_de_results.tsv"))
  sig_df <- df %>% filter(significant & !is.na(gene) & gene != "")
  
  # Volcano
  p_volcano <- ggplot(df %>% filter(is.finite(logFC) & !is.na(adj.P.Val)), aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Regulation), size = {pt_size}) +
    geom_hline(yintercept = log10_p_line, linetype = 2, color = "grey50") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = 2, color = "grey50") +
    scale_color_manual(values = c("Upregulated" = "{col_up_val}", "Downregulated" = "{col_down_val}", "Dropout (1-Condition)" = "grey40", "Not Significant" = "grey75"), drop = FALSE) +
    theme_bw(base_size = {txt_size}) +
    theme(panel.grid = element_blank()) +
    labs(title = paste("Volcano Plot:", curr_title), x = paste0("Log2 Fold Change (", exp_name, " / ", ref_name, ")"), y = "-Log10 Adjusted p-value")
  
  if (show_labels && nrow(sig_df) > 0) {{
    p_volcano <- p_volcano + geom_label_repel(data = sig_df, aes(label = gene), size = {lbl_size}, max.overlaps = {max_overlaps_val}, show.legend = FALSE)
  }}
  ggsave(paste0("{export_pfx}_", comp_name, "_volcano.png"), plot = p_volcano, width = 8, height = 6, dpi = 300)
  
  # MA
  # MA
  finite_fc <- df$logFC[is.finite(df$logFC)]
  y_cap <- if (length(finite_fc) > 0) max(abs(finite_fc), na.rm = TRUE) + 2 else 6
  plot_ma_df <- df %>%
    mutate(
      AveExpr = case_when(
        !is.na(AveExpr) & is.finite(AveExpr) ~ AveExpr,
        !is.na(ExpQuant) & is.finite(ExpQuant) ~ ExpQuant,
        !is.na(RefQuant) & is.finite(RefQuant) ~ RefQuant,
        TRUE ~ 0
      ),
      plot_logFC = case_when(logFC == Inf ~ y_cap, logFC == -Inf ~ -y_cap, TRUE ~ logFC)
    ) %>%
    filter(!is.na(AveExpr) & !is.na(plot_logFC))
  
  p_ma <- ggplot(plot_ma_df, aes(x = AveExpr, y = plot_logFC)) +
    geom_point(aes(color = Regulation), size = {pt_size}) +
    geom_hline(yintercept = 0, color = "grey30") +
    geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = 2, color = "grey50") +
    geom_hline(yintercept = c(-y_cap, y_cap), linetype = 3, color = "grey70") +
    coord_cartesian(ylim = c(-y_cap - 0.5, y_cap + 0.5)) +
    scale_color_manual(values = c("Upregulated" = "{col_up_val}", "Downregulated" = "{col_down_val}", "Dropout (1-Condition)" = "grey40", "Not Significant" = "grey75"), drop = FALSE) +
    theme_bw(base_size = {txt_size}) +
    theme(panel.grid = element_blank()) +
    labs(title = paste("MA Plot:", curr_title), x = "Log2 Average Expression", y = paste0("Log2 Fold Change (", exp_name, " / ", ref_name, ")"))
  if (show_labels && nrow(sig_df) > 0) {{
    p_ma <- p_ma + geom_label_repel(data = sig_df, aes(x = AveExpr, y = logFC, label = gene), size = {lbl_size}, max.overlaps = {max_overlaps_val}, show.legend = FALSE)
  }}
  ggsave(paste0("{export_pfx}_", comp_name, "_ma.png"), plot = p_ma, width = 8, height = 6, dpi = 300)
  
  # Scatter
  all_quants <- c(df$RefQuant[is.finite(df$RefQuant)], df$ExpQuant[is.finite(df$ExpQuant)])
  axis_floor <- if (length(all_quants) > 0) min(all_quants, na.rm = TRUE) - 1 else 0
  plot_sc_df <- df %>%
    mutate(
      plot_Ref = ifelse(is.na(RefQuant) | !is.finite(RefQuant), axis_floor, RefQuant),
      plot_Exp = ifelse(is.na(ExpQuant) | !is.finite(ExpQuant), axis_floor, ExpQuant)
    )
  
  p_scatter <- ggplot(plot_sc_df, aes(x = plot_Ref, y = plot_Exp)) +
    geom_point(aes(color = Regulation), size = {pt_size}) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey50") +
    scale_color_manual(values = c("Upregulated" = "{col_up_val}", "Downregulated" = "{col_down_val}", "Dropout (1-Condition)" = "grey40", "Not Significant" = "grey75"), drop = FALSE) +
    theme_bw(base_size = {txt_size}) +
    theme(panel.grid = element_blank()) +
    labs(title = paste("Scatter Plot:", curr_title), x = paste("Log2", ref_name, "Average Abundance"), y = paste("Log2", exp_name, "Average Abundance"))
  
  if (show_labels && nrow(sig_df) > 0) {{
    p_scatter <- p_scatter + geom_label_repel(data = sig_df, aes(x = RefQuant, y = ExpQuant, label = gene), size = {lbl_size}, max.overlaps = {max_overlaps_val}, show.legend = FALSE)
  }}
  ggsave(paste0("{export_pfx}_", comp_name, "_scatter.png"), plot = p_scatter, width = 8, height = 6, dpi = 300)
}}

cat("Pipeline completed successfully!\n")
')
      
      writeLines(paste0(header_code, ingestion_code, pipeline_code), file)
    }
  )
}