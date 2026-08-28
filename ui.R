# ui.R
ui <- page_navbar(
  title = "Proteomics Explorer Dashboard",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "1. Data Import & Setup",
    sidebarLayout(
      sidebarPanel(
        # ui.R (Snippet inside Nav Panel 1)
        selectInput("data_source", "Select Source Format:", 
                    choices = c("Spectronaut", "MSFragger", "Generic Table (Long or Wide)", "Re-import Exported RDS/State")),
        fileInput("file_upload", "Upload Protein File", accept = c(".tsv", ".csv", ".rds", ".txt")),
        hr(),
        h5("Sample Metadata Mapping"),
        helpText("Assign the Condition and Replicate number for each sample column:"),
        uiOutput("sample_mapping_ui"), # Dynamic mapping interface
        hr(),
        actionButton("btn_process_import", "Process & Map Samples", class = "btn-primary w-100")
      ),
      mainPanel(
        h4("Raw Data Preview"),
        DTOutput("import_preview_table")
      )
    )
  ),
  
  nav_panel(
    "2. Normalization & Imputation",
    sidebarLayout(
      sidebarPanel(
        selectInput("norm_method", "Normalization Method:",
                    choices = c("none", "center.median", "center.mean", "quantiles", "quantiles.robust")),
        selectInput(
          inputId = "impute_method",
          label   = "Select Imputation Method:",
          choices = c(
            "Hybrid (MAR: KNN / MNAR: MinDet)",
            "No Imputation (Show 1-Condition Dropouts on Margins)",
            "MinDet",
            "knn",
            "nbavg",
            "None"
          ), selected = "None"),
        actionButton("btn_apply_norm", "Apply Transformation", class = "btn-success w-100")
      ),
      mainPanel(
        plotOutput("norm_qc_plot", height = "450px"),
        hr(),
        h5("Imputation Statistics"),
        verbatimTextOutput("imputation_stats_text"),
        hr(),
        h5("Transformed Data Preview"),
        DTOutput("transformed_data_table")
      )
    )
  ),
  
  nav_panel(
    title = "3. Differential Expression & Visualization",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        
        # Comparison & Contrast Selection
        tags$h4("Contrast Settings"),
        uiOutput("comparison_selector"),
        checkboxInput("flip_direction", "Invert Contrast Direction", value = FALSE),
        
        hr(style = "margin: 10px 0;"),
        
        # Plot Type & Cutoff Thresholds
        tags$h4("Cutoffs & Plot Style"),
        selectInput(
          inputId = "plot_type",
          label   = "Select Plot Layout:",
          choices = c("Volcano", "MA", "Scatter"),
          selected = "Volcano"
        ),
        sliderInput(
          inputId = "adj_p_cutoff", 
          label   = "Adjusted p-value Cutoff (FDR):", 
          min     = 0.0001, 
          max     = 0.20, 
          value   = 0.05, 
          step    = 0.005
        ),
        numericInput(
          inputId = "fc_cutoff",
          label   = "Log2 Fold Change Cutoff:",
          value   = 1.0,
          min     = 0,
          step    = 0.1
        ),
        
        hr(style = "margin: 10px 0;"),
        
        # Standard DE Styling
        tags$h4("Standard DE Highlighting & Labels"),
        checkboxInput("show_de_colors", "Color Differentially Expressed Proteins", value = TRUE),
        checkboxInput("show_de_labels", "Label Differentially Expressed Proteins", value = FALSE),
        
        fluidRow(
          column(6, colourpicker::colourInput("col_up", "Upregulated:", value = "#66CCFE")),
          column(6, colourpicker::colourInput("col_down", "Downregulated:", value = "#FF0066"))
        ),
        
        hr(style = "margin: 10px 0;"),
        
        # Custom Protein Set Highlighting & Dynamic Inputs
        tags$h4("Custom Protein Sets / Subsets"),
        uiOutput("custom_protein_sets_ui"),
        div(
          style = "margin-top: 8px; margin-bottom: 15px;",
          actionButton("btn_add_custom_set", "+ Add Protein Set", class = "btn-sm btn-primary"),
          actionButton("btn_remove_custom_set", "- Remove Last Set", class = "btn-sm btn-outline-danger")
        ),
        
        hr(style = "margin: 10px 0;"),
        
        # Plot Aesthetics & Export Sizing
        tags$h4("Display & Label Tuning"),
        fluidRow(
          column(6, numericInput("point_size", "Point Size:", value = 2, min = 0.5, max = 10, step = 0.5)),
          column(6, numericInput("text_size", "Axis Text Size:", value = 12, min = 6, max = 24, step = 1))
        ),
        fluidRow(
          column(6, numericInput("label_size", "Label Text Size:", value = 3.5, min = 1, max = 10, step = 0.5)),
          column(6, numericInput("max_overlaps", "Max Label Overlaps:", value = 10, min = 0, max = 100, step = 1))
        ),
        fluidRow(
          column(6, numericInput("plot_width", "Export Width (in):", value = 8, min = 3, max = 20, step = 0.5)),
          column(6, numericInput("plot_height", "Export Height (in):", value = 6, min = 3, max = 20, step = 0.5))
        ),
        
        div(
          style = "margin-top: 10px;",
          downloadButton("download_plot_png", "Download PNG", class = "btn-sm btn-secondary"),
          downloadButton("download_plot_pdf", "Download PDF", class = "btn-sm btn-secondary")
        )
      ),
      
      mainPanel(
        width = 8,
        accordion(
          open = c("Interactive Plotly Explorer", "Publication ggplot"),
          accordion_panel(
            title = "Interactive Plotly Explorer",
            plotlyOutput("plotly_view", height = "550px")
          ),
          accordion_panel(
            title = "Publication ggplot",
            plotOutput("ggplot_view", height = "550px")
          )
        )
      )
    )
  ),
  
  nav_panel(
    "4. Export & Audit Trail",
    fluidRow(
      column(12,
             wellPanel(
               h4("Experiment Summary Statistics"),
               DTOutput("experiment_summary_table")
             ),
             hr()
      )
    ),
    fluidRow(
      column(6,
             wellPanel(
               h4("Export Analysis Bundle"),
               textInput("export_prefix", "Project Name Prefix:", value = "Proteomics_Project"),
               radioButtons(
                 "script_source_type", 
                 "Script Ingestion Mode:",
                 choices = c("Raw Input File (End-to-End)" = "raw", "Saved State (.rds)" = "rds"),
                 selected = "raw"
               ),
               downloadButton("download_excel", "Download Data (Excel)", class = "btn-primary w-100 mb-2"),
               downloadButton("download_script", "Download Reproducible R Script", class = "btn-outline-success w-100 mb-2"),
               downloadButton("download_audit", "Download Audit Trail (JSON)", class = "btn-secondary w-100 mb-2"),
               downloadButton("download_state", "Export Full State (RDS for Re-upload)", class = "btn-info w-100")
             )
      ),
      column(6,
             wellPanel(
               h4("Audit Log Preview"),
               verbatimTextOutput("audit_preview")
             )
      )
    )
  )
)