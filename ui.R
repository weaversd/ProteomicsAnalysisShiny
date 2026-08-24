# ui.R
ui <- page_navbar(
  title = "Proteomics Explorer Dashboard",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "1. Data Import & Setup",
    sidebarLayout(
      sidebarPanel(
        selectInput("data_source", "Select Source Format:", 
                    choices = c("Spectronaut", "MSFragger", "Re-import Exported RDS/State")),
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
        selectInput("impute_method", "Imputation Method:",
                    choices = c("None", "Hybrid (MAR: KNN / MNAR: MinDet)", "knn", "MinDet", "MinProb", "zero")),
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
    "3. Interactive Differential Analysis",
    sidebarLayout(
      sidebarPanel(
        width = 3, # Keeps sidebar concise so main panel gets maximum width
        uiOutput("comparison_selector"),
        checkboxInput("flip_direction", "Invert Contrast Direction (Flip Sign)", value = FALSE),
        selectInput("plot_type", "Plot Type:", choices = c("Volcano", "MA", "Scatter")),
        hr(),
        h5("Plot Styling Controls"),
        sliderInput("point_size", "Point Size:", min = 0.5, max = 5, value = 2, step = 0.5),
        sliderInput("text_size", "Axis Label Size:", min = 8, max = 20, value = 12),
        hr(),
        hr(),
        h5("Gene Label Settings"),
        checkboxInput("show_labels", "Display Gene Labels", value = TRUE),
        
        conditionalPanel(
          condition = "input.show_labels == true",
          sliderInput("label_size", "Label Text Size:", min = 1, max = 8, value = 3),
          sliderInput("max_overlaps", "Max Allowed Overlaps:", min = 1, max = 100, value = 15, step = 1)
        ),
        colourpicker::colourInput("col_up", "Upregulated Color:", value = "#66CCFE"),
        colourpicker::colourInput("col_down", "Downregulated Color:", value = "#FF0066"),
        sliderInput(
          inputId = "adj_p_cutoff", 
          label   = "Adjusted p-value cutoff (FDR):", 
          min     = 0.0001, 
          max     = 0.20, 
          value   = 0.05, 
          step    = 0.005
        ),
        sliderInput("fc_cutoff", "|log2 FC| cutoff:", min = 0, max = 4, value = 1, step = 0.25),
        hr(),
        h5("Export Current Plot"),
        fluidRow(
          column(6, numericInput("plot_width", "Width (in):", value = 8, min = 2, max = 20)),
          column(6, numericInput("plot_height", "Height (in):", value = 6, min = 2, max = 20))
        ),
        downloadButton("download_plot_png", "Download PNG", class = "btn-outline-primary w-100 mb-2"),
        downloadButton("download_plot_pdf", "Download PDF", class = "btn-outline-danger w-100")
      ),
      
      # Main panel taking up full width of remaining window space
      mainPanel(
        width = 9,
        accordion(
          id = "plots_accordion",
          multiple = TRUE, # Allows BOTH panels to be open simultaneously
          
          # Panel 1: Interactive Plotly (On top)
          accordion_panel(
            title = "Interactive Plotly View",
            value = "plotly_panel",
            plotlyOutput("plotly_view", height = "600px", width = "100%")
          ),
          
          # Panel 2: Publication ggplot (Below)
          accordion_panel(
            title = "Publication-Ready ggplot View",
            value = "ggplot_panel",
            plotOutput("ggplot_view", height = "600px", width = "100%")
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