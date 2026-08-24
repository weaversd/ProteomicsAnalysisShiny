parse_spectronaut <- function(path) {
  df <- read.delim(path, sep = "\t", check.names = FALSE)
  
  # ----------------------------------------------------------------------------
  # Robust Accession Column Selection
  # ----------------------------------------------------------------------------
  acc_col <- intersect(c("PG.ProteinAccessions", "Protein.ID", "ProteinAccessions"), names(df))[1]
  if (is.na(acc_col)) {
    stop("Could not find a valid Protein Accession column in the uploaded Spectronaut file.")
  }
  
  # ----------------------------------------------------------------------------
  # Robust Gene & Description Column Extraction
  # ----------------------------------------------------------------------------
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
      # Replace empty strings or NAs with Accession
      gene = ifelse(is.na(gene) | gene == "", Accession, gene),
      description = ifelse(is.na(description) | description == "", Accession, description)
    ) %>%
    distinct(Accession, .keep_all = TRUE)
  
  # ----------------------------------------------------------------------------
  # Process Intensity Data
  # ----------------------------------------------------------------------------
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
}


# R/parsers.R (Snippet for parse_msfragger)

parse_msfragger <- function(path) {
  df <- read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  acc_col <- intersect(c("Protein.ID", "Protein", "Protein ID", "Accession"), names(df))[1]
  if (is.na(acc_col)) {
    stop("Could not find a valid Protein ID column in the uploaded MSFragger file.")
  }
  
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
  
  # Regex accommodates either dot or space separators and explicitly excludes MaxLFQ
  int_cols <- names(df)[which(grepl("[\\.\\s]Intensity$", names(df), ignore.case = TRUE) & !grepl("MaxLFQ", names(df), ignore.case = TRUE))]
  
  if (length(int_cols) == 0) {
    int_cols <- names(df)[which(grepl("[\\.\\s]MaxLFQ[\\.\\s]Intensity$", names(df), ignore.case = TRUE))]
  }
  
  if (length(int_cols) == 0) {
    stop("No valid intensity columns found in MSFragger file.")
  }
  
  df_proc <- df %>%
    mutate(Protein.ID = .data[[acc_col]]) %>%
    select(Protein.ID, all_of(int_cols)) %>%
    pivot_longer(cols = -Protein.ID, names_to = "ID", values_to = "Intensity") %>%
    mutate(
      Intensity = as.numeric(Intensity),
      Intensity = ifelse(Intensity == 0, NA, Intensity),
      ID        = str_remove(ID, "[\\.\\s](MaxLFQ[\\.\\s])?Intensity$"),
      condition = str_remove_all(str_extract(ID, "x[0-9]+|[A-Za-z]+"), "x"),
      BR        = str_extract(ID, "\\d+$"),
      LogInt    = log2(Intensity)
    ) %>%
    select(Protein.ID, ID, Intensity, condition, BR, LogInt)
  
  list(data = df_proc, geneDict = geneDict)
}

parse_generic <- function(path) {
  # Read either CSV or TSV
  delim <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  df <- read.delim(path, sep = delim, check.names = FALSE, stringsAsFactors = FALSE)
  
  # 1. Identify Accession / Protein ID column
  acc_col <- intersect(c("Accession", "Protein.ID", "Protein_ID", "Protein", "ProteinID"), names(df))[1]
  if (is.na(acc_col)) {
    stop("Generic upload requires a protein ID column named 'Accession', 'Protein.ID', or 'Protein'.")
  }
  
  # 2. Extract Gene and Description dictionaries
  gene_col <- intersect(c("Gene", "Gene.Name", "Gene_Name", "Symbol", "GeneName"), names(df))[1]
  desc_col <- intersect(c("Description", "ProteinDescription", "Protein.Description"), names(df))[1]
  
  geneDict <- df %>%
    mutate(
      Accession   = as.character(.data[[acc_col]]),
      gene        = if (!is.na(gene_col) && gene_col %in% names(df)) as.character(.data[[gene_col]]) else as.character(.data[[acc_col]]),
      description = if (!is.na(desc_col) && desc_col %in% names(df)) as.character(.data[[desc_col]]) else as.character(.data[[acc_col]])
    ) %>%
    mutate(
      gene = ifelse(is.na(gene) | gene == "", Accession, gene),
      description = ifelse(is.na(description) | description == "", Accession, description)
    ) %>%
    select(Accession, gene, description) %>%
    distinct(Accession, .keep_all = TRUE)
  
  # 3. Detect Long vs. Wide format
  sample_col <- intersect(c("Sample", "Sample.ID", "Sample_ID", "ID", "Run", "File"), names(df))[1]
  intensity_col <- intersect(c("Intensity", "Quantity", "Abundance", "LogInt", "Value"), names(df))[1]
  
  is_long <- !is.na(sample_col) && !is.na(intensity_col)
  
  if (is_long) {
    # --------------------------------------------------------------------------
    # LONG FORMAT PARSING
    # --------------------------------------------------------------------------
    cond_col <- intersect(c("Condition", "Group", "Treatment"), names(df))[1]
    br_col   <- intersect(c("Replicate", "BR", "BioRep", "TechRep"), names(df))[1]
    
    df_proc <- df %>%
      mutate(
        Protein.ID = as.character(.data[[acc_col]]),
        ID         = as.character(.data[[sample_col]]),
        Intensity  = as.numeric(.data[[intensity_col]]),
        Intensity  = ifelse(Intensity == 0, NA, Intensity),
        condition  = if (!is.na(cond_col) && cond_col %in% names(df)) as.character(.data[[cond_col]]) else str_remove_all(str_extract(ID, "x[0-9]+|[A-Za-z]+"), "x"),
        BR         = if (!is.na(br_col) && br_col %in% names(df)) as.character(.data[[br_col]]) else str_extract(ID, "\\d+$"),
        LogInt     = log2(Intensity)
      ) %>%
      mutate(
        condition = ifelse(is.na(condition) | condition == "", "Cond1", condition),
        BR        = ifelse(is.na(BR) | BR == "", "1", BR)
      ) %>%
      select(Protein.ID, ID, Intensity, condition, BR, LogInt)
    
  } else {
    # --------------------------------------------------------------------------
    # WIDE FORMAT PARSING
    # --------------------------------------------------------------------------
    # Exclude metadata columns to isolate quantitative sample columns
    meta_cols <- c(acc_col, gene_col, desc_col, "Organism", "Length", "Coverage")
    candidate_cols <- setdiff(names(df), meta_cols)
    
    # Target intensity-labeled columns if present, otherwise take numeric columns
    int_pattern_cols <- candidate_cols[grepl("\\.Intensity$|^Intensity_|_Intensity$", candidate_cols)]
    
    if (length(int_pattern_cols) > 0) {
      quant_cols <- int_pattern_cols
    } else {
      # Fallback: select numeric columns
      quant_cols <- candidate_cols[sapply(df[, candidate_cols, drop = FALSE], is.numeric)]
    }
    
    if (length(quant_cols) == 0) {
      stop("Could not identify quantitative sample columns. Use '.Intensity' suffix or numeric columns.")
    }
    
    df_proc <- df %>%
      mutate(Protein.ID = as.character(.data[[acc_col]])) %>%
      select(Protein.ID, all_of(quant_cols)) %>%
      pivot_longer(cols = -Protein.ID, names_to = "ID", values_to = "Intensity") %>%
      mutate(
        Intensity = as.numeric(Intensity),
        Intensity = ifelse(Intensity == 0, NA, Intensity),
        ID        = str_remove_all(ID, "\\.Intensity$|^Intensity_|_Intensity$"),
        condition = str_remove_all(str_extract(ID, "x[0-9]+|[A-Za-z]+"), "x"),
        BR        = str_extract(ID, "\\d+$"),
        LogInt    = log2(Intensity)
      ) %>%
      mutate(
        condition = ifelse(is.na(condition) | condition == "", "Cond1", condition),
        BR        = ifelse(is.na(BR) | BR == "", "1", BR)
      ) %>%
      select(Protein.ID, ID, Intensity, condition, BR, LogInt)
  }
  
  list(data = df_proc, geneDict = geneDict)
}
