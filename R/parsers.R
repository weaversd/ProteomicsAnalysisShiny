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
  
  # 1. Accession, Gene, Description Extraction
  acc_col <- intersect(c("Protein.ID", "Protein", "Protein ID", "Accession"), names(df))[1]
  if (is.na(acc_col)) {
    stop("Could not find a valid Protein ID column in the uploaded MSFragger file.")
  }
  
  gene_col <- intersect(c("Gene", "PG.Genes", "Gene Name", "Gene.Name", "Gene Symbol"), names(df))[1]
  desc_col <- intersect(c("Description", "PG.ProteinDescriptions", "Protein Description"), names(df))[1]
  
  geneDict <- df %>% 
    mutate(
      Accession   = as.character(.data[[acc_col]]),
      gene        = if (!is.na(gene_col) && gene_col %in% names(df)) as.character(.data[[gene_col]]) else as.character(.data[[acc_col]]),
      description = if (!is.na(desc_col) && desc_col %in% names(df)) as.character(.data[[desc_col]]) else as.character(.data[[acc_col]])
    ) %>% 
    select(Accession, gene, description) %>% 
    mutate(
      gene = ifelse(is.na(gene) | gene == "", Accession, gene),
      description = ifelse(is.na(description) | description == "", Accession, description)
    ) %>%
    distinct(Accession, .keep_all = TRUE)
  
  # 2. Strict Selection of Quantitative Columns (Select ONE type only)
  # Look for standard Intensity (excluding MaxLFQ)
  std_int_cols <- names(df)[which(grepl("[ .]Intensity$", names(df)) & !grepl("MaxLFQ", names(df)))]
  
  # Look for MaxLFQ Intensity
  maxlfq_cols <- names(df)[which(grepl("[ .]MaxLFQ[ .]Intensity$", names(df)))]
  
  # Default to Standard Intensity if found; otherwise fallback to MaxLFQ
  if (length(std_int_cols) > 0) {
    int_cols <- std_int_cols
  } else if (length(maxlfq_cols) > 0) {
    int_cols <- maxlfq_cols
  } else {
    stop("No valid intensity columns found in MSFragger file.")
  }
  
  # 3. Process into Long Format
  df_proc <- df %>%
    mutate(Protein.ID = as.character(.data[[acc_col]])) %>%
    select(Protein.ID, all_of(int_cols)) %>%
    pivot_longer(cols = -Protein.ID, names_to = "Raw_ID", values_to = "Intensity") %>%
    mutate(
      Intensity = suppressWarnings(as.numeric(Intensity)),
      Intensity = ifelse(Intensity == 0, NA_real_, Intensity),
      # Clean ID to just the sample name (e.g. "C_1", "N_1")
      ID = str_remove(Raw_ID, "[ .](MaxLFQ[ .])?Intensity$"),
      condition = str_remove_all(str_extract(ID, "x[0-9]+|[A-Za-z]+"), "x"),
      BR        = str_extract(ID, "\\d+$"),
      condition = ifelse(is.na(condition) | condition == "", "Sample", condition),
      BR        = ifelse(is.na(BR) | BR == "", "1", BR),
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

# Helper: Extract clean, concise sample names from full file paths
clean_sample_names <- function(raw_cols) {
  # 1. Extract base filename (handles both Windows \ and Unix / slashes)
  basenames <- gsub(".*[\\\\/]", "", raw_cols)
  
  # 2. Strip standard mass-spectrometry raw data extensions
  names_no_ext <- gsub("(?i)\\.(wiff|raw|mzml|d|dia|tsv|txt)$", "", basenames, perl = TRUE)
  
  if (length(names_no_ext) <= 1) return(names_no_ext)
  
  # 3. Detect longest common prefix
  s_min <- min(names_no_ext)
  s_max <- max(names_no_ext)
  chars_min <- strsplit(s_min, "")[[1]]
  chars_max <- strsplit(s_max, "")[[1]]
  len <- min(length(chars_min), length(chars_max))
  p_len <- 0
  while (p_len < len && chars_min[p_len + 1] == chars_max[p_len + 1]) {
    p_len <- p_len + 1
  }
  
  # 4. Detect longest common suffix
  rev_names <- vapply(names_no_ext, function(x) paste(rev(strsplit(x, "")[[1]]), collapse = ""), character(1))
  r_min <- min(rev_names)
  r_max <- max(rev_names)
  r_chars_min <- strsplit(r_min, "")[[1]]
  r_chars_max <- strsplit(r_max, "")[[1]]
  r_len <- min(length(r_chars_min), length(r_chars_max))
  s_len <- 0
  while (s_len < r_len && r_chars_min[s_len + 1] == r_chars_max[s_len + 1]) {
    s_len <- s_len + 1
  }
  
  # 5. Extract the variable core tokens
  cand <- vapply(names_no_ext, function(x) {
    total_len <- nchar(x)
    core_start <- p_len + 1
    core_end <- total_len - s_len
    if (core_start <= core_end) {
      substr(x, core_start, core_end)
    } else {
      x
    }
  }, character(1), USE.NAMES = FALSE)
  
  # Strip leftover leading/trailing delimiter artifacts
  cand_clean <- gsub("^[_.-]+|[_.-]+$", "", cand)
  
  # Ensure the simplified tokens are non-empty and unique
  if (all(nchar(cand_clean) > 0) && length(unique(cand_clean)) == length(names_no_ext)) {
    return(cand_clean)
  }
  
  return(names_no_ext)
}

# Parser for DIA-NN Matrix Output (report.pg_matrix.tsv)
parse_diann <- function(path) {
  df <- read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  # 1. Identify Protein Group / Accession column
  acc_col <- intersect(c("Protein.Group", "Protein.Ids", "Protein.Names", "Protein.ID", "Accession"), names(df))[1]
  if (is.na(acc_col)) {
    stop("Could not find a valid Protein Group / Accession column in the DIA-NN matrix file.")
  }
  
  # 2. Extract Gene and Description dictionaries
  gene_col <- intersect(c("Genes", "Gene", "Gene.Names", "Gene Name"), names(df))[1]
  desc_col <- intersect(c("First.Protein.Description", "Protein.Description", "Description"), names(df))[1]
  
  geneDict <- df %>% 
    mutate(
      Accession   = as.character(.data[[acc_col]]),
      gene        = if (!is.na(gene_col) && gene_col %in% names(df)) as.character(.data[[gene_col]]) else as.character(.data[[acc_col]]),
      description = if (!is.na(desc_col) && desc_col %in% names(df)) as.character(.data[[desc_col]]) else as.character(.data[[acc_col]])
    ) %>% 
    select(Accession, gene, description) %>% 
    mutate(
      gene = ifelse(is.na(gene) | gene == "", Accession, gene),
      description = ifelse(is.na(description) | description == "", Accession, description)
    ) %>%
    distinct(Accession, .keep_all = TRUE)
  
  # 3. Identify Quantitative Columns (exclude known DIA-NN metadata columns)
  meta_cols <- c(
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "Gene.Names",
    "First.Protein.Description", "Protein.Description", "Description",
    "N.Sequences", "N.Proteotypic.Sequences", "Global.Q.Value", "Global.PG.Q.Value",
    "PG.Q.Value", "Q.Value", "Precursor.Id", "Modified.Sequence", "Stripped.Sequence"
  )
  raw_sample_cols <- setdiff(names(df), meta_cols)
  
  if (length(raw_sample_cols) == 0) {
    stop("No quantitative sample columns found in DIA-NN matrix file.")
  }
  
  # 4. Generate clean sample IDs and map
  clean_ids <- clean_sample_names(raw_sample_cols)
  col_map   <- setNames(clean_ids, raw_sample_cols)
  
  # 5. Pivot long, convert linear intensity to Log2, and infer conditions
  df_proc <- df %>%
    mutate(Protein.ID = as.character(.data[[acc_col]])) %>%
    select(Protein.ID, all_of(raw_sample_cols)) %>%
    pivot_longer(cols = -Protein.ID, names_to = "Raw_Col", values_to = "Intensity") %>%
    mutate(
      Intensity = suppressWarnings(as.numeric(Intensity)),
      Intensity = ifelse(Intensity <= 0 | is.na(Intensity), NA_real_, Intensity),
      ID        = unname(col_map[Raw_Col]),
      condition = str_remove_all(str_extract(ID, "x[0-9]+|[A-Za-z]+"), "x"),
      BR        = str_extract(ID, "\\d+$"),
      condition = ifelse(is.na(condition) | condition == "", "Sample", condition),
      BR        = ifelse(is.na(BR) | BR == "", "1", BR),
      LogInt    = log2(Intensity)
    ) %>%
    select(Protein.ID, ID, Intensity, condition, BR, LogInt)
  
  list(data = df_proc, geneDict = geneDict)
}

# R/parsers.R (PEAKS Parser)

parse_peaks <- function(path) {
  # PEAKS outputs are typically comma-separated CSV files
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # 1. Identify Accession column
  acc_col <- intersect(c("Accession", "Protein Accession", "Protein.Accession", "Protein"), names(df))[1]
  if (is.na(acc_col)) {
    stop("Could not find a valid Accession column in the PEAKS file.")
  }
  
  # 2. Extract Gene from GN= in Description, and Clean Accession (strip |...)
  desc_col <- intersect(c("Description", "Protein Description", "Protein.Description"), names(df))[1]
  
  # Clean accession: discard |gene_ID part
  clean_accessions <- gsub("\\|.*$", "", as.character(df[[acc_col]]))
  
  # Extract gene symbol from GN= in description
  if (!is.na(desc_col) && desc_col %in% names(df)) {
    descriptions <- as.character(df[[desc_col]])
    extracted_genes <- stringr::str_match(descriptions, "\\bGN=([^\\s]+)")[, 2]
    # Fallback to cleaned accession if GN= is missing or empty
    gene_symbols <- ifelse(is.na(extracted_genes) | extracted_genes == "", clean_accessions, extracted_genes)
  } else {
    descriptions <- clean_accessions
    gene_symbols <- clean_accessions
  }
  
  geneDict <- data.frame(
    Accession   = clean_accessions,
    gene        = gene_symbols,
    description = descriptions,
    stringsAsFactors = FALSE
  ) %>%
    distinct(Accession, .keep_all = TRUE)
  
  # 3. Locate Quantitative Area Columns
  # Handles DB search ('Area Sample 1') and LFQ ('Sample 1 Area')
  # Excludes 'Group X Area' summary columns
  all_area_cols <- names(df)[grepl("(?i)\\barea\\b", names(df), perl = TRUE)]
  area_cols <- all_area_cols[!grepl("(?i)\\bgroup\\b|\\bprofile\\b|\\bratio\\b", all_area_cols, perl = TRUE)]
  
  if (length(area_cols) == 0) {
    stop("No quantitative sample Area columns found in PEAKS file.")
  }
  
  # 4. Standardize Sample Names (e.g., 'Area Sample 1' or 'Sample 1 Area' -> 'Sample 1')
  sample_clean_ids <- gsub("(?i)\\barea\\b", "", area_cols, perl = TRUE)
  sample_clean_ids <- trimws(gsub("\\s+", " ", sample_clean_ids))
  col_map <- setNames(sample_clean_ids, area_cols)
  
  # 5. Process into Long Format and Log2 Transform
  df_proc <- df %>%
    mutate(Protein.ID = clean_accessions) %>%
    select(Protein.ID, all_of(area_cols)) %>%
    pivot_longer(cols = -Protein.ID, names_to = "Raw_Col", values_to = "Area") %>%
    mutate(
      Intensity = suppressWarnings(as.numeric(Area)),
      Intensity = ifelse(Intensity <= 0 | is.na(Intensity), NA_real_, Intensity),
      ID        = unname(col_map[Raw_Col]),
      condition = str_remove_all(str_extract(ID, "x[0-9]+|[A-Za-z]+"), "x"),
      BR        = str_extract(ID, "\\d+$"),
      condition = ifelse(is.na(condition) | condition == "", "Sample", condition),
      BR        = ifelse(is.na(BR) | BR == "", "1", BR),
      LogInt    = log2(Intensity)
    ) %>%
    select(Protein.ID, ID, Intensity, condition, BR, LogInt)
  
  list(data = df_proc, geneDict = geneDict)
}
