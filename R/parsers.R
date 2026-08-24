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


parse_msfragger <- function(path) {
  df <- read.delim(path, sep = "\t", check.names = FALSE)
  
  # ----------------------------------------------------------------------------
  # 1. Accession Column
  # ----------------------------------------------------------------------------
  acc_col <- intersect(c("Protein.ID", "Protein", "Protein ID", "Accession"), names(df))[1]
  if (is.na(acc_col)) {
    stop("Could not find a valid Protein ID column in the uploaded MSFragger file.")
  }
  
  # ----------------------------------------------------------------------------
  # 2. Gene & Description Dictionary
  # ----------------------------------------------------------------------------
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
  
  # ----------------------------------------------------------------------------
  # 3. Target Standard Intensity Columns ONLY
  # ----------------------------------------------------------------------------
  # Regex explanation: Matches columns ending in '.Intensity' while excluding '.MaxLFQ.Intensity'
  int_cols <- names(df)[which(grepl(" Intensity$", names(df)) & !grepl("MaxLFQ", names(df)))]
  print(names(df))
  
  # Fallback: If standard Intensity isn't present, check for MaxLFQ
  if (length(int_cols) == 0) {
    int_cols <- names(df)[which(grepl(" MaxLFQ Intensity$", names(df)))]
  }
  
  if (length(int_cols) == 0) {
    stop("No valid intensity columns found in MSFragger file.")
  }
  
  # ----------------------------------------------------------------------------
  # 4. Reshape & Clean Metadata
  # ----------------------------------------------------------------------------
  df_proc <- df %>%
    mutate(Protein.ID = .data[[acc_col]]) %>%
    select(Protein.ID, all_of(int_cols)) %>%
    pivot_longer(cols = -Protein.ID, names_to = "ID", values_to = "Intensity") %>%
    mutate(
      Intensity = ifelse(Intensity == 0, NA, Intensity),
      # Clean '.Intensity' or '.MaxLFQ.Intensity' off sample IDs
      ID        = str_remove(ID, "\\.(MaxLFQ\\.)?Intensity$"),
      condition = str_remove_all(str_extract(ID, "x[0-9]+|[A-Za-z]+"), "x"),
      BR        = str_extract(ID, "\\d+$"),
      LogInt    = log2(Intensity)
    ) %>%
    select(Protein.ID, ID, Intensity, condition, BR, LogInt)
  
  list(data = df_proc, geneDict = geneDict)
}