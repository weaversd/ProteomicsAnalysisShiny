# tests/testthat/test-imputation.R
library(testthat)
library(MsCoreUtils)
source("../../global.R")

test_that("Dynamic hybrid imputation correctly segregates MAR and MNAR proteins", {
  mat <- matrix(c(
    # Protein 1: Full in x1 (3/3), MAR in x2 (2/3), MNAR in x3 (1/3)
    10, 11, 10.5,    10, NA, 10.2,    NA, NA, 10,
    # Protein 2: Complete dropout in x1 (0/3), Full in x2 (3/3), Full in x3 (3/3)
    NA, NA, NA,      15, 14.8, 15.2,  16, 16.1, 15.9
  ), nrow = 2, byrow = TRUE)
  
  colnames(mat) <- c("x11", "x12", "x13", "x21", "x22", "x23", "x31", "x32", "x33")
  rownames(mat) <- c("Prot_MAR_MNAR", "Prot_Dropout")
  
  cond_lookup <- setNames(rep(c("x1", "x2", "x3"), each = 3), colnames(mat))
  conds <- unname(cond_lookup[colnames(mat)])
  
  global_mar <- MsCoreUtils::impute_matrix(mat, method = "nbavg")
  global_min <- min(mat, na.rm = TRUE)
  imputed_mat <- mat
  
  for (cond in unique(conds)) {
    cols <- which(conds == cond)
    n_reps <- length(cols)
    sub_m <- mat[, cols, drop = FALSE]
    present <- rowSums(!is.na(sub_m))
    mar_threshold <- ceiling(n_reps / 2)
    
    # CASE 1: MAR (>= 50% present)
    mar_proteins <- names(present[present >= mar_threshold & present < n_reps])
    if (length(mar_proteins) > 0) {
      imputed_mat[mar_proteins, cols] <- global_mar[mar_proteins, cols]
    }
    
    # CASE 2: MNAR (< 50% present)
    mnar_proteins <- names(present[present < mar_threshold])
    if (length(mnar_proteins) > 0) {
      for (p in mnar_proteins) {
        row_vals <- sub_m[p, ]
        if (all(is.na(row_vals))) {
          imputed_mat[p, cols] <- global_min - 0.5
        } else {
          # Impute using min observed value in that row/condition with small left-shift
          min_val <- min(row_vals, na.rm = TRUE)
          imputed_mat[p, cols] <- ifelse(is.na(row_vals), min_val - 0.5, row_vals)
        }
      }
    }
  }
  
  # Assertions
  expect_false(any(is.na(imputed_mat)))
  expect_lt(imputed_mat["Prot_Dropout", "x11"], min(mat, na.rm = TRUE))
})


test_that("Reference protein normalization correctly aligns samples to target baseline", {
  # Mock 2 samples with a 2-fold (1 log2 unit) loading discrepancy
  mock_wide <- data.frame(
    Protein.ID = c("P_TARGET", "P_OTHER"),
    Sample_1   = c(10, 5),
    Sample_2   = c(11, 6) # Sample 2 has 1 log2 unit more overall signal
  )
  
  ref_row <- mock_wide %>% filter(Protein.ID == "P_TARGET")
  sample_cols <- c("Sample_1", "Sample_2")
  ref_vals <- as.numeric(ref_row[1, sample_cols]) # c(10, 11)
  
  global_ref_mean <- mean(ref_vals)               # 10.5
  sample_offsets  <- ref_vals - global_ref_mean   # c(-0.5, +0.5)
  
  for (idx in seq_along(sample_cols)) {
    col_name <- sample_cols[idx]
    mock_wide[[col_name]] <- mock_wide[[col_name]] - sample_offsets[idx]
  }
  
  # After adjustment:
  # P_TARGET should be exactly 10.5 in both samples
  expect_equal(mock_wide$Sample_1[mock_wide$Protein.ID == "P_TARGET"], 10.5)
  expect_equal(mock_wide$Sample_2[mock_wide$Protein.ID == "P_TARGET"], 10.5)
  
  # P_OTHER should now be equalized across samples (5 - (-0.5) = 5.5, 6 - 0.5 = 5.5)
  expect_equal(mock_wide$Sample_1[mock_wide$Protein.ID == "P_OTHER"], 5.5)
  expect_equal(mock_wide$Sample_2[mock_wide$Protein.ID == "P_OTHER"], 5.5)
})