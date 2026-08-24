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