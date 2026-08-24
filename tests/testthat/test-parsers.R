# tests/testthat/test-parsers.R
library(testthat)
source("../../global.R")
source("../../R/parsers.R")

test_that("parse_generic handles wide-format tables correctly", {
  # 1. Mock wide format data
  wide_mock <- data.frame(
    Protein.ID = c("P001", "P002"),
    Gene = c("GENEA", "GENEB"),
    Description = c("Protein A", "Protein B"),
    x11.Intensity = c(1000, 2000),
    x12.Intensity = c(1500, 2500),
    x21.Intensity = c(3000, 4000),
    x22.Intensity = c(3500, 4500)
  )
  
  tmp_file <- tempfile(fileext = ".tsv")
  write.table(wide_mock, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  res <- parse_generic(tmp_file)
  
  expect_type(res, "list")
  expect_named(res, c("data", "geneDict"))
  expect_equal(nrow(res$geneDict), 2)
  expect_equal(nrow(res$data), 8) # 2 proteins * 4 samples
  expect_true(all(c("Protein.ID", "ID", "Intensity", "condition", "BR", "LogInt") %in% names(res$data)))
  expect_equal(res$data$LogInt[res$data$Protein.ID == "P001" & res$data$ID == "x11"], log2(1000))
  
  unlink(tmp_file)
})

test_that("parse_generic handles long-format tables and missing gene fallbacks", {
  # 2. Mock long format data without Gene/Description columns
  long_mock <- data.frame(
    Accession = rep(c("PROT_A", "PROT_B"), each = 3),
    Sample = rep(c("Ctrl_1", "Ctrl_2", "Trt_1"), 2),
    Intensity = c(100, 200, 300, 400, 500, 0) # Contains a 0 intensity
  )
  
  tmp_file <- tempfile(fileext = ".csv")
  write.csv(long_mock, tmp_file, row.names = FALSE)
  
  res <- parse_generic(tmp_file)
  
  # Gene names should fall back to Accession IDs
  expect_equal(res$geneDict$gene, c("PROT_A", "PROT_B"))
  
  # Zero intensities must be converted to NA
  na_entry <- res$data[res$data$Protein.ID == "PROT_B" & res$data$ID == "Trt_1", ]
  expect_true(is.na(na_entry$Intensity))
  expect_true(is.na(na_entry$LogInt))
  
  unlink(tmp_file)
})

# tests/testthat/test-parsers.R (Snippet around line 71)

test_that("parse_msfragger excludes MaxLFQ columns when raw Intensity is present", {
  fragger_mock <- data.frame(
    Protein.ID = c("P1", "P2"),
    Gene = c("G1", "G2"),
    Description = c("D1", "D2"),
    x11.Intensity = c(10, 20),
    x12.Intensity = c(15, 25),
    x11.MaxLFQ.Intensity = c(100, 200),
    x12.MaxLFQ.Intensity = c(150, 250)
  )
  
  tmp_file <- tempfile(fileext = ".tsv")
  # Explicitly pass sep = "\t"
  write.table(fragger_mock, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  res <- parse_msfragger(tmp_file)
  
  expect_equal(unique(res$data$ID), c("x11", "x12"))
  expect_equal(nrow(res$data), 4)
  
  unlink(tmp_file)
})