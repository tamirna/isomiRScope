#!/usr/bin/env Rscript
#
# Regression test: classify the bundled example mrd and compare the result
# against the stored expected output.
#
# Usage: Rscript tests/test_example.R

suppressPackageStartupMessages(library(tidyverse))

script_dir <- function() {
  call_args <- base::commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", call_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  getwd()
}

repo_root <- dirname(script_dir())
example_mrd <- file.path(repo_root, "inst", "extdata", "example.mrd")
expected_tsv <- file.path(repo_root, "inst", "extdata", "example.expected.tsv")

for (f in c(example_mrd, expected_tsv)) {
  if (!file.exists(f)) {
    stop("missing test fixture: ", f, call. = FALSE)
  }
}

observed_tsv <- tempfile(fileext = ".tsv")
on.exit(unlink(observed_tsv), add = TRUE)

status <- system2("Rscript", c(
  file.path(repo_root, "bin", "quantify_isomirs.R"),
  paste0("--readsFile=", example_mrd),
  paste0("--outFile=", observed_tsv),
  "--cores=1"
))

if (status != 0) {
  stop("quantify_isomirs.R exited with status ", status, call. = FALSE)
}

observed <- read_delim(observed_tsv, delim = "\t", show_col_types = FALSE)
expected <- read_delim(expected_tsv, delim = "\t", show_col_types = FALSE)

if (!identical(dim(observed), dim(expected))) {
  stop("dimension mismatch: got ", nrow(observed), "x", ncol(observed),
       ", expected ", nrow(expected), "x", ncol(expected), call. = FALSE)
}

if (!isTRUE(all.equal(observed, expected))) {
  print(all.equal(observed, expected))
  stop("output does not match ", basename(expected_tsv), call. = FALSE)
}

message("OK: ", nrow(observed), " classified reads match ", basename(expected_tsv))
