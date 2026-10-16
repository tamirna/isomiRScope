#!/usr/bin/env Rscript
#
# isomiRScope - classification of isomiRs from miRDeep2 mrd files
#
# Copyright 2026, TAmiRNA GmbH
# Released under the GNU General Public License v3. See LICENSE and COPYRIGHT.
#
# Reads one miRBase.mrd file as produced by the miRDeep2 quantifier.pl and
# writes a tab-separated table with one row per collapsed read, annotated with
# the canonical arm it belongs to, its isomiR type and whether it carries a
# polymorphism. See docs/algorithm.md for the classification rules.
#
# Usage:
#   Rscript bin/quantify_isomirs.R --readsFile=<miRBase.mrd> --outFile=<out.tsv> [--cores=4]
#
# Extracted from the miND small RNA sequencing pipeline.

suppressPackageStartupMessages({
  library(tidyverse)
  library(doParallel)
  library(slider)
  library(R.utils)
})

usage <- function() {
  cat(
    "isomiRScope - classify isomiRs from a miRDeep2 mrd file\n\n",
    "Usage:\n",
    "  Rscript bin/quantify_isomirs.R --readsFile=<miRBase.mrd> --outFile=<out.tsv> [--cores=N]\n\n",
    "Arguments:\n",
    "  --readsFile  miRBase.mrd written by miRDeep2 quantifier.pl\n",
    "  --outFile    output path for the tab-separated isomiR table\n",
    "  --cores      number of cores used for classification (default: 4)\n",
    sep = ""
  )
}

# Locate the repository root relative to this script so that the functions are
# found regardless of the working directory the script is called from.
script_dir <- function() {
  call_args <- base::commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", call_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  getwd()
}

source(file.path(dirname(script_dir()), "R", "isomir_functions.R"))

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE, defaults = list(
  cores = 4
))

if (!is.null(args$help) || is.null(args$readsFile) || is.null(args$outFile)) {
  usage()
  quit(status = if (is.null(args$help)) 1 else 0)
}

readsFile <- args$readsFile
outFile <- args$outFile
cores <- as.integer(args$cores)

if (!file.exists(readsFile)) {
  stop("readsFile does not exist: ", readsFile, call. = FALSE)
}
if (is.na(cores) || cores < 1) {
  stop("cores must be a positive integer, got: ", args$cores, call. = FALSE)
}

dir.create(dirname(outFile), showWarnings = FALSE, recursive = TRUE)

# read in mrd file line by line
# skip empty lines
data_raw <- tibble(value = read_lines(readsFile, skip_empty_rows = TRUE))

# split data based on > at the row beginning
split_data <- lapply(split(data_raw, cumsum(grepl('^>', data_raw$value))), function(x) x[-1, ])

# specify list names
names(split_data) <- data_raw$value[grepl('^>', data_raw$value)]

# use list mrd miRNA list to prepare a named mrd tibble list
# imap is used to extract .x tibble and .y list name
# lists containing NULL are discarded
tb_list_out <- split_data %>%
  imap(., prepare_tb_list) %>%
  compact()

if (length(tb_list_out) == 0) {
  stop("no precursors with mapped reads found in ", readsFile, call. = FALSE)
}

# use the named mrd tibble list to prepare annotation tibble
# bind rows collapses the tibble list to one tibble
# per row annotation information for 1 miRNA
tb_annotated <- tb_list_out %>%
  map(., prepare_annotation_tb)

tb_annotated <- tb_annotated %>% bind_rows()

# remove annotation information from named mrd tibble list
# only keep miRbase id and read information
tb_list_out_2 <- map(tb_list_out, ~ filter(.x, grepl("^seq", id) | id == "mirbase_id"))

# add the postiion information to each mrd tibble
tb_list_position <- map(tb_list_out_2, prepare_position_tb)

# register parallel processing
registerDoParallel(cores = cores)

# parallel looping
# over rows of annotation_tb
# classify isomirs
tb_list_isomir <- foreach(mir_annotation = slide(tb_annotated, ~.x),
                          .combine = "append") %dopar% annotate_isomirs(mir_annotation,
                                                                        tb_list_position)

# add the mirbase id in an extra mir_name column and remove respective row
tb_isomir <- map(tb_list_isomir, ~ mutate(.x, mir_name = .x %>%
                                            filter(id == "mirbase_id") %>%
                                            pull(value)) %>%
                   filter(id != "mirbase_id")) %>%
  bind_rows()

# remove > from mir name
# add the information of the canonical form to the mir name
# if we only have one mature version remove the -M information from the miR name end
# unite the information of polymorphic and 3' or 5' info
tb_isomir <- tb_isomir %>%
  mutate(mir_name = str_replace(mir_name, '^>', '')) %>%
  unite("mir_name", all_of(c("mir_name", "canonical_type")), sep = "-", na.rm = TRUE) %>%
  mutate(mir_name = str_replace(mir_name, '-M$', '')) %>%
  unite("isomir_name", all_of(c("polymorphic", "isomir_type")), remove = FALSE, na.rm = TRUE)

# isomir type wrong is assigned if no canonical form was found
# remove rows with no canonical form
# seperate the collapsed read column based on _x to get a count column
# set count to numeric
tb_isomir_counts <- tb_isomir %>%
  filter((is.na(isomir_type)) | isomir_type != "wrong") %>%
  separate(id, c("read_id", "count"), sep = "_x") %>%
  mutate(count = as.numeric(count))

tb_isomir_counts %>% write_delim(outFile, delim = "\t")

message("Wrote ", nrow(tb_isomir_counts), " classified reads to ", outFile)
