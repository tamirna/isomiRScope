# isomiRScope - classification of isomiRs from miRDeep2 mrd files
#
# Copyright 2026, TAmiRNA GmbH
# Released under the GNU General Public License v3. See LICENSE and COPYRIGHT.
#
# Extracted from the miND small RNA sequencing pipeline.

#' mrd miRNA list to named mrd tibbles
#'
#' @param tb tibble containing one column with mrd info
#' @param name name of the miR/list
#'
#' @return a tibble with an id column and a value column storing the name of the miRNA inside the tibble
#' @export
#'
#' @examples
prepare_tb_list <- function(tb, name) {
  
  # seperate the tibble column based on unknown number of white spaces
  # add a row based on the list name storing the miR name info
  # remove line endings and mutation information
  tb_list_out <- tb %>% 
    separate(value, into = c("id", "value"), sep = "\\s\\s+") %>% 
    add_row(id = "mirbase_id" , value = name, .before = 1) %>% 
    mutate(value = str_replace(value, '\t.*$', ''))
  
  # return only a tibble if reads (seq) were found for a miRNA
  if(tb_list_out %>% 
     pull(id) %>% 
     grepl(pattern = "^seq") %>% 
     any()) {
    return(tb_list_out)
  } else {}
  
}

#' Prepare mir precursor information
#'
#' @param tb_list_in a named tibble list containing mrd tibbles
#'
#' @return a tibble containing precursor information
#' @export
#'
#' @examples
prepare_annotation_tb <- function(tb_list_in) { 
  
  # columns storing annotation information of the precursor
  # mirbase id ... miRNA name
  # exp ... information of the location of the canonical miRNAs
  # pri_seq ... sequence of the pre miR
  # pri_struct ... structure/folding information
  
  annotation_info <- c("mirbase_id", 
                       "exp", 
                       "pri_seq", 
                       "pri_struct")
  
  # only selection rows with annotation info
  # bring in wide format
  # add columns for starting and end position of the canonical sequence
  tb_annotated <- tb_list_in %>% 
    dplyr::filter(id %in% annotation_info) %>%
    pivot_wider(names_from = id) 
    
  tb_annotated <- tb_annotated %>%
    mutate("5_start" = case_when(grepl("5", exp) ~ (gregexpr("5", exp) %>%
                                   unlist() %>%
                                   head(n=1) %>%
                                   as.numeric()),
                                 TRUE ~ NA),
           "5_end" = case_when(grepl("5", exp) ~ (gregexpr("5", exp) %>%
                                                    unlist() %>%
                                                    tail(n=1) %>%
                                                    as.numeric()),
                               TRUE ~ NA),
           "3_start" = case_when(grepl("3", exp) ~ (gregexpr("3", exp) %>%
                                                      unlist() %>%
                                                      head(n=1) %>%
                                                      as.numeric()),
                                 TRUE ~ NA),
           "3_end" = case_when(grepl("3", exp) ~ (gregexpr("3", exp) %>%
                                                    unlist() %>%
                                                    tail(n=1) %>%
                                                    as.numeric()),
                               TRUE ~ NA),
           "M_start" = case_when(grepl("M", exp) ~ (gregexpr("M", exp) %>%
                                                      unlist() %>%
                                                      head(n=1) %>%
                                                      as.numeric()),
                                 TRUE ~ NA),
           "M_end" = case_when(grepl("M", exp) ~ (gregexpr("M", exp) %>%
                                                    unlist() %>%
                                                    tail(n=1) %>%
                                                    as.numeric()),
                               TRUE ~ NA))
  return(tb_annotated)
}

#' Prepare start and end positions of mapped reads
#'
#' @param tb_list_in list of named mrd tibbles
#'
#' @return named mrd tibble containing position information
#' @export
#'
#' @examples
prepare_position_tb <- function(tb_list_in) {  
  
  # specify row wise operation
  # add the character position of the first and last occurrence of a character value
  # seq_start ... character position of first character (not .)
  # seq_end ... character position of last character (before . continuous)
  # ungroup to remove row wise information 
  tb_position <- tb_list_in %>% 
    rowwise() %>%
    mutate("seq_start" = case_when(grepl('\\.', value) ~ (gregexpr('[a-zA-Z]', value) %>%
                                                            unlist() %>%
                                                            head(n=1)),
                                   TRUE ~ NA),
           "seq_end" = case_when(grepl('\\.', value) ~ (gregexpr('[a-zA-Z]', value) %>%
                                                          unlist() %>%
                                                          tail(n=1)),
                                 TRUE ~ NA)) %>%
    ungroup()
  
  return(tb_position)
}

#' Classify isoMirs
#'
#' @param mir_annotation tibble row containing precursor information (start and end postions of the canonical sequences)
#' @param tb_list_position named mrd tibble list with position information of each sequence
#'
#' @return named mrd tibble list with position and isomir info
#' @export
#'
#' @examples
annotate_isomirs <- function(mir_annotation, tb_list_position) {
  
  # extract the miRNA name of the annotation
  mir_name <- mir_annotation %>% 
    pull(mirbase_id)
  
  start_5 <- mir_annotation %>% 
    pull(`5_start`)
  
  end_5 <- mir_annotation %>% 
    pull(`5_end`)
  
  start_3 <- mir_annotation %>% 
    pull(`3_start`)
  
  end_3 <- mir_annotation %>% 
    pull(`3_end`)
  
  start_M <- mir_annotation %>% 
    pull(`M_start`)
  
  end_M <- mir_annotation %>%
    pull(`M_end`)
  
  
  # perform row wise operation
  # add the information of the canonical sequence
  # mirDeep2 - The number of reads falling into an interval 2nt upstream and 5nt downstream of the mature/star sequence is determined.
  # start position has to be between -2nt from the canonical start and + 50% of the sequence length e.g. 20 nt 5'_start + 10
  # end position has to be between  50% of the sequence length and +5nt from the canonical end e.g. 20 nt 5'_end - 10
  # isomir type is added based on the classification rules described in docs/algorithm.md
  # iso-3 / iso-5 / iso-3-5
  # check for polymorphisms - upper case letter in sequence indicates a mutation
  # ungroup to remove row wise condition
  tb_list_isomir <- tb_list_position[[mir_name]] %>%
    rowwise() %>%
    mutate(canonical_type = case_when(
      between(
        seq_start,
        start_5 - 2,
        start_5 + ((seq_end - seq_start) / 2)
      ) & between(
        seq_end,
        end_5 - ((seq_end - seq_start) / 2),
        end_5 + 5
      ) ~ "5p",
      between(
        seq_start,
        start_3 - 2,
        start_3 + ((seq_end - seq_start) / 2)
      ) & between(
        seq_end,
        end_3 - ((seq_end - seq_start) / 2),
        end_3 + 5
      ) ~ "3p",
      between(
        seq_start,
        start_M - 2,
        start_M + ((seq_end - seq_start) / 2)
      ) & between(
        seq_end,
        end_M - ((seq_end - seq_start) / 2),
        end_M + 5
      ) ~ "M",
      TRUE ~ NA
    )) %>%
    mutate(isomir_type = case_when(
      canonical_type == "5p" & seq_start == start_5 & seq_end == end_5 ~ "canonical",
      canonical_type == "5p" & seq_start != start_5 & seq_end != end_5 ~ "iso-3-5",
      canonical_type == "5p" & seq_start != start_5 ~ "iso-5",
      canonical_type == "5p" & seq_end != end_5 ~ "iso-3",
      canonical_type == "3p" & seq_start == start_3 & seq_end == end_3 ~ "canonical",
      canonical_type == "3p" & seq_start != start_3 & seq_end != end_3 ~ "iso-3-5",
      canonical_type == "3p" & seq_start != start_3 ~ "iso-5",
      canonical_type == "3p" & seq_end != end_3 ~ "iso-3",
      canonical_type == "M" & seq_start == start_M & seq_end == end_M ~ "canonical",
      canonical_type == "M" & seq_start != start_M & seq_end != end_M ~ "iso-3-5",
      canonical_type == "M" & seq_start != start_M ~ "iso-5",
      canonical_type == "M" & seq_end != end_M ~ "iso-3",
      TRUE ~ "wrong"
    )) %>%
    mutate(polymorphic = case_when(
      grepl("[A-Z]", value) ~ "poly",
      TRUE ~ NA
    )) %>%
    ungroup()
  
  # return a named list
  tmplist <- list(tb_list_isomir)
  names(tmplist) <- mir_name
  return(tmplist)
}

