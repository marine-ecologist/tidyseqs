#' Function to extract sequences from symportal output
#'
#'
#' same as extract_seqs but in long format
#'
#' @param folder location of the root Symportal output
#' @param type returns either "relative" or "absolute"
#' @param metadata location of a metadata file in csv, must contain a column called sample_name with matches
#' @param factors option to include a vector list of column names, e.g: c("location", "latitude") to not import every column in a metadata file, if left NULL then will import all columns
#' @export
#' @return A data.frame of seq.ID (columns) and sample_name (rows) with either relative or absolute abundance of sequences.


extract_seqs <-  function(folder, metadata=NULL, type = "absolute", factors=NULL, onlyprofile=FALSE, remove_zero=TRUE, ...) {

  # read absolute abundances:
  file_list <- list.files(path = folder, pattern = "seqs.absolute.abund_and_meta.txt", include.dirs = TRUE, recursive = TRUE)
  full_data <- read.delim(paste0(folder, "/", file_list)) %>%
    dplyr::select(sample_name, 40:ncol(.)) %>% # select just the symbiodinium columns
    dplyr::slice(-dplyr::n()) # remove the last row, summary data

  # read profiles
  its2_profile <- extract_its2_profile(folder)

  ### absolute abundance
  absolute <- full_data




  #columns matching "onlyprofile"
  if (isTRUE(onlyprofile)) {
    absolute <- absolute %>%
      dplyr::filter(sample_name %in% its2_profile$sample_name)
  }


  # tidy
  absolute <- absolute %>%
    tibble::column_to_rownames("sample_name") %>% # sample_name column to rowname
    dplyr::filter(rowSums(dplyr::select(., dplyr::where(is.numeric))) != 0) %>% # drop zero sum rows
    dplyr::select(dplyr::where(~ sum(. != 0) > 0)) %>% # drop zero sum columns
    dplyr::select(dplyr::where(~ any(!is.na(.)))) # drop blank columns


  relative <- absolute %>%
    dplyr::mutate(row_sum = rowSums(dplyr::select(., dplyr::where(is.numeric)))) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ . / row_sum)) %>%
    dplyr::select(-row_sum)

  # pivot

  absolute <- absolute %>%
    tibble::rownames_to_column("sample_name") %>%
    tidyr::pivot_longer(cols = -sample_name, names_to = "seq.ID", values_to = "abundance") %>%
    dplyr::filter(abundance>0.0001) %>%
    dplyr::mutate(seq.ID = stringr::str_replace(seq.ID, "^X", "")) %>% # drop X if first in seq.ID
    dplyr::group_by(sample_name) |>
    dplyr::arrange(desc(abundance))

  relative <- relative %>%
    tibble::rownames_to_column("sample_name") %>%
    tidyr::pivot_longer(cols = -sample_name, names_to = "seq.ID", values_to = "abundance") %>%
    dplyr::filter(abundance>0.0001) %>%
    dplyr::mutate(seq.ID = stringr::str_replace(seq.ID, "^X", "")) %>% # drop X if first in seq.ID
    dplyr::arrange(sample_name, desc(abundance))

  if (!is.null(metadata)) {
    absolute <- left_join(absolute, metadata, by="sample_name")
  }


  # return functions:
  if (type == "absolute") {
    return(absolute)
  } else if (type == "relative") {
    return(relative)
  }
}
