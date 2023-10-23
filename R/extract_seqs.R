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
#' @return A data.frame of seq.ID (columns) and sample.ID (rows) with either relative or absolute abundance of sequences.


extract_seqs <-  function(folder, metadata=NULL, type = "absolute", factors=NULL, onlyprofile=FALSE, remove_zero=TRUE, ...) {

  # read absolute abundances:
  file_list <- list.files(path = folder, pattern = "seqs.absolute.abund_and_meta.txt", include.dirs = TRUE, recursive = TRUE)
  absolute <- utils::read.delim(paste0(folder, "/", file_list)) %>%
    dplyr::select(sample_name, 40:ncol(.)) %>% # select just the symbiodinium columns
    dplyr::slice(-dplyr::n()) # remove the last row, summary data

  seqlist <- absolute |> dplyr::select(-sample_name) |> colnames()


  # remove zero rows and columns
  # if (!isTRUE(remove_zero)){
  # absolute <- absolute %>%
  # #  tibble::column_to_rownames("sample_name") %>% # sample_name column to rowname
  #   dplyr::filter(rowSums(dplyr::select(., dplyr::where(is.numeric))) != 0) %>% # drop zero sum rows
  #   dplyr::select(dplyr::where(~ sum(. != 0) > 0)) %>% # drop zero sum columns
  #   dplyr::select(dplyr::where(~ any(!is.na(.)))) # drop blank columns
  # }

  # create relative df
  relative <- absolute %>%
    dplyr::mutate(row_sum = rowSums(dplyr::select(., dplyr::where(is.numeric)))) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ . / row_sum)) %>%
    dplyr::select(-row_sum)



  # pivot_longer
  absolute <- absolute %>%
    #tibble::rownames_to_column("sample_name") %>%
    tidyr::pivot_longer(cols = -sample_name, names_to = "seq.ID", values_to = "abundance") %>%
    dplyr::filter(abundance>0.0001) %>%
    dplyr::mutate(seq.ID = stringr::str_replace(seq.ID, "^X", "")) %>% # drop X if first in seq.ID
    dplyr::group_by(sample_name) |>
    dplyr::arrange(sample_name,desc(abundance))

  relative <- relative %>%
   # tibble::rownames_to_column("sample.ID") %>%
    tidyr::pivot_longer(cols = -sample_name, names_to = "seq.ID", values_to = "abundance") %>%
    dplyr::filter(abundance>0.0001) %>%
    dplyr::mutate(seq.ID = stringr::str_replace(seq.ID, "^X", "")) %>% # drop X if first in seq.ID
    dplyr::arrange(sample_name, desc(abundance))


  #columns matching "onlyProfiles"
  if (isTRUE(onlyprofile)) {
    its2_profile <- extract_its2_profile(folder)

    options(warn=-1)
    # extract unique seq.ID from matching its2.profile
    its_seqs <- its2_profile %>%
      dplyr::rowwise() %>%
      dplyr::mutate(collapsed_column = paste(na.omit(dplyr::c_across(dplyr::starts_with("ITS2.profile"))), collapse = "-")) %>%
      dplyr::select(collapsed_column)
    its_seqs <- its_seqs |>
      tidyr::separate(collapsed_column, into = paste0("col", 1:max(stringr::str_count(its_seqs$collapsed_column, "[/-]"))), sep = "/|-") |>
      tidyr::pivot_longer(dplyr::everything(), names_to = "col", values_to = "value") %>%
      dplyr::select(-col) %>%
      dplyr::distinct() |> na.omit() |> dplyr::pull(value) |> suppressMessages()
    options(warn=0)
#      print(head(its_seqs))

    absolute <- absolute %>%
      dplyr::filter(seq.ID %in% its_seqs)
  }


  # read metadata
  if (!is.null(metadata)) {
    if (!is.null(factors)) {
      metadf <- metadata |>
        dplyr::select(sample_name, all_of(factors)) |>
        dplyr::rename(sample_name)
    } else {
      metadf <- metadata
    }
    absolute <- dplyr::left_join(absolute, metadf, by="sample_name") |> dplyr::ungroup()
    relative <- dplyr::left_join(relative, metadf, by="sample_name") |> dplyr::ungroup()
  }

  # if (isTRUE(remove_zero)) {
  #   absolute <- absolute %>% dplyr::filter(abundance != 0)
  #   relative <- relative %>% dplyr::filter(abundance != 0)
  # }

  # return functions:
  if (type == "absolute") {
    return(absolute)
  } else if (type == "relative") {
    return(relative)
  }
}
