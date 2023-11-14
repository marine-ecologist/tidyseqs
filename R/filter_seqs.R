#' Function to extract sequences from symportal output
#'
#'
#' same as extract_seqs but in long format
#'
#' @param input location of input file
#' @param clade filter by single "C" or multiple clades c("C", "D") to filter sequences by
#' @param drop_samples drop samples by named vector, e.g. c("H00B07", "H00B06"), or by one or more partial matches, e.g. c("07","B06")
#' @param keep_samples drop samples by named vector, e.g. c("H00B07", "H00B06"), or by one or more partial matches, e.g. c("07","B06")
#' @param drop_seqs drop seqs by named vector, e.g. c("X2777817_G", "X2777816_G"), or by one or more partial matches, e.g. c("X2","OT")
#' @param keep_seqs drop seqs by named vector, e.g. c("X2777817_G", "X2777816_G"), or by one or more partial matches, e.g. c("X2","OT")
#' @param keep_profiles keep only seq.ID in matching ITS2 profile
#' @param drop_profiles keep only seq.ID in matching ITS2 profile
#' @param threshold Set threshold to remove samples if less than the threshold (defaults to 0)
#' @export
#' @return A data.frame of seq.ID (columns) and sample.ID (rows) with either relative or absolute abundance of sequences.



filter_seqs <- function(input, folder=NULL, type="relative",
                        drop_samples=NULL, keep_samples=NULL,
                        drop_seqs=NULL, keep_seqs=NULL,
                        keep_profiles=FALSE, drop_profiles=FALSE,
                        clade = NULL, threshold=0, ...){


  #-------- sample_name --------#
  # Drop
  if (!is.null(drop_samples)) {
    input <- input %>% as.data.frame() %>%
      dplyr::filter(!grepl(paste(drop_samples, collapse = "|"), sample_name))
  }

  # Keep
  if (!is.null(keep_samples)) {
    input <- input %>% as.data.frame() %>%
      dplyr::filter(grepl(paste(drop_samples, collapse = "|"), sample_name))
  }

  #-------- seq.ID --------#
  # drop
  if (!is.null(drop_seqs)) {
    input <- input %>%
      dplyr::filter(!seq.ID %in% drop_seqs)
  }

  # Keep
  if (!is.null(keep_seqs)) {
    input <- input %>%
      dplyr::filter(seq.ID %in% keep_seqs)
  }

  #-------- profiles --------#

  # Drop profiles
  if (isTRUE(drop_profiles)) {

    itsprofiles <- extract_its2_profile(folder_path) %>%
      mutate(seqids = str_split(ITS2.profile.1, pattern = "(?<!\\.)[-/]", n = Inf, simplify = FALSE)) |>
      pull(seqids) |>
      unlist() |>
      unique()

    input <- input %>%
      dplyr::filter(!seq.ID %in% itsprofiles)
  }

  # Keep profiles
  if (isTRUE(keep_profiles)) {

    itsprofiles <- extract_its2_profile(folder_path) %>%
      mutate(seqids = str_split(ITS2.profile.1, pattern = "(?<!\\.)[-/]", n = Inf, simplify = FALSE)) |>
      pull(seqids) |>
      unlist() |>
      unique()

    input <- input %>%
      dplyr::filter(seq.ID %in% itsprofiles)
  }

  #-------- clade --------#

  if (!is.null(clade)) {
    input <- input %>%
      dplyr::filter(grepl(clade, seq.ID)) |>
      as.data.frame()
  }

  #-------- threshold --------#

   if (!threshold==0) {

    precheck <- input |>
      dplyr::group_by(sample_name) |>
      dplyr::summarise(total=sum(abundance, na.rm=TRUE))

    if (sum(precheck$total == 1) >= 5) {
      warning("\n
              Filtering failed:\n
              Filtering by an abundance threshold will only work with abundance input data, e.g. extract_seqs(type=\"relative\").\n
              Re-extract with extract_seqs(type=\"absolute\") to filter by a threshold of sequence abundances \n \n")
      stop()
    }


    threshold_list <- input %>%
      dplyr::group_by(sample_name) %>%
      dplyr::summarize(total_abundance = sum(abundance)) %>%
      dplyr::filter(total_abundance < threshold)

    input <- input |>
      dplyr::filter(!sample_name %in% threshold_list$sample_name)

  }

  # convert to relative / absolute abundance
  absolute <- input

  relative <- input %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(total=sum(abundance, na.rm = TRUE)) %>%
    dplyr::mutate(abundance=abundance/total) %>%
    dplyr::select(-total)


  #return functions:
  if (type == "absolute") {
    return(absolute)
    } else if (type == "relative") {
    return(relative)
  }

}
