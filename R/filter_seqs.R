#' Function to extract sequences from symportal output
#'
#'
#' same as extract_seqs but in long format
#'
#' @param input location of input file
#' @param matchprofiles keep only seq.ID in matching profile?
#' @param clade filter by single "C" or multiple clades c("C", "D") to filter sequences by
#' @param drop_samples drop samples by named vector, e.g. c("H00B07", "H00B06"), or by one or more partial matches, e.g. c("07","B06")
#' @param keep_samples drop samples by named vector, e.g. c("H00B07", "H00B06"), or by one or more partial matches, e.g. c("07","B06")
#' @param drop_seqs drop seqs by named vector, e.g. c("X2777817_G", "X2777816_G"), or by one or more partial matches, e.g. c("X2","OT")
#' @param keep_seqs drop seqs by named vector, e.g. c("X2777817_G", "X2777816_G"), or by one or more partial matches, e.g. c("X2","OT")
#' @param threshold Set threshold to remove samples if less than the threshold (defaults to 1000)
#' @export
#' @return A data.frame of seq.ID (columns) and sample.ID (rows) with either relative or absolute abundance of sequences.


filter_seqs <- function(input, type="relative", drop_samples=NULL, clade = NULL, keep_samples=NULL, drop_seqs=NULL, keep_seqs=NULL, threshold=0, ...){


  # Drop rows by sample name
  if (!is.null(drop_samples)) {
    input <- input %>% as.data.frame() %>%
      dplyr::filter(!grepl(paste(drop_samples, collapse = "|"), sample_name))
  }

  # Drop rows by sample name
  if (!is.null(keep_samples)) {
    input <- input %>% as.data.frame() %>%
      dplyr::filter(grepl(paste(keep_samples, collapse = "|"), sample_name))
  }

  # Keep rows only with the sample names
  if (!is.null(drop_samples)) {
    input <- input %>% as.data.frame() %>%
      dplyr::filter(!grepl(paste(drop_samples, collapse = "|"), sample_name))
  }

  # Drop seq.ID by full match
  if (!is.null(drop_seqs)) {
    input <- input %>%
      dplyr::filter(!seq.ID %in% drop_seqs)
  }

  # Keep seq.ID by full match
  if (!is.null(keep_seqs)) {
    input <- input %>%
      dplyr::filter(seq.ID %in% keep_seqs)
  }

  # Keep seq.ID by full match
  if (!is.null(clade)) {
    input <- input %>%
      dplyr::filter(grepl(clade, seq.ID)) |>
      as.data.frame()
  }

  # filter less than 0
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
      dplyr::group_by(seq.ID) %>%
      dplyr::summarize(total_abundance = sum(abundance)) %>%
      dplyr::filter(total_abundance < threshold)

    input <- input |>
      dplyr::filter(!seq.ID %in% threshold_list$seq.ID)

  }

    # make relative
    relative <- input |>
      dplyr::group_by(sample_name) |>
      dplyr::mutate(total=sum(abundance, na.rm = TRUE))

    relative <- relative |>
      dplyr::mutate(abundance=abundance/total) |>
      dplyr::select(-total)


    # return functions:
    if (type == "absolute") {
      return(input)
    } else if (type == "relative") {
      return(relative)
    }



return(input)

}
