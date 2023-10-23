#' compare_samples
#'
#'
#'
#'
#' @param input1 first output
#' @param input2 second output
#' @export
#' @return A data.frame of seq.ID (columns) and sample.ID (rows) with either relative or absolute abundance of sequences.


compare_seqs <- function(input1, input2){


  output <- setdiff(unique(input1$seq.ID), unique(input2$seq.ID))

  return(output)

}
