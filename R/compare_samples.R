#' compare_samples
#'
#'
#'
#'
#' @param input1 first output
#' @param input2 second output
#' @export
#' @return A data.frame of seq.ID (columns) and sample.ID (rows) with either relative or absolute abundance of sequences.


compare_samples <- function(input1, input2){


  output <- setdiff(unique(input1$sample_name), unique(input2$sample_name))

  return(output)

}
