folder_path <- ("/Volumes/Macintosh HD/Users/rof011/symbiodinium/20230120T102936_esampayo")

# get sequences
sequences <- extract_seqs(folder_path)

# get sequences and drop sample_names
sequences <- extract_seqs(folder_path) |>
             filter_seqs(drop_samples=c("AU18_0062","AU18_0063","AU18_0064","AU18_0065"))

#or
remove_these_seqs <- c("AU18_0062","AU18_0063","AU18_0064","AU18_0065")

sequences <- extract_seqs(folder_path) |>
  filter_seqs(drop_samples=remove_these_seqs)

