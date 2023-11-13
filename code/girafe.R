extract_seqs(folder_path, type="absolute") |>
  filter_seqs(threshold=2000, type="relative") |>
  plot_seqs(type="ggplot")

extract_seqs(folder_path, type="absolute") |>
  filter_seqs(type="relative") |>
  plot_seqs(type="ggplot")


extract_seqs(folder_path, type="absolute") |>
  filter_seqs(keep_profiles=TRUE) |>
  plot_seqs(type="ggplot")

extract_seqs(folder_path, type="absolute") |>
  filter_seqs(drop_profiles=TRUE, type="relative") |>
  filter_seqs(keep_seqs=c("1611788_C", "37253_C", "1147954_C", "1645568_C"), type="absolute") |>
  plot_seqs(type="ggplot")


itsprofiles <- extract_its2_profile(folder_path) %>%
  mutate(seqids = str_split(ITS2.profile.1, pattern = "(?<!\\.)[-/]", n = Inf, simplify = FALSE)) |>
  pull(seqids) |>
  unlist() |>
  unique()


tmp <- extract_seqs(folder_path, type="absolute") |>
  filter_seqs(keep_profiles=TRUE)

library(ggiraph)
tmpplot <- ggplot(data = tmp,aes(x = sample_name, y = abundance, fill = seq.ID)) +
  #scale_y_reverse() +
  ggplot2::theme_bw() +
  ggiraph::geom_bar_interactive(aes(tooltip = seq.ID), color = "black", linewidth = 0.1, show.legend = FALSE, stat = "identity") +
  ggplot2::scale_fill_manual(values = colour.seqs_new) +
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#  geom_col_interactive() +
  theme(legend.position = 'none')

ggiraph::girafe(code = print(tmpplot), width = 6, point = 4, width_svg = 8, height_svg = 4)

