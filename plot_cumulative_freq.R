# Generate a cumulative-frequency plot for categorical cancer-type levels.
#
# Reads cancer_terms.tsv, classifies each term into a cancer-type category
# (carcinoma, sarcoma, lymphoma, etc.), then plots cumulative proportions
# for each category level at six sampling depths:
#   20, 50, 100, 200, 500, and the complete data frame.
#
# Output: cumulative_freq_plot.png

library(ggplot2)

# ---------------------------------------------------------------------------- #
# Configuration
# ---------------------------------------------------------------------------- #

INPUT_FILE  <- "cancer_terms.tsv"
OUTPUT_FILE <- "cumulative_freq_plot.png"

SAMPLE_SIZES <- c(20, 50, 100, 200, 500)   # plus complete

# ---------------------------------------------------------------------------- #
# Load data
# ---------------------------------------------------------------------------- #

message("Reading ", INPUT_FILE, " ...")
cancer_df <- read.table(INPUT_FILE, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)

message("Loaded ", nrow(cancer_df), " cancer terms.")

# ---------------------------------------------------------------------------- #
# Classify descriptions into cancer-type categories
# ---------------------------------------------------------------------------- #

# Classify a cancer description string into one of 11 mutually exclusive
# categories.  Patterns are checked from most-specific to least-specific so
# that compound terms are captured by the most appropriate bucket:
#   - "carcinoma" is checked BEFORE "adenoma" so that "adenocarcinoma"
#     (a malignant tumour) is correctly placed in Carcinoma rather than Adenoma.
#   - "glioblastoma" is included in the glioma branch and checked BEFORE the
#     generic "blastoma" branch, so glioblastomas remain in Glioma while
#     other blastomas (medulloblastoma, neuroblastoma, etc.) fall into Blastoma.
# Returns a single character value.
classify_type <- function(desc) {
  d <- tolower(desc)
  if (grepl("carcinoma",                                          d)) return("Carcinoma")
  if (grepl("sarcoma",                                            d)) return("Sarcoma")
  if (grepl("lymphoma",                                           d)) return("Lymphoma")
  if (grepl("leukemia|leukaemia",                                 d)) return("Leukemia")
  if (grepl("glioma|glioblastoma|astrocytoma|oligodendroglioma",  d)) return("Glioma")
  if (grepl("blastoma",                                           d)) return("Blastoma")
  if (grepl("melanoma",                                           d)) return("Melanoma")
  if (grepl("myeloma",                                            d)) return("Myeloma")
  if (grepl("adenoma",                                            d)) return("Adenoma")
  if (grepl("tumor|tumour",                                       d)) return("Tumor")
  return("Other")
}

cancer_df$cancer_type <- vapply(cancer_df$description, classify_type,
                                character(1))

# ---------------------------------------------------------------------------- #
# Determine category order from the complete data set (descending frequency)
# ---------------------------------------------------------------------------- #

full_freq   <- sort(table(cancer_df$cancer_type), decreasing = TRUE)
cat_levels  <- names(full_freq)           # fixed order for all plots

# ---------------------------------------------------------------------------- #
# Build cumulative-frequency data frame for each sampling depth
# ---------------------------------------------------------------------------- #

# Use a fixed random seed so results are reproducible across runs.
# Rows are sampled without replacement; the complete data frame uses all rows.
set.seed(42)
row_sample <- sample(nrow(cancer_df))

all_sizes  <- c(SAMPLE_SIZES, nrow(cancer_df))
size_labels <- c(as.character(SAMPLE_SIZES), "complete")

build_cum_freq <- function(n, label) {
  idx          <- row_sample[seq_len(n)]
  subset_types <- cancer_df$cancer_type[idx]
  freq_vec     <- table(factor(subset_types, levels = cat_levels))
  cum_prop     <- cumsum(as.numeric(freq_vec)) / sum(freq_vec)
  data.frame(
    category     = factor(cat_levels, levels = cat_levels),
    cum_freq     = cum_prop,
    sample_size  = factor(label, levels = size_labels)
  )
}

plot_df <- do.call(rbind, Map(build_cum_freq, all_sizes, size_labels))

# ---------------------------------------------------------------------------- #
# Plot
# ---------------------------------------------------------------------------- #

p <- ggplot(plot_df,
            aes(x = category, y = cum_freq,
                colour = sample_size, group = sample_size)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_brewer(palette = "Dark2",
                      name    = "Sample size") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "Cumulative frequencies of cancer-type categories",
    subtitle = "At sample sizes 20, 50, 100, 200, 500, and complete data frame",
    x        = "Cancer-type category (ordered by descending frequency in full data)",
    y        = "Cumulative frequency"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )

ggsave(OUTPUT_FILE, plot = p, width = 10, height = 6, dpi = 150)
message("Plot written to ", OUTPUT_FILE)
