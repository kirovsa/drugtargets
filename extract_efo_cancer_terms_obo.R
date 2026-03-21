# Extract cancer-related terms from the EFO (Experimental Factor Ontology) OBO file.
#
# Reads efo.obo (OBO format), finds the top-level "cancer" term (EFO:0000311),
# and traverses all descendant branches via is_a relationships using
# breadth-first search.
#
# Output: a tab-delimited file with two columns:
#   efo_id      - EFO identifier in EFO:XXXXXXX format
#   description - human-readable term label

# Dependencies:
#   - obo: for reading OBO format into R
# Install:
#   install.packages("obo")

library(obo)

# ---------------------------------------------------------------------------- #
# Configuration
# ---------------------------------------------------------------------------- #

INPUT_FILE  <- "efo.obo"
OUTPUT_FILE <- "efo_cancer_terms.tsv"

# Top-level "cancer" class
CANCER_ROOT_ID <- "EFO:0000311"

# ---------------------------------------------------------------------------- #
# Helpers
# ---------------------------------------------------------------------------- #

# Build a children map from parent_id -> character vector of child_ids
build_children_map <- function(obo_obj) {
  children_map <- list()

  # obo_obj$stanzas is typically a named list keyed by IDs (e.g. EFO:0000311)
  stanzas <- obo_obj$stanzas
  if (is.null(stanzas) || length(stanzas) == 0) {
    stop("No stanzas found in parsed OBO.")
  }

  for (id in names(stanzas)) {
    stanza <- stanzas[[id]]

    # 'is_a' is usually a character vector like:
    #   "EFO:0000311 ! cancer"
    parents <- stanza[["is_a"]]
    if (is.null(parents)) next

    # extract parent IDs before any annotation after space
    parent_ids <- sub("[[:space:]].*$", "", parents)
    parent_ids <- parent_ids[nzchar(parent_ids)]

    for (p in parent_ids) {
      children_map[[p]] <- c(children_map[[p]], id)
    }
  }

  children_map
}

# Get label/name for a given term ID
get_term_label <- function(obo_obj, id) {
  stanza <- obo_obj$stanzas[[id]]
  if (is.null(stanza)) return(NA_character_)

  nm <- stanza[["name"]]
  if (is.null(nm) || length(nm) == 0) return(NA_character_)
  nm[[1]]
}

# ---------------------------------------------------------------------------- #
# Load ontology
# ---------------------------------------------------------------------------- #

message("Reading ", INPUT_FILE, " ...")
test <- read.obo(INPUT_FILE)

# ---------------------------------------------------------------------------- #
# Build lookup structures
# ---------------------------------------------------------------------------- #

children_map <- build_children_map(test)

# ---------------------------------------------------------------------------- #
# Breadth-first traversal from the cancer root
# ---------------------------------------------------------------------------- #

if (is.null(test$stanzas[[CANCER_ROOT_ID]])) {
  stop("Could not find cancer root term ", CANCER_ROOT_ID, " in ", INPUT_FILE)
}

visited <- character(0)
queue   <- CANCER_ROOT_ID

while (length(queue) > 0) {
  current <- queue[[1]]
  queue <- queue[-1]

  if (current %in% visited) next

  visited <- c(visited, current)

  child_ids <- children_map[[current]]
  if (!is.null(child_ids)) {
    new_children <- child_ids[!child_ids %in% visited]
    queue <- c(queue, new_children)
  }
}

message("Found ", length(visited), " cancer-related terms (including root).")

# ---------------------------------------------------------------------------- #
# Assemble output
# ---------------------------------------------------------------------------- #

descriptions <- vapply(visited, function(id) get_term_label(test, id), character(1))

result <- data.frame(
  efo_id      = visited,
  description = descriptions,
  stringsAsFactors = FALSE
)

# Remove any terms without a label (e.g. obsolete / anonymous nodes)
result <- result[!is.na(result$description) & nzchar(result$description), ]

# Sort by EFO number for reproducible output
# (keep only EFO:XXXXXXX IDs; anything else will sort last)
get_num <- function(x) {
  m <- regexpr("EFO:[0-9]{7}", x)
  if (m[[1]] == -1) return(Inf)
  as.numeric(sub("EFO:", "", regmatches(x, m)))
}

result <- result[order(vapply(result$efo_id, get_num, numeric(1))), ]
rownames(result) <- NULL

# ---------------------------------------------------------------------------- #
# Write output
# ---------------------------------------------------------------------------- #

write.table(result, file = OUTPUT_FILE, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

message("Wrote ", nrow(result), " rows to ", OUTPUT_FILE)