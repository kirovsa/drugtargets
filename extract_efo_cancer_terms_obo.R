# Extract cancer-related terms from the EFO (Experimental Factor Ontology) OBO file.
#
# Reads efo.obo (OBO format), finds the top-level "cancer" term (EFO:0000311),
# and traverses all descendant branches via is_a relationships (and all other
# relationship tags present in the OBO) using breadth-first search.
#
# Output: a tab-delimited file with two columns:
#   efo_id      - EFO identifier in EFO:XXXXXXX format
#   description - human-readable term label
#
# Dependencies:
#   - ontologyIndex: for reading OBO format into R
# Install:
#   install.packages("ontologyIndex")

library(ontologyIndex)

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

# Normalize a vector of relationship values to parent IDs.
#
# ontologyIndex may store relationship values either as plain IDs ("EFO:...")
# or as strings with trailing annotations. This helper extracts the first token
# and keeps only non-empty strings.
extract_parent_ids <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  x <- as.character(x)
  x <- sub("[[:space:]].*$", "", x)  # keep first token only
  x <- x[nzchar(x)]
  unique(x)
}

# Build a children map from parent_id -> character vector of child_ids.
#
# This uses:
#   - ont$is_a
#   - plus *all other relationship-like tags* that ontologyIndex exposed as
#     list-valued fields aligned to ont$id.
#
# We intentionally skip core structural fields (parents/children/ancestors etc.)
# to avoid double-counting or introducing computed links.
build_children_map <- function(ont) {
  children_map <- list()

  if (is.null(ont$id) || length(ont$id) == 0) {
    stop("No term IDs found in parsed OBO.")
  }

  ids <- ont$id

  # Fields we do not want to treat as relationship tags.
  skip_fields <- c(
    "id", "name",
    "parents", "children", "ancestors",
    "obsolete",
    # metadata fields frequently present when extract_tags = "everything"
    "format-version", "data-version", "subsetdef", "synonymtypedef",
    "idspace", "remark", "ontology", "owl-axioms",
    "domain", "range", "inverse_of",
    "is_transitive", "is_metadata_tag", "is_class_level",
    "is_functional", "is_inverse_functional", "is_symmetric",
    "holds_over_chain", "transitive_over",
    "expand_assertion_to", "expand_expression_to",
    "union_of", "intersection_of", "disjoint_from",
    "equivalent_to"
  )

  # Candidate relationship fields are list-valued and length == length(ids)
  # (i.e., aligned by term index).
  rel_fields <- setdiff(names(ont), skip_fields)
  rel_fields <- rel_fields[vapply(rel_fields, function(nm) {
    v <- ont[[nm]]
    is.list(v) && length(v) == length(ids)
  }, logical(1))]

  # Ensure is_a is always included (even if user provided it in skip_fields).
  rel_fields <- unique(c("is_a", rel_fields))

  for (i in seq_along(ids)) {
    child_id <- ids[[i]]

    for (field in rel_fields) {
      parents_raw <- ont[[field]][[i]]
      parent_ids <- extract_parent_ids(parents_raw)
      if (length(parent_ids) == 0) next

      for (p in parent_ids) {
        children_map[[p]] <- c(children_map[[p]], child_id)
      }
    }
  }

  # Deduplicate children lists
  for (p in names(children_map)) {
    children_map[[p]] <- unique(children_map[[p]])
  }

  children_map
}

# Get label/name for a given term ID
get_term_label <- function(ont, id) {
  idx <- match(id, ont$id)
  if (is.na(idx)) return(NA_character_)

  nm <- ont$name[[idx]]
  if (is.null(nm) || length(nm) == 0) return(NA_character_)
  nm[[1]]
}

# ---------------------------------------------------------------------------- #
# Load ontology
# ---------------------------------------------------------------------------- #

message("Reading ", INPUT_FILE, " ...")
toont <- get_ontology(INPUT_FILE, extract_tags = "everything")

# ---------------------------------------------------------------------------- #
# Build lookup structures
# ---------------------------------------------------------------------------- #

children_map <- build_children_map(ont)

# ---------------------------------------------------------------------------- #
# Breadth-first traversal from the cancer root
# ---------------------------------------------------------------------------- #

if (!(CANCER_ROOT_ID %in% ont$id)) {
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

message("Found \