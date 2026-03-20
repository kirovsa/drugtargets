# Extract cancer-related terms from the Mondo Disease Ontology JSON file.
#
# Reads mondo.json, finds the top-level "cancer" term (MONDO:0004992), and
# traverses all descendant branches via is_a relationships using breadth-first
# search. Only MONDO-prefixed nodes are included in the output.
# Writes the results as a tab-delimited file with two columns:
#   mondo_id    - Mondo Disease Ontology identifier in MONDO:XXXXXXX format
#   description - human-readable term label

library(jsonlite)

# ---------------------------------------------------------------------------- #
# Configuration
# ---------------------------------------------------------------------------- #

INPUT_FILE  <- "mondo.json"
OUTPUT_FILE <- "mondo_cancer_terms.tsv"

# OBO URI prefix and the MONDO number for the top-level "cancer" class
MONDO_URI_PREFIX  <- "http://purl.obolibrary.org/obo/MONDO_"
CANCER_ROOT_MONDO <- "0004992"

# ---------------------------------------------------------------------------- #
# Helper: convert OBO URI to MONDO:XXXXXXX identifier
# ---------------------------------------------------------------------------- #

uri_to_mondo_id <- function(uri) {
  mondo_num <- sub(MONDO_URI_PREFIX, "", uri, fixed = TRUE)
  paste0("MONDO:", mondo_num)
}

# ---------------------------------------------------------------------------- #
# Load ontology
# ---------------------------------------------------------------------------- #

message("Reading ", INPUT_FILE, " ...")
ont <- fromJSON(INPUT_FILE, simplifyVector = FALSE)
graph <- ont[["graphs"]][[1]]

nodes <- graph[["nodes"]]
edges <- graph[["edges"]]

# ---------------------------------------------------------------------------- #
# Build lookup structures
# ---------------------------------------------------------------------------- #

# Map from URI -> label for all MONDO nodes that have a label
node_labels <- list()
for (node in nodes) {
  id <- node[["id"]]
  if (!startsWith(id, MONDO_URI_PREFIX)) next
  lbl <- node[["lbl"]]
  if (!is.null(lbl) && nchar(lbl) > 0) {
    node_labels[[id]] <- lbl
  }
}

# Build children map: parent URI -> vector of child URIs (is_a edges only,
# restricted to MONDO-prefixed nodes on both ends)
children_map <- list()
for (edge in edges) {
  if (edge[["pred"]] != "is_a") next
  parent <- edge[["obj"]]
  child  <- edge[["sub"]]
  if (!startsWith(parent, MONDO_URI_PREFIX)) next
  if (!startsWith(child,  MONDO_URI_PREFIX)) next
  children_map[[parent]] <- c(children_map[[parent]], child)
}

# ---------------------------------------------------------------------------- #
# Breadth-first traversal from the cancer root
# ---------------------------------------------------------------------------- #

cancer_root_uri <- paste0(MONDO_URI_PREFIX, CANCER_ROOT_MONDO)

visited <- character(0)
queue   <- cancer_root_uri

while (length(queue) > 0) {
  current <- queue[[1]]
  queue   <- queue[-1]

  if (current %in% visited) next

  visited <- c(visited, current)

  child_uris <- children_map[[current]]
  if (!is.null(child_uris)) {
    new_children <- child_uris[!child_uris %in% visited]
    queue <- c(queue, new_children)
  }
}

message("Found ", length(visited), " cancer-related terms (including root).")

# ---------------------------------------------------------------------------- #
# Assemble output data frame
# ---------------------------------------------------------------------------- #

mondo_ids <- uri_to_mondo_id(visited)
get_label <- function(uri) {
  lbl <- node_labels[[uri]]
  if (is.null(lbl)) NA_character_ else lbl
}

descriptions <- vapply(visited, get_label, character(1))

result <- data.frame(
  mondo_id    = mondo_ids,
  description = descriptions,
  stringsAsFactors = FALSE
)

# Remove any terms without a label (e.g. obsolete / anonymous nodes)
result <- result[!is.na(result$description), ]

# Sort by MONDO number for reproducible output
mondo_nums     <- as.numeric(sub("MONDO:", "", result$mondo_id))
result         <- result[order(mondo_nums), ]
rownames(result) <- NULL

# ---------------------------------------------------------------------------- #
# Write output
# ---------------------------------------------------------------------------- #

write.table(result, file = OUTPUT_FILE, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

message("Wrote ", nrow(result), " rows to ", OUTPUT_FILE)
