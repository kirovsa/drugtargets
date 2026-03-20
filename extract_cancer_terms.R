# Extract cancer-related terms from the Disease Ontology (DO) JSON file.
#
# Reads doid-base.json, finds the top-level "cancer" term (DOID:162), and
# traverses all descendant branches via is_a relationships using breadth-first
# search. Writes the results as a tab-delimited file with two columns:
#   do_id       - Disease Ontology identifier in DOID:XXXXXX format
#   description - human-readable term label

library(jsonlite)

# ---------------------------------------------------------------------------- #
# Configuration
# ---------------------------------------------------------------------------- #

INPUT_FILE  <- "doid-base.json"
OUTPUT_FILE <- "cancer_terms.tsv"

# OBO URI prefix and the DOID number for the top-level "cancer" class
DOID_URI_PREFIX  <- "http://purl.obolibrary.org/obo/DOID_"
CANCER_ROOT_DOID <- "162"

# ---------------------------------------------------------------------------- #
# Helper: convert OBO URI to DOID:XXXXXX identifier
# ---------------------------------------------------------------------------- #

uri_to_do_id <- function(uri) {
  doid_num <- sub(DOID_URI_PREFIX, "", uri, fixed = TRUE)
  paste0("DOID:", doid_num)
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

# Map from URI -> label for all nodes that have a label
node_labels <- list()
for (node in nodes) {
  lbl <- node[["lbl"]]
  if (!is.null(lbl) && nchar(lbl) > 0) {
    node_labels[[node[["id"]]]] <- lbl
  }
}

# Build children map: parent URI -> vector of child URIs (is_a edges only)
children_map <- list()
for (edge in edges) {
  if (edge[["pred"]] == "is_a") {
    parent <- edge[["obj"]]
    child  <- edge[["sub"]]
    children_map[[parent]] <- c(children_map[[parent]], child)
  }
}

# ---------------------------------------------------------------------------- #
# Breadth-first traversal from the cancer root
# ---------------------------------------------------------------------------- #

cancer_root_uri <- paste0(DOID_URI_PREFIX, CANCER_ROOT_DOID)

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

do_ids       <- uri_to_do_id(visited)
get_label <- function(uri) {
  lbl <- node_labels[[uri]]
  if (is.null(lbl)) NA_character_ else lbl
}

descriptions <- vapply(visited, get_label, character(1))

result <- data.frame(
  do_id       = do_ids,
  description = descriptions,
  stringsAsFactors = FALSE
)

# Remove any terms without a label (e.g. obsolete / anonymous nodes)
result <- result[!is.na(result$description), ]

# Sort by DOID number for reproducible output
doid_nums      <- as.numeric(sub("DOID:", "", result$do_id))
result         <- result[order(doid_nums), ]
rownames(result) <- NULL

# ---------------------------------------------------------------------------- #
# Write output
# ---------------------------------------------------------------------------- #

write.table(result, file = OUTPUT_FILE, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

message("Wrote ", nrow(result), " rows to ", OUTPUT_FILE)
