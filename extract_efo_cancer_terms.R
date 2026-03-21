# Extract cancer-related terms from the EFO (Experimental Factor Ontology) OWL file.
#
# Reads efo.owl, finds the top-level "cancer" term (EFO:0000311), and
# traverses all descendant branches via subClassOf relationships using breadth-first
# search. Only EFO-prefixed nodes are included in the output.
# Writes the results as a tab-delimited file with two columns:
#   efo_id      - EFO identifier in EFO:XXXXXXX format
#   description - human-readable term label

library(xml2)

# ---------------------------------------------------------------------------- #
# Configuration
# ---------------------------------------------------------------------------- #

INPUT_FILE  <- "efo.owl"
OUTPUT_FILE <- "efo_cancer_terms.tsv"

# EFO URI prefixes and the EFO number for the top-level "cancer" class
EFO_URI_PREFIXES <- c(
  "http://www.ebi.ac.uk/efo/EFO_",
  "http://purl.obolibrary.org/obo/EFO_"
)
CANCER_ROOT_EFO <- "0000311"

# ---------------------------------------------------------------------------- #
# Helpers
# ---------------------------------------------------------------------------- #

is_efo_uri <- function(uri) {
  any(startsWith(uri, EFO_URI_PREFIXES))
}

# Convert any EFO URI to an EFO:XXXXXXX identifier
uri_to_efo_id <- function(uri) {
  m <- regexpr("EFO_[0-9]{7}", uri)
  if (m[[1]] == -1) return(NA_character_)
  efo_token <- regmatches(uri, m)
  sub("^EFO_", "EFO:", efo_token)
}

# ---------------------------------------------------------------------------- #
# Load ontology
# ---------------------------------------------------------------------------- #

message("Reading ", INPUT_FILE, " ...")
doc <- read_xml(INPUT_FILE)

ns <- c(
  owl  = "http://www.w3.org/2002/07/owl#",
  rdf  = "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
  rdfs = "http://www.w3.org/2000/01/rdf-schema#"
)

# ---------------------------------------------------------------------------- #
# Build lookup structures
# ---------------------------------------------------------------------------- #

# Find all owl:Class elements
classes <- xml_find_all(doc, "//owl:Class", ns)

# Map from URI -> label for all EFO nodes that have a label
node_labels <- list()

# Build children map: parent URI -> vector of child URIs (subClassOf only,
# restricted to EFO-prefixed nodes on both ends)
children_map <- list()

for (cls in classes) {
  attrs <- xml_attrs(cls)
  if (!"rdf:about" %in% names(attrs)) next
  uri <- attrs[["rdf:about"]]
  if (!is_efo_uri(uri)) next

  # Get rdfs:label
  lbl_nodes <- xml_find_all(cls, "rdfs:label", ns)
  if (length(lbl_nodes) > 0) {
    lbl <- xml_text(lbl_nodes[[1]])
    if (!is.na(lbl) && nchar(lbl) > 0) {
      node_labels[[uri]] <- lbl
    }
  }

  # Get named superclasses (subClassOf with rdf:resource, not anonymous restrictions)
  sc_nodes <- xml_find_all(cls, "rdfs:subClassOf", ns)
  for (sc in sc_nodes) {
    sc_attrs <- xml_attrs(sc)
    if (!"rdf:resource" %in% names(sc_attrs)) next
    parent_uri <- sc_attrs[["rdf:resource"]]
    if (!is_efo_uri(parent_uri)) next
    children_map[[parent_uri]] <- c(children_map[[parent_uri]], uri)
  }
}

# ---------------------------------------------------------------------------- #
# Breadth-first traversal from the cancer root
# ---------------------------------------------------------------------------- #

# Resolve the cancer root URI robustly: pick whichever IRI form (EBI or OBO
# PURL) actually appears in the parsed ontology (has a label or has children).
root_candidates <- c(
  paste0("http://www.ebi.ac.uk/efo/EFO_", CANCER_ROOT_EFO),
  paste0("http://purl.obolibrary.org/obo/EFO_", CANCER_ROOT_EFO)
)

cancer_root_uri <- root_candidates[
  which(vapply(root_candidates, function(u) {
    (!is.null(node_labels[[u]])) || (!is.null(children_map[[u]]))
  }, logical(1)))[1]
]

if (is.null(cancer_root_uri) || length(cancer_root_uri) == 0 || is.na(cancer_root_uri)) {
  stop("Could not resolve cancer root URI for EFO_", CANCER_ROOT_EFO,
       " in this OWL file (neither label nor children found for either candidate IRI).")
}

message("Using cancer root URI: ", cancer_root_uri)

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

efo_ids <- uri_to_efo_id(visited)
get_label <- function(uri) {
  lbl <- node_labels[[uri]]
  if (is.null(lbl)) NA_character_ else lbl
}

descriptions <- vapply(visited, get_label, character(1))

result <- data.frame(
  efo_id      = efo_ids,
  description = descriptions,
  stringsAsFactors = FALSE
)

# Remove any terms without a label (e.g. obsolete / anonymous nodes)
result <- result[!is.na(result$description), ]

# Sort by EFO number for reproducible output
efo_nums       <- as.numeric(sub("EFO:", "", result$efo_id))
result         <- result[order(efo_nums), ]
rownames(result) <- NULL

# ---------------------------------------------------------------------------- #
# Write output
# ---------------------------------------------------------------------------- #

write.table(result, file = OUTPUT_FILE, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

message("Wrote ", nrow(result), " rows to ", OUTPUT_FILE)