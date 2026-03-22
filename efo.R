library(ontologyIndex)

# Function to extract parent IDs
extract_parent_ids <- function(term_id, ont) {
  parents <- ont$obo$id[ont$obo$relationship=="is_a" & ont$id==term_id]
  return(parents)
}

# Function to build a map of children terms
build_children_map <- function(ont) {
  children_map <- list()
  for (term in ont$id) {
    children <- ont$obo$id[ont$obo$relationship=="is_a" & ont$obo$parent==term]
    children_map[[term]] <- children
  }
  return(children_map)
}

# Function to get the term label
get_term_label <- function(term_id, ont) {
  label <- ont$name[ont$id == term_id]
  return(label)
}

# Main processing function
INPUT_FILE <- "path/to/ontology/file.obo"
CANCER_ROOT_ID <- "EFO:0000311"

ont <- get_ontology(INPUT_FILE, extract_tags="everything")
children_map <- build_children_map(ont)

# Perform BFS to find all cancer terms
visited <- c()
queue <- c(CANCER_ROOT_ID)

while (length(queue) > 0) {
  term_id <- queue[1]
  queue <- queue[-1]

  if (term_id %in% visited) next
  visited <- c(visited, term_id)

  children <- children_map[[term_id]]
  for (child in children) {
    if (!(child %in% visited)) {
      queue <- c(queue, child)
    }
  }
}

# Write results to efo_cancer_terms.tsv
write.table(data.frame(term_id = visited, label = sapply(visited, get_term_label, ont = ont)),
            "efo_cancer_terms.tsv", sep="\t", row.names=FALSE, quote=FALSE)

message("Found ", length(visited), " cancer terms.")