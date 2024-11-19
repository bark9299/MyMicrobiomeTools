#' Find Taxonomic Information
#' @title Find specific hash id in a phyloseq object by searching the name of a taxa
#'
#' @description Searches the taxonomic table of a phyloseq object for a specified taxonomic level and name.
#'
#' @param physeq Your phyloseq object.
#' @param tax_level Taxonomic level you are searching at (genus, species, etc.)
#' @param tax_name The taxonomic name you are searching for.
#'
#' @import phyloseq
#' @import microbiome
#' @import dplyr
#' @import devtools
#'
#' @author Elyse Barker
#'
#' @param physeq A phyloseq object containing the taxonomic table.
#' @param tax_level The taxonomic level to search (e.g., "Genus").
#' @param tax_name The name within the specified taxonomic level to find.
#'
#' @return A data frame with taxonomic information for the matching taxa.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'physeq' is a phyloseq object loaded in your environment
#' find_taxonomy_info(physeq, "Genus", "Acinetobacter")
#' }
#'
#' @name find_taxonomy_info
#' @export
find_taxonomy_info <- function(physeq, tax_level, tax_name) {
  # Check if physeq is a valid phyloseq object
  if (!inherits(physeq, "phyloseq")) {
    stop("Input must be a phyloseq object.")
  }

  # Extract the taxonomic table from the phyloseq object
  tax_table <- tax_table(physeq)
  tax_table_df <- data.frame(tax_table)

  # Ensure tax_level and tax_name are character strings
  if (!is.character(tax_level) || !is.character(tax_name)) {
    stop("Taxonomic level and name must be character strings.")
  }

  # Check if the specified taxonomic level exists in the taxonomic table
  if (!(tax_level %in% colnames(tax_table_df))) {
    stop("Specified taxonomic level not found in the taxonomic table.")
  }

  # Find the rows where the specified taxonomic level matches the input tax_name
  tax_matches <- which(tax_table_df[, tax_level] == tax_name, arr.ind = TRUE)

  # Check if any matches are found
  if (length(tax_matches) == 0) {
    stop("No matches found for the specified taxonomic name.")
  }

  # Extract the information for the matching taxa

  tax_info <- tax_table_df[tax_matches, ]

  # Include the hash IDs as rownames
  tax_info <- data.frame(HashID = rownames(tax_info), tax_info, row.names = NULL)

  return(tax_info)
}
