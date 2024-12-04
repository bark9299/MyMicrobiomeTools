#' Get Taxonomic Information for a Given Hash ID
#'
#' This function retrieves taxonomic information for a specific hash ID from a given phyloseq object.
#'
#' @param physeq A phyloseq object containing taxonomic information.
#' @param hash_id The hash ID of the taxon of interest.
#' @return A data frame containing the taxonomic ranks for the specified hash ID.
#' @export
hash_to_taxa <- function(physeq, hash_id) {
  # Check if input is a valid phyloseq object
  if (!inherits(physeq, "phyloseq")) {
    stop("The input must be a phyloseq object.")
  }

  # Extract the taxonomy table
  tax_table <- tax_table(physeq)

  # Check if the hash ID exists in the taxonomy table
  if (!hash_id %in% rownames(tax_table)) {
    stop("The specified hash ID does not exist in the phyloseq object.")
  }

  # Retrieve the taxonomic information for the given hash ID
  tax_info <- as.data.frame(tax_table[hash_id, , drop = FALSE])

  return(tax_info)
}

# Example usage (assuming you have a phyloseq object `physeq`):
# tax_info <- hash_to_taxa(physeq, "3cc0594c91659e5a54927e2f16deb354")
# print(tax_info)
