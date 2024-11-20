#' Find unique taxa in a specific sample type that is not present in any other sample type
#'
#'@title A tool to search for taxa that is unique to a specific sample type and not found in any other sample type in the phyloseq table
#'
#' @description
#' This function calculates the mean relative abundance and prevalence of each unique taxon.
#' automatically calculates mean relative abundance within the function so input a phyloseq object with raw counts
#' It includes both the hash ID and associated taxonomic ranks.
#'
#' @usage find_unique_taxa(physeq_obj, sample_type_column, target_sample_type)
#'
#' @param physeq_obj A phyloseq object.
#' @param sample_type_column The column title that has the sample type information (e.g "body.site").
#' @param target_sample_type The sample type of interest (e.g., "Gill").
#' @return A data frame with columns for taxon, taxonomic ranks, mean relative abundance, and prevalence that are all unique to the specific sample type.
#' @export
find_unique_taxa <- function(physeq_obj, sample_type_column, target_sample_type) {
  # Extract metadata
  metadata <- as.data.frame(sample_data(physeq_obj))

  # Extract sample names for the target sample type
  target_samples <- rownames(metadata[metadata[[sample_type_column]] == target_sample_type, , drop = FALSE])

  # Extract the OTU table and adjust orientation if necessary
  otu_table <- otu_table(physeq_obj)
  if (taxa_are_rows(physeq_obj)) {
    otu_table <- t(otu_table)
  }

  # Check if target samples exist in the OTU table
  valid_samples <- intersect(target_samples, rownames(otu_table))
  if (length(valid_samples) == 0) {
    stop("No matching samples found for the target sample type.")
  }

  # Convert to relative abundance
  relative_abundance <- otu_table / rowSums(otu_table)

  # Subset the relative abundance table for the target samples
  target_abundance <- relative_abundance[valid_samples, , drop = FALSE]

  # Subset the OTU table for all other samples
  other_samples <- setdiff(rownames(otu_table), valid_samples)
  other_abundance <- relative_abundance[other_samples, , drop = FALSE]

  # Identify taxa present only in the target sample type
  target_taxa <- colnames(target_abundance)[colSums(target_abundance > 0) > 0]
  other_taxa <- colnames(other_abundance)[colSums(other_abundance > 0) > 0]
  unique_taxa <- setdiff(target_taxa, other_taxa)

  # Extract taxonomic information
  tax_table <- tax_table(physeq_obj)
  tax_info <- as.data.frame(tax_table[unique_taxa, , drop = FALSE])

  # Calculate prevalence and mean relative abundance for the unique taxa
  if (length(unique_taxa) > 0) {
    prevalence <- colSums(target_abundance[, unique_taxa, drop = FALSE] > 0) / nrow(target_abundance)
    mean_relative_abundance <- colMeans(target_abundance[, unique_taxa, drop = FALSE])

    # Create the output dataframe
    result <- data.frame(
      Taxa = unique_taxa,
      Prevalence = prevalence,
      MeanRelativeAbundance = mean_relative_abundance
    )

    # Merge with taxonomic information
    result <- cbind(result, tax_info[match(result$Taxa, rownames(tax_info)), ])

    # Sort by prevalence (highest to lowest) and then by mean relative abundance
    result <- result[order(-result$Prevalence, -result$MeanRelativeAbundance), ]
  } else {
    result <- data.frame(
      Taxa = character(),
      Prevalence = numeric(),
      MeanRelativeAbundance = numeric()
    )
    warning("No unique taxa found for the target sample type.")
  }

  return(result)
}
