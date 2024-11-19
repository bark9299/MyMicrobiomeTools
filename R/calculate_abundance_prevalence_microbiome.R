#' Calculate Abundance Prevalence
#'
#' @title A tool to search for a specific taxa's prevalence and mean relative abundance from your phyloseq object
#'
#' @description
#' Calculates the mean relative abundance and prevalence from a phyloseq object.
#' You input the hash id of the taxa of interest, your phyloseq object, and the category you want to look at.
#' This function automatically calculates relative abundances, so the input phyloseq object needs to be raw counts.
#' You can filter your phyloseq object beforehand to just have one site of interest.
#'
#' @param physeq Your phyloseq object.
#' @param taxon_id The hash id of interest from your taxa table.
#' @param category The taxonomic level of interest (e.g., genus or species).
#'
#' @return A data frame with mean relative abundance and prevalence by category.
#'
#' @usage calculate_abundance_prevalence(physeq, taxon_id, category)
#'
#' @examples
#' \dontrun{
#'   # Assuming 'physeq' is a phyloseq object loaded in your environment
#'   calculate_abundance_prevalence(physeq, taxon_id, category)
#' }
#'
#' @import phyloseq
#' @import microbiome
#' @import dplyr
#' @import devtools
#' @import stats
#'
#' @author Elyse Barker
#'
#' @name calculate_abundance_prevalence
#' @export
calculate_abundance_prevalence <- function(physeq, taxon_id, category) {
  # Check if the taxon exists in the phyloseq object
  if (!(taxon_id %in% taxa_names(physeq))) {
    stop("Taxon ID not found in the phyloseq object.")
  }

  # Calculate relative abundance for the entire phyloseq object
  relative_abundance_physeq <- transform_sample_counts(physeq, function(x) x / sum(x))

  # Extract the relative abundance vector for the specific taxon
  relative_abundance_values <- as.numeric(otu_table(relative_abundance_physeq)[taxon_id, ])

  # Extract sample data
  sample_data <- sample_data(physeq)

  # Create a data frame with category and relative abundance
  results <- data.frame(
    Category = sample_data[[category]],
    RelativeAbundance = relative_abundance_values
  )

  # Calculate mean relative abundance by category
  mean_relative_abundance <- aggregate(RelativeAbundance ~ Category, data = results, FUN = mean)

  # Calculate prevalence (proportion of samples with non-zero abundance) by category
  results$Presence <- relative_abundance_values > 0
  prevalence_by_category <- aggregate(Presence ~ Category, data = results, FUN = mean)

  # Combine the two metrics
  aggregated_results <- merge(mean_relative_abundance, prevalence_by_category, by = "Category")
  colnames(aggregated_results) <- c("Category", "MeanRelativeAbundance", "Prevalence")

  return(aggregated_results)
}

