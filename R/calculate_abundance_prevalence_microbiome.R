#' Calculate Abundance and Prevalence by Metadata Category
#'
#' @title A tool to search for a specific taxon's prevalence and mean relative abundance from your phyloseq object
#'
#' @description
#' Calculates the mean relative abundance and prevalence of a specific taxon across all levels of a specified category in your phyloseq object.
#' You input the hash ID of the taxon of interest, your phyloseq object, and the category you want to look at.
#' This function offers an option to automatically calculate relative abundances if the input phyloseq object contains raw counts.
#'
#' @param physeq Your phyloseq object.
#' @param taxon_id The hash ID of interest from your taxa table.
#' @param category The column name of your metadata (sample_data) of interest (e.g., "body.site").
#' @param convert_to_relative Logical; whether to convert the data to relative abundance. Default is TRUE.
#'
#' @return A data frame with mean relative abundance and prevalence by category.
#'
#' @usage calculate_abundance_prevalence(physeq, taxon_id, category, convert_to_relative = TRUE)
#'
#' @examples
#' \dontrun{
#'   # Assuming 'physeq' is a phyloseq object loaded in your environment
#'   calculate_abundance_prevalence(physeq, taxon_id, category, convert_to_relative = TRUE)
#' }
#'
#' @import phyloseq
#' @import dplyr
#'
#' @name calculate_abundance_prevalence
#' @export
calculate_abundance_prevalence <- function(physeq, taxon_id, category, convert_to_relative = TRUE) {
  # Check if the taxon exists in the phyloseq object
  if (!(taxon_id %in% taxa_names(physeq))) {
    stop("Taxon ID not found in the phyloseq object.")
  }

  # Convert to relative abundance if specified
  if (convert_to_relative) {
    physeq <- transform_sample_counts(physeq, function(x) x / sum(x))
  }

  # Determine the orientation of the OTU table
  otu_table_data <- otu_table(physeq)

  # Extract the relative abundance vector for the specific taxon
  if (taxa_are_rows(physeq)) {
    relative_abundance_values <- as.numeric(otu_table_data[taxon_id, ])
  } else {
    relative_abundance_values <- as.numeric(otu_table_data[, taxon_id])
  }

  # Extract sample data
  sample_data_df <- as.data.frame(sample_data(physeq))

  # Create a data frame with category and relative abundance
  results <- data.frame(
    Category = sample_data_df[[category]],
    RelativeAbundance = relative_abundance_values
  )

  # Calculate mean relative abundance and prevalence by category
  aggregated_results <- results %>%
    group_by(Category) %>%
    summarize(
      MeanRelativeAbundance = mean(RelativeAbundance, na.rm = TRUE),
      Prevalence = mean(RelativeAbundance > 0, na.rm = TRUE)
    )

  return(aggregated_results)
}
