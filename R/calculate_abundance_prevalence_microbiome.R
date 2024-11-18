#' calculate_abundance_prevalence
#' @name MyMicrobiomeTools
#'
#' @title "A tool to search for a specific taxa's prevalence and mean relative abundance from your phyloseq object"
#'
#' @author Elyse Barker
#'
#' @description
#' calculates the mean relative abundance and prevalence from a phyloseq object
#' You input the hash id of the taxa of interest, your phyloseq object, and what category you want to look at
#' this function automatically calculates to relative abundances so input phyloseq object needs to be raw counts
#' need to have one column in the metadata be only one body site/environment of interest
#'
#' @usage calculate_abundance_prevalence(physeq, taxon_id, category)
#' @examples
#' # example code
#' #"body.site" column of the metadata only has the body site "gill". You can filter your phyloseq object beforehand to just have one site of interest.
#' result <- calculate_abundance_prevalence(gill_not_rel, "8d9fa267695027600ad8cf2eaca55c9c", "body.site")

# Load necessary libraries
library(phyloseq)
library(microbiome)
#this function automatically calculates to relative abundances so input phyloseq object needs to be raw counts
#need to have one column in the metadata be only one body site/environment of interest

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

