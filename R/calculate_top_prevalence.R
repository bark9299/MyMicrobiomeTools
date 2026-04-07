#' Calculate Mean Relative Abundance and Prevalence by Group
#'
#'@title A tool to search for taxa's prevalence and mean relative abundance from your phyloseq object by the top n most prevalent taxa
#'
#' @description
#' This function calculates the mean relative abundance and prevalence of each taxon
#' from a phyloseq object, grouped by a specified metadata column and sorted by prevalence.
#' It includes both the hash ID and associated taxonomic ranks.
#'
#' @usage calculate_top_prevalence(physeq, n = 20, group_var = NULL)
#'
#' @param physeq A phyloseq object.
#' @param n The number of top taxa to show (sorted by prevalence).
#' @param group_var The metadata column name to group by (e.g., "body.site").
#' @return A data frame with columns for group, taxon, taxonomic ranks, mean relative abundance, and prevalence.
#' @export
calculate_top_prevalence <- function(physeq, n = 20, group_var = NULL) {
  # Check if phyloseq object is valid
  if (!inherits(physeq, "phyloseq")) {
    stop("The input must be a phyloseq object.")
  }

  # Get the OTU table and taxonomic table
  otu_table <- as(otu_table(physeq), "matrix")
  otu_table <- apply(otu_table, 2, function(x) x / sum(x)) # Convert to relative abundance
  tax_table <- as(tax_table(physeq), "matrix")

  metadata <- sample_data(physeq)

  # Filter by group if specified
  if (!is.null(group_var)) {
    if (!group_var %in% colnames(metadata)) {
      stop("The specified group variable does not exist in the metadata.")
    }
    groups <- unique(metadata[[group_var]])
    results_list <- list()

    for (group in groups) {
      samples <- rownames(metadata[metadata[[group_var]] == group, ])
      otu_sub <- otu_table[, samples, drop = FALSE]

      mean_abundance <- rowMeans(otu_sub)
      prevalence <- rowSums(otu_sub > 0) / length(samples)

      df <- data.frame(
        TaxonID = rownames(otu_sub),
        Group = group,
        MeanRelativeAbundance = mean_abundance,
        Prevalence = prevalence
      )

      # Add taxonomic ranks to the data frame
      tax_ranks <- tax_table[rownames(df), , drop = FALSE]
      df <- cbind(df, tax_ranks)

      df <- df[order(-df$Prevalence), ] # Sort by prevalence
      df <- df[1:min(n, nrow(df)), ] # Select top n taxa
      results_list[[group]] <- df
    }

    results <- do.call(rbind, results_list)
  } else {
    mean_abundance <- rowMeans(otu_table)
    prevalence <- rowSums(otu_table > 0) / ncol(otu_table)

    results <- data.frame(
      TaxonID = rownames(otu_table),
      MeanRelativeAbundance = mean_abundance,
      Prevalence = prevalence
    )

    # Add taxonomic ranks to the data frame
    tax_ranks <- tax_table[rownames(results), , drop = FALSE]
    results <- cbind(results, tax_ranks)

    results <- results[order(-results$Prevalence), ] # Sort by prevalence
    results <- results[1:min(n, nrow(results)), ] # Select top n taxa
  }

  return(results)
}
