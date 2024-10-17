# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

library(GenomicRanges)
library(dplyr)
library(tibble)
library(Seurat)
library(tidyr)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(BiocParallel)
library(ggplot2)

# re <- readRDS("/homes/users/gfuentes/scratch/projects/cyt_singlecell/data/_targets/objects/scATAC_CLUSTERED_LESS.ID")



# Step 0: Convert Seurat object to SingleCellExperiment
##' Function that converts a Seurat object to a SingleCellExperiment object.
##' The function allows specifying which assay from the Seurat object
##' should be used for the conversion.
##' @param SO A Seurat object to be converted.
##' @param assay A character string specifying the assay to be used
##' for the conversion.
##' @return A SingleCellExperiment object containing the data from the
##' specified assay of the input Seurat object.
SeuratToSCE <- function(SO, assay) {
  obj <- Seurat::as.SingleCellExperiment(SO, assay = assay)
  return(obj)
}




# Step 1: Normalize counts in a SingleCellExperiment object
##' Function that normalizes raw count data in a
##' SingleCellExperiment object using a scaling factor.
##' The scaling factor is generated from a random normal distribution
##' and is applied to each cell to adjust the raw counts.
##' @param scSE A SingleCellExperiment object containing raw count data.
##' @return A SingleCellExperiment object with normalized counts,
##' where the raw counts have been adjusted based on the scaling factor.
##' The normalized counts can be accessed via the `normcounts` function.
getNormCounts <- function(scSE) {
  sf <- 2^rnorm(ncol(scSE))
  sf <- sf/mean(sf)
  normcounts(scSE) <- t(t(counts(scSE))/sf)
  return(scSE)
}




# Step 2: Extrat mean of normalized count data for each region and cell type
##' Function that computes the mean of normalized counts for each
##' regulatory region, grouped by cell type from a
##' SingleCellExperiment object.
##' @param scSE A SingleCellExperiment object containing normalized count data,
##' cell type identities, and cell IDs.
##' @return A data.frame with the mean of normalized counts for each
##' regulatory region, grouped by cell type and cell ID. The data frame
##' includes the following columns:
##' - `regulatory_element`: The regulatory region identifier.
##' - `mean_value`: The mean normalized count for the corresponding regulatory element.
##' - `cell_type`: The cell type associated with the mean value.
##' - `cell_id`: The ID of one cell belonging to the cell type (can be modified as needed).
TissueSpecAvg <- function(scSE) {
  counts_matrix <- normcounts(scSE)  # Normalized counts matrix
  cell_types <- colData(scSE)$ident  # Cell types
  cell_ids <- colnames(scSE)  # Cell IDs

  # Get unique cell types
  unique_cell_types <- unique(cell_types)

  # Register parallel processing
  num_workers <- BiocParallel::bpparam()$workers   # Automatically detect the number of workers
  register(MulticoreParam(workers = num_workers))  # Register parallel backend

  # Parallelize the operation to calculate mean for each unique cell type
  mean_list <- bplapply(unique_cell_types, function(current_cell_type) {
    # Get indices of cells for the current cell type
    cell_indices <- which(cell_types == current_cell_type)

    # Calculate row-wise means for the current cell type
    mean_values <- rowMeans(counts_matrix[, cell_indices, drop = FALSE], na.rm = TRUE)

    # Create a result data frame for this cell type
    data.frame(
      regulatory_element = rownames(counts_matrix),
      mean_value = mean_values,
      cell_type = current_cell_type,
      stringsAsFactors = FALSE
    )
  })

  # Combine results in single data frame
  mean_df <- do.call(rbind, mean_list)
  rownames(mean_df) <- NULL

  return(mean_df)
}





# Step 3: Extract sum of means for each regulatory region
##' Function that calculates the total sum of mean normalized counts
##' for each regulatory region, grouped by `regulatory_element`.
##' @param mean_df A data frame containing mean normalized counts,
##' grouped by `regulatory_element` and `cell_type`.
##' @return A data.frame with the total sum of mean normalized counts
##' for each regulatory element, including the following columns:
##' - `regulatory_element`: The regulatory element.
##' - `sum_mean_value`: The total sum of mean values for each regulatory element.
TotalTissueAvg <- function(mean_df) {
  sum_means <- mean_df %>%
    group_by(regulatory_element) %>%
    summarise(sum_mean_value = sum(mean_value, na.rm = TRUE), .groups = 'drop')

  return(sum_means)
}



# Step 4: Compute Regulatory Element Tissue Specific Index (RETSI)
##' Function that calculates the Regulatory Element Tissue Specific Index
##' for each regulatory region, using mean normalized counts and the total
##' sum of mean values for those regulatory elements.
##' @param mean_df A data frame containing mean normalized counts,
##' grouped by `regulatory_element` and `cell_type`.
##' @return A data.frame in wide format where each row corresponds to a
##' regulatory element and columns represent each `cell_type` with their
##' respective RETSI values, including the following columns:
##' - `regulatory_element`: The regulatory element (set as row names).
##' - Columns for each `cell_type` with their respective RETSI values.
computeRETSI <- function(mean_df, sum_means) {
  retsi_df <- mean_df %>%
    left_join(sum_means, by = "regulatory_element") %>%
    mutate(RETSI = ifelse(sum_mean_value > 0, mean_value / sum_mean_value, 0)) %>%
    select(cell_type, regulatory_element, RETSI) %>%
    distinct() %>%
    pivot_wider(names_from = cell_type, values_from = RETSI, names_prefix = "RETSI_") %>%
    column_to_rownames("regulatory_element")
  return(retsi_df)
}


# Step 5: Add RETSI data to Single Cell Experiment
##' Function that adds the RETSI data to the row data of a
##' Single Cell Experiment (SCE) object.
##' @param scSE A SingleCellExperiment object
##' @param retsi_df A data frame containing the RETSI values for
##' each regulatory element. It should be in wide format where row names
##' correspond to `regulatory_element`.
##' @return A SingleCellExperiment object with the updated row data
##' including the RETSI values.
addRETSIToSCE <- function(scSE, retsi_df) {
  rowData(scSE) <- retsi_df
  return(scSE)
}






# Main function to calculate RETSI and add it to the SingleCellExperiment
computeAndAddRETSI <- function(SO, assay) {

  scSE <- SeuratToSCE(SO=re, assay = "ATAC")
  scSE <- getNormCounts(scSE)
  mean_df <- TissueSpecAvg(scSE)
  sum_means <- TotalTissueAvg(mean_df)
  retsi_df <- computeRETSI(mean_df, sum_means)
  scSE <- addRETSIToSCE(scSE, retsi_df)

  return(scSE)
}




# ## PLOT distribution values
# library(ggplot2)
# final_retsi_df2 %>%
#   as.data.frame() %>%
#   mutate(cell_type = factor(cell_type)) %>%
#   ggplot(aes(x = RETSI,
#              fill = cell_type, group = cell_type)) +
#   geom_density(alpha = 0.5) +
#   geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed") +
#   theme(legend.position = "top") +
#   scale_fill_manual(values = c("Alpha" = "#f5bc00",
#                                "Ductal" = "#7CAE00",
#                                "Beta-H" = "#f6483c",
#                                "Beta-L" = "#fcbbb6",
#                                "Acinar" = "#00BE67",
#                                "Stellate" = "#C77CFF")) +
#   scale_color_manual(values = c("Alpha" = "#f5bc00",
#                                 "Ductal" = "#7CAE00",
#                                 "Beta-H" = "#f6483c",
#                                 "Beta-L" = "#fcbbb6",
#                                 "Acinar" = "#00BE67",
#                                 "Stellate" = "#C77CFF")) +
#   ggtitle("RETSI")
#
#
# ## applying
# # RE_bulk <- readRDS("/homes/users/gfuentes/scratch/projects/coculture/data/_targets/objects/NEW_REGULATORY_ELEMENTS_IRF1_BIND")
# # ols <- findOverlaps(RETSI_gr, RE_bulk, maxgap=200)
# # re <- RETSI_gr[queryHits(ols),]
# # mcols(RE_bulk) <- mcols(RE_bulk)[,ncol(mcols(RE_bulk))]
# # colnames(mcols(RE_bulk)) <- paste0("BULK.", colnames(mcols(RE_bulk)))
# # mcols(re) <- cbind(mcols(re), BULK_RE=mcols(RE_bulk)[subjectHits(ols),])
#
