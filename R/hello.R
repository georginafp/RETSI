# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
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

re <- readRDS("/homes/users/gfuentes/scratch/projects/cyt_singlecell/data/_targets/objects/scATAC_CLUSTERED_LESS.ID")


## TESTING --------------------------------------------------------------------
# Randomly sample 500 cells and 500 peaks

## For testing
set.seed(123)  # Set seed for reproducibility
sampled_cells <- sample(colnames(re), size = 500, replace = FALSE)
sampled_peaks <- sample(rownames(re), size = 500, replace = FALSE)
subset_seurat_object <- subset(re, cells = sampled_cells, features = sampled_peaks)



# Subset the Seurat object
counts_matrix <- Seurat::GetAssayData(subset_seurat_object,
                                      assay = "ATAC",
                                      layer = "data") # normalized data

# Convert counts matrix to a data frame
counts_df <- as.data.frame(as.matrix(counts_matrix))


# Step 1: Get cell types for each cell ID
cell_types <- sapply(colnames(counts_df), function(cell_id)
  as.character(subset_seurat_object@active.ident[[cell_id]]))


# Step 2: Add cell type information to dataframe
long_df <- as.data.frame(t(counts_df)) %>%
  mutate(cell_id = names(cell_types),          # Assign cell IDs from names of cell_types
         cell_type = cell_types)               # Assign cell types from values of cell_types



# Step 2: Calculate mean values for each regulatory element by cell type
mean_df <- long_df %>%
  # filter(cell_type != "Unknown") %>%
  tidyr::pivot_longer(cols = starts_with("chr"),
               names_to = "regulatory_element",
               values_to = "expression_value") %>%
  group_by(cell_type, regulatory_element) %>%
  summarise(mean_value = mean(expression_value, na.rm = TRUE), .groups = 'drop')


# Step 3: Compute sum of means for each regulatory element
sum_means <- mean_df %>%
  group_by(regulatory_element) %>%
  summarise(sum_mean_value = sum(mean_value, na.rm = TRUE), .groups = 'drop')


# Step 4: Calculate RETSI
final_retsi_df <- distinct(mean_df, cell_type) %>%
  inner_join(mean_df %>%
               left_join(sum_means, by = "regulatory_element") %>%
               mutate(RETSI = ifelse(sum_mean_value > 0, mean_value / sum_mean_value, 0)) %>%
               select(cell_type, regulatory_element, RETSI),
             by = "cell_type")



RETSI_gr <- makeGRangesFromDataFrame(final_retsi_df %>%
                                       mutate(
                                         chromosome = sub("^(chr[0-9XY]+)-.*$", "\\1", regulatory_element),
                                         start = as.numeric(sub("^chr[0-9XY]+-(\\d+)-(\\d+)$", "\\1", regulatory_element)),
                                         end = as.numeric(sub("^chr[0-9XY]+-(\\d+)-(\\d+)$", "\\2", regulatory_element))
                                       ),
                                     keep.extra.columns = T)



library(ggplot2)
final_retsi_df %>%
  as.data.frame() %>%
  mutate(cell_type = factor(cell_type)) %>%
  ggplot(aes(x = RETSI,
             fill = cell_type, group = cell_type)) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed") +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("Alpha" = "#f5bc00",
                               "Ductal" = "#7CAE00",
                               "Beta-H" = "#f6483c",
                               "Beta-L" = "#fcbbb6",
                               "Acinar" = "#00BE67",
                               "Stellate" = "#C77CFF")) +
  scale_color_manual(values = c("Alpha" = "#f5bc00",
                               "Ductal" = "#7CAE00",
                               "Beta-H" = "#f6483c",
                               "Beta-L" = "#fcbbb6",
                               "Acinar" = "#00BE67",
                               "Stellate" = "#C77CFF")) +
  ggtitle("RETSI")





## applying
# RE_bulk <- readRDS("/homes/users/gfuentes/scratch/projects/coculture/data/_targets/objects/NEW_REGULATORY_ELEMENTS_IRF1_BIND")
# ols <- findOverlaps(RETSI_gr, RE_bulk, maxgap=200)
# re <- RETSI_gr[queryHits(ols),]
# mcols(RE_bulk) <- mcols(RE_bulk)[,ncol(mcols(RE_bulk))]
# colnames(mcols(RE_bulk)) <- paste0("BULK.", colnames(mcols(RE_bulk)))
# mcols(re) <- cbind(mcols(re), BULK_RE=mcols(RE_bulk)[subjectHits(ols),])



# Step 1: Get the normalized sparse matrix from Seurat object
counts_matrix <- Seurat::GetAssayData(re, assay = "ATAC", layer = "data")  # Retain sparse matrix format

# To avoid crushing work with the sparse format directly
# Step 2: Get cell types for each cell ID
cell_types <- as.character(re@active.ident)
cell_ids <- colnames(counts_matrix)

# Step 3: Instead of converting the entire matrix, we'll process it chunk-wise
# Add cell type information
# We'll handle the matrix in a chunk-wise manner to avoid memory overload

# Prepare to store results in chunks
mean_df_list <- list()

chunk_size <- 500  # Adjust this based on your system's available memory

# Process the matrix in chunks
for (i in seq(1, nrow(counts_matrix), by = chunk_size)) {
  # Define the range of rows (peaks) to process in the current chunk
  peak_range <- i:min(i + chunk_size - 1, nrow(counts_matrix))

  # Subset the sparse matrix for current chunk
  chunk_matrix <- counts_matrix[peak_range, ]

  # Convert sparse matrix chunk to dense format (only for the current chunk)
  chunk_df <- as.data.frame(as.matrix(chunk_matrix))

  # Transpose the chunk matrix and add cell type info
  long_chunk_df <- as.data.frame(t(chunk_df)) %>%
    mutate(cell_id = cell_ids,          # Assign cell IDs
           cell_type = cell_types)      # Assign cell types

  # Filter out unknown cell types and pivot longer for regulatory elements
  chunk_mean_df <- long_chunk_df %>%
    # filter(cell_type != "Unknown") %>%
    tidyr::pivot_longer(cols = starts_with("chr"),
                 names_to = "regulatory_element",
                 values_to = "expression_value") %>%
    group_by(cell_type, regulatory_element) %>%
    summarise(mean_value = mean(expression_value, na.rm = TRUE), .groups = 'drop')

  # Append result to list
  mean_df_list[[length(mean_df_list) + 1]] <- chunk_mean_df

  # Clean up memory manually
  rm(chunk_matrix, chunk_df, long_chunk_df, chunk_mean_df)
  gc()
}

# Combine all chunks into a single data frame after processing
mean_df <- bind_rows(mean_df_list)

# Step 4: Compute sum of means for each regulatory element
sum_means <- mean_df %>%
  group_by(regulatory_element) %>%
  summarise(sum_mean_value = sum(mean_value, na.rm = TRUE), .groups = 'drop')

# Step 5: Calculate RETSI
final_retsi_df2 <- distinct(mean_df, cell_type) %>%
  inner_join(mean_df %>%
               left_join(sum_means, by = "regulatory_element") %>%
               mutate(RETSI = ifelse(sum_mean_value > 0, mean_value / sum_mean_value, 0)) %>%
               select(cell_type, regulatory_element, RETSI),
             by = "cell_type")



final_retsi_df2 %>%
  as.data.frame() %>%
  mutate(cell_type = factor(cell_type)) %>%
  ggplot(aes(x = RETSI,
             fill = cell_type, group = cell_type)) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed") +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("Alpha" = "#f5bc00",
                               "Ductal" = "#7CAE00",
                               "Beta-H" = "#f6483c",
                               "Beta-L" = "#fcbbb6",
                               "Acinar" = "#00BE67",
                               "Stellate" = "#C77CFF")) +
  scale_color_manual(values = c("Alpha" = "#f5bc00",
                                "Ductal" = "#7CAE00",
                                "Beta-H" = "#f6483c",
                                "Beta-L" = "#fcbbb6",
                                "Acinar" = "#00BE67",
                                "Stellate" = "#C77CFF")) +
  ggtitle("RETSI")
