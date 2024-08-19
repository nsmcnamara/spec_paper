### This script is part of the Leaf Optical Properties Paper from the ACORN project.
### This script imports the raw data of spectral measurements, adds metadata,
### checks for outliers, calculates uncertainties and reflectance.
### Established 2024-08-15
### Last Update 2024-08-19
### Author: Simone McNamara

#### SETUP ####
# install.packages("devtools")
# library("devtools")
# install_github("meireles/spectrolab")
# install.packages("tidyverse")
# install.packages("conflicted")

# load libraries
library(spectrolab)
library(tidyverse)
library(conflicted)

# check and resolve conflicts
conflict_scout()
conflicts_prefer(
  dplyr::filter,
  dplyr::lag,
  spectrolab::smooth,
  spectrolab::combine
)

# define directory paths
res.path1 <- "/output/"
dat.path1 <- "/data/raw/"

# load functions
# source("source/rdadapt.R")

# select data subset
ss <- "AT_pubescens_2"
# ss <- "AT_robur_2"
# ss <- "CH_pubescens_2"
# ss <- "CH_robur_2"
# ss <- "CH_robur_3"

# select species
species <- if_else(str_detect(ss, "pubescens"), "Q.pubescens", "Q.robur")


#### DATA IMPORT ####
# import metadata
metadata <- read_csv("data/raw/metadata.csv")

# get file paths for data subset
file_paths <- list.dirs(paste0(getwd(), dat.path1, ss), recursive = FALSE)

# get file paths for metadata for each data subset
meta_files <- list.files(paste0(getwd(), dat.path1, ss), pattern = ".csv", include.dirs = FALSE)

# Create df to store data in
spec_df <- data.frame()

# import spectral measurement and merge with metadata
# Loop through each folder within the file paths
for (i in seq_along(file_paths)) {
  # Load spectral measurements
  spectra <- read_spectra(path = file_paths[i])

  # Load metadata
  meta <- read.csv(paste0(getwd(), dat.path1, ss, "/", meta_files[i]))

  # Create a table for merging spectra with measurement type
  mmt <- tibble(
    type = rep(c("WR", "WRL", "BR", "BRL"), each = 5, times = nrow(meta)),
    planting_location = rep(meta$planting_location, each = 20)
  )

  # Merge with spectra
  spectra_mmt <- cbind(mmt, spectra)

  # Merge with metadata, make sure species matches, because planting location in Austria is not unique
  spectra_mmt_meta <- right_join(metadata[metadata$species == species, ], spectra_mmt, by = "planting_location")

  # Store data in df
  spec_df <- rbind(spec_df, spectra_mmt_meta) # Append the new data
}


#### OUTLIER DETECTION 1: VISUAL INSPECTION ####
# visually inspect the plots of different measurement types

# WR
plot(as_spectra(filter(spec_df, type == "WR")[, 15:2165])) # check no single line looks "substantially different"

# check where mean is below 0.9
outlier_plants <- spec_df |>
  filter(type == "WR") |>
  rowwise() |>
  mutate(mean = mean(c_across(15:2165))) |>
  filter(mean < 0.9) |>
  ungroup() |>
  select(planting_location, sample_name, mean)

length(outlier_plants$planting_location)

# fine for AT_pubescens_2




# Create an empty matrix for Calculated Reflectance
# CR_current <- matrix(NA, nrow = nrow(current_merged) / 20, ncol = ncol(current_merged) - 3, dimnames = NULL)

# n = 1

# for (j in seq(1, nrow(current_merged), 20)) {
#   for (k in 4:ncol(current_merged)) {
#     WR = mean(current_merged[(j + 1):(j + 4), k])
#     WRL = mean(current_merged[(j + 6):(j + 9), k])
#     BR = mean(current_merged[(j + 11):(j + 14), k])
#     BRL = mean(current_merged[(j + 16):(j + 19), k])
#     CR_current[n, k - 3] = (BRL * WR - WRL * BR) / (WR - BR)
#   }
#   n = n + 1
# }


# # Rename column names
# colnames(CR_current_meta)[14:2164] <- 350:2500

# # Store processed data frame in the list
# processed_data_list[[i]] <- CR_current_meta

### NEXT TIME ###
# outlier functions: for all mmt types
# add LOF
# run on all data subsets.
# then calculate reflectance and uncertainties
