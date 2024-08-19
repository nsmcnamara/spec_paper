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
# install.packages("DescTools")

# load libraries
library(spectrolab)
library(tidyverse)
library(conflicted)
library(DescTools)

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

inspect_by_type <- function(type, plot_range_min, plot_range_max, check_range_min, check_range_max, mean_threshold, filter_condition = "below") {
  # Calculate the mean and filter based on the threshold condition
  if (filter_condition == "below") {
    outlier_plants <- spec_df |>
      filter(type == !!type) |>
      rowwise() |>
      mutate(mean = mean(c_across(all_of(check_range_min:check_range_max)))) |>
      filter(mean < mean_threshold) |>
      ungroup() |>
      select(planting_location, sample_name, mean)
  } else {
    outlier_plants <- spec_df |>
      filter(type == !!type) |>
      rowwise() |>
      mutate(mean = mean(c_across(all_of(check_range_min:check_range_max)))) |>
      filter(mean > mean_threshold) |>
      ungroup() |>
      select(planting_location, sample_name, mean)
  }

  # Plot the spectra for the specific type (WR, WRL, BR, BRL)
  plot(as_spectra(filter(spec_df, type == !!type)[, 15:2165]),
       main = paste0(ss, " ", type))
if (nrow(outlier_plants) > 0) {
    plot(as_spectra(filter(spec_df, type == !!type & planting_location %in% outlier_plants$planting_location)[, plot_range_min:plot_range_max]), 
         col = "red", add = TRUE)
  }
  
  # Print the number of unique planting locations
  num_unique_locations <- length(unique(outlier_plants$planting_location))

  if (num_unique_locations > 0) {
    print(paste0("You have outliers for ", ss, " ", type))
    print(unique(outlier_plants$planting_location))
  } else {
    print(paste0("No outliers for ", ss, " ", type))
  }
}

# Call the function for each type
inspect_by_type("WR", 15, 2165, 15, 2165, 0.9, "below") # WR: Check where mean is below 0.9
inspect_by_type("WRL", 15, 2165, 145, 185, 0.4, "above") # WRL: Check where mean is above 0.4 (480 nm - 520 nm) - I added this myself
inspect_by_type("BR", 15, 2165, 15, 2165, 0.1, "above") # BR: Check where mean is above 0.1
inspect_by_type("BRL", 15, 2165, 145, 185, 0.2, "above") # BRL: Check where mean is above 0.2 (480 nm - 520 nm) - I added this myself




#### OUTLIER DETECTION 2: LOF ####
## k: The kth-distance to be used to calculate the LOFs.
## k 5 makes sense: always 5 measurements that should make a local cluster
## this identifies if one or more scans are abnormally far apart from each other
## irrespective of whether a single plant is just different from other plants


# LOF outliers detector
df <- spec_df
scan_type <- "BRL"

detect_lof_outliers <- function(df, lof_threshold = 2, k = 5) {
  # Define types
  scan_types <- levels(as.factor(df$type))
  
  for (scan_type in scan_types) {
    print(paste0("Processing: ", ss, " ", scan_type))
    
    # Filter by type
    spec_by_type <- df |>
      filter(type == scan_type)
    
    # Calculate LOF scores and insert before column of 350nm
    lof_by_type <- spec_by_type |>
      mutate(lof_score = LOF(subset(spec_by_type, select = `350`:`2500`), k = k), .before = "350")
    
    # Identify outliers
    lof_outliers <- lof_by_type |>
      filter(lof_score > lof_threshold)
    
    # if only one scan per plant is outlier: set all nm values of that single scan to NA
    single_outliers <- lof_outliers |>
      group_by(planting_location) |>
      filter(n() == 1)
    
    if (nrow(single_outliers) > 0) {
      print("You have single scan outliers:")
      print(paste0(single_outliers$planting_location, " scan: ", single_outliers$sample_name))
      
      # set all nm values of that single scan to NA
      df <- df |>
        mutate(across(`350`:`2500`,
                      ~ ifelse(sample_name %in% single_outliers$sample_name, NA, .)
        ))
      
    } else {
      print("You have no single scan outliers.")
    }
    
    
    # if more than one scan per plant is outlier: set all nm values for all scans for this plant to NA (eg B7, A20)
    mult_outliers <- lof_outliers |>
      group_by(planting_location) |>
      filter(n() > 1)
    
    if (nrow(mult_outliers) > 0) {
      print("You have multiple scan outliers:")
      print(paste0(mult_outliers$planting_location, " scan: ", mult_outliers$sample_name))
      
      df <- df |>
        mutate(across(`350`:`2500`,
                      ~ ifelse(planting_location %in% mult_outliers$planting_location, NA, .)
        ))
      
    } else {
      print("You have no multiple scan outliers.")
    }
  }
  
  return(df)
}

# Call the function
outliers_rm <- detect_lof_outliers(spec_df, lof_threshold = 1.2, 5)


plot(as_spectra(filter(spec_df, planting_location == "B_1")[, 15:2165]),
     main = paste0(ss, " ", "B_1"))
plot(as_spectra(filter(spec_df, sample_name == "at_pub_3_00035")[, 15:2165]), 
       col = "red", add = TRUE)
plot(as_spectra(filter(spec_df, sample_name %in% out)[, 15:2165]), 
     col = "red", add = TRUE)

out <- c("at_pub_1_00116", "at_pub_1_00117", "at_pub_1_00118", "at_pub_1_00119", "at_pub_1_00120")






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
# add LOF
# run on all data subsets.
# then calculate reflectance and uncertainties
