### This script is part of the Leaf Optical Properties Paper from the ACORN project.
### This script imports the raw data of spectral measurements, adds metadata,
### checks for outliers, calculates uncertainties and reflectance.
### Established 2024-08-15
### Last Update 2024-08-20
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
ss <- "AT_robur_2"
ss <- "CH_pubescens_2"
ss <- "CH_robur_2"
ss <- "CH_robur_3"

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
inspect_by_type <- function(df) {
  # Define scan types
  scan_types <- levels(as.factor(df$type))

  # Conditions for each scan type
  conditions <- data.frame(
    scan_types = scan_types,
    threshold = c(0.1, 0.2, 0.9, 0.4),
    operator = c("greater", "greater", "smaller", "greater"),
    mean_range_min = c(350, 480, 350, 480),
    mean_range_max = c(2500, 520, 2500, 520)
  )

  # Loop through scan types
  for (scan_type in scan_types) {
    print(paste0("Processing: ", ss, " ", scan_type))

    # Retrieve conditions for the current scan type
    threshold <- conditions$threshold[conditions$scan_types == scan_type]
    operator <- conditions$operator[conditions$scan_types == scan_type]
    mean_range_min <- as.character(conditions$mean_range_min[conditions$scan_types == scan_type])
    mean_range_max <- as.character(conditions$mean_range_max[conditions$scan_types == scan_type])

    # Filter the data by type
    spec_by_type <- df |>
      filter(type == scan_type)

    # Calculate mean for the specified range
    mean_by_type <- spec_by_type |>
      rowwise() |>
      mutate(mean = mean(c_across(all_of(mean_range_min):all_of(mean_range_max))), .before = `350`) |>
      ungroup()

    # Identify outliers based on the operator condition
    if (operator == "greater") {
      outliers <- mean_by_type |>
        filter(mean > threshold)
    } else {
      outliers <- mean_by_type |>
        filter(mean < threshold)
    }

    # Plot the spectra for the specific type
    plot(as_spectra(subset(spec_by_type, select = `350`:`2500`)),
      main = paste0(ss, " ", scan_type)
    )
    # If outliers are detected, plot them in red
    if (nrow(outliers) > 0) {
      plot(as_spectra(subset(outliers, select = `350`:`2500`)),
        col = "red", add = TRUE
      )
      legend("topright", legend = outliers$planting_location, title = "Outliers")

      print(paste0("You have outliers for ", ss, " ", scan_type))
      print(unique(outliers$planting_location))

      # Set outliers to NA in the main df
      df <- df |>
        mutate(across(
          `350`:`2500`,
          ~ ifelse(planting_location %in% outliers$planting_location, NA, .)
        ))
    } else {
      print(paste0("No outliers for ", ss, " ", scan_type))
    }
  }

  # Return modified dataframe
  return(df)
}

# call
vis_outlier_rm <- inspect_by_type(spec_df)

# at pub: rm B20
# at rob: all clear
# ch pub: 8E3 too low around 1000 in BRL, manual rm:

# ch rob 2: all clear
# ch rob 3: check BRL high in 1000 but seems otherwise ok so leave it

# manual check:
# tst <- spec_df |>
#   filter(type == "BRL")
# plot(as_spectra(subset(tst, select = `350`:`2500`)),
#      main = paste0(ss, " BRL ")
# )
# tst2 <- tst |>
#   filter(`1000` > 0.55)
# plot(as_spectra(subset(tst2, select = `350`:`2500`)),
#      col = "red", add = TRUE
# )
# manual rm:
# vis_outlier_rm <- vis_outlier_rm |>
#   mutate(across(
#     `350`:`2500`,
#     ~ ifelse(planting_location == "8_E_3", NA, .)))

#### OUTLIER DETECTION 2: LOF ####
## k: The kth-distance to be used to calculate the LOFs.
## k 5 makes sense: always 5 measurements that should make a local cluster
## this identifies if one or more scans are abnormally far apart from each other
## irrespective of whether a single plant is just different from other plants


# LOF outliers detector
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

      # plot
      plot(as_spectra(subset(spec_by_type, select = `350`:`2500`)),
        main = paste0(ss, " ", scan_type)
      )
      plot(as_spectra(subset(single_outliers, select = `350`:`2500`)),
        col = "red", add = TRUE
      )
      legend(x = "topright", legend = single_outliers$planting_location, title = "Single Scan Outliers")

      # set all nm values of that single scan to NA
      df <- df |>
        mutate(across(
          `350`:`2500`,
          ~ ifelse(sample_name %in% single_outliers$sample_name, NA, .)
        ))
    } else {
      print("You have no single scan outliers.")
    }


    # if more than one scan per plant is outlier: set all nm values for all scans for this plant to NA
    mult_outliers <- lof_outliers |>
      group_by(planting_location) |>
      filter(n() > 1)

    if (nrow(mult_outliers) > 0) {
      print("You have multiple scan outliers:")
      print(paste0(mult_outliers$planting_location, " scan: ", mult_outliers$sample_name))

      # plot
      plot(as_spectra(subset(spec_by_type, select = `350`:`2500`)),
        main = paste0(ss, " ", scan_type)
      )
      plot(as_spectra(subset(mult_outliers, select = `350`:`2500`)),
        col = "red", add = TRUE
      )
      legend(
        x = "topright", legend = paste0(mult_outliers$planting_location, " ", mult_outliers$sample_name),
        title = "Multi Scan Outliers"
      )

      # Set all 20 scans for that plant to NA
      df <- df |>
        mutate(across(
          `350`:`2500`,
          ~ ifelse(planting_location %in% mult_outliers$planting_location, NA, .)
        ))
    } else {
      print("You have no multiple scan outliers.")
    }
  }

  return(df)
}

# Call the function
lof_outliers_rm <- detect_lof_outliers(vis_outlier_rm, lof_threshold = 2, 5)

# man plotting 
# tst <- vis_outlier_rm |>
#   filter(type == "BR")
# plot(as_spectra(subset(tst, select = `350`:`2500`)),
#      main = paste0(ss, " BR ")
# )
# tst2 <- tst |>
#   filter(planting_location == "38_G_2")
# plot(as_spectra(subset(tst2, select = `350`:`2500`)),
#      col = "red", add = TRUE
# )

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
# then calculate reflectance and uncertainties
