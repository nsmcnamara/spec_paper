### This script is part of the Leaf Optical Properties Paper from the ACORN project.
### This script imports the raw data of spectral measurements, adds metadata,
### checks for outliers, calculates uncertainties and reflectance.
### Established 2024-08-15
### Last Update 2024-09-09
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
source("source/inspect_by_type.R")
source("source/detect_lof_outliers.R")


# select data subset
ss <- "AT_pubescens_2"
# ss <- "AT_robur_2"
# ss <- "CH_pubescens_2"
# ss <- "CH_robur_2"
# ss <- "CH_robur_3"

# select species
if (str_detect(ss, "pubescens")) {
  species <- "Q.pubescens"
} else if (str_detect(ss, "robur")) {
  species <- "Q.robur"
} else {
  print("Could not detect species")
}


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

# call function
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
#     ~ ifelse(planting_location == "B_20", NA, .)))

#### OUTLIER DETECTION 2: LOF ####
## k: The kth-distance to be used to calculate the LOFs.
## k 5 makes sense: always 5 measurements that should make a local cluster
## this identifies if one or more scans are abnormally far apart from each other
## irrespective of whether a single plant is just different from other plants

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

#### CALCULATE REFLECTANCE ####

# enter here the metadata columns for grouping WITHOUT type and sample name
metadata_cols <- colnames(lof_outliers_rm[1:12])

# calculate means of each measurement type for each plant (i.e., mean of 5 scans)
mean_by_type <- lof_outliers_rm |>
  # group by all metadata columns
  group_by(across(c(metadata_cols, type))) |>
  # caluculate mean for each of the wavelengths
  summarise(across(`350`:`2500`, ~ mean(., na.rm = TRUE)))

# calculate reflectance
CR <- mean_by_type |>
  # pivot data frame so each measurement has its own row
  pivot_longer(
    cols = `350`:`2500`,
    names_to = "nm",
    values_to = "val"
  ) |>
  # make nm as numeric so it will be ordered in descending
  mutate(nm = as.numeric(nm)) |>
  # pivot again so we have one row for each nm and can calculate reflectance per nm
  pivot_wider(
    names_from = type,
    values_from = val
  ) |>
  # group by all metadata columns
  group_by(across(c(metadata_cols, nm))) |>
  # calculate reflectance
  summarise(CR = (BRL * WR - WRL * BR) / (WR - BR)) |>
  # pivot again so each plant has one row
  pivot_wider(
    names_from = nm,
    values_from = CR
  ) |>
  ungroup() |>
  # replace NaN with NA
  mutate_all(~ifelse(is.nan(.), NA, .))
  
 
#### TRIM ####
CR_trim <- CR |>
  select(-c(`350`:`400`))

#### OUTLIER DETECTION 3: VISUAL INSPECTION OF CR ####

inspect_CR <- function(df) {
  
    # Calculate mean across the trimmed spectrum
    CR_trim_mean <- CR_trim |>
      rowwise() |>
      mutate(mean = mean(c_across(`401`:`2500`)), .before = `401`)
    
    # Identify outliers 
      outliers_above <- CR_trim_mean |>
        filter(mean > 0.31)
      outliers_below <- CR_trim_mean |>
        filter(mean < 0.24)
    
    outliers <- rbind(outliers_above, outliers_below)
    
    # Plot the spectra for the specific type
    plot(as_spectra(subset(df, select = `401`:`2500`)),
         main = paste0(ss)
    )
   
    # If outliers are detected, plot them in red
    if (nrow(outliers) > 0) {
      plot(as_spectra(subset(outliers, select = `401`:`2500`)),
           col = "red", add = TRUE
      )
      legend("topright", legend = outliers$planting_location, title = "Outliers")

      print(paste0("You have potential outliers for ", ss))
      print(unique(outliers$planting_location))


   }
  
  # Return modified dataframe
  return(df)
}

inspect_CR(CR_trim)


#   # Set outliers to NA in the main df
#   df <- df |>
#     mutate(across(
#       `350`:`2500`,
#       ~ ifelse(planting_location %in% outliers$planting_location, NA, .)
#     ))
# } else {
#   print(paste0("No outliers for ", ss, " ", scan_type))
# }




### NEXT TIME ####
# Uncertainties: https://www.sciencedirect.com/science/article/pii/S0034425721003217
# 
# 
# The relative measurement uncertainty was defined as the absolute uncertainty divided by 
# the mean reflectance of the material (or across the dataset, for leaf measurements) with units 
# of per cent (Eq. (4)).
# 
# Absolute measurement uncertainty:(3)
# 
# where UR,abs is the absolute uncertainty associated with the target reflectance.
# R is the spectral reflectance of the target.xi are the readings (leaf clip; Rw, Tw, Rb, Tb; 
# integrating sphere: Ir, Is, Id).Uxi is the standard uncertainty associated with the reading xi.
# STDxi is the standard deviation among scans of each reading xi.N is the number of scans per reading.
# 
# Relative measurement uncertainty:(4)
# 
# where UR,rel is the relative uncertainty associated with the target reflectance.UR,abs is the absolute uncertainty associated with the target reflectance.R is the spectral reflectance of the target.
# 
# The standard uncertainty of each reading corresponds to the probability distributions associated with all different sources of uncertainty, including the instrument characteristics and experimental conditions. The contribution of individual sources of uncertainty was not considered in our uncertainty calculation, except for the ambient temperature (see 2.6).

# Step 1: Calculate standard uncertainty for each reading Uxi
# with Uxi = STDxi/rad(N), 
# where STDxi is the standard deviation among scans of each reading 
# and N is the number of scans per reading



# enter here the metadata columns for grouping WITHOUT type and sample name
metadata_cols <- colnames(lof_outliers_rm[1:12])

# calculate means of each measurement type for each plant (i.e., mean of 5 scans)

std_by_type <- lof_outliers_rm |>
  # pivot data frame so each measurement has its own row
  pivot_longer(
    cols = `350`:`2500`,
    names_to = "nm",
    values_to = "val"
  ) |>
  # make nm as numeric so it will be ordered in descending
  mutate(nm = as.numeric(nm)) |>
  # group by all metadata columns
  group_by(across(c(metadata_cols, type, nm))) |>
  # add n_scans = number of scans per reading kept after outliers and calculate sd and standard uncertainty (uxi)
  summarise(n_scans = n(),
              mean = mean(val, na.rm = TRUE),
              sd = sd(val, na.rm = TRUE),
            uxi = sd/n_scans,
  ) |>
  # pivot again so we have one row for each nm and can calculate AU per nm
  pivot_wider(
  names_from = type,
  values_from = c(mean, sd, uxi)
  ) |>
  group_by(across(c(metadata_cols, nm))) |>
  # calculate CR
  mutate(CR = (mean_BRL * mean_WR - mean_WRL * mean_BR) / (mean_WR - mean_BR)) |>
  # calculate AU
  mutate(AU = (sqrt(CR/mean_BRL) * uxi_BRL^2) + 
  (sqrt(CR/mean_BR) * uxi_BR^2) + 
  (sqrt(CR/mean_WRL) * uxi_WRL^2) + 
  (sqrt(CR/mean_WR) * uxi_WR^2)) |>
  # select 
  select(c(metadata_cols, nm, CR, AU)) |>
  ungroup() |>
  # replace NaN with NA
  mutate_all(~ifelse(is.nan(.), NA, .))

CR <- std_by_type |>
  select(c(metadata_cols, nm, CR)) |>
  # pivot
  pivot_wider(
    names_from = nm,
    values_from = CR
  )
  
AU <- std_by_type |>
  select(c(metadata_cols, nm, AU)) |>
  pivot_wider(
    names_from = nm,
    values_from = AU
  )






# generate bib file
knitr::write_bib(c(.packages()), "temp/packages.bib")
