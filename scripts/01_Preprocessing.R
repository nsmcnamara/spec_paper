### This script is part of the Leaf Optical Properties Paper from the ACORN project.
### This script imports the raw data of spectral measurements, adds metadata,
### checks for outliers, calculates uncertainties and reflectance.
### Established 2024-08-15
### Last Update 2024-09-24
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
library(data.table)

# check and resolve conflicts
conflict_scout()
conflicts_prefer(
  dplyr::filter
)

# load functions
source("source/inspect_by_type.R")
source("source/detect_lof_outliers.R")
source("source/inspect_CR.R")

## Define Variables
# define directory paths
res.path1 <- "/output/"
dat.path1 <- "/data/raw/"

# select folder w/ data subset
ss <- "AT_pubescens_2"
ss <- "AT_robur_2"
ss <- "CH_pubescens_2"
ss <- "CH_robur_2"
ss <- "CH_robur_3"

# select species
if (str_detect(ss, "pubescens")) {
  species <- "Q.pubescens"
} else if (str_detect(ss, "robur")) {
  species <- "Q.robur"
} else {
  print("Could not detect species")
}

# define number of leaf scans per scan type (e.g. 5 white refs --> 5)
n_leaf_scans <- 5

#### DATA IMPORT ####
# import metadata
metadata <- read_csv("data/raw/metadata.csv")

# base path
base_path <- paste0(getwd(), dat.path1, ss)

# get file paths for data subset
file_paths <- list.dirs(base_path, recursive = FALSE)

# get file paths for metadata for each data subset
meta_files <- list.files(base_path, pattern = ".csv", include.dirs = FALSE)

# Create data frame to store data in
spec_df <- data.frame()

# import spectral measurement and merge with metadata
# Loop through each folder within the file paths
for (i in seq_along(file_paths)) {
  # Load spectral measurements
  spectra <- read_spectra(path = file_paths[i])

  # Load metadata
  meta <- read.csv(paste0(base_path, "/", meta_files[i]))  

  # Create a table for merging spectra with measurement type
  mmt <- data.frame(
    type = rep(c("WR", "WRL", "BR", "BRL"), each = n_leaf_scans, times = nrow(meta)),
    planting_location = rep(meta$planting_location, each = 4*n_leaf_scans)
  )

  # Merge with spectra
  spectra_mmt <- cbind(mmt, spectra)

  # Merge with metadata, make sure species matches, because planting location in Austria is not unique
  spectra_mmt_meta <- right_join(metadata[metadata$species == species, ], 
                                 spectra_mmt, by = "planting_location")
  
  # Store data in df
  spec_df <- rbind(spec_df, spectra_mmt_meta) # Append the new data
}


#### OUTLIER DETECTION 1: VISUAL INSPECTION ####
# This function plots the spectra by scan type.
# It also calculates the means for each scan type for certain wavelengths
# and colors these red as an additional warning.

# open png device and specify file name 
png(filename = paste0(getwd(), res.path1, species, "_", ss, "_", "outliers_vis_1.png"))

# print 4 plots in window so all scan types can be seen together
par(mfrow = c(2,2))

# call function and store the outliers in a data frame
outliers_vis_1 <- inspect_by_type(spec_df)

# close device
dev.off()

# manual check: plot by scan type, e.g. "BRL"
# tst <- spec_df |>
#   filter(type == "BRL")
# plot(as_spectra(subset(tst, select = `350`:`2500`)),
#      main = paste0(ss, " BRL ")
# )

# manual check: plot where a specific value is above a threshold, e.g. at 1000 nm is above 0.55
# tst2 <- tst |>
#   filter(`1000` > 0.55)
# plot(as_spectra(subset(tst2, select = `350`:`2500`)),
#      col = "red", add = TRUE
# )

# manual addition of outlier to outlier df: e.g. for planting location "8_E_3")
# outliers_vis_1_man <- spec_df |>
#   filter(planting_location == "8_E_3") |>
#   mutate(outlier_type = "vis_man", .before = 1)
# outliers_vis_1 <- rbind(outliers_vis_1_man, outliers_vis_1)


#### OUTLIER DETECTION 2: LOF ####
# This function is based on this paper: https://dl.acm.org/doi/pdf/10.1145/335191.335388
# It calculates the local outlier factor (LOF), i.e. if something is an outlier given its local neighbourhood.

# k: The kth-distance to be used to calculate the LOFs.
# if k is set to the number of scans per scan type,
# this identifies if one or more scans are abnormally far apart from each other
# irrespective of whether a single plant is just different from other plants.
# i.e. it aims to identify measurement error.
# LOF threshold 2 is somewhat arbitrary but reasonable based on the paper.


# open png device and specify file name
png(filename = paste0(getwd(), res.path1, species, "_", ss, "_", "outliers_lof.png"))

# print 4 plots in window so all scan types can be seen together
par(mfrow = c(2,2))

# Call the function
outliers_lof <- detect_lof_outliers(spec_df, lof_threshold = 2, k = n_leaf_scans)

# close device
dev.off()

# for AT_robur: remove D_11
# spec_df <- subset(spec_df, !(sample_name %in% outliers_lof$sample_name))
# for CH_pubescens: remove 9_C_2
# spec_df <- subset(spec_df, !(sample_name %in% outliers_lof$sample_name))
# for CH robur 3: remove 27_G_3 all, 38_G_2 single scan
# spec_df <- subset(spec_df, !(sample_name %in% outliers_lof$sample_name))
# spec_df <- spec_df |>
#   filter(planting_location != "27_G_3")

#### CALCULATE REFLECTANCE AND ABSOLUTE UNCERTAINTY ####
# Calculation of Reflectance follows this paper: https://doi.org/10.1080/01431169208904118
# Calculation of Absolute Uncertainty follows this paper: https://doi.org/10.1016/j.rse.2021.112601

# enter here the metadata columns for grouping WITHOUT type and sample name
metadata_cols <- colnames(spec_df[1:12])

# calculate reflectance and absolute uncertainty
CR_AU_comb <- spec_df |>
  # pivot data frame so each measurement has its own row
  pivot_longer(
    cols = `350`:`2500`,
    names_to = "nm",
    values_to = "val"
  ) |>
  # make nm as numeric so it will be ordered in descending
  mutate(nm = as.numeric(nm)) |>
  # group by all metadata columns
  group_by(across(c((all_of(metadata_cols)), type, nm))) |>
  # add n_scans = number of scans per reading and calculate sd and standard uncertainty (uxi)
  summarise(
    n_scans = n(),
    mean = mean(val, na.rm = TRUE),
    sd = sd(val, na.rm = TRUE),
    uxi = sd / sqrt(n_scans),
  ) |>
  # drop n_scans: if we remove a scan based on outlier analysis, it will produce different values and result in NA
  select(-n_scans) |>
  # pivot again so we have one row for each nm and can calculate AU per nm
  pivot_wider(
    names_from = type,
    values_from = c(mean, sd, uxi)
  ) |>
  group_by(across(c(metadata_cols, nm))) |>
  # calculate Reflectanc
  mutate(CR = (mean_BRL * mean_WR - mean_WRL * mean_BR) / (mean_WR - mean_BR)) |>
  # calculate Absolute Uncertainty
  mutate(AU = sqrt(
    (mean_BR * (mean_WRL - mean_BRL) / ((mean_WR-mean_BR)^2))^2 * uxi_WR^2 +
      (mean_BR/(mean_WR-mean_BR))^2 * uxi_WRL^2 +
      (mean_WR * (mean_WRL-mean_BRL)/((mean_WR-mean_BR)^2))^2 * uxi_BR^2 +
      (mean_WR / (mean_WR-mean_BR))^2 * uxi_BRL^2 
    
  ))  |>
  # select relevant columns
  select(c(metadata_cols, nm, CR, AU)) |>
  ungroup() |>
  # replace NaN with NA
  mutate_all(~ ifelse(is.nan(.), NA, .))

# make df w Reflectance
CR <- CR_AU_comb |>
  select(c(metadata_cols, nm, CR)) |>
  # pivot
  pivot_wider(
    names_from = nm,
    values_from = CR
  )

# make df w Absolute Uncertainties
AU <- CR_AU_comb |>
  select(c(metadata_cols, nm, AU)) |>
  pivot_wider(
    names_from = nm,
    values_from = AU
  )

# Plot
# specify file name
png(filename = paste0(getwd(), res.path1, species, "_", ss, "_", "CR_AU.png"))

# use entire window to plot
par(mfrow = c(1,1))

# add margins
par(mar = c(5, 4, 4, 6))

# plot reflectance
plot(as_spectra(subset(CR, select = `350`:`2500`)),
     main = paste0(ss, " untrimmed"),
     ylab = "Reflectance [no unit]",
     xlab = "Wavelength [nm]",
     las = 1,
     ylim = c(0,1))

# add line where data will be trimmed
abline(v = 400, col = "blue", lty = 3, lwd = 2)
text(x = 500, y = 0.9, labels = "trim", col = "blue")

# plot absolute uncertainty
par(new = TRUE)
plot(as_spectra(subset(AU, select = `350`:`2500`)),
     ylim = c(0, 0.1),
     col = "red", axes = FALSE, xlab = "", ylab = "")

axis(side = 4, col = "red", col.axis = "red",
     at = seq(0, 0.1, by = 0.02),
     labels = format(seq(0, 0.1, by = 0.02), scientific = FALSE),
     las = 1)

mtext("Absolute Uncertainty [no unit]", side = 4, line = 4, col = "red")

dev.off()


# trim the noisy section of 350 - 400 nm

CR_trim <- CR |>
  select(-c(`350`:`400`)) |>
  mutate(feat = "CR", .before = `401`)

AU_trim <- AU |>
  select(-c(`350`:`400`)) |>
  mutate(feat = "AU", .before = `401`)


# save trimmed data
write_csv(CR_trim, 
          file = paste0(getwd(), "/data/processed/", species, "_", ss, "_CR_trimmed.csv"))
write_csv(AU_trim, 
          file = paste0(getwd(), "/data/processed/", species, "_", ss, "_AU_trimmed.csv"))


#### OUTLIER DETECTION 3: VISUAL INSPECTION OF CR ####
# Plot
# specify file name
png(filename = paste0(getwd(), res.path1, species, "_", ss, "_", "vis_2.png"))

# use entire window to plot
par(mfrow = c(1,1))

# call function
outliers_vis_2 <- inspect_CR(CR_trim, AU_trim)

dev.off()


#### REMOVING OUTLIERS IF NECESSARY ####
# combine outlier dfs
all_outliers <- bind_rows(
  outliers_vis_1 |> select(c("outlier_type", metadata_cols)),
  outliers_lof |> select(c("outlier_type", metadata_cols)),
  outliers_vis_2 |> select(c("outlier_type", metadata_cols))
)

# for AT pubescens: B20, additionally inspect E_10, H_29: both seem ok, do not remove
# all_outliers <- all_outliers |>
#   filter(planting_location == "B_20")

# for AT robur: remove single scan for D_11, additionally inspect C_7, D_26, E_10: all seem ok, do not remove
all_outliers <- all_outliers |>
  filter(planting_location == "D_11")

# for CH pubescens: remove single scan for 9_C_2, additionally inspect 10_D_1, 3_A_3, 8_E_3: remove 10_D_1, 8_E_3
all_outliers <- all_outliers |>
  filter(planting_location == "10_D_1" | planting_location == "8_E_3" | planting_location == "9_C_2")

# for CH robur 2: remove single scans for 19_B_2, 21_I_2

# for CH robur 3: remove single scans for 38_G_2, all for 27_G_3, 27_F_1
all_outliers <- all_outliers |>
  filter(planting_location == "38_G_2" | planting_location == "27_G_3" | planting_location == "27_F_1")

plot(as_spectra(subset(CR_trim, select = `401`:`2500`)),
     main = paste0(ss)
)
plot(as_spectra(subset(CR_trim, select = `401`:`2500`, planting_location == "27_F_1")),
     col = "red", add = TRUE
)



# save combined outlier df
write_csv(all_outliers, 
          file = paste0(getwd(), "/data/processed/", species, "_", ss, "_outliers.csv"))

# remove outliers from data
all_outliers <- all_outliers |>
  filter(planting_location != "38_G_2")

CR_outliers_rm <- anti_join(CR_trim, all_outliers, by = "planting_location")
AU_outliers_rm <- anti_join(AU_trim, all_outliers, by = "planting_location")



# save new and processed data
write_csv(CR_outliers_rm, 
          file = paste0(getwd(), "/data/processed/", species, "_", ss, "_CR_outliers_rm.csv"))

write_csv(AU_outliers_rm, 
          file = paste0(getwd(), "/data/processed/", species, "_", ss, "_AU_outliers_rm.csv"))

### Combining Species ###
at_pub <- read_csv("data/processed/Q.pubescens_AT_pubescens_2_CR_outliers_rm.csv")
ch_pub <- read_csv("data/processed/Q.pubescens_CH_pubescens_2_CR_outliers_rm.csv")

at_rob <- read_csv("data/processed/Q.robur_AT_robur_2_CR_outliers_rm.csv")
ch_rob <- read_csv("data/processed/Q.robur_CH_robur_2_CR_outliers_rm.csv")

pub_spectra <- rbind(at_pub, ch_pub)
write_csv(pub_spectra, 
          file = paste0(getwd(), "/data/processed/", "pub_spectra.csv"))

rob_spectra <- rbind(at_rob, ch_rob)
write_csv(rob_spectra, 
          file = paste0(getwd(), "/data/processed/", "rob_spectra.csv"))


# generate bib file
knitr::write_bib(c(.packages()), "temp/packages.bib")

