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