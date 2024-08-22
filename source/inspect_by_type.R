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
