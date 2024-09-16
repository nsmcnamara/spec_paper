#### OUTLIER DETECTION 1: VISUAL INSPECTION ####
# This function plots the spectra by scan type.
# It also calculates the means for each scan type for certain wavelengths
# and colors these red as an additional warning. 


inspect_by_type <- function(df) {
  # Define scan types
  scan_types <- levels(as.factor(df$type))
  
  # Define the conditions for each scan type to plot in red and throw warning
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
    
    
    # Calculate mean for the specified range
    mean_by_type <- df |>
      filter(type == scan_type) |>
      mutate(mean = rowMeans(pick(all_of(mean_range_min):all_of(mean_range_max))), .before = `350`)
    
    # Identify outliers based on the conditions
    if (operator == "greater") {
      outliers <- mean_by_type |>
        filter(mean > threshold)
    } else {
      outliers <- mean_by_type |>
        filter(mean < threshold)
    }
    
    # Plot the spectra for the specific type
    plot(as_spectra(subset(mean_by_type, select = `350`:`2500`)),
         main = paste0(ss, " ", scan_type)
    )
    # If outliers are detected, plot them in red
    if (nrow(outliers) > 0) {
      plot(as_spectra(subset(outliers, select = `350`:`2500`)),
           col = "red", add = TRUE
      )
      legend("topright", legend = outliers$planting_location, title = "Outliers")
      
      print(paste0("You have potential outliers for ", ss, " ", scan_type))
      print(unique(outliers$planting_location))
      
    } else {
      print(paste0("No outliers for ", ss, " ", scan_type))
    }
  }
  
  # Return modified dataframe
  tot_outliers <- data.frame()
  tot_outliers <- rbind(tot_outliers, outliers)
  return(tot_outliers)
}
