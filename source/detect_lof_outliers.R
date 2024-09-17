detect_lof_outliers <- function(df, lof_threshold = 2, k = n_leaf_scans) {
  # Define types
  scan_types <- levels(as.factor(df$type))
  
  # Define df for returning outliers
  tot_outliers <- data.frame()
  
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
    
    # plot
    plot(as_spectra(subset(spec_by_type, select = `350`:`2500`)),
         main = paste0(ss, " ", scan_type)
    )
    
    if (nrow(lof_outliers) > 0) {
      print("You have LOF outliers:")
      print(paste0(lof_outliers$planting_location))
      
      # plot outliers
      plot(as_spectra(subset(lof_outliers, select = `350`:`2500`)),
           col = "red", add = TRUE
      )
      legend(x = "topright", legend = lof_outliers$planting_location, title = "LOF Outliers")
      
      # Add outliers to df to be returned
      tot_outliers <- rbind(tot_outliers, df[df$sample_name %in% lof_outliers$sample_name, ])
      
    } else {
      print("You have no LOF outliers.")
    }
    
  }
  
  # add column to identify which outlier type it is
  tot_outliers <- tot_outliers |>
    mutate(outlier_type = "lof") |>
    relocate(outlier_type)
  
  return(tot_outliers)
}