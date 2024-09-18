inspect_CR <- function(df1, df2) {
  
  # Calculate mean across the trimmed spectrum
  CR_trim_mean <- df1 |>
    mutate(mean = rowMeans(pick(`401`:`2500`)), .before = `401`)
  
  # Identify potential outliers
  # calculate lower and upper bounds using IQR
  lower_bound <- quantile(CR_trim_mean$mean, 0.25) - 1.5 * IQR(CR_trim_mean$mean)
  upper_bound <- quantile(CR_trim_mean$mean, 0.75) + 1.5 * IQR(CR_trim_mean$mean)
  
  iqr_outliers <- CR_trim_mean |>
    filter(mean < lower_bound | mean > upper_bound) |>
    mutate(outlier_type = "IQR") |>
    relocate(outlier_type) |>
    select(-mean)
  
  # identify outliers where AU is too high
  au_outliers <- AU_trim |>
    filter(any(c_across(`401`:`2500`) > 0.01)) |>
    mutate(outlier_type = "AU") |>
    relocate(outlier_type)
  
  
  outliers <- rbind(iqr_outliers, au_outliers)
  
  # Plot the spectra for the specific type
  plot(as_spectra(subset(CR_trim_mean, select = `401`:`2500`)),
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
  
  # Return outliers dataframe
  return(outliers)
}