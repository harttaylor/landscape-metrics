# Load required packages
library(sf)
library(landscapemetrics)
library(tidyverse)
library(raster)
library(fasterize)
library(tools)

# Function to calculate metrics for a single buffer (keeping original function)
calculate_buffer_metrics <- function(habitat_polygons, point, buffer_radius) {
  # Create buffer
  point_buffer <- st_buffer(point, buffer_radius)
  
  # Clip habitat to buffer
  habitat_in_buffer <- st_intersection(habitat_polygons, point_buffer)
  
  # Check if there's any habitat in the buffer
  if (nrow(habitat_in_buffer) == 0) {
    return(list(
      prop_habitat = 0,
      edge_density = 0,
      aggregation_index = 0,
      cohesion = 0,
      n_patches = 0,
      mean_patch_area = 0,
      enn_mn = NA,
      shape_mn = NA,
      frac_mn = NA
    ))
  }
  
  # Calculate basic habitat amount
  buffer_area <- pi * buffer_radius^2
  habitat_area <- sum(st_area(habitat_in_buffer))
  prop_habitat <- as.numeric(habitat_area / buffer_area)
  
  # Create raster for landscapemetrics
  habitat_in_buffer$value <- 1
  
  # Create empty raster based on buffer extent
  bbox <- st_bbox(point_buffer)
  r <- raster(
    xmn = bbox["xmin"], xmx = bbox["xmax"],
    ymn = bbox["ymin"], ymx = bbox["ymax"],
    res = 5,  # 10m resolution
    crs = st_crs(habitat_polygons)
  )
  
  # Rasterize and calculate metrics
  tryCatch({
    r <- fasterize(habitat_in_buffer, r, field = "value", background = 0)
    
    metrics <- list(
      prop_habitat = prop_habitat,
      edge_density = lsm_l_ed(r)$value,
      aggregation_index = lsm_l_ai(r)$value,
      cohesion = lsm_l_cohesion(r)$value,
      n_patches = lsm_l_np(r)$value,
      mean_patch_area = lsm_l_area_mn(r)$value,
      enn_mn = lsm_l_enn_mn(r)$value,
      shape_mn = lsm_l_shape_mn(r)$value,
      frac_mn = lsm_l_frac_mn(r)$value
    )
    return(metrics)
    
  }, error = function(e) {
    warning(paste("Error calculating metrics for point buffer:", e$message))
    return(list(
      prop_habitat = prop_habitat,
      edge_density = NA,
      aggregation_index = NA,
      cohesion = NA,
      n_patches = NA,
      mean_patch_area = NA,
      enn_mn = NA,
      shape_mn = NA,
      frac_mn = NA
    ))
  })
}

# Modified main processing function for hypothesis 1
# Modified main processing function for hypothesis 1
process_landscape_metrics <- function(
    input_dir,           # Directory containing input files
    output_dir,          # Directory for output files
    buffer_sizes = c(150, 500, 1000),
    raster_resolution = 5  # Default 5m resolution
) {
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get list of all hypothesis 1 habitat shapefiles
  habitat_files <- list.files(input_dir, 
                              pattern = "hyp1.*\\.shp$", # You will have to adjust this for each hypo
                              full.names = TRUE)
  
  # Initialize master results dataframe
  all_results <- data.frame()
  
  # Process each habitat file (year)
  for(habitat_file in habitat_files) {
    # Extract year from filename
    year <- as.numeric(stringr::str_extract(basename(habitat_file), "\\d{4}"))
    
    if(is.na(year)) {
      warning(sprintf("Could not extract year from filename: %s", basename(habitat_file)))
      next
    }
    
    # Read habitat data
    habitat <- st_read(habitat_file, quiet = TRUE)
    
    # Look for corresponding points CSV
    points_csv <- file.path(input_dir, 
                            paste0("patch_metrics_hyp1_", year, ".csv"))
    
    if(!file.exists(points_csv)) {
      warning(sprintf("No matching points CSV found for year %d", year))
      next
    }
    
    # Read points data
    points_data <- read.csv(points_csv)
    
    # Create SF object from points data using Easting/Northing
    points_sf <- st_as_sf(points_data, 
                          coords = c("Easting", "Northing"),
                          crs = st_crs(habitat))
    
    # Process each point
    for (i in 1:nrow(points_sf)) {
      point <- points_sf[i,]
      cat(sprintf("Processing year %d, point %d of %d\n", 
                  year, i, nrow(points_sf)))
      
      # Process each buffer size
      for (size in buffer_sizes) {
        cat(sprintf("  Buffer size: %d m\n", size))
        
        metrics <- calculate_buffer_metrics(habitat, point, size)
        
        # Create a row for the results, including relevant metadata
        row_data <- data.frame(
          year = year,
          site_id = points_data$site_id[i],
          location_n = points_data$location_n[i],
          survey_year = points_data$survey_yea[i],
          veg_type = points_data$Veg_Type[i],
          origin_year = points_data$Origin_Yea[i],
          moisture_regime = points_data$Moisture_R[i],
          buffer_size = size,
          as.data.frame(t(unlist(metrics)))
        )
        
        # Add to results
        all_results <- rbind(all_results, row_data)
      }
    }
  }
  
  # Calculate correlations between metrics for each buffer size
  metric_cols <- c("prop_habitat", "edge_density", "aggregation_index", 
                   "cohesion", "n_patches", "mean_patch_area", 
                   "enn_mn", "shape_mn", "frac_mn")
  
  # Create a list to store correlation matrices
  cor_matrices <- list()
  
  # Open correlation log file
  cor_log <- file(file.path(output_dir, "hyp1_correlation_warnings.txt"))
  writeLines("Correlation Analysis Results:", cor_log)
  
  # Calculate correlations for each buffer size
  for(buffer in unique(all_results$buffer_size)) {
    # Subset data for this buffer size
    buffer_data <- all_results[all_results$buffer_size == buffer, ]
    
    # Calculate correlation matrix
    cor_matrix <- cor(buffer_data[, metric_cols], 
                      use = "pairwise.complete.obs")
    
    # Store correlation matrix
    cor_matrices[[paste0("buffer_", buffer)]] <- cor_matrix
    
    # Save individual correlation matrix
    write.csv(cor_matrix, 
              file.path(output_dir, 
                        sprintf("hyp1_metric_correlations_buffer_%dm.csv", buffer)))
    
    # Log highly correlated pairs for this buffer size
    writeLines(sprintf("\nBuffer size: %d m", buffer), cor_log)
    
    high_cors <- which(abs(cor_matrix) > 0.7 & cor_matrix < 1, 
                       arr.ind = TRUE)
    
    if (nrow(high_cors) > 0) {
      writeLines("Highly correlated metric pairs (r > 0.7):", cor_log)
      
      for (i in 1:nrow(high_cors)) {
        writeLines(sprintf("%s and %s: r = %.3f",
                           rownames(cor_matrix)[high_cors[i,1]],
                           colnames(cor_matrix)[high_cors[i,2]],
                           cor_matrix[high_cors[i,1], high_cors[i,2]]),
                   cor_log)
      }
    } else {
      writeLines("No highly correlated pairs found.", cor_log)
    }
  }
  
  # Calculate overall correlations across all buffer sizes
  cor_matrix_all <- cor(all_results[, metric_cols], 
                        use = "pairwise.complete.obs")
  
  # Save overall correlation matrix
  write.csv(cor_matrix_all, 
            file.path(output_dir, "hyp1_metric_correlations_all_buffers.csv"))
  
  # Log overall correlations
  writeLines("\nOverall correlations (across all buffer sizes):", cor_log)
  high_cors_all <- which(abs(cor_matrix_all) > 0.7 & cor_matrix_all < 1, 
                         arr.ind = TRUE)
  
  if (nrow(high_cors_all) > 0) {
    for (i in 1:nrow(high_cors_all)) {
      writeLines(sprintf("%s and %s: r = %.3f",
                         rownames(cor_matrix_all)[high_cors_all[i,1]],
                         colnames(cor_matrix_all)[high_cors_all[i,2]],
                         cor_matrix_all[high_cors_all[i,1], high_cors_all[i,2]]),
                 cor_log)
    }
  } else {
    writeLines("No highly correlated pairs found.", cor_log)
  }
  
  close(cor_log)

# Save results
write.csv(all_results, 
          file.path(output_dir, "hyp1_landscape_metrics_all_years.csv"), 
          row.names = FALSE)

return(all_results)
}

# Usage
process_landscape_metrics(
  input_dir = "C:/Users/hartt/Documents/CAWA Habitat Analysis/2_GIS/Outputs",
  output_dir = "C:/Users/hartt/Documents/CAWA Habitat Analysis/0_data/covariates"
  )

