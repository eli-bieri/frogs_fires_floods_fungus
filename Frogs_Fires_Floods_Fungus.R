
# This script divides australia into 10km grid cells, and classifies
# cells as either burnt or unburnt and flooded or unflooded based on fire 
# and flood extent data. Then it assigns frogid records to grid cells and 
# classifies records according to spatial and temporal parameters: Fire - 
# (prefire_burnt, prefire_unburnt, postfire_burnt, postfire_unburnt) and 
# flood - (preflood_flooded, preflood_unflooded, postflood_flooded, postflood_unflooded)
# unburnt and unflooded grid cells represent controls, making this a
# Before-After-Control-Impact (BACI) study design. The script compares species occurences
# and species richness between the categories

library(tidyverse)
library(magrittr)
library(lubridate)
library(sf)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(rgeos)
library(ggpubr)
library(readr)
library(tidyr)
library(ozmaps)
library(scico)
library(gam)
library(mgcv)
library(yaml)
library(ggpubr)
library(ggh4x)
library(viridis)
library(writexl)
library(tidyverse)
library(gridExtra)
library(lwgeom)

clean_data <- readRDS("FrogID_clean_data.RDS")

#Distribution change analysis bushfire

clean_data <- clean_data %>% filter(state == 'New South Wales')

# read in australia outline
aus <- rnaturalearth::ne_countries(country="Australia", returnclass="sf")
#st_crs(aus)
#make it the right data projection
aus.2 <- aus %>%
  st_transform(32754)
#make the grids
grids_10 <- st_make_grid(aus.2, cellsize=10000) %>%
  st_transform(crs=st_crs(aus)) %>%
  st_as_sf() %>%
  rename(geometry=x) %>%
  mutate(grid_id=1:nrow(.))
grids_10_2<- aus %>%
  st_intersects(grids_10) %>%
  as.data.frame()
aus_grids <- grids_10 %>%
  dplyr::filter(grid_id %in% grids_10_2$col.id) %>%
  mutate(grid_id=1:nrow(.)) %>%
  mutate(lon=st_coordinates(st_centroid(.$geometry))[,1]) %>%
  mutate(lat=st_coordinates(st_centroid(.$geometry))[,2]) %>%
  dplyr::select(grid_id, lon, lat, geometry)
#saved so can just load this in the future aus_grids <- readRDS("grids_10km.RDS")
saveRDS(aus_grids, "grids_10km.RDS")
#assign data to FrogID
grids_metadata <- data.frame(grid_id=aus_grids[[1]]) %>%
  mutate(col.id=1:nrow(.))
# Define the date ranges for the prefire and postfire seasons
prefire_range <- as.Date(c("2018-08-01", "2019-05-01"))
postfire_range <- as.Date(c("2020-08-01", "2021-05-01"))
# Filter the data for prefire and postfire periods
prefire_data <- clean_data %>% filter(date >= prefire_range[1] & date <= prefire_range[2])
postfire_data <- clean_data %>% filter(date >= postfire_range[1] & date <= postfire_range[2])
# Convert frog occurrences to spatial objects and assign them to grid cells
prefire_data_sf <- st_as_sf(prefire_data, coords = c("lng", "lat"), crs = 4326)
prefire_data_sf <- st_set_crs(prefire_data_sf, 4326) # Set the CRS explicitly
prefire_data_sf <- st_transform(prefire_data_sf, st_crs(aus_grids))
prefire_data_grid <- st_join(prefire_data_sf, aus_grids, left = FALSE)
postfire_data_sf <- st_as_sf(postfire_data, coords = c("lng", "lat"), crs = 4326)
postfire_data_sf <- st_set_crs(postfire_data_sf, 4326) # Set the CRS explicitly
postfire_data_sf <- st_transform(postfire_data_sf, st_crs(aus_grids))
postfire_data_grid <- st_join(postfire_data_sf, aus_grids, left = FALSE)
# Transform prefire_data_grid and postfire_data_grid to the desired CRS (EPSG:32754)
prefire_data_grid_proj <- st_transform(prefire_data_grid, 32754)
postfire_data_grid_proj <- st_transform(postfire_data_grid, 32754)
# Read in the fire layer
fire <- st_read("FESM_area_84.shp")
# Assign CRS to the fire layer
if (is.na(st_crs(fire))) {
  st_crs(fire) <- 4326
}
# Transform fire to the same CRS as aus_grids
fire_transformed <- st_transform(fire, st_crs(aus_grids))
# Classify the grid cells as burnt or unburnt
burnt_cells <- st_intersection(aus_grids, fire_transformed)
unburnt_cells <- st_difference(aus_grids, fire_transformed)
# Filter frog occurrences to only include those within burnt and unburnt grid cells
prefire_burnt <- prefire_data_grid[which(lengths(st_intersects(prefire_data_grid, burnt_cells)) > 0), ]
prefire_unburnt <- prefire_data_grid[which(lengths(st_intersects(prefire_data_grid, unburnt_cells)) > 0), ]
postfire_burnt <- postfire_data_grid[which(lengths(st_intersects(postfire_data_grid, burnt_cells)) > 0), ]
postfire_unburnt <- postfire_data_grid[which(lengths(st_intersects(postfire_data_grid, unburnt_cells)) > 0), ]
# Transform the filtered data back to WGS 84 (EPSG:4326)
prefire_burnt <- st_transform(prefire_burnt, 4326)
prefire_unburnt <- st_transform(prefire_unburnt, 4326)
postfire_burnt <- st_transform(postfire_burnt, 4326)
postfire_unburnt <- st_transform(postfire_unburnt, 4326)
# Randomly subsample prefire and postfire data to have an equal number of observations
min_n <- min(nrow(prefire_burnt), nrow(postfire_burnt), nrow(prefire_unburnt), nrow(postfire_unburnt))
subsampled_prefire_burnt <- prefire_burnt %>% slice_sample(n = min_n)
subsampled_postfire_burnt <- postfire_burnt %>% slice_sample(n = min_n)
subsampled_prefire_unburnt <- prefire_unburnt %>% slice_sample(n = min_n)
subsampled_postfire_unburnt <- postfire_unburnt %>% slice_sample(n = min_n)
# Combine the subsampled data
subsampled_data <- bind_rows(
  subsampled_prefire_burnt %>% mutate(status = "prefire_burnt"),
  subsampled_postfire_burnt %>% mutate(status = "postfire_burnt"),
  subsampled_prefire_unburnt %>% mutate(status = "prefire_unburnt"),
  subsampled_postfire_unburnt %>% mutate(status = "postfire_unburnt")
)

# Flood statistical analysis 

# Get the top 40 well-sampled species
top_species <- subsampled_data %>%
  count(species) %>%
  arrange(desc(n)) %>% # arrange by descending count
  head(40) %>%
  pull(species)

#view sample size for species (to adjust most well-sampled species threshold)
species_counts <- clean_data %>% 
  group_by(species) %>% 
  tally(sort = TRUE)

print(species_counts)


# Filter the data to only include the top 60 well-sampled species
subsampled_data_top_species <- subsampled_data %>%
  filter(species %in% top_species)

# Calculate the presence in grid cells for each species and status (prefire_burnt, postfire_burnt, prefire_unburnt, postfire_unburnt)
species_grid_presence <- subsampled_data_top_species %>%
  group_by(species, status) %>%
  summarize(n_grid_cells = n_distinct(st_coordinates(geometry))) %>%
  ungroup()

  # Print the results
  species_grid_presence
  
  # Statistical analysis
  
  # Create contingency table
  contingency_table <- table(subsampled_data_top_species$species, subsampled_data_top_species$status)
  
  # Chi-squared test
  chi_squared_test <- chisq.test(contingency_table)
  
  # Print the result
  chi_squared_test
  
  # Conduct an ANOVA for all statuses and all species (including controls)
  anova_result <- aov(n_grid_cells ~ status, data = species_grid_presence)
  summary(anova_result)
  
  
  #ANOVA for each species individually
  aov(n_grid_cells ~ species + status, data = species_grid_presence)
  
  model <- aov(n_grid_cells ~ species + status, data = species_grid_presence)
  summary(model)
  

  
# paired difference plot
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Reshape and summarize the data
  reshaped_data <- species_grid_presence %>%
    filter(status %in% c("postfire_burnt", "prefire_burnt")) %>% 
    spread(key = status, value = n_grid_cells, fill = 0) %>%
    group_by(species) %>%
    summarise(postfire_burnt = sum(postfire_burnt, na.rm = TRUE),
              prefire_burnt = sum(prefire_burnt, na.rm = TRUE)) %>%
    mutate(difference = postfire_burnt - prefire_burnt)
  
  # Plot the data
  ggplot(reshaped_data, aes(x = reorder(species, -difference), y = difference)) + 
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(title = "Change in Grid Cells (Postfire - Prefire) in Burnt Areas",
         x = "Species",
         y = "Difference in Grid Cells") +
    theme_minimal()
  
  # Print the species names along with their differences
  print(reshaped_data[, c("species", "difference")])
  
  View(reshaped_data)
  
  
  # Save the differences to a CSV file
  write.csv(reshaped_data[, c("species", "difference")], "species_difference.csv", row.names = FALSE)
  

#ANOVA for burnt cells only 
  # Filter the data to only include burnt grid cells
  burnt_data <- species_grid_presence %>%
    filter(status %in% c("postfire_burnt", "prefire_burnt"))
  
  # Run the ANOVA
  anova_result <- aov(n_grid_cells ~ status, data = burnt_data)
  
  
  # Display the results
  summary(anova_result)
  
#T-Test for differences among individual species
  
  # Filtering and reshaping data
  paired_data <- species_grid_presence %>%
    filter(status %in% c("postfire_burnt", "prefire_burnt")) %>% 
    spread(key = status, value = n_grid_cells, fill = 0)
  
  # Initializing the results list
  results_list <- list()
  
  # Running paired t-tests
  for(spec in unique(paired_data$species)) {
    spec_data <- paired_data[paired_data$species == spec, ]
    test_result <- t.test(spec_data$postfire_burnt, spec_data$prefire_burnt, paired = TRUE)
    results_list[[spec]] <- test_result
  }
  
  # Print results
  results_list
  

#Species richness analysis
  
  library(dplyr)
  
  # Determine the minimum number of observations
  min_samples <- min(nrow(prefire_unburnt), nrow(prefire_burnt), nrow(postfire_burnt), nrow(postfire_unburnt))
  
  # Define the function to compute species richness
  compute_avg_species_richness <- function(data) {
    return(mean(data$number_of_sp_calling, na.rm = TRUE))
  }
  
  num_iterations <- 100  # Number of bootstrap iterations
  bootstrap_results <- data.frame(category=character(), avg_richness=numeric(), iteration=integer())
  
  for(i in 1:num_iterations) {
    # Randomly subsample the datasets
    subsampled_prefire_unburnt <- prefire_unburnt %>% slice_sample(n = min_samples)
    subsampled_prefire_burnt <- prefire_burnt %>% slice_sample(n = min_samples)
    subsampled_postfire_burnt <- postfire_burnt %>% slice_sample(n = min_samples)
    subsampled_postfire_unburnt <- postfire_unburnt %>% slice_sample(n = min_samples)
    
    categories <- list(prefire_unburnt = subsampled_prefire_unburnt, 
                       prefire_burnt = subsampled_prefire_burnt,
                       postfire_burnt = subsampled_postfire_burnt,
                       postfire_unburnt = subsampled_postfire_unburnt)
    
    for(cat_name in names(categories)) {
      avg_richness <- compute_avg_species_richness(categories[[cat_name]])
      bootstrap_results <- rbind(bootstrap_results, data.frame(category=cat_name, avg_richness=avg_richness, iteration=i))
    }
  }
  
  # Average richness across bootstrap iterations for each category
  final_avg_richness <- bootstrap_results %>%
    group_by(category) %>%
    summarise(avg_richness = mean(avg_richness))
  
  print(final_avg_richness)
  
  # ANOVA using bootstrapped data
  anova_result <- aov(avg_richness ~ category, data = bootstrap_results)
  print(summary(anova_result))
  
  post_hoc <- TukeyHSD(anova_result)
  print(post_hoc)

   #bar plot with error bars
  SE_data <- bootstrap_results %>%
    group_by(category) %>%
    summarise(
      mean_richness = mean(avg_richness),
      SE = sd(avg_richness) / sqrt(n())
    )
  
  print(SE_data)
  

  ggplot(SE_data, aes(x = category, y = mean_richness, fill = category)) + 
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_richness - SE, ymax = mean_richness + SE), width = 0.25) +
    labs(y = "Average Number of Species Calling per Recording", x = "Status") +
    theme_minimal() + 
    theme(legend.position = "none")
  
  library(ggplot2)
  
  # Custom color palette
  custom_colors <- c("postfire_unburnt" = scales::alpha("palegreen", 0.7),
                     "postfire_burnt" = scales::alpha("salmon", 0.7),
                     "prefire_unburnt" = scales::alpha("palegreen", 0.7),
                     "prefire_burnt" = scales::alpha("salmon", 0.7))
  
  ggplot(SE_data, aes(x = category, y = mean_richness, fill = category)) + 
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_richness - SE, ymax = mean_richness + SE), width = 0.2, size = .5) +
    scale_fill_manual(values = custom_colors) +
    labs(y = "Mean # of Species Calling per Recording", x = "Status") +
    theme_minimal() +
    labs(fill = "Status") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c("Prefire burnt", "Prefire unburnt", "Postfire burnt", "Postfire unburnt"))
  
  
  #Species trait analysis
  
  data <- read.csv("species_trait_matrix.csv")
  colnames(data)
  
  
  # Read the data
  data <- read.csv("species_trait_matrix.csv")
  
  # Replace spaces with underscores in column names
  names(data) <- gsub(" ", "_", names(data))
  
  library(tidyverse)
  library(car)
  
  # Read the data
  data <- read.csv("species_trait_matrix.csv")
  
  library(tidyverse)
  library(car)
  
  # Read the data
  data <- read.csv("species_trait_matrix.csv")
  
  # Convert the categorical variables to factors
  data$Ecological.Group <- as.factor(data$Ecological.Group)
  data$Breeding.strategy <- as.factor(data$Breeding.strategy)
  data$Phylogenetic.Group <- as.factor(data$Phylogenetic.Group)
  data$Size.Group <- as.factor(data$Size.Group)
  data$Color <- as.factor(data$Color)
  data$Conservation.Status <- as.factor(data$Conservation.Status)
  data$Disease.related.declines. <- as.factor(data$Disease.related.declines.)
  
  # Fit the GLM model
  glm_model <- glm(Net.grid.cell.change ~ Ecological.Group + Breeding.strategy + Phylogenetic.Group +
                     Size.Group + Color + Conservation.Status + Dependence.on.fire.adapted.vegetation +
                     Wildfire.Overlap + Disease.related.declines., data = data, family = gaussian(link="identity"))
  
  # Display the summary of the model
  summary(glm_model)
  
  # Check for model assumptions, residuals vs fitted
  plot(glm_model, which=1)
  
  # Check for normality of residuals
  plot(glm_model, which=2)
  
  # Check for influential outliers
  influence_plot(glm_model, main="Influence Plot", sub="Circle size is proportional to Cook's distance")
  
  
  # *** Flood analysis ***
  library(httr)
  library(sf)
  
  # Define the date ranges for the preflood and postflood seasons
  preflood_range_flood <- as.Date(c("2021-08-01", "2022-05-01"))
  postflood_range_flood <- as.Date(c("2023-08-01", "2024-05-01"))
  
  # Filter the data for preflood and postflood periods
  preflood_data_flood <- clean_data %>% filter(date >= preflood_range_flood[1] & date <= preflood_range_flood[2])
  postflood_data_flood <- clean_data %>% filter(date >= postflood_range_flood[1] & date <= postflood_range_flood[2])
  
  # Define the URL for the GeoJSON query and initialize variables for pagination
  query_url_flood <- "https://portal.data.nsw.gov.au/arcgis/rest/services/dataenablement/FloodExtent_2022/FeatureServer/0/query"
  offset_flood <- 0
  record_count_flood <- 15
  all_data_flood <- NULL
  
  repeat {
    # Query parameters for retrieving data in GeoJSON format
    params_flood <- list(
      f = "geojson",
      where = "1=1",
      outFields = "*",
      outSR = "4326",
      resultOffset = offset_flood,
      resultRecordCount = record_count_flood
    )
    
    # GET request
    response_flood <- GET(url = query_url_flood, query = params_flood)
    if (response_flood$status_code == 200) {
      geojson_data_flood <- content(response_flood, as = "text", encoding = "UTF-8")
      sf_data_flood <- st_read(geojson_data_flood, quiet = TRUE)
      all_data_flood <- rbind(all_data_flood, sf_data_flood)
      if (nrow(sf_data_flood) < record_count_flood) break
      offset_flood <- offset_flood + record_count_flood
    } else {
      stop("Failed to retrieve data: HTTP status code ", response_flood$status_code)
    }
  }
  print(all_data_flood)
  
  saveRDS(all_data_flood, "all_data_flood.RDS")
  
  # Using the already created grid cells
  aus_grids <- readRDS("grids_10km.RDS")
  
  # Transforming and simplifying flood data
  flood_data_transformed_flood <- st_transform(all_data_flood, st_crs(aus_grids))
  flood_data_transformed_flood <- st_make_valid(flood_data_transformed_flood)
  tolerance_flood <- 100
  flood_data_simplified_flood <- st_simplify(flood_data_transformed_flood, dTolerance = tolerance_flood)
  
  # Classifying grid cells as flooded or unflooded
  flooded_cells_flood <- st_intersects(aus_grids, flood_data_simplified_flood)
  unflooded_cells_flood <- aus_grids[!unlist(lapply(flooded_cells_flood, function(x) length(x) > 0)), ]
  
  # Extracting indices of flooded cells
  flooded_indices_flood <- which(unlist(lapply(flooded_cells_flood, length)) > 0)
  actual_flooded_cells_flood <- aus_grids[flooded_indices_flood, ]
  st_crs(actual_flooded_cells_flood) <- st_crs(aus_grids)
  
  # Assign occurrences to grid cells
  assign_to_grid_flood <- function(data, grid_cells) {
    data_sf_flood <- st_as_sf(data, coords = c("lng", "lat"), crs = 4326) %>%
      st_set_crs(4326) %>%
      st_transform(st_crs(grid_cells)) %>%
      st_join(grid_cells, left = FALSE)
    return(data_sf_flood)
  }
  
  preflood_flooded_flood <- assign_to_grid_flood(preflood_data_flood, actual_flooded_cells_flood)
  preflood_unflooded_flood <- assign_to_grid_flood(preflood_data_flood, unflooded_cells_flood)
  postflood_flooded_flood <- assign_to_grid_flood(postflood_data_flood, actual_flooded_cells_flood)
  postflood_unflooded_flood <- assign_to_grid_flood(postflood_data_flood, unflooded_cells_flood)
  
  # Subsampling the data
  subsample_groups_flood <- function(...) {
    args <- list(...)
    min_n_flood <- min(sapply(args, nrow))
    lapply(args, function(x) x %>% slice_sample(n = min_n_flood))
  }
  
  subsampled_data_list_flood <- subsample_groups_flood(preflood_flooded_flood, preflood_unflooded_flood,
                                                       postflood_flooded_flood, postflood_unflooded_flood)
  
  # Combining the subsampled data
  subsampled_data_flood <- bind_rows(
    subsampled_data_list_flood[[1]] %>% mutate(status = "preflood_flooded"),
    subsampled_data_list_flood[[2]] %>% mutate(status = "preflood_unflooded"),
    subsampled_data_list_flood[[3]] %>% mutate(status = "postflood_flooded"),
    subsampled_data_list_flood[[4]] %>% mutate(status = "postflood_unflooded")
  )
  

    # Statistical analysis (flood analysis)
  
  # Get the top 40 well-sampled species (for flood analysis)
  top_species_flood <- subsampled_data_flood %>%
    count(species) %>%
    arrange(desc(n)) %>%
    head(40) %>%
    pull(species)
  
  # View sample size for species (flood analysis)
  species_counts_flood <- clean_data %>% 
    group_by(species) %>% 
    tally(sort = TRUE)
  print(species_counts_flood)
  
  # Filter the data to include only the top 40 well-sampled species (flood analysis)
  subsampled_data_top_species_flood <- subsampled_data_flood %>%
    filter(species %in% top_species_flood)
  
  # Calculate the presence in grid cells for each species and status (flood analysis)
  species_grid_presence_flood <- subsampled_data_top_species_flood %>%
    group_by(species, status) %>%
    summarize(n_grid_cells = n_distinct(st_coordinates(geometry))) %>%
    ungroup()
  print(species_grid_presence_flood)
  
  # Create contingency table
  contingency_table_flood <- table(subsampled_data_top_species_flood$species, subsampled_data_top_species_flood$status)
  # Chi-squared test
  chi_squared_test_flood <- chisq.test(contingency_table_flood)
  print(chi_squared_test_flood)
  
  # ANOVA for all statuses and species (including controls) in flood analysis
  anova_result_flood <- aov(n_grid_cells ~ status, data = species_grid_presence_flood)
  summary(anova_result_flood)
  
  # ANOVA for each species individually in flood analysis
  model_flood <- aov(n_grid_cells ~ species + status, data = species_grid_presence_flood)
  summary(model_flood)
  
  # Paired difference plot (flood analysis)
  reshaped_data_flood <- species_grid_presence_flood %>%
    filter(status %in% c("postflood_flooded", "preflood_flooded")) %>%
    spread(key = status, value = n_grid_cells, fill = 0) %>%
    group_by(species) %>%
    summarise(postflood_flooded = sum(postflood_flooded, na.rm = TRUE),
              preflood_flooded = sum(preflood_flooded, na.rm = TRUE)) %>%
    mutate(difference = postflood_flooded - preflood_flooded)
  ggplot(reshaped_data_flood, aes(x = reorder(species, -difference), y = difference)) + 
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(title = "Change in Grid Cells (Postflood - Preflood) in Flooded Areas",
         x = "Species", y = "Difference in Grid Cells") +
    theme_minimal()
  print(reshaped_data_flood[, c("species", "difference")])
  View(reshaped_data_flood)
  write.csv(reshaped_data_flood[, c("species", "difference")], "species_difference_flood.csv", row.names = FALSE)
  
 
  
  
  
  # ANOVA for flooded cells only (flood analysis)
  flooded_data_flood <- species_grid_presence_flood %>%
    filter(status %in% c("postflood_flooded", "preflood_flooded"))
  anova_result_flooded_flood <- aov(n_grid_cells ~ status, data = flooded_data_flood)
  summary(anova_result_flooded_flood)
  
  # T-Test for differences among individual species (flood analysis)
  results_list_flood <- list()
  for(spec in unique(paired_data_flood$species)) {
    spec_data_flood <- paired_data_flood[paired_data_flood$species == spec, ]
    test_result_flood <- t.test(spec_data_flood$postflood_flooded, spec_data_flood$prefire_flooded, paired = TRUE)
    results_list_flood[[spec]] <- test_result_flood
  }
  print(results_list_flood)
  

  
  # Species richness analysis (flood analysis)
  bootstrap_results_flood <- data.frame(category=character(), avg_richness=numeric(), iteration=integer())
  
  for(i in 1:num_iterations) {
    categories_flood <- list(
      preflood_flooded = subsampled_data_list_flood[[1]], 
      preflood_unflooded = subsampled_data_list_flood[[2]],
      postflood_flooded = subsampled_data_list_flood[[3]],
      postflood_unflooded = subsampled_data_list_flood[[4]]
    )
    for(cat_name in names(categories_flood)) {
      avg_richness_flood <- compute_avg_species_richness(categories_flood[[cat_name]])
      bootstrap_results_flood <- rbind(bootstrap_results_flood, data.frame(category=cat_name, avg_richness=avg_richness_flood, iteration=i))
    }
  }
  
  final_avg_richness_flood <- bootstrap_results_flood %>%
    group_by(category) %>%
    summarise(avg_richness = mean(avg_richness))
  
  print(final_avg_richness_flood)
  
  anova_result_richness_flood <- aov(avg_richness ~ category, data = bootstrap_results_flood)
  summary(anova_result_richness_flood)
  
  post_hoc_flood <- TukeyHSD(anova_result_richness_flood)
  print(post_hoc_flood)
  
  # Bar plot with error bars (flood analysis)
  SE_data_flood <- bootstrap_results_flood %>%
    group_by(category) %>%
    summarise(mean_richness = mean(avg_richness), SE = sd(avg_richness) / sqrt(n()))
  ggplot(SE_data_flood, aes(x = category, y = mean_richness, fill = category)) + 
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_richness - SE, ymax = mean_richness + SE), width = 0.25) +
    labs(y = "Average Number of Species Calling per Recording", x = "Status") +
    theme_minimal() +
    theme(legend.position = "none")
  
  
  # Calculate the number of occurrences in each group
  group_counts <- subsampled_data_flood %>%
    group_by(status) %>%
    summarise(count = n())
  
  # Print the result
  print(group_counts)
  
  
  # Count the occurrences in each group before subsampling
  count_preflood_flooded <- nrow(preflood_flooded_flood)
  count_preflood_unflooded <- nrow(preflood_unflooded_flood)
  count_postflood_flooded <- nrow(postflood_flooded_flood)
  count_postflood_unflooded <- nrow(postflood_unflooded_flood)
  
  # Create a data frame to display the counts
  group_counts_before_subsampling <- data.frame(
    status = c("preflood_flooded", "preflood_unflooded", "postflood_flooded", "postflood_unflooded"),
    count = c(count_preflood_flooded, count_preflood_unflooded, count_postflood_flooded, count_postflood_unflooded)
  )
  
  # Print the data frame
  print(group_counts_before_subsampling)
  
  