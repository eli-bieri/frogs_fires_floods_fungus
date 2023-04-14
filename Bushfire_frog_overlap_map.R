library(tidyverse)
library(sf)
library(here)
library(concaveman)
library(ggplot2)
library(cowplot)

# Set the working directory
setwd("C:/Users/ebier/Downloads/R")

# Read FrogID_bushfire dataset
frogs <- read_rds("FrogID_clean_data.rds")

# Filter the dataset for the three species
selected_species <- c("Mixophyes fleayi", "Mixophyes balbus", "Mixophyes iteratus")

filtered_frogs <- frogs %>%
  filter(species %in% selected_species) %>%
  dplyr::select(species, lat, lng)

# Convert to sf object
filtered_frogs_sf <- st_as_sf(filtered_frogs, coords = c("lng", "lat"), crs = 4326)

# Read the fire layer
fire <- st_read("FESM_area_84.shp")

# Assign CRS to the fire layer if missing
if (is.na(st_crs(fire))) {
  st_crs(fire) <- 4326
}


# Transform filtered_frogs_sf and fire to a projected CRS
filtered_frogs_projected <- st_transform(filtered_frogs_sf, 3577)
fire_projected <- st_transform(fire, 3577)

# Calculate the overlap
overlap_sf <- st_intersection(filtered_frogs_projected, fire_projected)

# Join the original frog dataset with the overlap data
joined_data <- st_join(filtered_frogs_projected, overlap_sf, left = TRUE)

# Transform back to geographic CRS
joined_data <- st_transform(joined_data, 4326)

colnames(centroids_df)


# Calculate centroids
centroids <- st_centroid(st_geometry(joined_data))
centroids_df <- cbind(as.data.frame(joined_data), st_coordinates(centroids))

# Convert to a data frame for plotting
map_data <- centroids_df

# Plot the map
ggplot() +
  geom_sf(data = fire, fill = "#FEC3A6", colour = "#FEC3A6") +
  geom_point(data = map_data,
             aes(x = X, y = Y, color = species.x),
             alpha = 0.6, size = 3) +
  scale_color_manual(values = c("Mixophyes fleayi" = "palegoldenrod",
                                "Mixophyes balbus" = "lightblue",
                                "Mixophyes iteratus" = "lightgreen")) +
  labs(title = "Mixophyes Bushfire Overlap",
       subtitle = "M. fleayi, M. balbus, and M. iteratus",
       caption = "Data source: FrogID and FESM NSW",
       color = "Species") +
  theme_void() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        plot.caption = element_text(size = 10))
