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

setwd("C:/Users/ebier/Downloads/R")

clean_data <- readRDS("FrogID_clean_data.RDS")

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

# Get the top 10 well-sampled species
top_species <- subsampled_data %>%
  count(species) %>%
  top_n(10, wt = n) %>% # Rank by count (n)
  pull(species)


# Filter the data to only include the top 10 well-sampled species
subsampled_data_top_species <- subsampled_data %>%
  filter(species %in% top_species)

# Calculate the presence in grid cells for each species and status (prefire_burnt, postfire_burnt, prefire_unburnt, postfire_unburnt)
species_grid_presence <- subsampled_data_top_species %>%
  group_by(species, status) %>%
  summarize(n_grid_cells = n_distinct(st_coordinates(geometry))) %>%
  ungroup()

# Print the results
species_grid_presence
