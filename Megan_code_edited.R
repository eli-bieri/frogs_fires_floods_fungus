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

setwd("C:/Users/ebier/Downloads/R")

clean_data <- readRDS("FrogID_clean_data.RDS")

# read in australia outline
aus <- rnaturalearth::ne_countries(country="Australia", returnclass="sf")
#st_crs(aus)

#make it the right data projection
aus.2 <- aus %>%
  st_transform(32754)

#make the grids- you can change cell size and see what happens!
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

crs <- "+proj=longlat +datum=WGS84 +no_defs"

points_sf <- st_as_sf(clean_data, coords = c("lng", "lat"), 
                      crs = crs, agr = "constant")

#set the crs for the grids
aus_grids <- aus_grids %>%
  st_set_crs(crs)

assigned_points <- as.data.frame(st_intersects(points_sf, aus_grids))


#going to be confusing but frog_id_10km doesnt mean its 10km grid - cant be bothered to change everything
FrogID_clean_10km <- clean_data %>%
  mutate(row.id = 1:nrow(.)) %>%
  left_join(assigned_points, by="row.id") %>%
  right_join(., grids_metadata, by="col.id") %>%
  dplyr::select(-row.id, -col.id)

# remove any frog records which have location accuracy > than 3000 m 
FrogID_clean_10km <- FrogID_clean_10km %>%
  dplyr::filter(location_accuracy <=3000)
  
FrogID_clean_10km$summer[FrogID_clean_10km$date > "2018-08-01" & FrogID_clean_10km$date < "2019-05-01"] <- "Prefire"
FrogID_clean_10km$summer[FrogID_clean_10km$date > "2020-08-01" & FrogID_clean_10km$date < "2021-05-01"] <- "Postfire"

FrogID_10km_indep <- unique(FrogID_clean_10km, by = c("user_id", "lat", "lng"))

FrogID_10km_indep <- FrogID_10km_indep %>%
  dplyr::select(-states_name, -state, -number_of_sp_calling, -user_group_name, -user_group_id, -validated_by, -image_urls, -quality_call, -inappropriate_content, -has_people_activity, -water_body, -habitat, -lga_name)

FrogID_10km_crop <- FrogID_10km_indep%>%
  filter(between(lat, -40, -21) & between(lng, 138, 155))

#take out anything not in the time frame
FrogID_10km_omit <-FrogID_10km_crop[complete.cases(FrogID_10km_crop$summer), ]

#number of submissions per grid, Prefire- and Postfire-
FrogID_10km_step1 <- FrogID_10km_omit %>%
  dplyr::filter(summer=="Prefire") %>%
  group_by(grid_id) %>%
  summarise(n_submissions_before = length(unique(id))) %>%
  full_join(FrogID_10km_omit, by = "grid_id")

FrogID_10km_step2 <- FrogID_10km_omit %>%
  dplyr::filter(summer=="Postfire") %>%
  group_by(grid_id) %>%
  summarise(n_submissions_after = length(unique(id))) %>%
  full_join(FrogID_10km_step1, by = "grid_id")

length(unique(FrogID_10km_step2$grid_id))

#number of submissions each species has
species_wellsamp <- FrogID_10km_step2%>%
  group_by(species) %>%
  summarise(n_spp_wellsamp = length(unique(id))) %>%
  full_join(FrogID_10km_step2)

#filtering dataset
species_wellsamp <- species_wellsamp %>%
  filter(n_submissions_before>25 & n_submissions_after>25)

#selecting the most well sampled species
top_N <- species_wellsamp %>%
  dplyr::select(species, n_spp_wellsamp)

top_N <- unique(top_N, by=c("species"))

top_30 <- top_N%>%
  arrange(desc(n_spp_wellsamp))%>%
  slice(1:30)%>%
  inner_join(species_wellsamp, by = "species") %>%
  dplyr::select(-n_spp_wellsamp.y) %>% 
  dplyr::rename(n_spp_wellsamp = n_spp_wellsamp.x)

#length(unique(species_wellsamp$grid_id))


#redo n_submission column after excluding things

FrogID_10km_step3 <- top_30 %>%
  dplyr::filter(summer=="Prefire") %>%
  group_by(grid_id) %>%
  summarise(n_submissions_before_2 = length(unique(id))) %>%
  full_join(top_30, by = "grid_id")

FrogID_10km_step4 <- top_30 %>%
  dplyr::filter(summer=="Postfire") %>%
  group_by(grid_id) %>%
  summarise(n_submissions_after_2 = length(unique(id))) %>%
  full_join(FrogID_10km_step3, by = "grid_id")

clean_top_30 <- FrogID_10km_step4 %>%
  select(grid_id, species, id, lat, lng, date, summer, n_spp_wellsamp, n_submissions_after_2, n_submissions_before_2)


#for the equal grid iterations
it_step <- clean_top_30

current_grid <- it_step %>% subset(grid_id == 883) #just need a random grid to start off with to have the right columns
equal_grids_full <- current_grid[-c(1:nrow(current_grid)), ]
for(i in unique(it_step$grid_id)){
  current_grid <- it_step %>% subset(grid_id == i)
  current_grid_Prefire <- current_grid %>% filter(summer=="Prefire")
  current_grid_Postfire <- current_grid %>% filter(summer=="Postfire")
  
  a <- unique(current_grid$n_submissions_after_2)
  b <- unique(current_grid$n_submissions_before_2)
  
  unique_id_Prefire <- current_grid %>%
    filter(summer=="Prefire") %>%
    distinct(id)
  
  unique_id_Postfire <- current_grid %>%
    filter(summer=="Postfire") %>%
    distinct(id)
  
  if(a<b){
    equal_Prefire <- sample_n(unique_id_Prefire, a)
    equal_Prefire_full <- current_grid_Prefire[(current_grid_Prefire$id %in% equal_Prefire$id), ]
    equal_grids <- rbind(equal_Prefire_full, current_grid_Postfire)
  } else if (b<a){
    equal_Postfire <- sample_n(unique_id_Postfire, b)
    equal_Postfire_full <- current_grid_Postfire[(current_grid_Postfire$id %in% equal_Postfire$id), ]
    equal_grids <- rbind(equal_Postfire_full, current_grid_Prefire)
  } else if (b==a) {
    equal_grids <- rbind(current_grid_Prefire, current_grid_Postfire)
  }
  equal_grids_full <- rbind(equal_grids_full, equal_grids)
}


#test<- equal_grids_full
#length(unique(test$grid_id))

#check_Prefire <- equal_grids_full %>% filter(summer=="Prefire") 
#check_Postfire <- equal_grids_full %>% filter(summer=="Postfire") 
#length(unique(check_Postfire$species))


#making change variables on equal grid data set
EG1 <- equal_grids_full %>%
  dplyr::filter(summer=="Prefire")%>%
  group_by(species)%>%
  summarise(n_occurrence_before = length(unique(grid_id)))%>%
  full_join(equal_grids_full)

EG2 <- equal_grids_full %>%
  dplyr::filter(summer=="Postfire")%>%
  group_by(species)%>%
  summarise(n_occurrence_after = length(unique(grid_id)))%>%
  full_join(EG1)

sub_rep_10km_30 <- EG2%>%
  mutate(
    change = (n_occurrence_after-n_occurrence_before),
    percent_change = (((n_occurrence_after-n_occurrence_before)/n_occurrence_before)*100)
  )







#putting it together for anlaysis 
simp_1 <- sub_rep_10km_30 %>%
  mutate(percent_change1 = percent_change)%>%
  select(species, percent_change1) 
simp_1<- unique(simp_1, by= c("species"))

simp_2 <- sub_rep_10km_30 %>%
  mutate(percent_change2 = percent_change) %>%
  select(species, percent_change2)
simp_2<- unique(simp_2, by= c("species"))

simp_3 <- sub_rep_10km_30 %>%
  mutate(percent_change3 = percent_change)%>%
  select(species, percent_change3)
simp_3<- unique(simp_3, by= c("species"))

simp_4 <- sub_rep_10km_30 %>%
  mutate(percent_change4 = percent_change)%>%
  select(species, percent_change4)
simp_4<- unique(simp_4, by= c("species"))

simp_5 <- sub_rep_10km_30 %>%
  mutate(percent_change5 = percent_change)%>%
  select(species, percent_change5)
simp_5<- unique(simp_5, by= c("species"))

simp_6 <- sub_rep_10km_30 %>%
  mutate(percent_change6 = percent_change)%>%
  select(species, percent_change6)
simp_6<- unique(simp_6, by= c("species"))

simp_7 <- sub_rep_10km_30 %>%
  mutate(percent_change7 = percent_change)%>%
  select(species, percent_change7)
simp_7<- unique(simp_7, by= c("species"))

simp_8 <- sub_rep_10km_30 %>%
  mutate(percent_change8 = percent_change)%>%
  select(species, percent_change8)
simp_8<- unique(simp_8, by= c("species"))

simp_9 <- sub_rep_10km_30 %>%
  mutate(percent_change9 = percent_change)%>%
  select(species, percent_change9)
simp_9<- unique(simp_9, by= c("species"))

simp_10 <- sub_rep_10km_30 %>%
  mutate(percent_change10 = percent_change)%>%
  select(species, percent_change10)
simp_10<- unique(simp_10, by= c("species"))

simp_11 <- sub_rep_10km_30 %>%
  mutate(percent_change11 = percent_change)%>%
  select(species, percent_change11)
simp_11<- unique(simp_11, by= c("species"))

simp_12 <- sub_rep_10km_30 %>%
  mutate(percent_change12 = percent_change)%>%
  select(species, percent_change12)
simp_12<- unique(simp_12, by= c("species"))

simp_13 <- sub_rep_10km_30 %>%
  mutate(percent_change13 = percent_change)%>%
  select(species, percent_change13)
simp_13<- unique(simp_13, by= c("species"))

simp_14 <- sub_rep_10km_30 %>%
  mutate(percent_change14 = percent_change)%>%
  select(species, percent_change14)
simp_14<- unique(simp_14, by= c("species"))

simp_15 <- sub_rep_10km_30 %>%
  mutate(percent_change15 = percent_change)%>%
  select(species, percent_change15)
simp_15<- unique(simp_15, by= c("species"))

simp_16 <- sub_rep_10km_30 %>%
  mutate(percent_change16 = percent_change)%>%
  select(species, percent_change16)
simp_16<- unique(simp_16, by= c("species"))

simp_17 <- sub_rep_10km_30 %>%
  mutate(percent_change17 = percent_change)%>%
  select(species, percent_change17)
simp_17<- unique(simp_17, by= c("species"))

simp_18 <- sub_rep_10km_30 %>%
  mutate(percent_change18 = percent_change)%>%
  select(species, percent_change18)
simp_18<- unique(simp_18, by= c("species"))

simp_19 <- sub_rep_10km_30 %>%
  mutate(percent_change19 = percent_change)%>%
  select(species, percent_change19)
simp_19<- unique(simp_19, by= c("species"))

simp_20 <- sub_rep_10km_30 %>%
  mutate(percent_change20 = percent_change)%>%
  select(species, percent_change20)
simp_20<- unique(simp_20, by= c("species"))

simp_21 <- sub_rep_10km_30 %>%
  mutate(percent_change21 = percent_change)%>%
  select(species, percent_change21)
simp_21<- unique(simp_21, by= c("species"))

simp_22 <- sub_rep_10km_30 %>%
  mutate(percent_change22 = percent_change)%>%
  select(species, percent_change22)
simp_22<- unique(simp_22, by= c("species"))

simp_23 <- sub_rep_10km_30 %>%
  mutate(percent_change23 = percent_change)%>%
  select(species, percent_change23)
simp_23<- unique(simp_23, by= c("species"))

simp_24 <- sub_rep_10km_30 %>%
  mutate(percent_change24 = percent_change)%>%
  select(species, percent_change24)
simp_24<- unique(simp_24, by= c("species"))

simp_25 <- sub_rep_10km_30 %>%
  mutate(percent_change25 = percent_change)%>%
  select(species, percent_change25)
simp_25<- unique(simp_25, by= c("species"))

simp_26 <- sub_rep_10km_30 %>%
  mutate(percent_change26 = percent_change)%>%
  select(species, percent_change26)
simp_26<- unique(simp_26, by= c("species"))

simp_27 <- sub_rep_10km_30 %>%
  mutate(percent_change27 = percent_change)%>%
  select(species, percent_change27)
simp_27<- unique(simp_27, by= c("species"))

simp_28 <- sub_rep_10km_30 %>%
  mutate(percent_change28 = percent_change)%>%
  select(species, percent_change28)
simp_28<- unique(simp_28, by= c("species"))

simp_29 <- sub_rep_10km_30 %>%
  mutate(percent_change29 = percent_change)%>%
  select(species, percent_change29)
simp_29<- unique(simp_29, by= c("species"))

simp_30 <- sub_rep_10km_30 %>%
  mutate(percent_change30 = percent_change)%>%
  select(species, percent_change30)
simp_30<- unique(simp_30, by= c("species"))




#rep_x30_full <- full_join(simp_1, simp_2, simp_3, simp_4, simp_5, simp_6, simp_7, simp_8, simp_9, simp_10, simp_11, simp_12, simp_13, simp_14, simp_15, simp_16, simp_17, simp_18, simp_19, simp_20, simp_21, simp_22, simp_23, simp_24, simp_25, simp_26, simp_27, simp_28, simp_29, simp_30, by = "species")

rep_x30_full <- full_join(simp_1, simp_2, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_3, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_4, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_5, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_6, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_7, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_8, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_9, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_10, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_11, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_12, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_13, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_14, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_15, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_16, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_17, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_18, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_19, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_20, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_21, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_22, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_23, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_24, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_25, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_26, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_27, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_28, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_29, by = "species")
rep_x30_full <- full_join(rep_x30_full, simp_30, by = "species")

rep_x30_full_long <- rep_x30_full %>%
  gather(Iteration, percent_change, percent_change1:percent_change30)

rep_x30_full_long$Grid_size <- "10"


#make mean and standard deviation and error for each species
rep_x30_full_long <- rep_x30_full_long %>%
  group_by(species)%>%
  mutate(mean = mean(percent_change),
         sd = sd(percent_change),
         se = sd(percent_change) / sqrt(length(percent_change)))




rep_x30_full_long_10km <- rep_x30_full_long %>%
  arrange(mean, species)


#rep_x20_full_long$species = as.factor(rep_x20_full_long$species)

my_spp_order <- # specify the order to plot species based on duration trend
  rep_x30_full_long_10km %>%
  arrange(desc(mean)) %>%
  .$species

my_spp_order <- unique(my_spp_order, by= c("species"))

unique(rep_x30_full_long$species)

magic <- ggplot() +
  geom_hline(yintercept = 0,
             linetype = "dotted") + 
  geom_boxplot(data=rep_x30_full_long_10km, aes(x=factor(species, my_spp_order), y= percent_change), colour = "black", fill = "white", shape = 21)+
  labs(y = "Change in grid occurence (%)", x= "Species") +
  ylim(-25, 50)+
  coord_flip() +
  force_panelsizes(rows = unit(20, "cm"),
                   cols = unit(10, "cm")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))+
  #theme(axis.text.y = element_text(face = "italic"))
  #theme(axis.text.y = element_text(face = bold.labels))
  scale_x_discrete(labels=c("Litoria caerulea"=expression(bolditalic("Litoria caerulea")), "Litoria peronii"=expression(bolditalic("Litoria peronii")), 
                            "Limnodynastes peronii"=expression(bolditalic("Limnodynastes peronii")), "Rhinella marina"=expression(bolditalic("Rhinella marina")), 
                            "Litoria phyllochroa"=expression(bolditalic("Litoria phyllochroa")), "Adelotus brevis"=expression(italic("Adelotus brevis")), 
                            "Pseudophryne coriacea"=expression(italic("Pseudophryne coriacea")), "Uperoleia laevigata"=expression(italic("Uperoleia laevigata")),
                            "Litoria tyleri"=expression(italic("Litoria tyleri")), "Litoria latopalmata"=expression(italic("Litoria latopalmata")),
                            "Limnodynastes dumerilii"=expression(italic("Limnodynastes dumerilii")), "Litoria rubella"=expression(italic("Litoria rubella")),
                            "Crinia parinsignifera"=expression(italic("Crinia parinsignifera")), "Litoria fallax"=expression(italic("Litoria fallax")),
                            "Crinia signifera"=expression(italic("Crinia signifera")), "Litoria balatus"=expression(italic("Litoria balatus")),
                            "Litoria ewingii"=expression(italic("Litoria ewingii")), "Limnodynastes tasmaniensis"=expression(italic("Limnodynastes tasmaniensis")),
                            "Platyplectrum ornatum"=expression(italic("Platyplectrum ornatum")), "Litoria gracilenta"=expression(italic("Litoria gracilenta")),
                            "Pseudophryne australis"=expression(italic("Pseudophryne australis")), "Litoria quiritatus"=expression(italic("Litoria quiritatus")),
                            "Uperoleia fusca"=expression(italic("Uperoleia fusca")), "Mixophyes fasciolatus"=expression(italic("Mixophyes fasciolatus")),
                            "Litoria verreauxii"=expression(italic("Litoria verreauxii")), "Litoria chloris"=expression(italic("Litoria chloris")),
                            "Limnodynastes terraereginae"=expression(italic("Limnodynastes terraereginae")), "Litoria dentata"=expression(italic("Litoria dentata")),
                            "Litoria nasuta"=expression(italic("Litoria nasuta")), "Uperoleia rugosa"=expression(italic("Uperoleia rugosa")), parse=TRUE))

magic

install.packages("png")
capabilities("png")

ggsave("Output/rep30_boxplot_submission_2.pdf", width = 18, height = 22, units = "cm", dpi = 600)

ggsave("Output/rep30_boxplot_submission_2.png", width = 18, height = 22, units = "cm", dpi = 600)



length(unique(FrogID_10km_omit$id))

blah <- clean_top_30 %>%
  dplyr::select(species, n_spp_wellsamp)
blah <- unique(blah)







#models
#average between 0 and 1

current_it <- sub_rep_10km_30

joiner <- current_it%>%
  group_by(grid_id) %>%
  mutate(
    mean_lat = mean(lat),
    mean_lng = mean(lng)
  )%>%
  dplyr::select(grid_id, mean_lat, mean_lng)
joiner <- unique(joiner, by = c("grid_id, mean_lat, mean_lng"))

all_grids <- current_it %>%
  dplyr::select(grid_id) %>%
  distinct(grid_id)

change_gam_export <- mod_cut[-c(1:nrow(mod_cut)), ]

for(i in unique(current_it$species)){
  
  current_species <- current_it %>%
    filter(species == i)
  
  grids_Prefire <- current_species %>%
    dplyr::filter(summer == "Prefire") %>%
    dplyr::select(grid_id) %>%
    distinct(grid_id)
  grids_Prefire$i_grid_id <- grids_Prefire$grid_id
  mod <- full_join(grids_Prefire, all_grids, by = "grid_id")
  mod$i_grid_id[is.na(mod$i_grid_id)] <- 0
  mod_bi_Prefire <- mod %>%
    mutate(binary = ifelse(i_grid_id>1, 1, 0)) %>%
    dplyr::select(-i_grid_id)
  mod_bi_Prefire$fire <- "Prefire"
  
  
  grids_Postfire <- current_species %>%
    dplyr::filter(summer == "Postfire") %>%
    dplyr::select(grid_id) %>%
    distinct(grid_id)
  grids_Postfire$i_grid_id <- grids_Postfire$grid_id
  mod <- full_join(grids_Postfire, all_grids, by = "grid_id")
  mod$i_grid_id[is.na(mod$i_grid_id)] <- 0
  mod_bi_Postfire <- mod %>%
    mutate(binary = ifelse(i_grid_id>1, 1, 0)) %>%
    dplyr::select(-i_grid_id)
  mod_bi_Postfire$fire <- "Postfire"
  
  
  mod_bi <- rbind(mod_bi_Prefire, mod_bi_Postfire)
  mod_all <- full_join(mod_bi, joiner, by = "grid_id")
  mod_all$species <- i
  
  mod_cut <- unique(mod_all)
  
  change_gam_export <- rbind(change_gam_export, mod_all)
  
}

sub_gam_rep_30 <- change_gam_export


sub_gam_rep_30$Iteration <- "1"
sub_gam_rep_30$Iteration <- "2"
sub_gam_rep_30$Iteration <- "3"
sub_gam_rep_30$Iteration <- "4"
sub_gam_rep_30$Iteration <- "5"
sub_gam_rep_30$Iteration <- "6"
sub_gam_rep_30$Iteration <- "7"
sub_gam_rep_30$Iteration <- "8"
sub_gam_rep_30$Iteration <- "9"
sub_gam_rep_30$Iteration <- "10"
sub_gam_rep_30$Iteration <- "11"
sub_gam_rep_30$Iteration <- "12"
sub_gam_rep_30$Iteration <- "13"
sub_gam_rep_30$Iteration <- "14"
sub_gam_rep_30$Iteration <- "15"
sub_gam_rep_30$Iteration <- "16"
sub_gam_rep_30$Iteration <- "17"
sub_gam_rep_30$Iteration <- "18"
sub_gam_rep_30$Iteration <- "19"
sub_gam_rep_30$Iteration <- "20"
sub_gam_rep_30$Iteration <- "21"
sub_gam_rep_30$Iteration <- "22"
sub_gam_rep_30$Iteration <- "23"
sub_gam_rep_30$Iteration <- "24"
sub_gam_rep_30$Iteration <- "25"
sub_gam_rep_30$Iteration <- "26"
sub_gam_rep_30$Iteration <- "27"
sub_gam_rep_30$Iteration <- "28"
sub_gam_rep_30$Iteration <- "29"
sub_gam_rep_30$Iteration <- "30"


all_gam_mean <- rbind(sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30, sub_gam_rep_30)

mean_loc <- FrogID_10km_omit %>%
  group_by(grid_id)%>%
  mutate(mean_lat_2 = mean(lat),
         mean_lng_2 = mean(lng)) %>%
  dplyr::select(grid_id, mean_lat_2, mean_lng_2)
mean_loc <- unique(mean_loc)

all_gam_loc <- inner_join(all_gam_mean, mean_loc, by = "grid_id")
all_gam_loc <- all_gam_loc%>%
  dplyr::select(-mean_lat, -mean_lng)

Prefire_gam_mean<- all_gam_loc %>%
  filter(fire == "Prefire") %>%
  group_by(species, grid_id) %>%
  mutate(mean_bi = mean(binary)) %>%
  dplyr::select(-Iteration, -binary)
Prefire_gam_mean <- unique(Prefire_gam_mean, by = c("species", "grid_id", "mean_bi"))

Postfire_gam_mean<- all_gam_loc%>%
  filter(fire == "Postfire") %>%
  group_by(species, grid_id) %>%
  mutate(mean_bi = mean(binary))%>%
  dplyr::select(-Iteration, -binary)
Postfire_gam_mean <- unique(Postfire_gam_mean, by = c("species", "grid_id", "mean_bi", "fire"))

all_bi_mean <- rbind(Postfire_gam_mean, Prefire_gam_mean)


#gams with binary mean change
for(i in unique(all_bi_mean$species)){
  
  current_gam_spp <- all_bi_mean %>%
    filter(species == i)
  
  #current_gam_spp <- unique(current_gam_spp, by = c("grid_id, summer, mean_change"))
  
  gam_cut <- mgcv::gam(mean_bi ~ fire + s(mean_lng_2, mean_lat_2), data = current_gam_spp, method = "REML")
  
  tester <- capture.output(summary(gam_cut))
  
  write_yaml(tester, file = paste0(i, "_gam", ".yml"))
  
}




saveRDS(all_it_bi_mean, "Data/mean_bi_gams_sub.RDS")

summary(gam_cut)

gam.check(gam_cut)




#BH p-value correction

pvaluesTrial2 <- c(0.665, 0.579, 0.451, 0.441, 0.781, 0.829, 0.235, 0.898, 0.373, 0.558, 0.00728, 0.804, 0.168, 0.111, 0.558, 0.00901, 0.904, 0.502, 0.255, 0.559, 0.169, 0.0594, 0.375, 0.311, 0.687, 0.502, 0.471, 0.302, 0.267, 0.528)

p.adjust(pvaluesTrial2, method="BH")


