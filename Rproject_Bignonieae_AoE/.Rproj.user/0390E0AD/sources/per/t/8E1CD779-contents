## Processing NDM/VNDM outputs

# This script contains a set of of instructions and R functions to import 
# the outputs of VNDM/NDM into R. 

# The R project directory has the followin structure:

# Scripts are organized in separate files as follows:
# 1) Functions: Source files with the functions design for different tasks.
# 2) Workflow (This script): General pipeline of the analysis.
# 3) Additional scripts detailing the exploration of the consensus areas.

# This R project is organized in the following system of directories:
# ./Working_directory/ -> Working directory
#      data/           -> Data files
#      output/         -> Outputs from scripts
#      figs/           -> Figures
#      scr/            -> R Scripts with functions and general workflow

# The content in ./figs and ./output can be completely reproduced from ./scr.

## Nomenclature of the VNDM-NDM outputs.
# Results of the analysis of endemicity were saved and named using
# a specific nomenclature. This Area Name Code considers the spatial scale 
# of the analysis (grid size), the NDM code of the consensus area, 
# a relative geographical location, and the type of GIS file (i.e grids or points).

## Example: 1dgD_CA0_AF-grid 
#    1dg = grid size; 
#    D = Default analysis; 
#   CA0 = Consensus area number 0;
#    AF = A name for the geographical region where the area is located.
#    grid/point = ndm grid output file type.

# This Area Name Code allow area tracking along the exploration of areas of endemism.



# 0. Loading the functions to explore the VNDM-NDM results and R libraries -----

# Functions
source("./src/NDM_inputfiles.R")
source("./src/NDM_intoR.R")
source("./src/NDM_visualization.R")
source("./src/NDM_tables.R")
source("./src/NDM_acomparison.R")

# Libraries

# The following libraries are used in this R project. 
# Please install them first before executing the code below.

# install.packages(c("readxl", "tidyverse", "ggpubr", "gtable", "ggthemes", 
#                    "ggalt", "viridis",  "raster", "rgdal", "sp", "sf", 
#                    "rmapshaper", "UpSetR", "VennDiagram", "grid", "gridExtra",
#                    "patchwork", "ggrepel", "lwgeom"))

library(tidyverse)
library(readxl)
library(purrr)
library(raster)
library(ggpubr)
library(gridExtra)
library(viridis)



#------------------------- I) PREPARAING DATA FOR NDM/VNDM ------------------------- 
# 1. Select working directory or R project ----

# setwd("./foo") ## This is not necessary if you are working inside the RStudio project

# 2. Loading occurrence database into R ----

Total_DB <- read.csv("./data/Bignonieae_Database_sample.csv") 
glimpse(Total_DB)

# NOTE_0: This is a sample with the genera Pachyptera and Tanaecium. My database 
# had the following columns: 
# NAME1: for species names.
# XCOOR: latitude coordinate.
# YCOOR: longitude coordinate.
# All the functions used in this script requires to use these column names 

# 3. Selecting only unique records ----

Unique_DB <- Total_DB %>% distinct(NAME1, YCOOR, XCOOR, .keep_all = TRUE)

# Filling decimals to standardize up to 6 decimal digits.

Unique_DB$XCOOR <- format(round(Unique_DB$XCOOR, 6), nsmall = 6)
Unique_DB$YCOOR <- format(round(Unique_DB$YCOOR, 6), nsmall = 6)

# 4. Creating a list of coordinates in the .xyd format ----

# The function ndm_xyd() creates two files: NDM_DB.xyd, which contains species 
# occurrences in the format .xyd; and the file Spnames_list.txt which contains 
# the list of species and their respective NDM code.

ndm_xyd(Unique_DB)

# 5. Loading the list of species and NDM code into a data frame ----

spplist <- ndm_spplist("./output/Spnames_list.txt")



#---------------------- II) IMPORTING VNDM-NDM OUTPUTS INTO R ------------------------- 
# 1. Producing polygons from NDM grid output file at each spatial scale ----
# 1.1. Creating list of files and consensus areas names ----

# NOTE_1: These functions requires that original VNDM files have the area code name 
# specified above.

ndm_listgrid.files(path = "./data/", degree = 1)
ndm_listgrid.files(path = "./data/", degree = 2)
ndm_listgrid.files(path = "./data/", degree = 3)

# 1.2 Creating polygons from grid file in batch ####

## The function ndm_gridtopoly() creates the polygon shapefiles of consensus 
## areas from the grid coordinates ouput of NDM, and saves them in 
## "./output/GIS/". A copy of the ouput of this function can be found in
## "./data/GIS/Consensus_Areas_Shapefiles/".

pmap(list(file = files1dg,
          degree = 1,
          save.area.shp = TRUE,
          save.directory = "./output/GIS/",
          area.shp.name = anames1dg), ndm_gridtopoly)

pmap(list(file = files2dg,
          degree = 2,
          save.area.shp = TRUE,
          save.directory = "./output/GIS/",
          area.shp.name = anames2dg), ndm_gridtopoly)

pmap(list(file = files3dg,
          degree = 3,
          save.area.shp = TRUE,
          save.directory = "./output/GIS/",
          area.shp.name = anames3dg), ndm_gridtopoly)

# 2. Consensus areas data frames and NDM data matrices from NDM output file ----

# ndm_dataparsing function: Parsing the output text file into R
# ndmfile = c("./1dgD_total.out", "./2dgD_total.out", "./3dgD_total.out")
# grid_size = c(1, 2, 3)

e_ndm <- raster::extent(c(-111, -31, -40, 40))

df_1D <- ndm_dataparsing(ndmfile = "./data/NDM_1dg.out", 
                         grid_size = 1, ana_param = "D", 
                         create_matrix = TRUE, matrix_extent = e_ndm) 

df_2D <- ndm_dataparsing(ndmfile = "./data/NDM_2dg.out", 
                         grid_size = 2, ana_param = "D", 
                         create_matrix = TRUE, matrix_extent = e_ndm) 

df_3D <- ndm_dataparsing(ndmfile = "./data/NDM_3dg.out", 
                         grid_size = 3, ana_param = "D", 
                         create_matrix = TRUE, matrix_extent = e_ndm) 

ca_db <- rbind(df_1D, df_2D, df_3D)

ca_scores_files <- list.files("./data", pattern = ".acf", full.names = TRUE)


# Extracting consensus area scores from .acf files
ca_scores <- ndm_parsing_cascores(x = ca_scores_files)

ca_db$Cons_Area_Max_Score <- as.numeric(ca_db$Cons_Area_Max_Score)
ca_db$Cons_Area_Min_Score <- as.numeric(ca_db$Cons_Area_Min_Score)

ca_db <- left_join(ca_db, ca_scores, by = c("Spatial_Scale", "Consensus_Area")) %>% 
  mutate(Cons_Area_Min_Score.x = Cons_Area_Min_Score.y,
         Cons_Area_Max_Score.x = Cons_Area_Max_Score.y,
         Area_name = paste0(Analysis, Spatial_Scale, Consensus_Area) %>% 
           str_replace_all(pattern = " ", replacement = "")) %>% 
  dplyr::select(-c(Cons_Area_Min_Score.y, Cons_Area_Max_Score.y)) %>% 
  rename(Cons_Area_Max_Score = Cons_Area_Max_Score.x, 
         Cons_Area_Min_Score = Cons_Area_Min_Score.x) %>% 
  mutate(
    Species = replace(Species, Species == "Dolichandra unguis",
                      "Dolichandra unguis-cati"),
    Species = replace(Species, Species == "Bignonia sanctae", 
                      "Bignonia sanctae-crucis"))

# Producing area names for manuscript

geoname_df  <- ca_db %>%
  distinct(Area_name, Approx_geo, .keep_all = TRUE) %>% 
  group_by(Approx_geo) %>%
  add_tally()  %>%
  mutate(aoen = 1:n,
         Geo_name = paste(Approx_geo, aoen, sep = "_")) %>%
  dplyr::select(-aoen) %>%
  arrange(Approx_geo, Spatial_Scale) %>%
  dplyr::select(Approx_geo, Area_name, Geo_name)

ca_db <-left_join(ca_db, geoname_df)

# Save table
# write.csv(ca_db, "./output/CA_DB.csv")

ca_db <- ca_db %>%
  mutate(Geo_sector = case_when(
    Approx_geo %in% c("NEO", "SA") ~ "Continental patterns",
    Approx_geo %in% c("YUC", "MESO") ~ "Central America - Mesoamerica",
    Approx_geo %in% c("NWSA", "Col", "ColVen") ~ "North Western South America",
    Approx_geo %in% c("GU") ~ "Guiana Shield",
    Approx_geo %in% c("AM") ~ "Amazonia",
    Approx_geo %in% c("DD", "DDb", "THDD") ~ "Dry and Open Vegetation Diagonal",
    Approx_geo %in% c("ESA", "AF") ~ "Eastern South America and AtlanticForest")
  )

ca_db$Geo_sector <- factor(
  ca_db$Geo_sector, 
  levels = c("Continental patterns", "Amazonia", "Guiana Shield",
             "Eastern South America and AtlanticForest",  
             "Dry and Open Vegetation Diagonal",  "North Western South America",
             "Central America - Mesoamerica"), ordered = TRUE)


ca_db <- ca_db %>% 
  arrange(Geo_sector)

# Save table
# write.csv(ca_db, "./output/FINAL_NDMTable.csv")



#------------------------- III) EXPLORING CONSENSUS AREAS A POSTERIORI -------------------------
# 1. Summarizing results ----

# Loading consensus areas information data frame. 
ca_df <- read.csv("./output/CA_DB.csv", stringsAsFactors = FALSE)[-1]
glimpse(ca_df)

spp_list <- ndm_spplist(filepath = "./data/Spnames_list.txt")

# 1.1. Summary of information about consensus areas ----
(Summ_info <- ndm_casum(ca_df))


# 2. Maps of consensus areas ----

# Filtering the ca_df to have only the AoEs
ca_df <- read_csv("./output/CA_DB.csv")[-1] %>% 
  distinct(Area_name, Approx_geo, .keep_all = TRUE) 

# Vector of AoEs names original code 
aoe_original_code <- ca_df$Area_name

# Vector of geo_names
geo_names_vec <- ca_df$Geo_name

names(geo_names_vec) <- aoe_original_code

# vector of layers: sorting the layers to match all the names in the ca_df
ca_layer_names <- list.files(path = "./data/GIS/Consensus_Areas_Shapefiles/", 
                             pattern = "[1-3]dg.*grid.shp") %>%
  gsub(pattern = "([1-3]dg.*grid).shp", replacement = "\\1") 

names(ca_layer_names) <- gsub(ca_layer_names, pattern = "^([1-3])dgD_(.*)_.*_grid", replacement = "D\\1\\2")

ca_layer_names <- ca_layer_names[order(factor(names(ca_layer_names), 
                                              levels = aoe_original_code ))]

names(ca_layer_names) == names(geo_names_vec) # ALL TRUE

# List of maps with areas to plot
ALL_Aoe_list <- as.list(aoe_original_code)

# Color palette
aoe_pal <- c(rep("#CAB2D6", 2), rep("#A6CEE3", 3), rep("#1F78B4", 7), 
             rep("#B2DF8A", 8), rep("#33A02C", 22), rep("#FB9A99", 6),
             rep("#E31A1C", 16), rep("#FFFF99", 6))

all_aoe_plots_list <- purrr::pmap(.l = list(map_list_element = ALL_Aoe_list, 
                                            area_name = geo_names_vec,
                                            area_col = aoe_pal),
                                  .f = ndm_plotca4,
                                  dsn_aoes = "./data/GIS/Consensus_Areas_Shapefiles/",
                                  aoe_layers = ca_layer_names,
                                  alpha_aoe = 0.5)

names(all_aoe_plots_list) <- aoe_original_code

lay <- rbind(c(1, 2, 3, 4),
             c(5, 6, 7, 8),
             c(9, 10 ,11, 12),
             c(13, 14, 15, 16),
             c(17, 18, 19, 20)
)

ggsave(
  path = "./figs/",
  filename = "All_AoEs.pdf", 
  device = "pdf",
  plot = marrangeGrob(all_aoe_plots_list, layout_matrix = lay), 
  width = 21, height = 30, units = "cm"
)


# 3. Overlap plot ----

# 3.1. Rasterizing AoEs for general overlap ----

# Files for each spatial scale: 1, 2, 3
ndm_listgrid.files("./data/", degree = 3)
anames <- c(anames1dg, anames2dg, anames3dg)
afiles <- c(files1dg, files2dg, files3dg)

# Rasterize all AoEs
# Repeat for each scale: 
# 1, files1dg, anames1dg
# 2, files2dg, anames2dg
# 3, files3dg, anames3dg
purrr::pmap(list(file = files1dg,
                 degree = 1,
                 save.area.shp = FALSE,
                 save.raster = TRUE,
                 save.directory = "./output/AOE_raster/",
                 area.shp.name = anames1dg,
                 env_ras = "./data/AoE_raster/bio_1.tif",
                 region_shp_dns = "./data/GIS/Amer_land.shp",
                 layer_name = "Amer_land"), ndm_rasterize)

# Load AoEs raster layers
raster_files <- list.files(path = "./output/AOE_raster", full.names = TRUE, pattern = "grid")
ras_list <- lapply(raster_files, raster)

# Adjust resolution to 1 degree
ras_list <- lapply(ras_list, function(x) {aggregate(x, fact = 1/res(x))})

# Stack raster to compute the sum
AoE_stack <- stack(ras_list)

AoE_1dg_sum_ras <- sum(AoE_stack)  
plot(AoE_1dg_sum_ras)

# Plot of general overlap

# Preparing raster for ggplot2
# https://stackoverflow.com/questions/33227182/how-to-set-use-ggplot2-to-map-a-raster

test_spdf <- as(AoE_1dg_sum_ras, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

#test_df$value[which(test_df$value == 0)] <- "NA"
#(is.na(test_df$value))

# America croquis
Amer <- st_read(dsn = "./data/GIS/", layer = "Amer_land")
bigbbox <- extent(c(-112, -31, -40, 40))
Amer <- crop_shape(Amer, y = bigbbox)

# plot: Base map

overlap_plot <- ggplot() + 
  geom_tile(data = test_df, aes(x = x, y = y, fill = value)) +
  #scale_fill_viridis_b(breaks = c(1, 2, 5, 10, 15), alpha = 0.5) +
  scale_fill_viridis(alpha = 0.5) +
  geom_tile( data = subset(test_df, value %in% c(0)), aes(x = x, y = y, fill = value), fill = "white") +
  geom_sf(data = Amer, fill = NA, lwd = 0.5, color = "black") +
  scale_x_continuous(breaks = seq(-120, -30, 10)) +
  scale_y_continuous(breaks = seq(-40, 40, 10)) +
  theme_classic() +
  theme(panel.background = element_blank(),
        legend.position = c(0.17, 0.2),
        legend.background = element_rect(color = "black")) +
  labs(x = "Longitude", y = "Latitude",
       fill = "Overlapping\nAreas of Endemism\nper 1º cell")

# Save plot
ggsave(plot = overlap_plot, 
       filename = paste0("./figs/AoEOverlap_plot.pdf"), 
       device = "pdf",
       height = 6, width = 6)

# 3.2. Plotting overlap and cell aggregation ----

aggregation_list <- list(
  two = c(2, 4, 6, 8, 10, 12, 14, 16, 18),
  three = c(3, 6, 9, 12, 15, 18),
  four = c(4, 8, 12, 16),
  five = c(5, 10, 15),
  six = c(6, 12, 18),
  seven = c(7, 14),
  eight = c(8, 16),
  nine = c(9, 18),
  ten = c(10)
)

plot_overlap <- function(aggreg_vec) {
  
  ggplot() +
    geom_tile(data = test_df, aes(x = x, y = y, fill = value)) +
    scale_fill_viridis_b(breaks = aggreg_vec, alpha = 0.5) +
    geom_tile( data = subset(test_df, value %in% c(0)), aes(x = x, y = y, fill = value), fill = "white") +
    geom_sf(data = Amer, fill = NA, lwd = 0.5, color = "black") +
    scale_x_continuous(breaks = seq(-120, -30, 10)) +
    scale_y_continuous(breaks = seq(-40, 40, 10)) +
    theme_classic() +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          #legend.position = "none",
          legend.position = c(0.8, 0.8),
          #legend.direction = "horizontal",
          legend.background = element_rect(color = "black")) +
    labs(fill = "Overlapping\nAoEs")
  
  
}

overlap_plots_list <- lapply(aggregation_list, plot_overlap)

overlap_plots_list[["emptyplot"]]  <- ggplot() +
  geom_sf(data = Amer, fill = NA, lwd = 0.5, color = "black") +
  scale_x_continuous(breaks = seq(-120, -30, 10)) +
  scale_y_continuous(breaks = seq(-40, 40, 10)) +
  theme_classic()


for (i in seq_along(overlap_plots_list)) {
  
  ggsave(plot = overlap_plots_list[[i]], 
         filename = paste0("./figs/AoEOverlap_", names(overlap_plots_list[i]), ".tif"), 
         device = "tiff",
         height = 6, width = 6)
  
}

for (i in seq_along(overlap_plots_list)) {
  
  ggsave(plot = overlap_plots_list[[i]], 
         filename = paste0("./figs/AoEOverlap_", names(overlap_plots_list[i]), ".pdf"), 
         device = "pdf",
         height = 6, width = 6)
  
}


# 4. Comparing consensus areas against previous biogeographical regionalizations ----

# NOTE_14: Generally, areas of endemism are compared to the biogeographical  
# units of regionalization schemes by visual inspection of maps. The idea here
# is to do it quantitatively. This was inspired in the mapcurves method 
# implemented in the sabre R package to compare categorical maps (Nowosad and 
# Stepinski, 2018). However, our taks here is simpler as it implies to compare
# an individual consensus area against a biogeographical unit instead of complete
# regionalizations. The main diference is that mapcurves requires that the 
# spatial objects under comparison have the same extension, a condition 
# difficult to meet here. Therefore, I created the function ndm_contrast() 
# to adapt the general logic of the comparison made in mapcurves to individual
# areas and complete regionalizations with different extents. The idea is to 
# calculate how much of the consensus area is covered by the different
# biogeographical units in a regionalisation scheme, and conversely how much of 
# the biogeographical unit is covered by the consensus area. When two polygons
# are compared their different shapes make that the proportion of the area of 
# intersection between them be different for each one of them. So, considering 
# these proportional areas of intersection can give us an idea of how well the 
# two polygons fit into each other. Nowosad and Stepinski (2018) did this using
# more sophisticated methods and named these measures as Completeness and 
# Homogeneity. To develop the function ndm_contrast() I used those same terms, 
# but for the writing of the manuscript of Patterns of Endemism of Bignonieae I
# used the terms Uniformity of Consensus Area and Uniformity of the 
# Biogeographical unit to avoid unnecesary associations with their method. 
# Therefore, the dataframe will have the columns Completeness and Homogeneity 
# that stand for Uniformity of the Biogeographical unit (Ub in manuscript) and 
# Uniformity of Consensus area (Uc in manuscript), respectively. I took the average 
# of Uc and Ub to represent the spatial congruence (Sc in manuscript) between the 
# Consensus area and the biogeographical unit

# NOTE: The ndm_contrast() function requires to use a shapefile for the 
# biogeographical regionalisation, and the shapefile of the consensus area. 

# 4.1. Calculating the Uniformity of consensus areas and biogeographical units ----

# Consensus areas directory
conA_dsn <- "./data/GIS/Consensus_Areas_Shapefiles/"

# Consensus areas layer names for the Default analysis
# (70 consensus areas at 3 spatial scales)
ca_layer_names <- list.files(path = "./data/GIS/Consensus_Areas_Shapefiles/", 
                             pattern = "[1-3]dgD.*grid.shp") %>%
  gsub(pattern = "([1-3]dgD.*grid).shp", replacement = "\\1")

# Morrone Regionalization (Lowenberg-Neto, 2014)
Morrone <- "./data/GIS/Lowenberg_Neto_2014.shp"

# Gentry Regionalization (JPNG)
Gentry <- "./data/GIS/Gentry_Phytogeography.shp"

# Index to identify the consensus area to be compared (Change it!)
i <- 11 # 11 out of 70

# Example with Morrone
(ex_Morrone <- ndm_contrast(ca_dsn = conA_dsn, ca_layer = ca_layer_names[i], reg_file = Morrone))

morr_df <- map(.x = ca_layer_names,
               .f = ndm_contrast,
               ca_dsn = conA_dsn,
               reg_file = Morrone)

morr_df <- as.data.frame(do.call("rbind", morr_df))
colnames(morr_df) <- c("Consensus_area", "Bioregion", "Ub", "Uc")
#write.csv(morr_df, "./output/COMP_Morr.csv")

# Example with Gentry
(ex_Gentry <- ndm_contrast(ca_dsn = conA_dsn, ca_layer = ca_layer_names[i], reg_file = Gentry))

gen_df <- map(.x = ca_layer_names,
              .f = ndm_contrast,
              ca_dsn = conA_dsn,
              reg_file = Gentry)

gen_df <- as.data.frame(do.call("rbind", gen_df))
colnames(gen_df) <- c("Consensus_area", "Bioregion", "Ub", "Uc")
#write.csv(gen_df, "./output/COMP_Gen.csv")

# 4.2. Visualizing Uniformity and Spatial congruence ----
# 4.2.1. Calculating spatial congruence. ----

# Prepearing dataframes: Morrone
# Loading the dataframe produced by ndm_contrast()

MOR_CONTRAST <- read.csv("./data/COMP_Morr.csv", 
                         stringsAsFactors = FALSE) %>%
  mutate(Scale = ifelse(grepl("D1.*", Consensus_area), 1, 
                        ifelse(grepl("D2.*", Consensus_area), 2, 3)),
         Hierarchy = ifelse(grepl(".*province", Bioregion), "Province", 
                            ifelse(grepl(".*subregion", Bioregion), 
                                   "Subregion", "Dominion"))) %>% 
  group_by(Consensus_area) %>%
  mutate(Intersections = n(),
         Congruence = (as.numeric(Ub) + as.numeric(Uc)) / 2) %>%
  group_by(Consensus_area, Hierarchy) %>%
  mutate(Intersect_by = n())

head(MOR_CONTRAST)

# Prepearing dataframes: Gentry
GEN_CONTRAST <- read.csv("./data/COMP_Gen.csv", 
                         stringsAsFactors = FALSE) %>%
  mutate(Scale = ifelse(grepl("D1.*", Consensus_area), 1, 
                        ifelse(grepl("D2.*", Consensus_area), 2, 3)),
         Hierarchy = "Phytogeographical region") %>% 
  group_by(Consensus_area) %>%
  mutate(Intersections = n(),
         Congruence = (Ub + Uc)/2)


# Example of Comparison vs Congruence plot segregated by biogeographical units
{
  # GEN_CONTRAST and MOR_CONTRAST were created in 4.2.1
  
  # Joining dataframes of Morrone and Gentry
  GENMOR <- rbind(GEN_CONTRAST, MOR_CONTRAST) %>%
    group_by(Hierarchy) %>%
    mutate(med = median(Congruence),
           Hier_name = case_when(
             Hierarchy == "Dominion" ~ "Dominion (Morrone 2014)",
             Hierarchy == "Province" ~ "Province (Morrone 2014)",
             Hierarchy == "Subregion" ~ "Subregion (Morrone 2014)",
             Hierarchy == "Phytogeographical region" ~ "Phytogeographical region (Gentry 1982)"
           ))
  
  GENMOR$Hier_name <- factor(
    GENMOR$Hier_name, 
    levels = c("Province (Morrone 2014)", "Dominion (Morrone 2014)",
               "Subregion (Morrone 2014)", "Phytogeographical region (Gentry 1982)")
  )
  
  compcong_plot <-  GENMOR %>%
    ggplot(aes(x = Congruence), group = Hier_name) +
    geom_histogram(color = "white", alpha = 0.6, bins = 10) +
    scale_x_continuous(breaks = seq(0, 100, 10)) +
    scale_y_continuous(breaks = seq(0, 200, 25)) +
    coord_cartesian(xlim = c(0, 100), expand = TRUE) +
    facet_wrap(~Hier_name) +
    geom_vline(aes(xintercept = med, group = Hier_name), 
               colour = "red", lty = 2) +
    labs(y = "Pair-wise comparisons:\nAreas of Endemism vs Biogeographical units",
         x = "Spatial congruence") +
    # coord_flip() +
    theme_classic() +
    theme(strip.background = element_rect(fill = "#D4D4D4", colour = NA),
          strip.text = element_text(face = "bold", colour = "black"),
          text = element_text(size = 7))
  
  compcong_plot
  }
