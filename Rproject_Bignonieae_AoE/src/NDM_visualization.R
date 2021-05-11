### ndm_visualization: Script with functions to visualize ndm outputs.

################ Function ndm_plotca4 ################

# Funtion to make a basic plot of the area of endemism

# Arguments:
# 1) map_list_element: named list of consensus areas 
# 2) dsn_aoes; directory where the shapefile layers are saved 
# 4) aoe_layers: names of the shapefile layers
# 3) area_name: area name, 
# 5) area_col: a color to fill the area, 
# 6) alpha_aoe = 0.1: transparency of area
# 7) lwd_aoe = 0.2: linewidth of area polygon

ndm_plotca4 <- function(map_list_element, dsn_aoes, aoe_layers, 
                        area_name, area_col, alpha_aoe = 0.1, lwd_aoe = 0.2) {
  
  # libraries
  library(sf)
  library(ggplot2)
  library(stringr)
  library(tmaptools)
  library(raster)
  
  # America croquis
  Amer <- st_read(dsn = "./data/GIS/", layer = "Amer_land")
  bigbbox <- extent(c(-113, -34, -41, 39))
  Amer <- crop_shape(Amer, y = bigbbox)
  
  # Base map
  base_map <- ggplot() + 
    geom_sf(data = Amer, fill = NA, lwd = 0.1) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank())
  
  # Preparing aoe plots
  length_areas_vec <- length(map_list_element[[1]])
  areas_name_vec <- map_list_element[[1]]
  layers_vec <- aoe_layers[names(aoe_layers) %in% map_list_element[[1]]]
  
  jpng_palette <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF") # viridis(4)
  
  if (length_areas_vec == 1) {
    
    # Layer
    ca_grid_1 <- st_read(dsn = dsn_aoes, layer = layers_vec[1]) %>% st_union()
    st_crs(ca_grid_1) <- st_crs(Amer)
    
    
    # Names from layers
    name_from_layer <- paste(str_replace_all(layers_vec, 
                                             pattern = "^([1-3])dgD_(.*)_.*_grid", 
                                             replacement = "D\\1\\2"), 
                             collapse = "_")
    
    # Areas to plot
    areas_plot <- base_map +
      geom_sf(data = ca_grid_1, fill = area_col, color = area_col, 
              alpha = alpha_aoe, lwd = lwd_aoe) +
      theme(title = element_text(face = "bold", size = 12)) +
      labs(title = area_name)
    
    
  }
  
  if(paste(areas_name_vec, collapse = "_") != name_from_layer) {
    
    warning("Area names from AoE list and layers are not the same. Check vectors order. Areas in plot and file name might not coincide.")
    
  }
  
  areas_plot
  
}
