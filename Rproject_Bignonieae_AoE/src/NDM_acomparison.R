## Functions to compare consensus areas

############### Function ndm_contrast ##########

## This function compares consensus areas versus a specified biogeographical regionalization.
## Two measures are given:
## 1) Biogeographical region homogeneity (Hb): It refers to the proportion of the consensus area covering a particular bioregion. It is calculated by
##    calculating the intersecting area between de consensus area and the bioregion, divided by the area of the bioregion, and 
##    multiplied by 100. The value expresses the percentage of the bioregion covered by the specified consensus area.
##
## 2) Consensus area homogeneity (Hc): It refers to the proportion of the bioregion convering the consensus area. It expresses how much of the consensus area
##    is inside the bioregion. The greater the value, the greater the proportion covered by one bioregion; the lesser the value, the greater
##    the proportion covered by several bioregions. It is calculated by estimating the intersecting area between de consensus area and the
##    bioregion, divided by the area of the consensus area, divided by the area of the bioregion, and multiplied by 100.
##
## NOTE: This is inspired in the mapcurves method implemented in the sabre R package (Nowosad and Stepinski. 2018). However, our taks here is
## to compare an individual area against a regionalization and not comparing two regionalizations.

## Arguments:
## ca_dsn = directory in which the consensus area layer is saved
## ca_layer = name of consensus area layer
## reg_file = file path to the shape file with the regionalization. Attribute table must contain a column called "Regions" where the names
## of biogeographical areas area stores. If the shapefile for Morrone corresponds to Lowenberg-Neto et al. 2014, the named colums
## already included in this file are appropriate for this function. 

ndm_contrast <- function(ca_dsn, ca_layer, reg_file, reg_layer) {
  
  library(dplyr)
  library(rgdal)
  library(sf)
  library(rmapshaper)
  
  wgs84.proj4 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  ca_pol <- readOGR(dsn = ca_dsn, layer = ca_layer, p4s = wgs84.proj4)
  ca_pol <- st_as_sf(ca_pol)
  ca_pol <- st_union(ca_pol)
  
  # Selecting name
  if (grepl(ca_layer, pattern = ".*1dgD.*")) {
    x <- gsub(ca_layer, pattern = ".*1dgD_(.*)_.*_grid", replacement = "\\1")
    y <- paste("D1", x, sep = "")
  }
  
  if (grepl(ca_layer, pattern = ".*2dgD.*")) {
    x <- gsub(ca_layer, pattern = ".*2dgD_(.*)_.*_grid", replacement = "\\1")
    y <- paste("D2", x, sep = "")
  }
  
  if (grepl(ca_layer, pattern = ".*3dgD.*")) {
    x <- gsub(ca_layer, pattern = ".*3dgD_(.*)_.*_grid", replacement = "\\1")
    y <- paste("D3", x, sep = "")
  }
  
  # Creating dataframe with area name
  ca_pol <- st_sf(name = y, ca_pol)
  
  # Reading biogeographic regionalization shapefile
  Reg <- st_read(reg_file)
  Reg <- st_transform(Reg, crs = st_crs(ca_pol, asText = TRUE))
  
  # Obtaining the intersection of the consensus area with each biogeographical area
  CA_and_Reg <- st_intersection(ca_pol, Reg)
  
  # Extracting names of the biogeographical areas. Columns according to Morrone or using the column "Region" from attribute table in shapefile
  Neo_region <- unique(CA_and_Reg$REGION)
  Neo_region <- Neo_region[!is.na(Neo_region)]
  subreg <- unique(CA_and_Reg$Subregio_1) 
  subreg <- subreg[!is.na(subreg)]
  domini <- unique(CA_and_Reg$Dominions)
  domini <- domini[!is.na(domini)]
  provin <- unique(CA_and_Reg$Province_1)
  provin <- provin[!is.na(provin)]
  gen_region <- unique(CA_and_Reg$Region)
  gen_region <- gen_region[!is.na(gen_region)]
  # Calculating consensus area homogeneity and bioregion homogeneity 
  
  # Neotropical region in Morrone's proposal
  vec_hNReg <- vector(mode = "numeric", length = length(Neo_region))
  vec_cNReg <- vector(mode = "numeric", length = length(Neo_region))
  for (i in seq_along(Neo_region)) {
    
    a_inter_NReg <- CA_and_Reg %>%
      filter(REGION == Neo_region[i])  %>%
      st_union() %>%
      st_area()
    
    a_total_NReg <- Reg %>%
      filter(REGION == Neo_region[i]) %>%
      st_union() %>%
      st_area()  
    
    Homo_NReg <- round(a_inter_NReg / st_area(ca_pol) * 100, digits = 4) 
    Comp_NReg <- round(a_inter_NReg / a_total_NReg * 100, digits = 4) 
    
    vec_hNReg[i] <- Homo_NReg 
    vec_cNReg[i] <- Comp_NReg 
    
  }
  NReg_Comp_Hom_df <- cbind(Bioregion = as.character(Neo_region), Hb = vec_cNReg, Hc = vec_hNReg)
  
  # Subregions in Morrone's proposal
  vec_hsubr <- vector(mode = "numeric", length = length(subreg))
  vec_csubr <- vector(mode = "numeric", length = length(subreg))
  for (i in seq_along(subreg)) {
    
    a_inter_subr <- CA_and_Reg %>%
      filter(Subregio_1 == subreg[i])  %>%
      st_union() %>%
      st_area()
    
    a_total_subr <- Reg %>%
      filter(Subregio_1 == subreg[i]) %>%
      st_union() %>%
      st_area()  
    
    Homo_subr <- round(a_inter_subr / st_area(ca_pol) * 100, digits = 4) 
    Comp_subr <- round(a_inter_subr / a_total_subr * 100, digits = 4) 
    
    vec_hsubr[i] <- Homo_subr 
    vec_csubr[i] <- Comp_subr 
    
  }
  subr_Comp_Hom_df <- cbind(Bioregion = as.character(subreg), Hb = vec_csubr, Hc = vec_hsubr)
  
  # Dominions in Morrone's proposal
  vec_hdomi <- vector(mode = "numeric", length = length(domini))
  vec_cdomi <- vector(mode = "numeric", length = length(domini))
  for (i in seq_along(domini)) {
    
    a_inter_domi <- CA_and_Reg %>%
      filter(Dominions == domini[i])  %>%
      st_union() %>%
      st_area()
    
    a_total_domi <- Reg %>%
      filter(Dominions == domini[i]) %>%
      st_union() %>%
      st_area()  
    
    Homo_domi <- round(a_inter_domi / st_area(ca_pol) * 100, digits = 4) 
    Comp_domi <- round(a_inter_domi / a_total_domi * 100, digits = 4) 
    
    vec_hdomi[i] <- Homo_domi
    vec_cdomi[i] <- Comp_domi
    
  }
  domi_Comp_Hom_df <- cbind(Bioregion = as.character(domini), Hb = vec_cdomi, Hc = vec_hdomi)
  
  # Provinces in Morrone's proposal
  vec_hprov <- vector(mode = "numeric", length = length(provin))
  vec_cprov <- vector(mode = "numeric", length = length(provin))
  for (i in seq_along(provin)) {
    
    a_inter_prov <- CA_and_Reg %>%
      filter(Province_1 == provin[i])  %>%
      st_union() %>%
      st_area()
    
    a_total_prov <- Reg %>%
      filter(Province_1 == provin[i]) %>%
      st_union() %>%
      st_area()  
    
    Homo_prov <- round(a_inter_prov / st_area(ca_pol) * 100, digits = 4) 
    Comp_prov <- round(a_inter_prov / a_total_prov * 100, digits = 4) 
    
    vec_hprov[i] <- Homo_prov
    vec_cprov[i] <- Comp_prov
    
  }
  prov_Comp_Hom_df <- cbind(Bioregion = as.character(provin), Hb = vec_cprov, Hc = vec_hprov)
  
  # Region in other biogeographical proposals: Gentry in this case
  vec_hgentry <- vector(mode = "numeric", length = length(gen_region))
  vec_cgentry <- vector(mode = "numeric", length = length(gen_region))
  for (i in seq_along(gen_region)) {
    
    a_inter_gen <- CA_and_Reg %>%
      filter(Region == gen_region[i])  %>%
      st_union() %>%
      st_area()
    
    a_total_gen <- Reg %>%
      filter(Region == gen_region[i]) %>%
      st_union() %>%
      st_area()  
    
    Homo_gentry <- round(a_inter_gen / st_area(ca_pol) * 100, digits = 4) 
    Comp_gentry <- round(a_inter_gen / a_total_gen * 100, digits = 4) 
    
    vec_hgentry[i] <- Homo_gentry
    vec_cgentry[i] <- Comp_gentry
    
  }
  gentry_Comp_Hom_df <- cbind(Bioregion = as.character(gen_region), Hb = vec_cgentry, Hc = vec_hgentry)
  
  # Creating a dataframe for all the comparisons
  Comp_Hom_df <- rbind(NReg_Comp_Hom_df, subr_Comp_Hom_df, domi_Comp_Hom_df, prov_Comp_Hom_df, gentry_Comp_Hom_df)
  Comp_Hom_df <- cbind(Consensus_area = y, Comp_Hom_df)
  
  # Some fields in Morrone's shapefile attribute table are in blank for the columns specified above.
  # This is so because some regions are not classified under less inclusive bioeographical areas such as dominions.
  Comp_Hom_df <- Comp_Hom_df[which(complete.cases(Comp_Hom_df)), ]
  
  return(Comp_Hom_df)
  
}

