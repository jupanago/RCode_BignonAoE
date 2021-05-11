## NDM_tables

## Summarizing consensus areas information

################ ndm_casumm function ################

## This function creates a dataframe summarizing consensus area information for each analysis and scale.
## It gives the number of individual areas, species, groups, consensus areas, and ambiguous species per analysis and scale.

## Arguments:
## 1) x: A dataframe obtained from ndm_dataparsing function. 

ndm_casum <- function(x) {
  
  library(dplyr)
  
  ## Counting the number of areas formed by genera and phylogroups per analysis
  Grps <- x %>%
    group_by(Analysis, Spatial_Scale) %>%
    summarise(Groups_number = sum(Groups))
  
  Summary_info <- x %>%
    distinct(Analysis, Spatial_Scale, Consensus_Area, .keep_all = TRUE) %>%
    group_by(Analysis, Spatial_Scale) %>%
    # filter(Analysis %in% c("H", "S", "T")) %>%
    summarise(Indiv_Areas = sum(Individual_Areas),
              Spp = sum(Number_Spp),
              Groups = sum(Groups),
              Cons_Areas = n())
  
  Spp_uniq_number <- x %>%
    group_by(Analysis, Spatial_Scale) %>%
    summarise(Spp = length(unique(Species)))
  
  Summary_info$Spp <- Spp_uniq_number$Spp
  
  Summary_info$Groups[3] <- Grps$Groups_number[3]
  Summary_info$Spp[3] <- Summary_info$Spp[3] - Summary_info$Groups[3]
  
  # Summary_info$Groups[6] <- Grps$Groups_number[6]
  # Summary_info$Spp[6] <- Summary_info$Spp[6] - Summary_info$Groups[6]
  
  ### Ambiguous consensus area data
  
  CA_amb <- function(x, scl, ana) {
    
    ## Function to show consensus areas with ambiguous species for an specified analysis.
    ## Ambiguous species are those that appear as part of several consensus areas, therefore, they are duplicated records.
    ## x is a data.frame where the species part of consensus areas are in column named 'Species'
    ## spscale and analys are arguments for my own areas and analysis nomenclature. They indicate scale and analysis, respectively.
    
    require(dplyr)
    
    y <- x %>%
      filter(Spatial_Scale == scl, Analysis == ana) %>%
      dplyr::select(Species)
    
    index <- which(duplicated(y$Species))
    spp_vec <- unique(y$Species[index])
    
    
    Con_Areas_amb <- x %>% 
      filter(Spatial_Scale == scl, Analysis == ana) %>% 
      filter(Species %in% spp_vec) %>% 
      pull(Consensus_Area) %>% 
      unique() %>% 
      length()
    
    Con_Areas_amb
    
  }
  
  CAambD1 <- CA_amb(x, scl = "1", ana = "D")
  CAambD2 <- CA_amb(x, scl = "2", ana = "D")
  CAambD3 <- CA_amb(x, scl = "3", ana = "D")
  
  # CAambS1 <- CA_amb(x, scl = "1", ana = "S")
  # CAambS2 <- CA_amb(x, scl = "2", ana = "S")
  # CAambS3 <- CA_amb(x, scl = "3", ana = "S")
  # 
  # CAamb_vec <- c(CAambD1, CAambD2, CAambD3, CAambS1, CAambS2, CAambS3)
  CAamb_vec <- c(CAambD1, CAambD2, CAambD3)
  
  Summary_info$Cons_Areas_amb <- CAamb_vec
  
  ### Ambiguous species data
  
  spp_conflict <- function(x, scl, ana) {
    
    ## Function to show ambiguous species for an specified analysis.
    ## Ambiguous species are those that appear as part of several consensus areas, therefore, they are duplicated records.
    ## x is a data.frame where the species part of consensus areas are in column named 'Species'
    ## spscale and analys are arguments for my own areas and analysis nomenclature. They indicate scale and analysis, respectively.
    
    require(dplyr)
    
    y <- x %>%
      filter(Spatial_Scale == scl, Analysis == ana) %>%
      dplyr::select(Species)
    
    index <- which(duplicated(y$Species))
    unique(y$Species[index])
    
  }
  
  SppD1 <- spp_conflict(x, scl = "1", ana = "D")
  SppD2 <- spp_conflict(x, scl = "2", ana = "D")
  SppD3 <- spp_conflict(x, scl = "3", ana = "D")
  
  # SppS1 <- spp_conflict(x, scl = "1", ana = "S")
  # SppS2 <- spp_conflict(x, scl = "2", ana = "S")
  # SppS3 <- spp_conflict(x, scl = "3", ana = "S")
  
  # conflict_list <- list(SppD1, SppD2, SppD3, SppS1, SppS2, SppS3)
  conflict_list <- list(SppD1, SppD2, SppD3)
  
  Summary_info$spp_amb <- sapply(conflict_list, length, simplify = "array")
  
  return(Summary_info)
  
}