#### Habitat similarity code ####

library(tidyverse)
library(vegan)


# Get just substratum data -- imported in "1_data_curation.R"

UPC_sub<-UPC %>% filter(category == "SUBSTRATE",
                        campus == "UCSB") %>% 
  group_by(year,site,zone,transect,classcode) %>% 
  dplyr::summarize(cover = mean(pct_cov,na.rm = TRUE))



UPC_sub_wide<-UPC_sub %>% pivot_wider(names_from = classcode, 
                                      values_from = cover,
                                      values_fill = 0)

# site level
UPC_sub_wide_site_level <- UPC_sub_wide %>% 
  filter(site %in% site_table$site,
        year %in% pisco_troph_complete$year) %>%
  dplyr::group_by(site) %>% 
  dplyr::summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) 


UPC_sub_wide_site_level$transect<-NULL
UPC_sub_wide_site_level$year<-NULL


# run PCA
habitat_pca <- prcomp(UPC_sub_wide_site_level[,c(2:5)], center = TRUE,scale. = F)
summary(habitat_pca)

# extract loadings 
habitat_loadings <- as.data.frame(habitat_pca$x[,1:2])
habitat_loadings<-cbind(UPC_sub_wide_site_level[c(1)],habitat_loadings)

#View(habitat_loadings)

write.csv(habitat_loadings,"./Manuscript_scripts/Exported_data/hab_similarity.csv",row.names = F)






### Island level loadings ###
UPC_island_merge<-merge(UPC_sub_wide,site_table)

# summarize for MPAs and refs
UPC_sub_wide_island_level <- UPC_island_merge %>% 
  filter(site %in% site_table$site,
         year %in% pisco_troph_complete$year) %>%
  dplyr::group_by(Island,site_status) %>%
  dplyr::summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) 


UPC_sub_wide_island_level$transect<-NULL
UPC_sub_wide_island_level$year<-NULL

# summarize for all sites at an island
UPC_sub_wide_island_level_BOTH <- UPC_island_merge %>% 
  filter(site %in% site_table$site,
         year %in% pisco_troph_complete$year) %>%
  dplyr::group_by(Island) %>% 
  dplyr::summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

UPC_sub_wide_island_level_BOTH$site_status <- "all_sites"
UPC_sub_wide_island_level_BOTH$transect<-NULL
UPC_sub_wide_island_level_BOTH$year<-NULL


# bring back together
sub_island_level<-rbind(UPC_sub_wide_island_level,UPC_sub_wide_island_level_BOTH)



# run PCA
habitat_pca_island <- prcomp(sub_island_level[,c(3:6)], center = TRUE,scale. = F)
summary(habitat_pca_island)

# loadings
island_habitat_loadings <- as.data.frame(habitat_pca_island$x[,1:2])
island_habitat_loadings<-cbind(sub_island_level[c(1,2)],island_habitat_loadings)



write.csv(island_habitat_loadings,"./Manuscript_scripts/Exported_data/island_hab_similarity.csv",row.names = F)


# END # 
