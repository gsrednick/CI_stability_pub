# 5 - Network synchrony and stability calculations ###

# packages
library(tidyverse)
library(vegan)
library(phylin)
library(codyn)
library(pbapply)
library(parallel)

# 5A - MDS for pred-prey communities ####

pisco_troph_complete<-read.csv("./Manuscript_scripts/Exported_data/pisco_troph_complete.csv")

length(unique(pisco_troph_complete$site))

# do this only for links that have enough spatial replication
filter_links_df<-read.csv("./Data/PISCO/Trophic/trophic_network.csv")
networks_of_interest<-read.csv("./Manuscript_scripts/trophic/links_of_interest.csv")


filter_links_df_reduced<-filter_links_df %>% filter(network_ID %in% networks_of_interest$network_ID)


adonis_network_list<- split(filter_links_df_reduced, filter_links_df_reduced$network_ID)


core_number = detectCores()

pred_prey_mds<-pblapply(adonis_network_list,function(x){ #this function works and can be applied to list items with lapply
  #tryCatch({
  
  pred_com<-pisco_troph_complete %>% 
    dplyr::filter(species %in% x$consumer_ID) %>% 
    dplyr::select(year, species, site,cover) %>%
    spread(species, cover, fill = 0)
  
  prey_com<-pisco_troph_complete %>% 
    dplyr::filter(species %in% x$resource_ID) %>% 
    dplyr::select(year, species, site,cover) %>%
    spread(species, cover, fill = 0)
  
  pred_prey_com<-merge(pred_com,prey_com)
  
  pisco_env_net<-pred_prey_com[c(1:2)]
  pisco_com_net<-pred_prey_com[-c(1:2)]
  
  # make a dummy variable
  pisco_com_net$dummy<-1
  
  com_mds_net <- metaMDS(comm = pisco_com_net, distance = "bray", trace = F, autotransform = F, na.rm = FALSE)
  com_mds_net$stress # 
  
  com_mds_net_points<-data.frame(com_mds_net$points)
  mds_net<-merge(pisco_env_net,com_mds_net_points, by="row.names", all.x=TRUE)
  
  mds_net$network_ID<-unique(x$network_ID)
  
  
  return(mds_net)
  #}, error = function(e) return(NULL))
  
}#, mc.cores = core_number, mc.preschedule = FALSE)
)

pred_prey_mds_df<-do.call(rbind.data.frame, pred_prey_mds)

# write this for making figures later
write.csv(pred_prey_mds_df,"./Manuscript_scripts/Exported_data/mds_points_community.csv", row.names = F)


pred_prey_disper<-pblapply(adonis_network_list,function(x){ #this function works and can be applied to list items with lapply
  #tryCatch({
  
  pred_com<-pisco_troph_complete %>% 
    dplyr::filter(species %in% x$consumer_ID) %>% 
    dplyr::select(year, species, site,cover) %>%
    spread(species, cover, fill = 0)
  
  prey_com<-pisco_troph_complete %>% 
    dplyr::filter(species %in% x$resource_ID) %>% 
    dplyr::select(year, species, site,cover) %>%
    spread(species, cover, fill = 0)
  
  pred_prey_com<-merge(pred_com,prey_com)
  
  pisco_env_net<-pred_prey_com[c(1:2)]
  pisco_com_net<-pred_prey_com[-c(1:2)]
  
  # make a dummy variable
  pisco_com_net<-pisco_com_net
  pisco_com_net$dummy<-1
  
  # PERMANOVA -- results not used in formal analyses
  # com_mds_net <- metaMDS(comm = pisco_com_net, distance = "bray", trace = F, autotransform = F, na.rm = FALSE)
  # com_mds_net$stress # 
  
  vegdist_net<-vegdist(pisco_com_net, upper=F)
  
  disper_res<-betadisper(vegdist_net, pisco_env_net$site, type = c("centroid"))
  
  disper_res_df<-data.frame(distance = disper_res$distances, site = disper_res$group)
  
  disper_res_df_sum<-disper_res_df %>% group_by(site) %>% summarize_all(mean,na.rm = T)
  
  disper_res_df_sum$network_ID<-unique(x$network_ID)
  
  
  return(disper_res_df_sum)
  #}, error = function(e) return(NULL))
  
})

pred_prey_disper_df<-do.call(rbind.data.frame, pred_prey_disper)
dim(pred_prey_disper_df) # should be 1932 = 42 sites x 46 networks

disper_df<-merge(pred_prey_disper_df,site_table)

disper_df$distance_inv<- (disper_df$distance)*(-1) # Covert dispersion so that its more interpretable

ggplot(disper_df,aes(x = Island,y = distance_inv, color = site_status)) +
  geom_boxplot()

summary(aov(distance~site_status*Island,disper_df))


write.csv(disper_df,"./Manuscript_scripts/Exported_data/community_stability.csv",row.names = F)



# 5B - Metacommunity multivariate analyses ####

# aggregate to site_status x island

# make for MPA and ref metacommunities
DBall_meta_MPA_REF<-pisco_troph_complete %>% 
  #filter(species %in% trophic_species$species) %>% # this doesnt filter correctly
  dplyr::select(year, site, Island, site_status,species, cover) %>%
  spread(species, cover, fill = 0) %>%
  unite(island_status, c("Island", "site_status"), remove = F) %>% 
  group_by(year,island_status,Island) %>%
  summarize_if(is.numeric,sum)

# make for whole metacom
DBall_meta_both<-pisco_troph_complete %>% 
  #filter(species %in% trophic_species$species) %>% # this doesnt filter correctly
  dplyr::select(year, site, Island, site_status,species, cover) %>%
  spread(species, cover, fill = 0) %>%
  group_by(year,Island) %>%
  summarize_if(is.numeric,sum) %>%
  mutate(island_status = paste(Island, "all", sep = "_"))

# bring these together
DBall_meta<-rbind(DBall_meta_MPA_REF,DBall_meta_both)

# make long
DBall_meta_long <- DBall_meta %>% 
  pivot_longer(!c(year,island_status,Island), names_to = "species", values_to = "cover")


# split to list
network_list<- split(filter_links_df_reduced, filter_links_df_reduced$network_ID)

pred_prey_meta_mds<-pblapply(network_list,function(x){ #this function works and can be applied to list items with lapply
  #tryCatch({
  
  pred_com<-DBall_meta_long %>% 
    dplyr::filter(species %in% x$consumer_ID) %>% 
    dplyr::select(year, species, island_status,cover) %>%
    spread(species, cover, fill = 0)
  
  prey_com<-DBall_meta_long %>% 
    dplyr::filter(species %in% x$resource_ID) %>% 
    dplyr::select(year, species, island_status,cover) %>%
    spread(species, cover, fill = 0)
  
  pred_prey_com<-merge(pred_com,prey_com)
  
  pisco_env_net<-pred_prey_com[c(1:2)]
  pisco_com_net<-pred_prey_com[-c(1:2)]
  
  # make a dummy variable
  pisco_com_net<-pisco_com_net
  pisco_com_net$dummy<-1
  
  com_mds_net <- metaMDS(comm = pisco_com_net, distance = "bray", trace = F, autotransform = F, na.rm = FALSE)
  com_mds_net$stress # 
  
  com_mds_net_points<-data.frame(com_mds_net$points)
  mds_net<-merge(pisco_env_net,com_mds_net_points, by="row.names", all.x=TRUE)
  
  mds_net$network_ID<-unique(x$network_ID)
  
  
  return(mds_net)
  #}, error = function(e) return(NULL))
  
})

pred_prey_meta_mds_df<-do.call(rbind.data.frame, pred_prey_meta_mds)




pred_prey_meta_disper<-pblapply(network_list,function(x){ #this function works and can be applied to list items with lapply
  #tryCatch({
  
  pred_com<-DBall_meta_long %>% 
    dplyr::filter(species %in% x$consumer_ID) %>% 
    dplyr::select(year, species, island_status, cover) %>%
    spread(species, cover, fill = 0)
  
  prey_com<-DBall_meta_long %>% 
    dplyr::filter(species %in% x$resource_ID) %>% 
    dplyr::select(year, species, island_status, cover) %>%
    spread(species, cover, fill = 0)
  
  pred_prey_com<-merge(pred_com,prey_com)
  
  pisco_env_net<-pred_prey_com[c(1:2)]
  pisco_com_net<-pred_prey_com[-c(1:2)]
  
  # make a dummy variable
  pisco_com_net<-pisco_com_net
  pisco_com_net$dummy<-1
  
  # PERMANOVA -- results not used in formal analyses
  #com_mds_net <- metaMDS(comm = pisco_com_net, distance = "bray", trace = F, autotransform = F, na.rm = FALSE)
  #com_mds_net$stress # 
  
  vegdist_net<-vegdist(pisco_com_net, upper=F)
  
  disper_res<-betadisper(vegdist_net, pisco_env_net$island_status, type = c("centroid"))
  
  disper_res_df<-data.frame(distance = disper_res$distances, island_status = disper_res$group)
  
  disper_res_df_sum<-disper_res_df %>% group_by(island_status) %>% summarize_all(mean,na.rm = T)
  
  disper_res_df_sum$network_ID<-unique(x$network_ID)
  
  
  return(disper_res_df_sum)
  #}, error = function(e) return(NULL))
  
})

pred_prey_meta_disper_df<-do.call(rbind.data.frame, pred_prey_meta_disper)
dim(pred_prey_meta_disper_df) # should be 240 = 4 islands x 3 status x 20 networks (usable)

pred_prey_meta_disper_df<-pred_prey_meta_disper_df %>% 
  separate(island_status, c("Island", "site_status"), sep = "_")

pred_prey_meta_disper_df$distance_inv<- (pred_prey_meta_disper_df$distance)*(-1)

write.csv(pred_prey_meta_disper_df,"./Manuscript_scripts/Exported_data/metacommunity_stability.csv",row.names = F)



# 5C - Synchrony calculations ####

trophic_network<-filter_links_df_reduced
## 5C.1 - Community synchrony ####
cons_network_list<- split(trophic_network, trophic_network$network_ID)

sync_by_site<-lapply(cons_network_list,function(x){
  
  com_sync<-pisco_troph_complete %>% 
    filter(species %in% c(x$consumer_ID,x$resource_ID)) %>%
    synchrony(.,
              time.var = "year",
              species.var = "species",
              abundance.var = "cover",
              replicate.var = "site")
  
  com_sync$network_ID<-unique(x$network_ID)
  
  return(com_sync)
})

com_sync_df<-do.call(rbind.data.frame, sync_by_site)

write.csv(com_sync_df,"./Manuscript_scripts/Exported_data/community_synchrony.csv",row.names = F)




## 5C.2 - Metacommunity synchrony ####


sync_by_meta<-lapply(cons_network_list,function(x){
  
  meta_sync<-DBall_meta_long %>% 
    filter(species %in% c(x$consumer_ID,x$resource_ID)) %>%
    synchrony(.,
              time.var = "year",
              species.var = "species",
              abundance.var = "cover",
              replicate.var = "island_status")
  
  meta_sync$network_ID<-unique(x$network_ID)
  
  return(meta_sync)
})

meta_sync_df<-do.call(rbind.data.frame, sync_by_meta)

meta_sync_df<-meta_sync_df %>% 
  separate(island_status,
           into=c("Island", "site_status"),
           sep="_", convert = TRUE, extra = "merge")

write.csv(meta_sync_df,"./Manuscript_scripts/Exported_data/metacommunity_synchrony.csv",row.names = F)


## END ##

