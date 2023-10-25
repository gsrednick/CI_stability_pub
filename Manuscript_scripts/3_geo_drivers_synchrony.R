#### Drivers of geographic variation in community synchrony #####


# Written by: G. Srednick -- 4/14/2022; updated 10/9/2023

# packages
library(tidyverse)
library(sf)
library(ggExtra)
library(patchwork)

# try using older version of "sf" see if the results are different  -- so far its the same.... i have no idea why the values for these three sites have changed...
#library(devtools)
#install_version("sf", version = "1.0-12")

# unhash the following maybe?
#options("sp_evolution_status" = 2) # use sf instead of rgdal and rgeos in sp


# Spatial data

# for SST -- i have to bring this in as Raster and then get SST within 10 x 10 km (or smaller) rectangle centered as closely on a site as possible
# Temperature is MODIST
# SST monthly composite source:  https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1sstdmdayR20190.graph
# SST daily composite source: https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMWsstd1day_LonPM180.html 

#filter sites for ones of interest
sites_2<- sites %>% filter(site %in% site_table$site)

names(sites_2)<-c("site","latitude","longitude")



# work with sf
loc_sf <- sites_2 %>% st_as_sf(coords = c('longitude', 'latitude'), remove = F) # was true
#stop_sf <- california_cur %>% st_as_sf(coords = c('longitude', 'latitude'), remove = F) # was true

print("here's the heavy lifting!")


# attach geometries
stop_sf_d <- california_cur_2 %>% st_as_sf(coords = c('longitude', 'latitude'), remove = F) # daily
stop_sf_m <- california_cur_m_2 %>% st_as_sf(coords = c('longitude', 'latitude'), remove = F) # montly

# Use st_nearest_feature to cbind loc to stop by nearest points
# daily composite
joined_sf_d <- stop_sf_d %>% 
  cbind(
    loc_sf[st_nearest_feature(stop_sf_d, loc_sf),])

length(unique(joined_sf_d$site)) # 40

# monthly composite
joined_sf_m <- stop_sf_m %>% 
  cbind(
    loc_sf[st_nearest_feature(stop_sf_m, loc_sf),])

length(unique(joined_sf_m$site)) # 40

## mutate to add column showing distance between geometries
site_SST_d<-joined_sf_d %>%
  mutate(dist = st_distance(geometry, geometry.1, by_element = T))

site_SST_d_df<-as.data.frame(site_SST_d)

site_SST_m<-joined_sf_m %>%
  mutate(dist = st_distance(geometry, geometry.1, by_element = T))

site_SST_m_df<-as.data.frame(site_SST_m)

# format time; put it in PST
#site_SST_d_df$time_pst<-fastPOSIXct(site_SST_d_df$time)
#site_SST_m_df$time_pst<-fastPOSIXct(site_SST_m_df$time)

head(site_SST_d_df)
head(site_SST_m_df)

# Get cells that are closest to a site for a given date
site_SST_d_df_slice<-site_SST_d_df %>% dplyr::group_by(time_pst,site) %>% slice(which.min(dist))
site_SST_m_df_slice<-site_SST_m_df %>% dplyr::group_by(time_pst,site) %>% slice(which.min(dist))


# add a step -- filter to closest cell for a given site, then filter by cell from the larger dataset
#site_SST_m_df_site<-site_SST_m_df %>% dplyr::group_by(site) %>% slice(which.min(dist))
#site_SST_m_df_updated<-site_SST_m_df %>% filter(geometry %in% site_SST_m_df_site$geometry)

# Extract date

# Daily composite
site_SST_df_annual_MARSS <- site_SST_d_df_slice %>% 
  group_by(site,year) %>%
  dplyr::summarize(mean_SST = mean(sst,na.rm=T)) # mean for MARSS model and for Figure S1

# Monthly composite
site_SST_df_annual_CV <- site_SST_m_df_slice %>% 
  group_by(site,year) %>%
  dplyr::summarize(mean_SST_u = mean(sst,na.rm=T),
                   CV_SST = sd(sst,na.rm=T)/mean_SST_u) 


site_SST_df_site <- site_SST_m_df_slice %>% 
  filter(year %in% (2007:2020),
         !year == 2013) %>%
  group_by(site) %>%
  dplyr::summarize(mean_SST = mean(sst,na.rm=T),
                   CV_SST = sd(sst,na.rm=T)/mean_SST) # CV for SEM and for Figure S1


# plot for temporal coverage
ggplot(data = site_SST_m_df_slice) +
  geom_point(aes(x = year, y = site)) # full replication

# plot to check 
site_SST_df_annual_MARSS %>% 
  ggplot(aes(x= year,y=mean_SST, color = site)) +
  geom_point() +
  theme(legend.position = "none")

site_SST_df_annual_CV %>% 
  ggplot(aes(x= year,y=CV_SST, color = site)) +
  geom_point() +
  theme(legend.position = "none")



total_nas <- sum(is.na(site_SST_df_annual_MARSS)) # No NAs -- good
total_nas

# write annual for use in MARSS
write.csv(site_SST_df_annual_MARSS,"./Data/Env_data/SST_data/annual_SST_MARSS.csv", row.names = F)



### SST into site df ###
site_env_data<-site_SST_df_site

# to CSV for use in SEM
write.csv(site_env_data,"./Data/Env_data/site_env_data.csv",
          row.names = F)


#### Island Level ####

### Temperature 

# for annual
site_SST_d_df_metadata<-merge(site_SST_d_df_slice,site_table)

meta_SST_df_annual_pre<-site_SST_d_df_metadata %>%
  dplyr::filter(year %in% (2007:2019),
       !year == 2013) %>%
  group_by(Island,site_status,year) %>%
  dplyr::summarize(mean_SST = mean(sst,na.rm=T),
                   CV_SST = sd(sst,na.rm=T)/mean_SST)

meta_SST_df_annual_both <- site_SST_d_df_metadata %>% 
  filter(year %in% (2007:2019),
         !year == 2013) %>%
  group_by(Island,year) %>%
  dplyr::summarize(mean_SST = mean(sst,na.rm=T),
                   CV_SST = sd(sst,na.rm=T)/mean_SST)

meta_SST_df_annual_both$site_status <- "all"


meta_SST_df_annual<-rbind(meta_SST_df_annual_pre,meta_SST_df_annual_both) %>%
  unite(island_status, c("Island", "site_status"), remove = F)

write.csv(meta_SST_df_annual,"./Data/Env_data/SST_data/annual_meta_SST_MARSS.csv",row.names = F)





# For aggregate
site_SST_m_df_metadata<-merge(site_SST_m_df_slice,site_table)

SST_df_island <- site_SST_m_df_metadata %>% 
  dplyr::filter(year %in% (2007:2019),
         !year == 2013) %>%
  group_by(Island,site_status) %>%
  dplyr::summarize(mean_SST = mean(sst,na.rm=T),
                   CV_SST = sd(sst,na.rm=T)/mean_SST)

# sd(SST_df_island$CV_SST) # to test

# for MPAs and refs 
SST_df_island_BOTH <- site_SST_m_df_metadata %>% 
  filter(year %in% (2007:2019),
         !year == 2013) %>%
  group_by(Island) %>%
  dplyr::summarize(mean_SST = mean(sst,na.rm=T),
                   CV_SST = sd(sst,na.rm=T)/mean_SST)

SST_df_island_BOTH$site_status <- "all_sites"



SST_df_island<-rbind(SST_df_island_BOTH,SST_df_island)

### Into island level dataframe
island_env_data<-SST_df_island


write.csv(island_env_data,"./Manuscript_scripts/Exported_data/island_env_data.csv",
          row.names = F)









#### Plot mean temperature data for each MPAs and reference site at each island over time ####
# plot mean on one axis and CV on another

se <- function(x) sd(x) / sqrt(length(x)) # Create own se function
  
SST_annual_data<-merge(site_SST_df_annual_MARSS,site_SST_df_annual_CV)

temp_plot_df<-merge(SST_annual_data,site_table)

# THIS WAS MODIFIED 
temp_plot_df_summarized<-temp_plot_df %>% 
  group_by(site_status,Island,year) %>% 
  summarise(across(where(is.numeric), list(mean = ~mean(., na.rm = TRUE), se = ~se(.))))
            
temp_plot_df$Island <- factor(temp_plot_df$Island, 
                                     levels=c("Anacapa", "SCI", "SRI","SMI"))

temp_plot_df$year<-as.numeric(temp_plot_df$year)

# one for mean; one for CV
temp_plot_A<-temp_plot_df %>% 
  filter(!year == "2020") %>% 
  ggplot(aes(x = year, y = mean_SST,color = site_status, group = site_status)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(group = site),alpha = 0.2) +
 stat_summary(aes(color = site_status),
              fun.data = mean_se, 
              geom = "pointrange") +      
 stat_summary(aes(color = site_status), #color = "black", 
              geom = "point",
              fun = "mean", 
              size = 4) +
 stat_summary(geom = "line", 
              aes(color = site_status), 
              fun = "mean") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  removeGrid() +
  labs(x = "year", y = "annual mean SST (°C)", color = "MPA site status") +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"),
        #legend.position = c(0.15,0.2),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.background = element_rect(fill='transparent')) +
  scale_color_manual(values = c("MPA" = "red", "reference" = "blue")) +
  guides(color = "none") +
  facet_wrap(~Island, nrow =1) +
  coord_cartesian(ylim=c(12,20)) +
  scale_x_continuous(n.breaks = 10)
  

temp_plot_B<-temp_plot_df %>% 
  filter(!year == "2020") %>% 
  ggplot(aes(x = year, y = CV_SST,color = site_status, group = site_status)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(group = site),alpha = 0.2) +
  stat_summary(aes(color = site_status),
               fun.data = mean_se, 
               geom = "pointrange") +      
  stat_summary(aes(color = site_status), #color = "black", 
               geom = "point",
               fun = "mean", 
               size = 4) +
  stat_summary(geom = "line", 
               aes(color = site_status), 
               fun = "mean") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  removeGrid() +
  labs(x = "year", y = "annual CV SST (°C)", color = "MPA site status") +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_color_manual(values = c("MPA" = "red", "reference" = "blue")) +
  facet_wrap(~Island, nrow =1) +
  #coord_cartesian(ylim=c(12,20)) +
  scale_x_continuous(n.breaks = 10)
  
  

temp_plot_full<-temp_plot_A + temp_plot_B + plot_layout(ncol = 1, guides = 'collect') & theme(legend.position = "bottom") & plot_annotation(tag_levels = "A")

ggsave("./Manuscript_scripts/MS_figures/S1_temperature.pdf",
       plot = temp_plot_full,
       width = 10,
       height = 6.5)

ggsave("./Manuscript_scripts/MS_figures/Figure_S1.png",
       plot = temp_plot_full,
       width = 10,
       height = 6.5)

# If you're exhausting vector memory......unhash below to remove large SST files from environment to speed everything up following the running of this script
 #remove(california_cur)
 #remove(california_cur_nona)
 #remove(california_cur_updated)
 #remove(california_cur_2)

### END ###
