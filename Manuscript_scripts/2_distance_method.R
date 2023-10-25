## Calculate over-water distance between PISCO sampling sites ##
# Written by G. Srednick

# packages 
library(data.table)
library(marmap)
library(tidyverse)
library(devtools)
library(viridis)
library(ggspatial)
library(sf)
#library(sp)
library(stars)
#library(mapview)
library(raster)
library(fasttime)
sf_use_s2(T) # unhash maybe?


# Custom function to get location of closest point along a given isobath
# bathy: bathymetric object from marmap::getNOAA.bathy
# x, y: longitude and latitude of a single point
# isobath: negative integer. Minimum depth at which the point should be located
# depth_increment: positive integer to look for deeper cells. Big values allow
# the function to run faster but might lead to more imprecise results
find_closest_negative <- function(bathy, x, y, isobath = -1, depth_increment = 1) {
  # Duplicate point to avoid weird dist2isobath() error message
  point <- data.frame(x, y)
  a <- point[c(1,1), ]
  depth_a <- get.depth(bathy, a, locator = FALSE)
  depth <- unique(depth_a$depth)
  if (depth  < -9) {
    return(unique(depth_a))
  } else {
    while(depth >= -9 ) { # modified here to make minimum depth 14 m
      b <- dist2isobath(bathy, a, isobath = isobath)
      depth_b <- get.depth(bathy, b[,4:5], locator = FALSE)
      depth <- unique(depth_b$depth)
      isobath <- isobath - depth_increment
    }
    return(unique(depth_b))
  }
}



# Channel Islands PISCO spatial breadth
california <-getNOAA.bathy(lon1 = -120.5,lon2 = -119.1,lat1 = 34.2, lat2 =33.8, resolution = 1, keep = T)

# get it higher res from GEBCO
# https://download.gebco.net/
california_v2 <- readGEBCO.bathy("./Data/Spatial/GEBCO/GEBCO_08_Nov_2021_75a9b9262809/gebco_2021_n34.27384363487363_s33.818871416151524_w-120.70434931665659_e-119.25332218408585.nc", resolution = 1)
#summary(california_v2)


#plot(california_v2)

# get sites of interest
sites<-pisco_CI_full_level %>% dplyr::select(site,latitude,longitude)

sites_2<- sites %>% dplyr::filter(site %in% site_table$site)

names(sites_2)<-c("site","latitude","longitude")

print("This will take a bit of time...be patient.")


colnames(sites_2) <- c("site","y", "x")

points2 <- sites_2 %>% 
  mutate(map2_df(x,y, .f = ~find_closest_negative(california_v2, .x, .y,
                                                   isobath = -10,
                                                   depth_increment = 1)))
         
points2

# Some data wrangling to plot lines between pairs of points below
origin <- points2 %>% 
  dplyr::mutate(id = row_number()) %>% 
  dplyr::select(lon = x, lat = y, id)

destination <- points2 %>% 
  dplyr::mutate(id = row_number()) %>% 
  dplyr::select(lon, lat, id)

global <- bind_rows(origin, destination)

# Plot the original points and their new location
california_v2 %>% 
  ggplot(aes(x, y)) + 
  geom_tile(aes(fill = z)) +
  geom_contour(aes(z = z), breaks = 0, colour = 1, linewidth = 0.3) +
  #geom_text(aes(label = z), size = 2) +
  geom_point(data = points2, aes(x = lon, y = lat), color = "red", size = 2) +
  geom_point(data = points2, aes(x = x, y = y), color = "blue") +
  geom_line(data = global, aes(x = lon, y = lat, group = id), linewidth = 0.4) +
  scale_fill_gradient2(low = "skyblue", mid = "white", high = "red",
                       midpoint = 0) +
  coord_fixed()




## now bring in distances between points 
points_dist<-points2[c(1,4,5)]
points_dist_matrix<-points2[c(4,5)]
#points_dist_matrix<-points2[c(3,2)]

pal <- colorRampPalette(c("black","darkblue","blue","lightblue"))
plot(california_v2,image=TRUE,bpal=pal(100),asp=1,col="grey40",lwd=.7,
     main="Bathymetric map of Channel Islands")

plot(california_v2, n = 2, lwd = 0.9, add = TRUE)

points(points_dist_matrix$lon, points_dist_matrix$lat, pch = 21, col = "black",
       bg = "yellow", cex = 1.3)




# Compute transition object with no depth constraint
trans1 <- trans.mat(california_v2)

# Compute transition object with minimum depth constraint: 
# path impossible in waters shallower than -9 meters depth
trans2 <- trans.mat(california_v2,min.depth=-9)


# Computes least cost distances for both transition matrix and plots the results on the map
out1 <- lc.dist(trans1,points_dist_matrix,res="path")
out2 <- lc.dist(trans2,points_dist_matrix,res="path") 
lapply(out1,lines,col="yellow",lwd=4,lty=1) # No depth constraint (yellow paths)
lapply(out2,lines,col="red",lwd=1,lty=1) # Min depth set to -4 meters (red paths)

# Computes and display distance matrices for both situations
dist1 <- lc.dist(trans1,points_dist_matrix,res="dist",meters=TRUE) # in km
dist2 <- lc.dist(trans2,points_dist_matrix,res="dist",meters=TRUE) # in km
#head(dist1)
#head(dist2)





# name rows and columns of dataframe
site_name<-as.vector(points_dist[1])

site_distance_matrix<-as.matrix(dist2)
rownames(site_distance_matrix) <- site_name$site
colnames(site_distance_matrix) <- site_name$site


# dataframe for this.
site_distance_df<-as.data.frame.table(site_distance_matrix) %>%
  unite(Comp, Var1, Var2, sep='-')



site_distance_df<-as.data.frame(as.table(site_distance_matrix))

names(site_distance_df)<-c("SITE1","SITE2","geo_dist")


write.csv(site_distance_df,"./Data/PISCO/overwater_distance.csv")



# merge these with site table to get average distance among sites at the same island
site_table_merge_1<-site_table[c(1,4)]
site_table_merge_2<-site_table[c(1,4)]

names(site_table_merge_1)<-c("SITE1","ISLAND1")
names(site_table_merge_2)<-c("SITE2","ISLAND2")

site_distance_island_df<-merge(site_distance_df,site_table_merge_1)
site_distance_island_df<-merge(site_distance_island_df,site_table_merge_2)

# test
site_distance_island_df %>% 
  dplyr::filter(ISLAND1 == ISLAND2) %>% 
  dplyr::summarise_at(vars(geo_dist),mean)



### Manuscript Figure ###
sites_plot<-merge(sites_2,site_table)

# make this a spatial object
# get US coastline shapefile -- from https://maps.princeton.edu/catalog/stanford-xv279yj9196 
sf_object <- sf::st_read("./Data/Spatial/stanford-xv279yj9196-shapefile/xv279yj9196.shp")

print("This will also take a bit of time. You're almost there!")

# format time 
extractdate <- function(date) {
  day <- format(date, format="%d")
  month <- format(date, format="%m")
  year <- format(date, format="%Y")
  
  cbind(day, month, year)
}



#### SST ####
# This is used in the following script but some elements are used here 
california_cur<-read.csv("./Data/Env_data/SST_data/erdMWsstd1day_LonPM180_5662_c261_a92a.csv", header = T) # daily composite -- large file -- this takes quite some time to import and curate
california_cur_m<-read.csv("./Data/Env_data/SST_data/erdMWsstdmday_LonPM180_7868_a220_4b72.csv", header = T) # this is monthly composite

# Curate daily
california_cur$altitude<-NULL
california_cur<-california_cur[-1,]
california_cur$sst<-as.numeric(california_cur$sst)
california_cur$longitude<-as.numeric(california_cur$longitude)
california_cur$latitude<-as.numeric(california_cur$latitude)
california_cur_2 <- na.omit(california_cur)

# Curate monthly
california_cur_m$altitude<-NULL
california_cur_m<-california_cur_m[-1,]
california_cur_m$sst<-as.numeric(california_cur_m$sst)
california_cur_m$longitude<-as.numeric(california_cur_m$longitude)
california_cur_m$latitude<-as.numeric(california_cur_m$latitude)
california_cur_m_2 <- na.omit(california_cur_m)

# We use the data for the whole region for the figure
california_cur_2$time_pst<-fastPOSIXct(california_cur_2$time)
california_cur_2<- cbind(california_cur_2, extractdate(california_cur_2$time_pst))

california_cur_m_2$time_pst<-fastPOSIXct(california_cur_m_2$time)
california_cur_m_2<- cbind(california_cur_m_2, extractdate(california_cur_m_2$time_pst))

# CV SST - for entire grid
california_cur_m_3<-california_cur_m_2 %>% 
  dplyr::group_by(latitude,longitude) %>% 
  dplyr::mutate(group_id = cur_group_id())

SST_CV_per_point_m <- california_cur_m_3 %>% 
  dplyr::filter(year %in% (2007:2019),
                !year == 2013) %>%
  dplyr::group_by(group_id,latitude,longitude) %>%
  dplyr::summarize(mean_SST = mean(sst, na.rm = T),
                   CV_SST = sd(sst,na.rm=T)/mean_SST) 



SST_CV_per_point_nona<-na.omit(SST_CV_per_point_m)

SST_df_for_raster<-data.frame(x = SST_CV_per_point_nona$longitude,y = SST_CV_per_point_nona$latitude, z = SST_CV_per_point_nona$CV_SST)

# convert temperature to raster
sst_raster <- rasterFromXYZ(SST_df_for_raster, crs = crs(sf_object))
test_spdf <- as(sst_raster, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)

#Make sure they have the same CRS
#st_crs(sf_object_updated) = 4326

# crop these 
sf_object_updated <- st_crop(sf_object, c(xmin = -120.5, ymin = 33.6, xmax = -119.1, ymax = 36.08))
sf_object_small <- st_crop(sf_object, c(xmin = -120.55, ymin = 33.8, xmax = -119.2, ymax = 34.2))

# make polygon so i can fill
islands_poly<-sf_object_small
islands_poly$geometry <- sf::st_cast(islands_poly$geometry, "POLYGON")


bbox <- st_bbox(islands_poly)

# nudge locations
df_nudged <- sites_plot[c(1:3)]

# Add columns for the original x and y coordinates
df_nudged$x_end <- sites_plot$x + 0.01
df_nudged$y_end <- sites_plot$y + 0.01
df_nudged$site_status <-sites_plot$site_status

df_nudged <- df_nudged %>% mutate(x_jit = x + runif(40, min = -0.02, max = 0.02),
                    y_jit = y + runif(40, min = -0.02, max = 0.02))

df_nudged <- df_nudged %>% mutate(x_jit2 = jitter(x,1000), y_jit2 = jitter(y,1000))

#plot

df_nudged_spat <- st_as_sf(df_nudged, coords = c("x", "y"))
st_crs(df_nudged_spat) <- 4326

n=4
breaks = round(seq(min(SST_CV_per_point_nona$CV_SST),max(SST_CV_per_point_nona$CV_SST), length.out = n),3)
fill_breaks = c(0.05,0.1,.15,0.2,0.25,0.3)

cal_map_V2<-
  ggplot(data = islands_poly) +
  geom_sf(data = islands_poly, color = "white", fill = "white") +
  geom_tile(data = SST_CV_per_point_nona, aes(x = longitude, y = latitude, fill = CV_SST)) +
  scale_fill_viridis(breaks = fill_breaks) +
  geom_point(data = df_nudged, aes(x = x_jit2, y = y_jit2, color = site_status), size = 3, pch = 20) +
  geom_segment(data = df_nudged, aes(x = x_jit2, y = y_jit2, xend = x, yend = y)) +
  scale_color_manual(values = c(MPA = "red", reference = "blue")) +
  labs(x = "Longitude", y = "Latitude", fill = "CV SST", color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.justification = c("right","bottom"),
        legend.position = c(.98, .1),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        text = element_text(size=13, color = "black")) +
  ylim(c(33.8,34.2)) +
  xlim(c(-120.55,-119.3)) +
  ggsn::scalebar(islands_poly, 
                 dist = 10, st.size=3, st.dist = 0.05, location = "bottomleft",
                anchor = c(x = -120.5, y = 33.85),
                 height=0.03, transform = T, dist_unit = "km", model = 'WGS84') +
  ggspatial::annotation_north_arrow(location="tr", width = unit(1, "cm"),
                         height = unit(1, "cm")) +
  guides(fill = guide_colorbar(frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) +
  coord_sf(expand = 0)
  




ggsave("./Manuscript_scripts/MS_figures/site_map.pdf",
       plot=cal_map_V2,
       width = 10,
       height= 5)




# END #


