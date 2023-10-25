## MPA asynchrony -- Data curation and formatting ####

# Written by G. Srednick - 5/3/2022

# Purpose: Data import and processing for all downstream analyses.
# Schematic for processes included in file

# Data found in https://search.dataone.org/view/doi%3A10.6085%2FAA%2FPISCO_kelpforest.1.6
# Other curated data are provided in submitted code package and uploaded accordingly.



# packages
library(tidyverse)
library(zoo)
library(lubridate)
library(tidyr)

### data entry
# need a new filtering step here to filter based on taxon table
# every species should be included based on sampling
# notes will be in separate script for why this is.


# fish
pisco_fish<-read.csv("./Data/PISCO/PISCO_kelpforest_fish.1.3.csv")

ggplot(data = pisco_fish) +
  geom_line(aes(x = year, y = site))

replication_df<-pisco_fish %>% 
  dplyr::group_by(site,zone,transect) %>% 
  summarise()



# swath = algae + inverts
swath<-read.csv("./Data/PISCO/PISCO_kelpforest_swath.1.2.csv")

ggplot(data = swath) +
  geom_line(aes(x = year, y = site))

# algae + hab
UPC<-read.csv("./Data/PISCO/PISCO_kelpforest_upc.1.2.csv")

ggplot(data = UPC) +
  geom_line(aes(x = survey_year, y = site))



##  curate data to "per transect"
# sum across water column -- this is new on 4/12/2022
pisco_fish_reduced <- pisco_fish %>% 
  dplyr::group_by(year,site,zone,transect,classcode) %>% 
  dplyr::summarize(count = sum(count,na.rm = TRUE))

# Fish

pisco_fish_reduced_wide<-pisco_fish_reduced %>% 
  pivot_wider(names_from = classcode,
              values_from = count, 
              values_fill = NA)


# Swath
swath_reduced <- swath %>%  # mean across swaths
  dplyr::group_by(year,site,zone,transect,classcode) %>% 
  dplyr::summarize(count = mean(count,na.rm = TRUE))

# make wide
swath_reduced_wide<-swath_reduced %>% 
  pivot_wider(names_from = classcode,
              values_from = count, 
              values_fill = NA)


# UPC
UPC_reduced <- UPC %>% 
  dplyr::filter(category == "COVER") %>% 
  dplyr::group_by(year,site,zone,transect,classcode) %>% 
  dplyr::summarize(cover = mean(pct_cov,na.rm = TRUE))

UPC_reduced_wide<-UPC_reduced %>% pivot_wider(names_from = classcode, 
                                              values_from = cover,
                                              values_fill = NA)

dim(pisco_fish_reduced_wide)
dim(swath_reduced_wide)
dim(UPC_reduced_wide)

totalcols = (160-5) + (134-5) + (69-5) + 4

# merge data frames 
pisco_wide <- merge(pisco_fish_reduced_wide, swath_reduced_wide) %>%
  merge(UPC_reduced_wide)

dim(pisco_wide)

# get to single site level -- mean across zones and transects
pisco_wide_site_level <- pisco_wide %>% 
  dplyr::group_by(year,site) %>% 
  dplyr::summarise(across(where(is.numeric), ~mean(., na.rm = TRUE)))  # mean to get to site level; taking average across zones and transects



pisco_wide_site_level[sapply(pisco_wide_site_level, is.nan)] <- NA

write.csv(pisco_wide_site_level,"./Data/PISCO/PISCO_kelpforest_data_aggregated.csv", row.names = F)






# Cut to sites that we are interested in
pisco_sites<-read.csv("./Data/PISCO/PISCO_kelpforest_site_table.1.2_CI_sites.csv")

# filter to site level 
pisco_sites_level<- pisco_sites %>% 
  dplyr::group_by(site,MPA_Name,site_designation,site_status,Island,Side) %>% 
  dplyr::summarize()

pisco_sites_coord<- pisco_sites %>% 
  dplyr::group_by(site,latitude,longitude,site_status,Island,Side) %>% 
  dplyr::summarize()

pisco_CI <- merge(pisco_sites_level,pisco_wide_site_level) # filter for only the sites at CI


length(unique(pisco_CI$site))


## Duration: will have to filter out sites that have fewer years --- < 6 years
## 15 March 2022 -- will have to filter time series for years 2007 onward 

startyear = 2007



ggplot(data = pisco_CI) +
  geom_point(aes(x = year, y = site,color = site_status)) +
  facet_grid(~Island)

ggplot(data = pisco_CI) +
  geom_point(aes(x = year, y = site,color = site_status))

pisco_CI_full<- pisco_CI %>%  
  dplyr::group_by(site) %>% 
  dplyr::filter(!length(year)<= 6, # filter out sites with less than 6 years of sammpling
                year >= startyear,
                !site %in% c("SMI_TYLER_BIGHT_W", # remove these sites due to low abundance
                              "SMI_CROOK_POINT_W",
                             "SMI_CROOK_POINT_E", # added on 15 May 2023
                              "SMI_TYLER_BIGHT_E",
                             "SCI_PELICAN_CEN")) # removed due to missing temperature data
# JUST A TEST ---> see if not removing crook point does something

ggplot(data = pisco_CI_full) +
  geom_point(aes(x = year, y = site, color = site_status))

length(unique(pisco_CI_full$site)) # 50 sites; now 46 after dropping the four sites at SMI
dim(pisco_CI_full) # [1] 541 354


pisco_CI_full$transect <- NULL

pisco_CI_full_level<- pisco_CI_full %>% 
  dplyr::group_by(site,MPA_Name,site_designation,site_status,Island,Side) %>% 
  dplyr::summarize()

pisco_CI_full_level<-merge(pisco_CI_full_level,pisco_sites_coord)

length(unique(pisco_CI_full_level$site)) # 50 sites; now 46



### work in "spatial_sync.R" ###

# merge this back with distances here
#pisco_reduced_sites<-merge(pisco_CI_full_level,pisco_sites_dist) # distance between sites



# Make long for exploratory

pisco_CI_full_long<-pisco_CI_full %>% 
  gather(species, cover,-c(site,year,MPA_Name,site_designation,site_status,Island,Side))


psico_names<-data.frame(names(pisco_CI_full)) # names look okay, no duplicates. Although, may want to pick one MAPY variable


#ggplot(data = pisco_CI_full_long) +
#  geom_point(aes(x = year, y = species)) +
#  theme(legend.position = "none")



# filter out unnecessary species (from merge script) --- this isnt currently correct --- we should be using stuff from 'trophic_level_designation.R'
# Problem: this cuts species, and adds a few identifiers that arent important anymore

# I think we will cut this and just merge with 'trophic_level_designation.R' products later 
#pisco_CI_full_long_reduced<-merge(species_ready,pisco_CI_full_long, by = "species",fill=FALSE) 

pisco_CI_full_long_reduced<-pisco_CI_full_long

ggplot(data = pisco_CI_full_long_reduced) +
  geom_point(aes(x = year, y = site, color = site_status))


# Check number of species
length(unique(pisco_CI_full_long_reduced$species)) ### 274 species; updated following removal of "species_ready": 346

pisco_CI_full_long_reduced_V2<-na.omit(pisco_CI_full_long_reduced) # we want to do this beacuse it removes species that are all NAs

length(unique(pisco_CI_full_long_reduced_V2$species)) ### 174 species; now 207 spp




# filter out species that have fewer than 5 observations 
pisco_CI_full_long_reduced_V3<- pisco_CI_full_long_reduced_V2 %>% 
  dplyr::group_by(site,species) %>% 
  dplyr::filter(!length(species)<= 2) # observed more than one time 

ggplot(data = pisco_CI_full_long_reduced_V3) +
  geom_point(aes(x = year, y = site))

# no sites with consecutive years missing -- proceed with filling

length(unique(pisco_CI_full_long_reduced_V3$species)) ###  139 species; again up to 143




#################################### Data ##############################################

## Notes
# Only species with > 1 observations included
# Only years after 2006
# Spp. interpolated for missing survey years - this done in long format so not interpolated if they weren't present during a survey 
# NAs will be replaced with zeros in this -- "probability of encountering a given species is equal across sites" 

### Curation steps
# (1) remove sites with more than 2 years in a row missing -- DONE
# (2) approximate abundance for gap years


#Make df with full time series

temporal_res<-pisco_CI_full_long_reduced_V3 %>% 
  dplyr::group_by(site,year) %>% 
  dplyr::summarize() %>% 
  dplyr::select(site,year)
# ignore warnings 

dim(temporal_res)

DFfilled <- pisco_CI_full_long_reduced_V3 %>%
  tidyr::complete(year = startyear:2019, 
                  fill = list(cover = NA)) %>%
  as.data.frame() %>%
  dplyr::filter(!year == 2013)

temporal_res_afterfill<-DFfilled %>% 
  dplyr::group_by(site,year) %>% 
  dplyr::summarize() %>% 
  dplyr::select(site,year)

dim(temporal_res_afterfill)

missing_res<-anti_join(temporal_res_afterfill,temporal_res)
length(unique(missing_res$site)) # years missing at 18 sites

ggplot(data = DFfilled) +
  geom_point(aes(x = year, y = site))

dim(pisco_CI_full_long_reduced_V3) # BEFORE:  [1] 22157     9
dim(DFfilled) # AFTER: [1] 365960 rows.



# Fill NA's with linear interpolation
library(zoo)
pisco_data_actual<-DFfilled %>% 
  dplyr::group_by(site) %>%
  arrange(site,year) %>%
  dplyr::mutate(cover = na.approx(cover, maxgap = 2, rule = 2))


# Fill in missing site information
pisco_data_actual_filled<-pisco_data_actual %>%
  dplyr::group_by(site) %>%
  fill(site_designation,site_status,Island,Side,MPA_Name) %>%
  fill(site_designation,site_status,Island,Side,MPA_Name, .direction = "up")





# Make a table that shows the number of years per site
# unfilled data - what are the sites,years and ranges
survey_table_prelim<-pisco_CI_full_long_reduced_V3 %>% 
  dplyr::select(site,site_status,Island,year) %>%
  dplyr::group_by(site,site_status,Island,year) %>%
  dplyr::summarize_all(is.numeric) #%>%

survey_table_merge<-survey_table_prelim %>% 
  dplyr::group_by(site,site_status,Island) %>%
  dplyr::summarize(min_year = min(year),
            max_year = max(year),
            no_years = length(year))

# filled data - what years are missing at each site
filled_prelim<-pisco_data_actual_filled %>%
  dplyr::select(site,site_status,Island,year) %>%
  dplyr::group_by(site,site_status,Island,year) %>%
  dplyr::summarize_all(is.numeric)

missing_site_years<-anti_join(filled_prelim,survey_table_prelim)

# concatenate these for the table 
missing_site_concat<-missing_site_years %>% 
  dplyr::group_by(site) %>% 
  dplyr::reframe(missing_years = paste0(year, collapse = ", ")) 


# add them to the prelim table
final_table<-merge(survey_table_merge,missing_site_concat, all = T)
write.csv(final_table,"./Manuscript_scripts/Exported_data/Tables/Pubtables/temporal_res.csv",row.names = F)


# Trophic data ####


# Following from trophic_level_designation.R 


# merge with new trophic metadata
pisco_trophic<-read.csv("./Manuscript_scripts/trophic/pisco_trophic.csv")

pisco_troph_small<-pisco_trophic[c(7,8)]
names(pisco_troph_small)<-c("species","trophic_level")

missing_trophs<- pisco_data_actual_filled %>% 
  dplyr::filter(!species %in% pisco_troph_small$species)

unique(missing_trophs$species) # last check before merge -- some are removed still: CLAM EUGRUB
# these are all fine -- ignore

# merge proper 
pisco_troph_maindata<-merge(pisco_troph_small,pisco_data_actual_filled)  
length(unique(pisco_troph_maindata$species)) # number of species: 141; updated 149; currently 132

# Fill in missing species information
pisco_troph_maindata_filled<-pisco_troph_maindata %>%
  dplyr::group_by(species) %>%
  fill(trophic_level) %>%
  fill(trophic_level, .direction = "up")


# make wide
pisco_troph_wide<-pisco_troph_maindata_filled[-c(2,5:9)] %>% pivot_wider(names_from=species,values_from = cover,values_fill = list(Value = 0))
dim(pisco_troph_wide) # [1] 540  93 --- excellent; updated on 3/22/2022 -- now 540 151; updated 4/12/2022 -- now 540 134


# remove NAs
pisco_troph_wide[sapply(pisco_troph_wide, is.na)] <- 0


# make long again
pisco_troph_long<-pisco_troph_wide %>% 
  pivot_longer(!c(year,site),names_to = "species", values_to = "cover")

dim(pisco_troph_long) # [1] 80460     4 -- 71280 x 4 on 4/12/22


## make a site table and a species table
# site table 
site_table<-pisco_troph_maindata_filled %>% 
  dplyr::group_by(site,site_designation,site_status,Island,Side) %>% 
  dplyr::summarise() %>%
  dplyr::select(site,site_designation,site_status,Island,Side)
# ignore warnings

# Print this for later use 
write.csv(site_table,"./Manuscript_scripts/Exported_data/site_table.csv", row.names = F)


# species table
spp_table <- pisco_troph_maindata_filled %>% 
  dplyr::group_by(species,trophic_level) %>% 
  dplyr::summarise() %>%
  dplyr::select(species,trophic_level)
# ignore warnings

write.csv(spp_table,"./Manuscript_scripts/Exported_data/spp_table.csv", row.names = F)



# merge with site table
pisco_troph_maindata_complete_V1<-merge(pisco_troph_long,site_table) 
dim(pisco_troph_maindata_complete_V1) # [1] 80460     8  -- 71280 x 8 on 4/12/22

# merge with species table
pisco_troph_maindata_complete_V2<-merge(pisco_troph_maindata_complete_V1,spp_table) 
dim(pisco_troph_maindata_complete_V2) # [1] 80460     9  -- 71280 x 9 on 4/12/22

#ready to go
pisco_troph_complete<-pisco_troph_maindata_complete_V2


### Need to address differences in measurement -- scale this 
# Note: scaled cover not used in analyses --- there is an alternative normalization/standarization step in species interaction script
# Scaled cover is kept in code for any supplemental assessment

pisco_troph_complete$cover_scaled<-as.numeric(scale(pisco_troph_complete$cover)) # the scaled value equal to zero; constant added #0.1206656

length(unique(pisco_troph_complete$species)) # 132 

# wide
pisco_troph_complete_wide<-pisco_troph_complete[-c(4:9)] %>% pivot_wider(names_from=species,values_from = cover_scaled) # scaled
pisco_troph_complete_wide<-merge(pisco_troph_complete_wide,site_table)
dim(pisco_troph_complete_wide) # [1] 540     155 (138 on 4/12/22)


# wide, not scaled
pisco_troph_complete_wide_ns<-pisco_troph_complete[-c(5:10)] %>% pivot_wider(names_from=species,values_from = cover) # scaled
pisco_troph_complete_wide_ns<-merge(pisco_troph_complete_wide_ns,site_table)
dim(pisco_troph_complete_wide_ns) # [1] 540     155 (138 on 4/12/22)






# count number of species per trophic level

pisco_troph_complete %>% dplyr::group_by(trophic_level) %>% dplyr::summarise(spp_count = length(unique(species)))


write.csv(pisco_troph_complete,"./Manuscript_scripts/Exported_data/pisco_troph_complete.csv", row.names = F)

# Remove some objects to clear space
remove(pisco_fish)
remove(swath)
remove(UPC)
remove(UPC_reduced_wide)
remove(UPC_reduced)
#### End of data tidying ####

