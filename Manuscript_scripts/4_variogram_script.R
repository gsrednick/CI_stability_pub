######### Spatial synchrony and Moran's I #########


# use 'synchrony' package to calculate mantel correlations
# colors for MPA-MPA, ref-ref, MPA-ref -- presumably the biggest erosion with distance is MPA-ref?
# see how correlations erode with overwater distance
#library(synchrony)
library(ggExtra)
library(vegan)
library(phylin)


# also have to integrate side of island here.....
# so i would have to add another identifier on top of site_status that shows side of island. (e.g., same or diff)
# probably run this all through a loop for this to make it easier --- nah.


#### Fixing coordinates
# need to make a function in which i can run a variogram where the distances are already provided






# "fake" coordinates based on distances
fake_coord<-data.frame(cmd_test) # from "distance_method.R"
colnames(fake_coord)<-c("latitude","longitude")
fake_coord$site<-row.names(fake_coord)

# merge with site table 
new_coords<-merge(fake_coord,pisco_CI_full_level[c(1:4)])

# merge with species data

# real
#vario_data<-merge(pisco_CI_full_level[c(1:3,7:8)],pisco_troph_wide)
#vario_data_both<-vario_data[-c(1:3,6)]
#vario_data_MPA<-vario_data %>% filter(site_status == "MPA") %>% select(c(4:5,7:138))
#vario_data_ref<-vario_data %>% filter(site_status == "reference") %>% select(c(4:5,7:138))


# fake 
vario_data_fake<-merge(new_coords,pisco_troph_wide)


vario_data_both_both_fake<-vario_data_fake %>%
  dplyr::select(-c(1,4:7))

vario_data_both_N_fake<-vario_data_fake %>%
  filter(Side == "North") %>% 
  dplyr::select(-c(1,4:7))

vario_data_both_S_fake<-vario_data_fake %>%
  filter(Side == "South") %>% 
  dplyr::select(-c(1,4:7))

vario_data_MPA_both_fake<-vario_data_fake %>% 
  filter(site_status == "MPA") %>% 
  dplyr::select(-c(1,4:7))

vario_data_MPA_N_fake<-vario_data_fake %>% 
  filter(site_status == "MPA", Side == "North") %>% 
  dplyr::select(-c(1,4:7))

vario_data_MPA_S_fake<-vario_data_fake %>% 
  filter(site_status == "MPA", Side == "South") %>% 
  dplyr::select(-c(1,4:7))

vario_data_ref_both_fake<-vario_data_fake %>% 
  filter(site_status == "reference") %>% 
  dplyr::select(-c(1,4:7))

vario_data_ref_N_fake<-vario_data_fake %>% 
  filter(site_status == "reference", Side == "North") %>% 
  dplyr::select(-c(1,4:7))

vario_data_ref_S_fake<-vario_data_fake %>% 
  filter(site_status == "reference", Side == "South") %>% 
  dplyr::select(-c(1,4:7))


nrand_count = 100
nbin_count = 10

A_both <- vario(data = vario_data_both_both_fake, type = "Kendall", nrands = nrand_count, n.bins = nbin_count)
A_N <- vario(data = vario_data_both_N_fake, type = "Kendall", nrands = nrand_count, n.bins = nbin_count)
A_S <- vario(data = vario_data_both_S_fake, type = "Kendall", nrands = nrand_count, n.bins = nbin_count)
B_both <- vario(data = vario_data_MPA_both_fake, type = "Kendall", nrands = nrand_count, n.bins = nbin_count)
B_N <- vario(data = vario_data_MPA_N_fake, type = "Kendall", nrands = nrand_count, n.bins = nbin_count)
B_S <- vario(data = vario_data_MPA_S_fake, type = "Kendall", nrands = nrand_count, n.bins = nbin_count)
C_both <- vario(data = vario_data_ref_both_fake, type = "Kendall", nrands = nrand_count, n.bins = nbin_count)
C_N <- vario(data = vario_data_ref_N_fake, type = "Kendall", nrands = nrand_count, n.bins = nbin_count)
C_S <- vario(data = vario_data_ref_S_fake, type = "Kendall", nrands = nrand_count, n.bins = nbin_count)


A_both_DF<-data.frame(vario = A_both$vario,site_status="All sites",side = "both",bins = A_both$bins)
A_nort_DF<-data.frame(vario = A_N$vario,site_status="All sites",side = "North",bins = A_N$bins)
A_sout_DF<-data.frame(vario = A_S$vario,site_status="All sites",side = "South",bins = A_S$bins)
B_both_DF<-data.frame(vario = B_both$vario,site_status="MPA",side = "both",bins = B_both$bins)
B_nort_DF<-data.frame(vario = B_N$vario,site_status="MPA",side = "North",bins = B_N$bins)
B_sout_DF<-data.frame(vario = B_S$vario,site_status="MPA",side = "South",bins = B_S$bins)
C_both_DF<-data.frame(vario = C_both$vario,site_status="reference",side = "both",bins = C_both$bins)
C_nort_DF<-data.frame(vario = C_N$vario,site_status="reference",side = "North",bins = C_N$bins)
C_sout_DF<-data.frame(vario = C_S$vario,site_status="reference",side = "South",bins = C_S$bins)


vario_df_fake_new<-rbind(A_both_DF,
      A_nort_DF,
      A_sout_DF,
      B_both_DF,
      B_nort_DF,
      B_sout_DF,
      C_both_DF,
      C_nort_DF,
      C_sout_DF)


### could run this again but with it constrained to site

# fake
vario_overwater_plot<-vario_df_fake_new %>% 
  filter(side == "both") %>%
  ggplot(aes(x = bins/100, y = vario, color = site_status)) +
  geom_point(size = 4) +
  geom_line(aes(color = site_status)) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 0)+
  removeGrid() +
  labs(x = "Relative lag distances", y = "Kendall correlation", color = "site status") +
  theme(legend.position = c(0.8,0.8)) +
  #ylim(0.25,0.75) +
  scale_color_manual(values = c("All sites" = "black", "MPA" = "red", "reference" = "blue")) #+
  #ggtitle('Overwater coordinates')


# Takeaway: spatial synchrony greater among MPAs, reference sites desynchronize metacommunity.
# still need to fix distances....

# maybe bring in pvalues too and make the dots filled if significant -- based on rands -- add this to manuscript 

#varo_plot<-vario_real + vario_fake +plot_layout(ncol = 2, guides = "collect")

ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/Figures/vario_plot.pdf",
       plot=vario_overwater_plot,
       width = 6,
       height= 5)






# new #########
### Alternatively make a matrix and generate coordinates via MDS

# If i do this with "genetic distance" i can look at erosion of bray curtis dissimilarity with geographic distance


pisco_com_dummy_site<-data.frame(pisco_troph_wide[c(1:2)],pisco_com_dummy)

#sim_matrix<-pisco_com_dummy_site %>% group_by(site) %>% summarize_all(mean) # temporarily summarized to site; will have to make a distance for every year later.... maybe
sim_matrix<-pisco_com_dummy_site 


# bray curtis distance among sites
com_sim <- as.matrix(vegdist(sim_matrix[-1], distance = "bray", trace = FALSE, na.rm = F, autotransform = FALSE))

## need to do this for every year.... how....
# I want to show how multivariate similarity over time erodes with distance 

# (1) do vario for every year; show every point as a mean ± SE; we dont care about change over time that much -- care more about change over space

dim(com_sim)

# Rename matrix
sites_list<-sim_matrix$site
rownames(com_sim)<-sites
colnames(com_sim)<-sites
unique(sites_list)
dim(site_distance_matrix)

allsites<-data.frame(site = rownames(site_distance_matrix))


# filter out unused sites

realsite<-allsites %>% filter(site %in% pisco_troph_wide$site)
sites_vector<-realsite$site
site_distance_matrix_actual<-site_distance_matrix[rownames(site_distance_matrix)%in%sites_vector,colnames(site_distance_matrix)%in%sites_vector]



# make dfs for MPA and ref
MPA_sites<-site_table %>% filter(site_status == "MPA")
ref_sites<-site_table %>% filter(site_status == "reference")


# filter distance data for only refs and MPAs
site_distance_matrix_MPA<-site_distance_matrix_actual[rownames(site_distance_matrix_actual)%in%MPA_sites$site,
                                                      colnames(site_distance_matrix_actual)%in%MPA_sites$site]
dim(site_distance_matrix_MPA)
site_distance_matrix_ref<-site_distance_matrix_actual[rownames(site_distance_matrix_actual)%in%ref_sites$site,
                                                      colnames(site_distance_matrix_actual)%in%ref_sites$site]


mean(site_distance_matrix_actual)
# filter com data for only refs and MPAs
com_sim_MPA<-com_sim[rownames(com_sim)%in%MPA_sites$site,
                                         colnames(com_sim)%in%MPA_sites$site]

com_sim_ref<-com_sim[rownames(com_sim)%in%ref_sites$site,
                                         colnames(com_sim)%in%ref_sites$site]

## run variogram
# MPAs
vario_new_MPA<-gen.variogram(site_distance_matrix_MPA,com_sim_MPA, lag = 10000, bootstraps = 999)
plot(vario_new_MPA)

MPA_vario_data<-data.frame(lag = vario_new_MPA$lag, gamma = vario_new_MPA$gamma, site_status = "MPA")

# reference sites
vario_new_ref<-gen.variogram(site_distance_matrix_ref,com_sim_ref, lag = 10000, bootstraps = 999)
#plot(vario_new_ref)

ref_vario_data<-data.frame(lag = vario_new_ref$lag, gamma = vario_new_ref$gamma, site_status = "reference")

# All sites
vario_new_all<-gen.variogram(site_distance_matrix_actual,com_sim, lag = 10000, bootstraps = 999)
#plot(vario_new_all)

all_vario_data<-data.frame(lag = vario_new_all$lag, gamma = vario_new_all$gamma, site_status = "all sites")

# merge together
vario_data<-rbind(MPA_vario_data,ref_vario_data,all_vario_data)


plot(vario_new_all)
gv<-gv.model(vario_new_all)


# plot
vario_overwater_plot<-vario_data %>% 
  ggplot(aes(x = lag/1000, y = gamma, color = site_status)) + # convert to km for variogram
  geom_point(size = 4) +
  geom_line(aes(color = site_status)) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 0)+
  removeGrid() +
  labs(x = "Lag distances (km)", y = "Semivariance", color = "site status") +
  theme(legend.position = c(0.8,0.8)) +
  #ylim(-0.1,0.4) +
  scale_color_manual(values = c("all sites" = "black", "MPA" = "red", "reference" = "blue")) #+






# try in loop where every year is a different vario function -- take mean se for figure


# for MPAs and References
sim_year_list<- split(sim_matrix, sim_matrix$year)


vario_time<-lapply(sim_year_list,function(x){ #this function works and can be applied to list items with lapply
  
  
  com_sim <- as.matrix(vegdist(x[-c(1,2)], distance = "bray", trace = FALSE, na.rm = F, autotransform = FALSE))
  
  
  rownames(com_sim)<-x$site
  colnames(com_sim)<-x$site
  
  # filter out unused sites
  site_distance_matrix_actual<-site_distance_matrix[rownames(site_distance_matrix)%in%sites_vector,colnames(site_distance_matrix)%in%sites_vector]
  
  # filter distance data for only refs and MPAs
  site_distance_matrix_MPA<-site_distance_matrix_actual[rownames(site_distance_matrix_actual)%in%MPA_sites$site,
                                                        colnames(site_distance_matrix_actual)%in%MPA_sites$site]
  
  site_distance_matrix_ref<-site_distance_matrix_actual[rownames(site_distance_matrix_actual)%in%ref_sites$site,
                                                        colnames(site_distance_matrix_actual)%in%ref_sites$site]
  
  
  # filter com data for only refs and MPAs
  com_sim_MPA<-com_sim[rownames(com_sim)%in%MPA_sites$site,
                       colnames(com_sim)%in%MPA_sites$site]
  
  com_sim_ref<-com_sim[rownames(com_sim)%in%ref_sites$site,
                       colnames(com_sim)%in%ref_sites$site]
  
  ## run variogram
  # MPAs
  vario_new_MPA<-gen.variogram(site_distance_matrix_MPA,com_sim_MPA, lag = 10000, bootstraps = 999)
  #plot(vario_new_MPA)
  
  MPA_vario_data<-data.frame(lag = vario_new_MPA$lag, gamma = vario_new_MPA$gamma, site_status = "MPA")
  
  # reference sites
  vario_new_ref<-gen.variogram(site_distance_matrix_ref,com_sim_ref, lag = 10000, bootstraps = 999)
  #plot(vario_new_ref)
  
  ref_vario_data<-data.frame(lag = vario_new_ref$lag, gamma = vario_new_ref$gamma, site_status = "reference")
  
  # All sites
  vario_new_all<-gen.variogram(site_distance_matrix_actual,com_sim, lag = 10000, bootstraps = 999)
  #plot(vario_new_all)
  
  all_vario_data<-data.frame(lag = vario_new_all$lag, gamma = vario_new_all$gamma, site_status = "all sites")
  
  #merge together
  vario_data<-rbind(MPA_vario_data,ref_vario_data,all_vario_data)
  
  
  summary<-data.frame(vario_data,year = unique(x$year))
  
  return(summary)
  
})


vario_time_df<-do.call(rbind.data.frame, vario_time)


vario_overwater_plot_2<-vario_time_df %>% 
  ggplot(aes(x = lag/1000, y = gamma, color = site_status)) + # convert to km for variogram
  #geom_point(size = 4) +
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
  #geom_abline(intercept = 0, slope = 0)+
  removeGrid() +
  labs(x = "Lag distances (km)", y = "Semivariance", color = "MPA site status") +
  theme(legend.position = c(0.2,0.8)) +
  #ylim(-0.1,0.4) +
  scale_color_manual(values = c("all sites" = "black", "MPA" = "red", "reference" = "blue"))

# should add model to this to show the maximum spatial distance at which the variable is spatially autocorrelated


ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/Figures/vario_plot_new.pdf",
       plot=vario_overwater_plot_2,
       width = 6,
       height= 5)



# ANCOVA for this

vario_aov<-aov(gamma~lag * site_status, data = vario_time_df)
summary(vario_aov)
TukeyHSD(vario_aov)

# take away: MPAs that are far away from one another are more similar -- MPA effects are most important at local scale






#### Now just for targeted species ##
# for MPAs and References

pisco_targeted_wide<-pisco_targeted %>% 
  select(species,cover,site,year,Island,site_status) %>% 
  pivot_wider(names_from = species, values_from = cover,values_fill = 0)

pisco_targeted_wide$dummy<-1

sim_matrix_targ<-pisco_targeted_wide[-c(3,4)]


sim_targ_list<- split(sim_matrix_targ, sim_matrix_targ$year)

vario_time_targ<-pblapply(sim_targ_list,function(x){ #this function works and can be applied to list items with lapply
  
  
  com_sim <- as.matrix(vegdist(x[-c(1,2)], distance = "bray", trace = FALSE, na.rm = F, autotransform = FALSE))
  
  
  rownames(com_sim)<-x$site
  colnames(com_sim)<-x$site
  
  # filter out unused sites
  site_distance_matrix_actual<-site_distance_matrix[rownames(site_distance_matrix)%in%sites_vector,colnames(site_distance_matrix)%in%sites_vector]
  
  # filter distance data for only refs and MPAs
  site_distance_matrix_MPA<-site_distance_matrix_actual[rownames(site_distance_matrix_actual)%in%MPA_sites$site,
                                                        colnames(site_distance_matrix_actual)%in%MPA_sites$site]
  
  site_distance_matrix_ref<-site_distance_matrix_actual[rownames(site_distance_matrix_actual)%in%ref_sites$site,
                                                        colnames(site_distance_matrix_actual)%in%ref_sites$site]
  
  
  # filter com data for only refs and MPAs
  com_sim_MPA<-com_sim[rownames(com_sim)%in%MPA_sites$site,
                       colnames(com_sim)%in%MPA_sites$site]
  
  com_sim_ref<-com_sim[rownames(com_sim)%in%ref_sites$site,
                       colnames(com_sim)%in%ref_sites$site]
  
  ## run variogram
  # MPAs
  vario_new_MPA<-gen.variogram(site_distance_matrix_MPA,com_sim_MPA, lag = 10000, bootstraps = 999)
  #plot(vario_new_MPA)
  
  MPA_vario_data<-data.frame(lag = vario_new_MPA$lag, gamma = vario_new_MPA$gamma, site_status = "MPA")
  
  # reference sites
  vario_new_ref<-gen.variogram(site_distance_matrix_ref,com_sim_ref, lag = 10000, bootstraps = 999)
  #plot(vario_new_ref)
  
  ref_vario_data<-data.frame(lag = vario_new_ref$lag, gamma = vario_new_ref$gamma, site_status = "reference")
  
  # All sites
  vario_new_all<-gen.variogram(site_distance_matrix_actual,com_sim, lag = 10000, bootstraps = 999)
  #plot(vario_new_all)
  
  all_vario_data<-data.frame(lag = vario_new_all$lag, gamma = vario_new_all$gamma, site_status = "all sites")
  
  #merge together
  vario_data<-rbind(MPA_vario_data,ref_vario_data,all_vario_data)
  
  
  summary<-data.frame(vario_data,year = unique(x$year))
  
  return(summary)
  
})


vario_time_targ_df<-do.call(rbind.data.frame, vario_time_targ)


vario_overwater_targ_plot<-vario_time_targ_df %>% 
  ggplot(aes(x = lag/1000, y = gamma, color = site_status)) + # convert to km for variogram
  #geom_point(size = 4) +
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
  #geom_abline(intercept = 0, slope = 0)+
  removeGrid() +
  labs(x = "Lag distances (km)", y = "Semivariance", color = "MPA site status") +
  theme(legend.position = c(0.7,0.12),
        legend.background = element_blank()) +
  #ylim(-0.1,0.4) +
  scale_color_manual(values = c("all sites" = "black", "MPA" = "red", "reference" = "blue"))

# should add model to this to show the maximum spatial distance at which the variable is spatially autocorrelated


ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/Figures/vario_plot_targ.pdf",
       plot=vario_overwater_targ_plot,
       width = 6,
       height= 5)



vario_time_targ_df %>% 
  group_by(site_status) %>% 
  filter(lag < 25000) %>% 
  summarize_all(funs(mean,sd))

vario_time_targ_df %>%
  group_by(site_status) %>%
  filter(lag < 25000) %>% 
  dplyr::summarise(mean.sync = mean(gamma, na.rm = TRUE),
                   sd.sync = sd(gamma, na.rm = TRUE),
                   n.sync = n()) %>%
  mutate(se.sync = sd.sync / sqrt(n.sync),
         lower.ci = mean.sync - qt(1 - (0.05 / 2), n.sync - 1) * se.sync,
         upper.ci = mean.sync + qt(1 - (0.05 / 2), n.sync - 1) * se.sync)

#### Analytical End ####




### EXTRA ####


# there is also this ......

## Reconstruct worked example of Wagner (submitted):
X <- matrix(c(1, 2, 3, 2, 1, 0), 3, 2)
Y <- c(3, -1, -2)
tmat <- c(1:3)
## Canonical correspondence analysis (cca):
Example.cca <- cca(X, Y)
Example.cca <- mso(Example.cca, tmat)
msoplot(Example.cca)
Example.cca$vario
## Correspondence analysis (ca):
Example.ca <- mso(cca(X), tmat)
msoplot(Example.ca)

## Unconstrained ordination with test for autocorrelation
## using oribatid mite data set as in Wagner (2004)
data(mite)
data(mite.env)
data(mite.xy)
mite.cca <- cca(log(mite + 1))
mite.cca <- mso(mite.cca, mite.xy, grain = 1, permutations = 99)
msoplot(mite.cca)
mite.cca
## Constrained ordination with test for residual autocorrelation
## and scale-invariance of species-environment relationships
mite.cca <- cca(log(mite + 1) ~ SubsDens + WatrCont + Substrate + Shrub + Topo, mite.env)
mite.cca <- mso(mite.cca, mite.xy, permutations = 99)
msoplot(mite.cca)
mite.cca

env_reduced<-pisco_env %>% group_by(site) %>% summarize_all(mean) %>% select(1)

 
fake_coord_allsites<-fake_coord %>% filter(site %in% pisco_troph_wide$site)

mite.cca <- cca(sim_matrix[-c(1:2)],env_reduced)
mite.cca <- mso(mite.cca,site_distance_matrix_actual, grain = 10000, permutations = 999)
msoplot(mite.cca)



if(require(ade4)){
  data(oribatid)
  # Hellinger transformation
  fau <- sqrt(oribatid$fau / outer(apply(oribatid$fau, 1, sum), rep(1, ncol(oribatid$fau)), "*"))
  # Removing linear effect
  faudt <- resid(lm(as.matrix(fau) ~ as.matrix(oribatid$xy)))
  mvspec <- variogmultiv(faudt, oribatid$xy, nclass = 20)
  mvspec
  plot(mvspec$d, mvspec$var,type = 'b', pch = 20, xlab = "Distance", ylab = "C(distance)")
}


fau <- sqrt(sim_matrix[-c(1:2)] / outer(apply(sim_matrix[-c(1:2)], 1, sum), rep(1, ncol(sim_matrix[-c(1:2)])), "*"))
# Removing linear effect
faudt <- resid(lm(as.matrix(fau) ~ as.matrix(site_distance_matrix_actual)))
mvspec <- variogmultiv(faudt, site_distance_matrix_actual, nclass = 20)
mvspec
plot(mvspec$d, mvspec$var,type = 'b', pch = 20, xlab = "Distance", ylab = "C(distance)")





## END ##