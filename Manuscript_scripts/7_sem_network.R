# SEM for trophic networks ####

library(lavaan)
library(tidyverse)
#library(pbapply)
library(psych)
library(semTable)
library(dynamic)
library(performance)

# SEM Approach ####


# 8A - Community level SEM ####

# bring in data
#community_interactions<-read.csv("./Manuscript_scripts/Exported_data/community_interactions.csv")
community_interactions<-read.csv("./Manuscript_scripts/Exported_data/community_interactions_REV.csv")
community_synchrony<-read.csv("./Manuscript_scripts/Exported_data/community_synchrony.csv")
community_stability<-read.csv("./Manuscript_scripts/Exported_data/community_stability.csv")
habitat_loadings<-read.csv("./Manuscript_scripts/Exported_data/hab_similarity.csv")
site_env_data<-read.csv("./Data/Env_data/site_env_data.csv")
network_trophic_position<-read.csv("./Manuscript_scripts/Exported_data/network_trophic_position.csv")

# [For testing] This is the old CV SST file pre-rgdal depreciation - slightly different values for three SMI sites -- no major changes to story 
#site_env_data_old<-read.csv("/Users/icarus3/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/Data/Env_data/site_env_data.csv")


#qqnorm(log(island_env_data$CV_SST+1))

# merge together
B_sem_df <- Reduce(inner_join, list(site_env_data,
                                    community_interactions,
                                    community_synchrony,
                                    habitat_loadings,
                                    community_stability,
                                    network_trophic_position))

#B_sem_df <- Reduce(inner_join, list(site_env_data,Beta_and_sync,habitat_loadings)) # old
# transform interactions to improve clarity
B_sem_df$pred_B_3root_inv <- (B_sem_df$pred_B_3root * -1)
B_sem_df$prey_B_3root_inv <- (B_sem_df$prey_B_3root * -1)

length(unique(B_sem_df$network_ID)) # 20 networks 

B_sem_df_dc<-cbind(B_sem_df,dummy.code(B_sem_df$site_status))




path_network_specific<-
  '
prey_B_3root_inv ~ MPA + CV_SST + PC1 + count_comp_presnt
pred_B_3root_inv ~ MPA + CV_SST + PC1 + count_comp_presnt

#prey_B_3root_inv ~ reference
#prey_B_3root_inv ~ MPA

#pred_B_3root_inv ~ reference
#pred_B_3root_inv ~ MPA

synchrony ~ MPA + prey_B_3root_inv + pred_B_3root_inv + CV_SST + PC1 + count_comp_presnt

#synchrony ~ reference
# synchrony ~ MPA

distance_inv ~ prey_B_3root_inv + pred_B_3root_inv + synchrony

'



## (1) Community all - SEM ####

path_network_all<-
  '
prey_B_3root_inv ~ MPA
prey_B_3root_inv ~ CV_SST
prey_B_3root_inv ~ count_comp_presnt

pred_B_3root_inv ~ MPA
pred_B_3root_inv ~ CV_SST
pred_B_3root_inv ~ count_comp_presnt

synchrony ~ prey_B_3root_inv
synchrony ~ pred_B_3root_inv
synchrony ~ MPA
synchrony ~ CV_SST # dropping this leads to the best model performance (rmsea = 0.063)
synchrony ~ count_comp_presnt

count_comp_presnt ~ MPA

distance_inv ~ prey_B_3root_inv
distance_inv ~ pred_B_3root_inv
distance_inv ~ synchrony
distance_inv ~ MPA


# additional for test ----> marginally better.... 0.12 vs 0.15
pred_B_3root_inv ~~ prey_B_3root_inv # this definitely must stay 


distance_inv ~ CV_SST + count_comp_presnt #+ MPA + consumer_PA + resource_PA # have to figure out what to do here 

# seems like dropping the habitat PC1 is a good start

# this model is good; good rmsea, CFI, tli, and pvalue
'

# number of networks used for main SEM
length(unique(B_sem_df_dc$network_ID))
#networks_of_interest

fit_test_all_networks<-sem(path_network_all,
                  data = B_sem_df_dc,
                  ordered = names("site_status"))


summary(fit_test_all_networks)




# model indicators
varTable(fit_test_all_networks) # look at parameter variance
SEM_network_all_fit_df<-data.frame(model = lavInspect(fit_test_all_networks, "fit")[c("chisq","df", "pvalue", "rmsea","cfi","tli")])
SEM_network_all_fit_df # is model a good fit?
equivTest(fit_test_all_networks) # is model a good fit?
perform<-model_performance(fit_test_all_networks) # is model a good fit?
modificationindices(fit_test_all_networks, minimum.value = 20) #only print MIs > 20 

# look at covariance 
#plot_matrix(residuals(fit_test_all_networks, type="cor")$cov) # look for covariance


# Parameter estimates
param_summary_all<-parameterestimates(fit_test_all_networks) # non standardized
param_summary_all<-standardizedSolution(fit_test_all_networks, type="std.all") # standardized



### Spit to table ####
write.csv(SEM_network_all_fit_df,"./Manuscript_scripts/Exported_data/Tables/sem_all_fit_mod_2.csv")

# curate for plotting
SEM_all_plotting<-param_summary_all %>% filter(!rhs == lhs,
                                               !lhs == rhs)

SEM_all_plotting<-na.omit(SEM_all_plotting) # remove covariances that arent of interest

SEM_all_plotting$lhs<-factor(SEM_all_plotting$lhs, levels = c("CV_SST","PC1","MPA","count_comp_presnt",
                                                              "resource_PA","consumer_PA",
                                                              "prey_B_3root_inv","pred_B_3root_inv","synchrony",
                                                              "distance_inv"))

SEM_all_plotting$rhs<-factor(SEM_all_plotting$rhs, levels = c("CV_SST","PC1","MPA","count_comp_presnt",
                                                              "resource_PA","consumer_PA",
                                                              "prey_B_3root_inv","pred_B_3root_inv","synchrony",
                                                              "distance_inv"))


SEM_all_plotting<-SEM_all_plotting %>% mutate(direction = ifelse(est.std < 0, "-", "+"),
                                              sig = ifelse(pvalue < 0.05, "yes","no"))

# Write table for supplement
sem_com_all_table<-SEM_all_plotting[c(1,3,4,5,7)]
write.csv(sem_com_all_table,"./Manuscript_scripts/Exported_data/Tables/sem_com_all_coef_2.csv", row.names = F)



# 8B - Metacommunity level SEM ####

# bring in data
metacommunity_interactions<-read.csv("./Manuscript_scripts/Exported_data/metacommunity_interactions_REV.csv")
metacommunity_synchrony<-read.csv("./Manuscript_scripts/Exported_data/metacommunity_synchrony.csv")
metacommunity_stability<-read.csv("./Manuscript_scripts/Exported_data/metacommunity_stability.csv")
island_habitat_loadings<-read.csv("./Manuscript_scripts/Exported_data/island_hab_similarity.csv")
island_env_data<-read.csv("./Data/Env_data/island_env_data.csv")


island_env_data<-island_env_data %>% 
  mutate(site_status = fct_recode(site_status, "all" = "all_sites")) # recode


# merge these together 
B_meta_sem_df <- Reduce(inner_join, list(island_env_data,
                                         metacommunity_synchrony,
                                         metacommunity_interactions,
                                         metacommunity_stability,
                                         #island_habitat_loadings,
                                         network_trophic_position)) # revise disper_df for metacommunities

# add concatenated variable for interpretation
B_meta_sem_df<-B_meta_sem_df %>% 
  unite(island_status, c("Island", "site_status"), remove = F)
  
length(unique(B_meta_sem_df$network_ID))
length(unique(B_meta_sem_df$island_status))

# make coefficients inverse for clarity 
B_meta_sem_df$pred_B_3root_inv <- (B_meta_sem_df$pred_B_3root * -1)
B_meta_sem_df$prey_B_3root_inv <- (B_meta_sem_df$prey_B_3root * -1)

# dummy code for 3 levels of site status
B_meta_sem_df_dc<-cbind(B_meta_sem_df,dummy.code(B_meta_sem_df$site_status))


#count_comp_presnt
## Metacommunity scale ####

path_meta_network_specific<-
  '
prey_B_3root_inv ~ CV_SST + PC1 + count_comp_presnt
pred_B_3root_inv ~ CV_SST + PC1 + count_comp_presnt

prey_B_3root_inv ~ MPA
prey_B_3root_inv ~ all

pred_B_3root_inv ~ MPA
pred_B_3root_inv ~ all

synchrony ~ prey_B_3root_inv + pred_B_3root_inv + CV_SST + PC1 + count_comp_presnt

synchrony ~ MPA
synchrony ~ all

distance_inv ~ prey_B_3root_inv + pred_B_3root_inv + synchrony


'


## (2) Metacom all ####

### Across sites and networks ###
path_meta_network_all<-
  '
  # competition
prey_B_3root_inv ~ CV_SST
prey_B_3root_inv ~ count_comp_presnt
prey_B_3root_inv ~ MPA
prey_B_3root_inv ~ reference

# predation
pred_B_3root_inv ~ CV_SST
pred_B_3root_inv ~ count_comp_presnt
pred_B_3root_inv ~ MPA
pred_B_3root_inv ~ reference

# prey synchrony
synchrony ~ prey_B_3root_inv
synchrony ~ pred_B_3root_inv
synchrony ~ CV_SST
synchrony ~ count_comp_presnt
synchrony ~ MPA
synchrony ~ reference

# prey richness
count_comp_presnt ~ MPA
count_comp_presnt ~ reference

# stability
distance_inv ~ prey_B_3root_inv
distance_inv ~ pred_B_3root_inv
distance_inv ~ synchrony
distance_inv ~ count_comp_presnt
distance_inv ~ CV_SST
distance_inv ~ MPA # model is a much better fit with this path indicated
distance_inv ~ reference # model is a much better fit with this path indicated


# Additional - covariance between predation and competition
pred_B_3root_inv ~~ prey_B_3root_inv 

'

fit_test_all_meta_networks<-sem(path_meta_network_all,
                           data = B_meta_sem_df_dc)


summary(fit_test_all_meta_networks)
mean(abs(inspect(fit_test_all_meta_networks,what="std")$beta))


# model performance
SEM_network_meta_all_fit_df<-data.frame(model = lavInspect(fit_test_all_meta_networks, "fit")[c("chisq","df", "pvalue", "rmsea","srmr","cfi","tli","rmsea.pvalue")])
SEM_network_meta_all_fit_df
equivTest(fit_test_all_meta_networks)
varTable(fit_test_all_meta_networks)
modificationindices(fit_test_all_meta_networks, minimum.value = 20) #only print MIs > 20 

# check covariance table
#plot_matrix(residuals(fit_test_all_meta_networks, type="cor")$cov)

# standardized parameters
param_summary_meta_all<-standardizedSolution(fit_test_all_meta_networks, type="std.all") # might need to revise all with standardized coefficients.... probably even better for plotting



### Spit to table ####
write.csv(SEM_network_meta_all_fit_df,"./Manuscript_scripts/Exported_data/Tables/sem_meta_all_fit_mod_2.csv")

# curate for plotting
SEM_meta_all_plotting<-param_summary_meta_all %>% filter(!rhs == lhs,
                                               !lhs == rhs)

SEM_meta_all_plotting<-na.omit(SEM_meta_all_plotting) # remove covariances that arent of interest


SEM_meta_all_plotting$lhs<-factor(SEM_meta_all_plotting$lhs, levels = c("CV_SST","PC1","MPA","reference","count_comp_presnt",
                                                              "consumer_PA","resource_PA",
                                                              "prey_B_3root_inv","pred_B_3root_inv","synchrony",
                                                              "distance_inv"))

SEM_meta_all_plotting$rhs<-factor(SEM_meta_all_plotting$rhs, levels = c("CV_SST","PC1","MPA","reference","count_comp_presnt",
                                                                        "consumer_PA","resource_PA",
                                                                        "prey_B_3root_inv","pred_B_3root_inv","synchrony",
                                                                        "distance_inv"))


SEM_meta_all_plotting<-SEM_meta_all_plotting %>% mutate(direction = ifelse(est.std < 0, "-", "+"),
                                              sig = ifelse(pvalue < 0.05, "yes","no"))



# Write table for supplement
sem_meta_all_table<-SEM_meta_all_plotting[c(1,3,4,5,7)]


write.csv(sem_meta_all_table,"./Manuscript_scripts/Exported_data/Tables/sem_meta_all_coef_2.csv", row.names = F)

# table_dir<-dir("./Manuscript_scripts/Exported_data/Tables/")
# 
# semTable(fit_test_all_meta_networks,
#          columns = c("est","se","p"),
#          fits = c("chisq", "rmsea"),
#          type = "csv",
#          longtable = FALSE, 
#          table.float = TRUE,
#          file = "./Manuscript_scripts/Exported_data/Tables/all_meta_coef.csv")






### Metacommunity SEM reconfigured for comparing MPAs vs. ref metacoms

path_meta_network_all_MPAvRef<-
  '
  # competition
prey_B_3root_inv ~ CV_SST
prey_B_3root_inv ~ count_comp_presnt
prey_B_3root_inv ~ MPA
prey_B_3root_inv ~ all

# predation
pred_B_3root_inv ~ CV_SST
pred_B_3root_inv ~ count_comp_presnt
pred_B_3root_inv ~ MPA
pred_B_3root_inv ~ all

# prey synchrony
synchrony ~ prey_B_3root_inv
synchrony ~ pred_B_3root_inv
synchrony ~ CV_SST
synchrony ~ count_comp_presnt
synchrony ~ MPA
synchrony ~ all

# prey richness
count_comp_presnt ~ MPA
count_comp_presnt ~ all

# stability
distance_inv ~ prey_B_3root_inv
distance_inv ~ pred_B_3root_inv
distance_inv ~ synchrony
distance_inv ~ count_comp_presnt
distance_inv ~ CV_SST
distance_inv ~ MPA
distance_inv ~ all



# Additional - covariance between predation and competition
pred_B_3root_inv ~~ prey_B_3root_inv 

'

fit_test_all_meta_networks_MPAvRef<-sem(path_meta_network_all_MPAvRef,
                                data = B_meta_sem_df_dc)


summary(fit_test_all_meta_networks_MPAvRef)


# model performance
SEM_network_meta_all_fit_MPAvRef_df<-data.frame(model = lavInspect(fit_test_all_meta_networks_MPAvRef, "fit")[c("chisq","df", "pvalue", "rmsea","srmr","cfi","tli","rmsea.pvalue")])
SEM_network_meta_all_fit_MPAvRef_df
equivTest(fit_test_all_meta_networks_MPAvRef)
varTable(fit_test_all_meta_networks_MPAvRef)
modificationindices(fit_test_all_meta_networks_MPAvRef, minimum.value = 20) #only print MIs > 20 

# check covariance table
#plot_matrix(residuals(fit_test_all_meta_networks_MPAvRef, type="cor")$cov)

# standardized parameters
param_summary_meta_all_MPAvRef<-standardizedSolution(fit_test_all_meta_networks_MPAvRef, type="std.all") # might need to revise all with standardized coefficients.... probably even better for plotting



### Spit to table ####
write.csv(SEM_network_meta_all_fit_MPAvRef_df,"./Manuscript_scripts/Exported_data/Tables/sem_meta_all_fit_MPAvRef_mod.csv")

# curate for plotting
SEM_meta_all_plotting_MvR<-param_summary_meta_all_MPAvRef %>% filter(!rhs == lhs,
                                                         !lhs == rhs)

SEM_meta_all_plotting_MvR<-na.omit(SEM_meta_all_plotting_MvR) # remove covariances that arent of interest


SEM_meta_all_plotting_MvR$lhs<-factor(SEM_meta_all_plotting_MvR$lhs, levels = c("CV_SST","PC1","MPA","all","count_comp_presnt",
                                                                        "consumer_PA","resource_PA",
                                                                        "prey_B_3root_inv","pred_B_3root_inv","synchrony",
                                                                        "distance_inv"))

SEM_meta_all_plotting_MvR$rhs<-factor(SEM_meta_all_plotting_MvR$rhs, levels = c("CV_SST","PC1","MPA","all","count_comp_presnt",
                                                                        "consumer_PA","resource_PA",
                                                                        "prey_B_3root_inv","pred_B_3root_inv","synchrony",
                                                                        "distance_inv"))


SEM_meta_all_plotting_MvR<-SEM_meta_all_plotting_MvR %>% mutate(direction = ifelse(est.std < 0, "-", "+"),
                                                        sig = ifelse(pvalue < 0.05, "yes","no"))



# Write table for supplement
sem_meta_all_table_MvR<-SEM_meta_all_plotting_MvR[c(1,3,4,5,7)]


write.csv(sem_meta_all_table_MvR,"./Manuscript_scripts/Exported_data/Tables/sem_meta_all_MvR_coef.csv", row.names = F)







### Assess variation in species richness across MPAs

# ggplot(community_interactions,aes(x = site_status, y = log(count_comp_presnt+1))) +
#          geom_boxplot() +
#   facet_grid(~Island)
# 
# 
# div_model<-aov(log(count_comp_presnt + 1) ~ site_status * Island,community_interactions)
# summary(div_model)
# TukeyHSD(div_model)



## END ##