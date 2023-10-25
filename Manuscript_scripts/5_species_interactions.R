# 6 - Species interactions calculations ###

# for manuscript: 

# This script runs multivariate autoregressive models for each trophic network at community and metacommunity scale.
 
library(MARSS)
library(tidyverse)
library(parallel)
library(pbapply)


# define cores here!!
core_numbers = detectCores() - 2 # I use 10 cores for this.



# bring in env data
site_SST_df_annual<-read.csv("./Data/Env_data/SST_data/annual_SST_MARSS.csv")
meta_SST_df_annual<-read.csv("./Data/Env_data/SST_data/annual_meta_SST_MARSS.csv")

# competitor matrix

competitors<-trophic_network %>% select(-c(consumer_ID,ID))

competitor_species<-unique(competitors$resource_ID)

competitor_mat=matrix(list(0),nrow=length(competitor_species),
                      ncol=length(competitor_species),
                      dimnames=list(competitor_species,competitor_species))

dim(competitor_mat)

# split by network and assign competitive interactions
competitors_mat_list<-split(competitors, competitors$network_ID)


competitor_networks_list<-pblapply(competitors_mat_list,function(x){ #this function works and can be applied to list items with lapply
  
  comp_network_matrix<-competitor_mat %>% subset(rownames(competitor_mat) %in% x$resource_ID,
                                                 colnames(competitor_mat) %in% x$resource_ID)
  
  comp_network_matrix[]<-1
  
  L<-as.data.frame(as.table(comp_network_matrix)) %>% select(1:3)
  colnames(L)<-c("species_1","species_2","Freq")
  
  return(L)
  
})




competitor_networks_complete<-do.call(rbind.data.frame, competitor_networks_list)
rownames(competitor_networks_complete)<-NULL

links_for_merge<-links_for_matrix
colnames(links_for_merge) <- c("species_1","species_2","Freq")

consumers<-trophic_network %>% select(consumer_ID,network_ID)

# bind with predation matrix
complete_int_frequency<-rbind(links_for_merge,competitor_networks_complete)
dim(complete_int_frequency)

complete_int_frequency$species_1<-as.factor(complete_int_frequency$species_1)
complete_int_frequency$species_2<-as.factor(complete_int_frequency$species_2)


complete_int_frequency_reversed<-complete_int_frequency #  reverse order and renaming
names(complete_int_frequency_reversed)<-c("to","from","Freq")

complete_int_frequency_reversed<-complete_int_frequency_reversed[,c(2,1,3)]


complete_int_frequency_updated<-complete_int_frequency #%>% filter(!duplicated(cbind(species_1, species_2)))

int_reversed_updated<-complete_int_frequency_reversed %>% filter(!duplicated(cbind(to, from)))

dim(complete_int_frequency)
dim(complete_int_frequency_updated)


levs <- unique(unlist(complete_int_frequency_updated[c(1:2)], use.names = F))
ready_matrix<-table(lapply(complete_int_frequency_updated[c(1:2)], factor, levs))

sorted_matrix<-ready_matrix[, sort(colnames(ready_matrix))][sort(rownames(ready_matrix)), ]

int_matrix<-as.matrix(sorted_matrix)
lower.tri(sorted_matrix)
diag(int_matrix)<-1

# test -- remove this when done testing
int_matrix[sp, sp]


#Undirected Graph
#graph = graph_from_data_frame(complete_int_frequency_updated, directed = F) # original
graph = graph_from_data_frame(int_reversed_updated, directed = T) # this directs top down -- works for consumers and competitors

as_data_frame(graph, what="edges")

#Adjacency
int_matrix<-as_adjacency_matrix(graph, sparse = F)
int_matrix<-as_adjacency_matrix(graph, sparse = F)
diag(int_matrix) <- 1
dim(int_matrix)

int_matrix[int_matrix > .1] <- 1

## Write all data to csvs for analysis 
#write.table(int_matrix,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/MARS/Mars_for_spartan/species_interaction_matrix.txt")
#write.csv(int_matrix,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/MARS/Mars_for_spartan/species_interaction_matrix.csv")
#write.csv(pisco_troph_complete,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/MARS/Mars_for_spartan/time_series_data.csv", row.names = F)
#write.csv(trophic_species,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/MARS/Mars_for_spartan/trophic_species_data.csv", row.names = F)



# 6A - Community level interactions ####

# make new df for analyses
DBall<-pisco_troph_complete %>% 
  #filter(site %in% c("ANACAPA_EAST_ISLE_CEN")) %>%  # filtering stage for test -- remove this line later
  #filter(Island == "Anacapa") %>%  # filtering stage for test -- remove this line later
  filter(species %in% trophic_species$species) %>% # this doesnt filter correctly
  dplyr::select(year, site, species, cover) %>%
  spread(species, cover, fill = 0)

cons_network_list<- split(trophic_network, trophic_network$network_ID)

# for testing
#x = trophic_network %>% filter(network_ID == 43)


## Loop for calculation this ####
network_covariance<-mclapply(cons_network_list,function(x){ #this function works and can be applied to list items with lapply
  

  # get specific taxa for trophic network
  relevant_taxa<-x %>% 
    select(consumer_ID, resource_ID) %>% 
    gather() %>% 
    select(-key) %>% 
    rename(species = value)
  
  # filter time series for these species and reshape for analyses
  DBall<-pisco_troph_complete %>% 
    filter(species %in% relevant_taxa$species) %>%
    dplyr::select(year, site, species, cover) %>%
    spread(species, cover, fill = 0) %>% 
    group_by(site) %>% 
    mutate(site_ID =cur_group_id())
  
  dim(DBall)
  
  nsites=length(unique(DBall$site_ID)) #number of sites or repeats
  
  # for testing  
  #ksite = 1

  for (ksite in 1:nsites){

    DB=DBall[DBall$site_ID==ksite,] ## Select a site
    
    DB$site_ID<-NULL
    head(DB) 
    
    abundance_mat=as.matrix(DB[-c(1:2)]) ### Create matrix with time series of abundances
    
    sp=colnames(DB[-c(1:2)])
    nspecies=length(sp)
    
    
    dates=DB$year
    dates_bis=dates 
    tab_sp=DB[,sp]
    
    #Setting MARSS model -- known topology of the adjacency matrix
    ### This will define number of species as a by-product
    interaction_matrix<-int_matrix %>% subset(colnames(int_matrix) %in% sp, rownames(int_matrix) %in% sp)
    
    mat_ord2 <- function(mx) mx[sort(rownames(mx)), c(sort(rownames(mx)), sort(setdiff(colnames(mx), rownames(mx))))]
    
    n_species<-length(sp)
    
    # define all interactions in matrix as "1": all species interact; the context of interactions is defined later
    
    interaction_matrix=matrix(list(1),nrow=length(sp),
                              ncol=length(sp),
                              dimnames=list(sp,sp))

    B1=matrix(list(0),nrow=length(sp),ncol=length(sp),dimnames=list(rownames(interaction_matrix),colnames(interaction_matrix)))

    # organizing matrix
        
    for (i in 1:length(sp)){
      s=sp[i]
      for (j in 1:length(sp)){
        s2=sp[j]
        if(interaction_matrix[i,j]!=0){
          #B1[i,j]=paste("alpha_",i,",",j,sep="")
          B1[i,j]=paste(s2,s,sep="") ### beware, paste(s2,s,sep="") means we have the impacting species first  
        }else{
          if(s==s2){
            B1[i,j]=paste(s2,s,sep="")
          }
        }
      }
    }
    
   
    
    ### Defining other parameters
    iter_min=100 
    iter_estimate=10 
    iter_scenario=10
    iter_boot=1000 
    BFGS=T
    U1="zero"
    Q1="diagonal and unequal"
    Z1=diag(1,length(sp),length(sp))
    A1="zero"
    R1="zero"
    V1=diag(1,length(sp))
    pi1="zero"
    
    # Environmental effect is equal covariance 
    C <- "unconstrained"
    D <- "unconstrained"
    
    aalpha=0.05
    cntl.list=list(safe=TRUE, maxit=1500, allow.degen=TRUE) # from Dexter et al.
    
    tab_sp=log(tab_sp+1) # log+1 to remove NAs
    tab_sp=t(scale(tab_sp, scale = FALSE)) ### log and center the abundance data
    
  
    
    
    # Addition of covariates to model - 4 August 2023
    CV_temp_data<-site_SST_df_annual %>% 
      filter(site == unique(DB$site),
             year %in% DB$year) %>%
      select(mean_SST)
    
    covariates<-CV_temp_data$mean_SST
    covariate_matrix <- matrix(covariates, ncol = length(dates))
    
    tab_cov=t(scale(covariates, scale = FALSE)) # -- added on 4 Aug 2023
    # no NAs can be present in the vector
    
    #tab_cov=t(scale(tab_cov_bis, scale = FALSE)) ### same for the covariates // centered but not scaled here for comparison to Gregoire's code------ leaving abiotic out for right now!!!!
    rownames(tab_sp)=sp
    
    # Env. covariate parameters
    c1 = tab_cov # added on 4 Aug 2023
    d = tab_cov

    model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1,D=D,d=d) ### Specifying the parameters to MARSS - 

    # actual model run    
    fit_log=MARSS(tab_sp, method="BFGS",model=model.list) ### Fit of MAR model. -- this takes a very long time with the whole food web (> 15 mins)
    
    # resulting coefficients
    coef_table <- coef(fit_log, type = "matrix")$B[1:n_species, 1:n_species]
    rownames(coef_table)<-sp
    colnames(coef_table)<-sp
    
    B_table = coef(fit_log,type="matrix")$B[1:n_species, 1:n_species]
    Q_table = coef(fit_log,type="matrix")$Q[1:n_species, 1:n_species]
    

    # print results into external folder for analyses
    site <-unique(DB$site)
    network_ID <-unique(x$network_ID)
    coef_df_ready<-as.data.frame.table(coef_table)
    coef_df_ready$site <-unique(DB$site)
    coef_df_ready$network_ID <-unique(x$network_ID)
    
    write.csv(coef_table,file=paste("./MARS/Mars_model_B/Betas/Revised_AUG_23_autoC/tables/B_points_",site,"_NET-",network_ID,"_table.csv",sep=""))
    write.csv(coef_df_ready,file=paste("./MARS/Mars_model_B/Betas/Revised_AUG_23_autoC/dfs/B_points_",site,"_NET-",network_ID,"_df.csv",sep=""),row.names = F)
    
    print(site)
    print(network_ID)
    


  }

}, mc.cores = core_numbers, mc.preschedule = FALSE)





## 6A.2 - Import community level interactions ####

# bring back csvs by site and bind them
# label intraspecific B's as DD_B
# label interspecific B's as PB (predation) or CB (competition)
# make a label for this above

# list files
B_file_names <- list.files(path="./MARS/Mars_model_B/Betas/Revised_AUG_23_autoC/dfs", 
                           pattern=".csv", full.names = TRUE)

length(B_file_names)

# import, bind, and name them
B_dataframe <- do.call(rbind, lapply(B_file_names, function(x) 
  cbind(read.csv(x, stringsAsFactors = FALSE), filename = x)))


# clean them to show network ID; wont be necessary next time we run this
B_dataframe_mode<-B_dataframe %>% 
  separate(filename, sep = "-", into = c("filename","network_ID_temp")) %>%
  separate(network_ID_temp, sep = "_", into = c("network_ID","extra"))

# select columns of interest
B_dataframe_ready<-B_dataframe_mode %>% select(c(1,2,3,4,6)) # unhash this when ready


lvel_chck<-B_dataframe_ready %>% expand(site,network_ID) # all good




## Grouping ###

# group into density dependent - not used in these analyses 
B_df_DD<- B_dataframe_ready %>% filter(Var1 == Var2)
dim(B_df_DD)

# group into competitive
B_df_interspecific<- B_dataframe_ready %>% filter(!Var1 == Var2)


# Function to extract betas and reshape 

#x<-trophic_network %>% filter(network_ID == 38) # for testing

cons_network_list<- split(trophic_network, trophic_network$network_ID)

B_summarizing<-pblapply(cons_network_list,function(x){ #this function works and can be applied to list items with lapply
  
  # filter species for whether they were present in time series or not
  site_sums_species<-pisco_troph_complete %>% 
    filter(species %in% c(x$consumer_ID,x$resource_ID)) %>%
    group_by(site,species) %>% 
    summarise(sum_cover = sum(cover,na.rm = T))
  
  B_dataframe_ready_fornetwork<-B_dataframe_ready %>% 
    filter(network_ID %in% x$network_ID) # filter for network
  
  B_by_site_list<- split(B_dataframe_ready_fornetwork, B_dataframe_ready_fornetwork$site)
  
  consumer = unique(x$consumer_ID)
  

  B_by_site<-lapply(B_by_site_list,function(x){
    
    sums_for_site<-site_sums_species %>% filter(site %in% x$site)
    
    
    # i want to know if there is only one species; predation strength might be stronger for one vs > 1 species
    comp_richness_results <- sums_for_site %>%
      filter(!species == consumer) %>%
      summarise(count_comp = n_distinct(species),
                count_comp_presnt = n_distinct((species)[sum_cover > 0]))
    
    
    
    # for competitors; summarizing interaction strength
    competitors_results<-x %>% 
      filter(Var1 %in% sums_for_site$species, 
             Var2 %in% sums_for_site$species) %>%
      filter(!Var1 == Var2) %>% 
      filter(!Var1 == consumer,
             !Var2 == consumer) %>% 
      summarise(mean_comp_B = mean(Freq, na.rm =T),
                sd_comp_B = sd(Freq,na.rm = T)) 
    
    competitors_results_merged<-merge(competitors_results,comp_richness_results)
    
    
    # for consumers
    consumer_results<-x %>% 
      filter(!Var1 == Var2) %>%
      filter(Var1 == consumer) %>%
      summarise(mean_cons_B = mean(Freq, na.rm =T),
                sd_cons_B = sd(Freq,na.rm = T))
    
    # density dependent interactions
    B_df_interspecific<- x %>% 
      filter(Var1 == Var2) %>%
      summarise(mean_DD_B = mean(Freq, na.rm =T),
                sd_DD_B = sd(Freq,na.rm = T))
    
    B_results<-Reduce(function(x, y) merge(x, y, all=TRUE), list(competitors_results_merged, consumer_results, B_df_interspecific))
    
    B_results$site<-unique(x$site)
    
    cons_present<- sums_for_site %>%
      filter(species %in% consumer) %>% 
      mutate(cons_present = if_else(sum_cover > 0, TRUE,FALSE))
    
    
    B_results$cons_present <- cons_present$cons_present
    
    return(B_results)
  })
  
  B_by_site_df<-do.call(rbind.data.frame, B_by_site)
  B_by_site_df$network_ID<-unique(x$network_ID)
  
  return(B_by_site_df)
})


summarized_B_df<-do.call(rbind.data.frame, B_summarizing)
summarized_B_df<-merge(summarized_B_df,site_table)

# filter out networks with less than one prey -- temporary?
summarized_B_df<-summarized_B_df %>% filter(count_comp > 1)
dim(summarized_B_df) # all good

# only look at interactions where a predator was present
dim(summarized_B_df)
summarized_B_df_pred_present<-summarized_B_df %>% filter(cons_present==TRUE)
dim(summarized_B_df_pred_present)


# check for normality
qqnorm(summarized_B_df_pred_present$mean_comp_B)
qqnorm(summarized_B_df_pred_present$mean_cons_B)
qqnorm(summarized_B_df_pred_present$mean_DD_B)

summarized_B_df_pred_present$pred_B_3root<-sign(summarized_B_df_pred_present$mean_cons_B) * (abs(summarized_B_df_pred_present$mean_cons_B))^(1/3)
summarized_B_df_pred_present$prey_B_3root<-sign(summarized_B_df_pred_present$mean_comp_B) * (abs(summarized_B_df_pred_present$mean_comp_B))^(1/3)
summarized_B_df_pred_present$DD_B_3root<-sign(summarized_B_df_pred_present$mean_DD_B) * (abs(summarized_B_df_pred_present$mean_DD_B))^(1/3)

qqnorm(summarized_B_df_pred_present$pred_B_3root) # better
qqnorm(summarized_B_df_pred_present$prey_B_3root) # better



# write this as csv for importing in later scripts 

write.csv(summarized_B_df_pred_present,"./Manuscript_scripts/Exported_data/community_interactions_REV.csv",row.names = F)



# summarize to MPA/ref; island for every network - start without island
summarized_B_df_reduced<-summarized_B_df_pred_present %>% 
  group_by(site_status,network_ID) %>% 
  summarize_if(is.numeric,mean,na.rm = T)



# curate for MPA status
B_ratio_prep_MPA<-summarized_B_df_reduced %>% filter(site_status == "MPA")
B_ratio_prep_ref<-summarized_B_df_reduced %>% filter(site_status == "reference")

names(B_ratio_prep_MPA)
colnames(B_ratio_prep_MPA) <-c("site_status","network_ID","mean_comp_B_MPA","sd_comp_B_MPA","count_comp_MPA",       
                               "count_comp_presnt_MPA","mean_cons_B_MPA","sd_cons_B_MPA","mean_DD_B_MPA","sd_DD_B_MPA",
                               "pred_B_3root_MPA","prey_B_3root_MPA","DD_B_3root_MPA")   

colnames(B_ratio_prep_ref) <-c("site_status","network_ID","mean_comp_B_ref","sd_comp_B_ref","count_comp_ref",       
                               "count_comp_presnt_ref","mean_cons_B_ref","sd_cons_B_ref","mean_DD_B_ref","sd_DD_B_ref",
                               "pred_B_3root_ref","prey_B_3root_ref","DD_B_3root_ref")   


# simple
ratio_complete<-merge(B_ratio_prep_MPA[-1],B_ratio_prep_ref[-1])

# merge with trophic position -- put this in curation (1) later 
trophic_updated<-pisco_trophic
trophic_updated$consumer<-trophic_updated$classcode

trophic_updated_formerge<-trophic_updated[c("PA_trophic","consumer")]
names(trophic_updated_formerge)<-c("PA_trophic","consumer_ID")

consumer_w_ID<-trophic_network %>% 
  group_by(network_ID,consumer_ID) %>% 
  summarize_all(mean) %>% 
  select(network_ID,consumer_ID) %>% merge(., trophic_updated_formerge)


ratio_complete_withnetwork<-merge(ratio_complete,consumer_w_ID, by = "network_ID", all= F)
dim(ratio_complete_withnetwork)



## plot these for exploration
Beta_cons<-ggplot(ratio_complete_withnetwork,aes(x = pred_B_3root_ref, y = pred_B_3root_MPA, fill = PA_trophic)) +
  geom_pointrange(aes(ymax = pred_B_3root_MPA, ymin = pred_B_3root_ref)) +
  geom_point(pch = 21, size = 3)+
  geom_hline(yintercept = 0,linetype = "dashed") +
  geom_vline(xintercept = 0,linetype = "dashed") +
  geom_abline(intercept = 0 , slope = 1) +
  scale_fill_gradient(low = "white", high="darkgreen") +
  theme_bw() +
  theme(legend.position = "bottom") +
  removeGrid() +
  labs(x = expression("mean predation " * beta * " (cube root) at reference sites"),
       y = expression("mean predation " * beta * " (cube root) at MPAs"))


Beta_comp<-ggplot(ratio_complete,aes(x = prey_B_3root_ref, y = prey_B_3root_MPA)) +
  theme_bw() +
  geom_pointrange(aes(ymax = prey_B_3root_MPA, ymin = prey_B_3root_ref)) +
  geom_hline(yintercept = 0,linetype = "dashed") +
  geom_vline(xintercept = 0,linetype = "dashed") +
  geom_abline(intercept = 0 , slope = 1) +
  geom_point(color = "black", size = 3) +
  theme(legend.position = "none") +
  removeGrid() +
  labs(x = expression("mean competition " * beta * " (cube root) at reference sites"),
       y = expression("mean competition " * beta * " (cube root) at MPAs"))


Beta_DD<-ggplot(ratio_complete,aes(x = DD_B_3root_ref, y = DD_B_3root_MPA)) +
  theme_bw() +
  #geom_pointrange(aes(ymax = DD_B_3root_MPA, ymin = DD_B_3root_ref,color =  DD_B_3root_MPA > DD_B_3root_ref)) +
  geom_pointrange(aes(ymax = DD_B_3root_MPA, ymin = DD_B_3root_ref)) +
  geom_hline(yintercept = 0,linetype = "dashed") +
  geom_vline(xintercept = 0,linetype = "dashed") +
  geom_abline(intercept = 0 , slope = 1) +
  geom_point(color = "black", size = 3) +
  #scale_colour_manual(values = setNames(c('blue','red'),c(T, F))) +
  theme(legend.position = "none") +
  removeGrid() +
  lims(x = c(-3.5,1), y= c(-3.5,1)) +
  labs(x = expression("mean DD " * beta * " (cube root) at reference sites"),
       y = expression("mean DD " * beta * " (cube root) at MPAs"))




Beta_cons + Beta_comp + Beta_DD + plot_layout(ncol= 3)










# 6B - Metacommunity level interactions ####

# Reduced networks into list for function
cons_network_list<- split(filter_links_df_reduced, filter_links_df_reduced$network_ID)

# for testing model
#x = trophic_network %>% filter(network_ID == 7)

## Loop for calculation this ####
network_covariance_meta<-mclapply(cons_network_list,function(x){ 

  relevant_taxa<-x %>% 
    select(consumer_ID, resource_ID) %>% 
    gather() %>% 
    select(-key) %>% 
    rename(species = value)
  
  
  DBall<-DBall_meta_long %>% 
    filter(species %in% relevant_taxa$species) %>%
    group_by(island_status) %>% 
    mutate(site_ID = cur_group_id()) %>% 
    spread(species, cover, fill = 0)
  
  
  nsites=length(unique(DBall$site_ID)) #number of sites or repeats

  #ksite = 1 # for testing

  for (ksite in 1:nsites){

    DB=DBall[DBall$site_ID==ksite,] ## Select a site
    
    DB$site_ID<-NULL # remove site_ID; unnecessary after the first step
    head(DB)
    
    abundance_mat=as.matrix(DB[-c(1:3)]) ### Create matrix with time series of abundances
    sp=colnames(DB[-c(1:3)])
    nspecies=length(sp)
    
    ### Could use this with site instead of repeat
    dates=DB$year
    dates_bis=dates #already interpolated
    tab_sp=DB[,sp] # species abundance table
    
    #Setting MARSS model -- known topology of the adjacency matrix
    ### This will define number of species as a by-product
    interaction_matrix<-int_matrix %>% subset(colnames(int_matrix) %in% sp, rownames(int_matrix) %in% sp)
    int_matrix[sp, sp]
    
    mat_ord2 <- function(mx) mx[sort(rownames(mx)), c(sort(rownames(mx)), sort(setdiff(colnames(mx), rownames(mx))))] # old matrix
    n_species<-length(sp)
    
    # define all interactions in matrix as "1": all species interact; the context of interactions is defined later
    interaction_matrix=matrix(list(1),nrow=length(sp),
                              ncol=length(sp),
                              dimnames=list(sp,sp))

    B1=matrix(list(0),nrow=length(sp),ncol=length(sp),dimnames=list(rownames(interaction_matrix),colnames(interaction_matrix)))
    
    
# curate matrix
    for (i in 1:length(sp)){
      s=sp[i]
      for (j in 1:length(sp)){
        s2=sp[j]
        if(interaction_matrix[i,j]!=0){
          #B1[i,j]=paste("alpha_",i,",",j,sep="")
          B1[i,j]=paste(s2,s,sep="") ### beware, paste(s2,s,sep="") means we have the impacting species first  
        }else{
          if(s==s2){
            B1[i,j]=paste(s2,s,sep="")
          }
        }
      }
    }
    

    ### Defining other parameters
    iter_min=100 #
    iter_estimate=10 #not used for now
    iter_scenario=10
    iter_boot=1000
    BFGS=T
    
    U1="zero"
    Q1="diagonal and unequal" # changed from unequal 10.10.22
    Z1=diag(1,length(sp),length(sp))
    A1="zero"
    R1="zero"
    V1=diag(1,length(sp))
    pi1="zero"
    
    aalpha=0.05
    
    # Environmental effect is equal covariance 
    D <- "unconstrained"
    
    
    cntl.list=list(safe=TRUE, maxit=1500, allow.degen=TRUE) # from Dexter et al.
    
    tab_sp=log(tab_sp+1) # log +1 to deal with  NAs
    tab_sp=t(scale(tab_sp, scale = FALSE)) ### log and center the abundance data
    rownames(tab_sp)=sp
    
    
    # Addition of covariates to model - 4 August 2023
    temp_data<-meta_SST_df_annual %>% 
      filter(island_status %in% unique(DB$island_status),
             year %in% DB$year) %>%
      select(mean_SST)
    
    covariates<-temp_data$mean_SST
    covariate_matrix <- matrix(covariates, ncol = length(dates))
    tab_cov=t(scale(covariates, scale = FALSE)) # -- added on 4 Aug 2023

    d <- tab_cov
    
    model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1,D=D,d=d)
    
    # Actual model
    fit_log=MARSS(tab_sp, method="BFGS",model=model.list) ### Fit of MAR model. -- this takes a very long time with the whole food web (> 15 mins)
    
    coef_table <- coef(fit_log, type = "matrix")$B[1:n_species, 1:n_species]
    rownames(coef_table)<-sp
    colnames(coef_table)<-sp
    
    B_table = coef(fit_log,type="matrix")$B[1:n_species, 1:n_species]
    Q_table = coef(fit_log,type="matrix")$Q[1:n_species, 1:n_species]
    
    # Model results
    island_status <-unique(DB$island_status)
    network_ID <-unique(x$network_ID)
    coef_df_ready<-as.data.frame.table(coef_table)
    coef_df_ready$island_status <-unique(DB$island_status)
    coef_df_ready$Island <-unique(DB$Island)
    
    coef_df_ready$network_ID <-unique(x$network_ID)
    
    write.csv(coef_table,file=paste("./MARS/Mars_model_B/meta_Betas/Revised_AUG23_autoC/tables/B_points_",island_status,"_NET-",network_ID,"_table.csv",sep=""))
    write.csv(coef_df_ready,file=paste("./MARS/Mars_model_B/meta_Betas/Revised_AUG23_autoC/dfs/B_points_",island_status,"_NET-",network_ID,"_df.csv",sep=""),row.names = F)
    
    
  }
  }, mc.cores = 10, mc.preschedule = FALSE) # for mclapply




no_dfs<-length(unique(DBall_meta_long$island_status))*length(unique(trophic_network$network_ID))
no_dfs
# should be 552

## 6B.2 - Import metacommunity level interactions ####

# Bring files back in
B_meta_file_names <- list.files(path="./MARS/Mars_model_B/meta_Betas/Revised_AUG23_autoC/dfs", 
                                pattern=".csv", full.names = TRUE)

# import, bind, and name them
B_meta_dataframe <- do.call(rbind, lapply(B_meta_file_names, function(x) 
  cbind(read.csv(x, stringsAsFactors = FALSE), filename = x)))

# clean them to show network ID; wont be necessary next time we run this
B_meta_dataframe_mode<-B_meta_dataframe %>% 
  separate(filename, sep = "-", into = c("filename","network_ID_temp")) %>%
  separate(network_ID_temp, sep = "_", into = c("network_ID","extra"))

# select columns of interest
B_meta_dataframe_ready<-B_meta_dataframe_mode %>% select(1:5,7)


## Grouping ###

# group into density dependent 
B_meta_df_DD<- B_meta_dataframe_ready %>% filter(Var1 == Var2)
dim(B_meta_df_DD)

# group into competitive
B_meta_df_interspecific<- B_meta_dataframe_ready %>% filter(!Var1 == Var2)


# turn this into a function to pull values
#x<-trophic_network %>% filter(network_ID == 1) # for testing

cons_network_list<- split(trophic_network, trophic_network$network_ID)

B_meta_summarizing<-pblapply(cons_network_list,function(x){ #this function works and can be applied to list items with lapply
  
  # filter species for whether they were present in time series or not
  site_sums_species<-DBall_meta_long %>% 
    filter(species %in% c(x$consumer_ID,x$resource_ID)) %>%
    group_by(island_status,species) %>% 
    summarise(sum_cover = sum(cover,na.rm = T))
  
  B_meta_dataframe_ready_fornetwork<-B_meta_dataframe_ready %>% 
    filter(network_ID %in% x$network_ID) # filter for network
  
  B_by_site_list<- split(B_meta_dataframe_ready_fornetwork, B_meta_dataframe_ready_fornetwork$island_status)
  
  consumer = unique(x$consumer_ID)
  

  B_by_site<-lapply(B_by_site_list,function(x){
    
    sums_for_site<-site_sums_species %>% filter(island_status %in% x$island_status)
    
    # i want to know if there is only one species; predation strength might be stronger for one vs > 1 species
    comp_richness_results <- sums_for_site %>%
      filter(!species == consumer) %>%
      summarise(count_comp = n_distinct(species),
                count_comp_presnt = n_distinct((species)[sum_cover > 0]))
    
    
    # for competitors; summarizing interaction strength
    competitors_results<-x %>% 
      filter(Var1 %in% sums_for_site$species, 
             Var2 %in% sums_for_site$species) %>%
      filter(!Var1 == Var2) %>% 
      filter(!Var1 == consumer,
             !Var2 == consumer) %>% 
      summarise(mean_comp_B = mean(Freq, na.rm =T),
                sd_comp_B = sd(Freq,na.rm = T)) 
    
    competitors_results_merged<-merge(competitors_results,comp_richness_results)
    
    # for consumers; summarizing interaction strength
    consumer_results<-x %>% 
      filter(!Var1 == Var2) %>%
      filter(Var1 == consumer) %>%
      summarise(mean_cons_B = mean(Freq, na.rm =T),
                sd_cons_B = sd(Freq,na.rm = T))
    
    # density depedent interactions
    B_df_interspecific<- x %>% 
      filter(Var1 == Var2) %>%
      summarise(mean_DD_B = mean(Freq, na.rm =T),
                sd_DD_B = sd(Freq,na.rm = T))
    
    B_results<-Reduce(function(x, y) merge(x, y, all=TRUE), list(competitors_results_merged, consumer_results, B_df_interspecific))
    
    B_results$island_status<-unique(x$island_status)
    
    cons_present<- sums_for_site %>%
      filter(species %in% consumer) %>% 
      mutate(cons_present = if_else(sum_cover > 0, TRUE,FALSE))
    
    
    B_results$cons_present <- cons_present$cons_present
    
    B_results$network_ID<-unique(x$network_ID) 
    
    return(B_results)
  })
  
  B_by_site_df<-do.call(rbind.data.frame, B_by_site)
  
  return(B_by_site_df)
})


summarized_B_meta_df<-do.call(rbind.data.frame, B_meta_summarizing)

# filter out networks with less than one prey -- temporary?
summarized_B_meta_df<-summarized_B_meta_df %>% filter(count_comp > 1)
dim(summarized_B_meta_df)

# only look at interactions where a predator was present
dim(summarized_B_meta_df)
summarized_B_meta_df_pred_present<-summarized_B_meta_df %>% filter(cons_present==TRUE)
dim(summarized_B_meta_df_pred_present)



# check for normality
qqnorm(summarized_B_meta_df_pred_present$mean_comp_B)
qqnorm(summarized_B_meta_df_pred_present$mean_cons_B)
qqnorm(summarized_B_meta_df_pred_present$mean_DD_B)

summarized_B_meta_df_pred_present$pred_B_3root<-sign(summarized_B_meta_df_pred_present$mean_cons_B) * (abs(summarized_B_meta_df_pred_present$mean_cons_B))^(1/3)
summarized_B_meta_df_pred_present$prey_B_3root<-sign(summarized_B_meta_df_pred_present$mean_comp_B) * (abs(summarized_B_meta_df_pred_present$mean_comp_B))^(1/3)
summarized_B_meta_df_pred_present$DD_B_3root<-sign(summarized_B_meta_df_pred_present$mean_DD_B) * (abs(summarized_B_meta_df_pred_present$mean_DD_B))^(1/3)

qqnorm(summarized_B_meta_df_pred_present$pred_B_3root) # better
qqnorm(summarized_B_meta_df_pred_present$prey_B_3root) # better
#qqnorm(summarized_B_meta_df_pred_present$DD_B_3root) # better






# summarize to MPA/ref; island for every network - start without island

# this is likely unnecessary
summarized_B_meta_df_reduced<-summarized_B_meta_df_pred_present %>% 
  group_by(island_status,network_ID) %>% 
  summarize_if(is.numeric,mean,na.rm = T)



# curate for MPA status
B_ratio_prep_MPA<-summarized_B_df_reduced %>% filter(site_status == "MPA")
B_ratio_prep_ref<-summarized_B_df_reduced %>% filter(site_status == "reference")

names(B_ratio_prep_MPA)
colnames(B_ratio_prep_MPA) <-c("site_status","network_ID","mean_comp_B_MPA","sd_comp_B_MPA","count_comp_MPA",       
                               "count_comp_presnt_MPA","mean_cons_B_MPA","sd_cons_B_MPA","mean_DD_B_MPA","sd_DD_B_MPA",
                               "pred_B_3root_MPA","prey_B_3root_MPA","DD_B_3root_MPA")   

colnames(B_ratio_prep_ref) <-c("site_status","network_ID","mean_comp_B_ref","sd_comp_B_ref","count_comp_ref",       
                               "count_comp_presnt_ref","mean_cons_B_ref","sd_cons_B_ref","mean_DD_B_ref","sd_DD_B_ref",
                               "pred_B_3root_ref","prey_B_3root_ref","DD_B_3root_ref")   


# simple
ratio_complete<-merge(B_ratio_prep_MPA[-1],B_ratio_prep_ref[-1])

# merge with trophic position -- put this in curation (1) later 
trophic_updated<-pisco_trophic
trophic_updated$consumer<-trophic_updated$classcode

trophic_updated_formerge<-trophic_updated[c("PA_trophic","consumer")]
names(trophic_updated_formerge)<-c("PA_trophic","consumer_ID")

consumer_w_ID<-trophic_network %>% 
  group_by(network_ID,consumer_ID) %>% 
  summarize_all(mean) %>% 
  select(network_ID,consumer_ID) %>% merge(., trophic_updated_formerge)


ratio_complete_withnetwork<-merge(ratio_complete,consumer_w_ID, by = "network_ID", all= F)
dim(ratio_complete_withnetwork)






# merge synchrony data with interaction data
summarized_B_meta_df_reduced_updated<-summarized_B_meta_df_reduced %>% 
  separate(island_status,
           into=c("Island", "site_status"),
           sep="_", convert = TRUE, extra = "merge")

write.csv(summarized_B_meta_df_reduced_updated,"./Manuscript_scripts/Exported_data/metacommunity_interactions_REV.csv",row.names = F)




## 7- Trophic position ####

# Average trophic posiiton for prey; get trophic position for predators

pisco_trophic<-read.csv("./Data/PISCO/Trophic/pisco_trophic.csv")

# get just cons and res
consumer_network<-trophic_network %>% 
  dplyr::rename(classcode = consumer_ID) %>% 
  select(classcode, network_ID) %>% 
  group_by(classcode,network_ID) %>% 
  summarize_all(mean) 

resource_network<-trophic_network %>% 
  dplyr::mutate(classcode = resource_ID) %>% select(classcode, network_ID)

# merge with troph
consumer_network_troph<-merge(consumer_network,pisco_trophic) %>% 
  dplyr::rename(consumer_ID = classcode) %>%
  group_by(consumer_ID,network_ID) %>%
  dplyr::summarize(consumer_PA = mean(PA_trophic,na.rm=T))
  

resource_network_troph<-merge(resource_network,pisco_trophic) %>% 
  dplyr::rename(resource_ID = classcode)

# average res to network
resource_network_troph_sum<-resource_network_troph %>% 
  dplyr::group_by(network_ID) %>% 
  dplyr::summarize(resource_PA = mean(PA_trophic,na.rm=T))


# merge back together

network_trophic_position<-merge(consumer_network_troph,resource_network_troph_sum)

write.csv(network_trophic_position,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/Manuscript_scripts/Exported_data/network_trophic_position.csv", row.names = F)




# ready to bring into SEM script





# Extra for checking

predators_consumed<-filter_links_df_reduced %>% filter(species %in% networks_of_interest$consumer_ID)
not_predators<-filter_links_df_reduced %>% filter(!species %in% networks_of_interest$consumer_ID)

length(unique(predators_consumed$species))
length(unique(not_predators$species))

unique_competitors<-unique(competitors$resource_ID)

write.csv(unique_competitors,"./Manuscript_scripts/Exported_data/Tables/prey_eaten.csv", row.names = F)

### END ####