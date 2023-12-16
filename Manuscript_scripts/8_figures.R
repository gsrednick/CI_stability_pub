## Supplementary figures ####

library(emmeans)
library(lme4)
library(lmerTest)
library(tidyverse)
library(multcomp)
library(ggExtra)
library(patchwork)

# bring in site table
site_table<-read.csv("./Manuscript_scripts/Exported_data/site_table.csv")

#### Figure 3 - Community metrics  #### 

# Step 2 - Interaction strength within networks: MPA vs. ref
# currently figure 2 in stack
B_sem_df_filtered<- B_sem_df %>% filter(network_ID %in% networks_of_interest$network_ID)

## ANOVA for this
B_sem_df_filtered$Island<-factor(B_sem_df_filtered$Island, levels = c("Anacapa","SCI","SRI","SMI"))

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

# Stability (dispersion)
summary(aov(distance_inv~site_status*Island,B_sem_df_filtered))

com_stab_lmer<-lmer(distance_inv~site_status*Island + (1|network_ID),B_sem_df_filtered, REML = F)
summary(com_stab_lmer)
anova(com_stab_lmer)
ranova(com_stab_lmer)

#stab_emmeans <- lsmeans(com_stab_lmer, ~ site_status)
stab_com_PW<- emmeans(com_stab_lmer, list(pairwise ~ site_status*Island), adjust = "tukey")

com_stab_cld <- cld(stab_com_PW$`emmeans of site_status, Island`, Letters = letters)


# add replicates to x-axis
count_levels<-B_sem_df_filtered %>% 
  group_by(site_status,Island) %>% 
  summarize(count = length(distance_inv))

xlabs <- paste(unique(B_sem_df_filtered$site_status),
               "\n(N=",table(B_sem_df_filtered$site_status,B_sem_df_filtered$Island),")",sep="")


com_facet_labels <- B_sem_df_filtered %>%
  group_by(Island, site_status) %>%
  summarize(count = n()) %>%
  mutate(com_xlabs = paste(site_status, "\n(n=", count, ")", sep = ""))


B_sem_df_filtered_plotting <- B_sem_df_filtered %>%
  left_join(com_facet_labels, by = c("Island", "site_status"))

com_stab_cld<-merge(com_stab_cld,com_facet_labels)

ANOVA_stab<-B_sem_df_filtered_plotting %>%
  ggplot(aes(x = factor(com_xlabs), y = distance_inv, color = site_status)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  scale_color_manual(values = c(MPA = "red", reference = "blue")) +
  facet_grid(~Island, scales = "free_x") +
  theme_bw() +
  removeGrid() +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  labs(y = "network stability", x = "site status",  color = "site status") +
  geom_text(data = com_stab_cld,  # DataFrame with x, y, and letters
            aes(x = com_xlabs, y = -0.15, label = .group))


# Synchrony 
summary(aov(synchrony~site_status*Island,B_sem_df_filtered))

com_sync_lmer<-lmer(synchrony~site_status*Island + (1|network_ID),B_sem_df_filtered, REML = F)
anova(com_sync_lmer)
ranova(com_sync_lmer)

com_sync_PW<- emmeans(com_sync_lmer, list(pairwise ~ site_status*Island), adjust = "tukey")
com_sync_cld <- cld(com_sync_PW$`emmeans of site_status, Island`, Letters = letters)


ANOVA_sync<-ggplot(B_sem_df_filtered, aes(x = site_status, y = synchrony, color = site_status)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  scale_color_manual(values = c(MPA = "red", reference = "blue")) +
  facet_grid(~Island) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  removeGrid() +
  #labs(y = "mean predator-prey\n interaction strength (cube-root)") +
  labs(y = "prey synchrony", x = "site status",  color = "site status") +
  geom_text(data = com_sync_cld,  # DataFrame with x, y, and letters
            aes(x = site_status, y = 1.1, label = .group))




# Predation - no difference in B for predator-prey between MPAs and refs; yes islands
# Anacapa shows more positive predation coefficients
summary(aov(pred_B_3root~site_status*Island,B_sem_df_filtered))
summary(aov(mean_cons_B~site_status*Island,B_sem_df_filtered)) # untransformed is significant

com_pred_lmer<-lmer(pred_B_3root~site_status*Island + (1|network_ID),B_sem_df_filtered, REML = F)
summary(com_pred_lmer)
anova(com_pred_lmer)
ranova(com_pred_lmer)

com_pred_PW<- emmeans(com_pred_lmer, list(pairwise ~ site_status*Island), adjust = "tukey")
com_pred_cld <- cld(com_pred_PW$`emmeans of site_status, Island`, Letters = letters)


ANOVA_predation<-ggplot(B_sem_df_filtered, aes(x = site_status, y = pred_B_3root_inv, color = site_status)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  scale_color_manual(values = c(MPA = "red", reference = "blue")) +
  facet_grid(~Island) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  removeGrid() +
  #labs(y = "mean predator-prey\n interaction strength (cube-root)") +
  labs(y = expression("predation " * beta), x = "site status",  color = "site status") +
  geom_text(data = com_pred_cld,  # DataFrame with x, y, and letters
            aes(x = site_status, y = 2.4, label = .group))




# Competition -  difference in B for competition between MPAs and refs; also varies among islands
# difference in competitive dynamics stronger at SCI and SRI
summary(aov(prey_B_3root~site_status*Island,B_sem_df_filtered))
summary(aov(mean_comp_B~site_status*Island,B_sem_df_filtered)) # untransformed is significant
TukeyHSD(aov(prey_B_3root~site_status*Island,B_sem_df_filtered))

com_prey_lmer<-lmer(prey_B_3root~site_status*Island + (1|network_ID),B_sem_df_filtered, REML = F)
anova(com_prey_lmer)
ranova(com_prey_lmer)

com_prey_PW <- emmeans(com_prey_lmer, list(pairwise ~ site_status * Island), adjust = "tukey")
com_prey_cld <- cld(com_prey_PW$`emmeans of site_status, Island`, Letters = letters)

ANOVA_competition<-ggplot(B_sem_df_filtered, aes(x = site_status, y = prey_B_3root_inv, color = site_status)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c(MPA = "red", reference = "blue")) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  facet_grid(~Island) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  removeGrid() +
  #labs(y = "mean interspecific competition \n interaction strength (cube-root)") +
  labs(y = expression("competiton " * beta), x = "site status", color = "site status") +
  geom_text(data = com_prey_cld,  # DataFrame with x, y, and letters
            aes(x = site_status, y = 2, label = .group))




# Density dependence -  difference in B for DD between islands
# Negative DD weakest at Anacapa -- still all negative

summary(aov(DD_B_3root~site_status*Island,B_sem_df_filtered))
summary(aov(mean_comp_B~site_status*Island,B_sem_df_filtered))

TukeyHSD(aov(mean_comp_B~site_status*Island,B_sem_df_filtered))



ANOVA_DD<-ggplot(B_sem_df_filtered, aes(x = site_status, y = DD_B_3root, color = site_status)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c(MPA = "red", reference = "blue")) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  facet_grid(~Island) +
  theme_bw() +
  removeGrid() +
  #labs(y = "mean intraspecific density-depedent\n interaction strength (cube-root)") +
  labs(y = expression("mean density-dependent " * beta * "( cube root)"))





Figure_3_comp<-ANOVA_competition + ANOVA_predation + ANOVA_sync+ ANOVA_stab + plot_layout(ncol = 1, guides = "collect") & plot_annotation(tag_levels = "a") & theme(legend.position = "none")

ggsave("./Manuscript_scripts/MS_figures/Figure_3.pdf",
       plot = Figure_3_comp,
       width = 8,
       height = 8)

ggsave("./Manuscript_scripts/MS_figures/Figure_3.jpeg",
       plot = Figure_3_comp,
       width = 8,
       height = 8)

ggsave("./Manuscript_scripts/MS_figures/Figure_3.png",
       plot = Figure_3_comp,
       width = 8,
       height = 8,
       dpi = 300)


# community summary statistics

library(summarytools)


B_sem_df_filtered_summary <- B_sem_df_filtered %>%
  dplyr::group_by(Island, site_status) %>%
  dplyr::select(site,site_status,pred_B_3root,prey_B_3root,synchrony,distance_inv) %>%
  summarise(across(where(is.numeric),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE),
                        n = ~n())),
            .groups = "drop")


#B_sem_df_filtered_summary_CI <- B_sem_df_filtered_summary %>%
#  mutate(across(-c(Island, site_status),
#                list(se = ~sd / sqrt(n),
#                     lower.ci = ~mean - qt(1 - (0.05 / 2), n - 1) * (sd / sqrt(n)),
#                     upper.ci = ~mean + qt(1 - (0.05 / 2), n - 1) * (sd / sqrt(n)))))




B_sem_df_filtered_summary <- B_sem_df_filtered %>%
  dplyr::select(site,site_status,pred_B_3root,prey_B_3root,synchrony,distance_inv) %>%
  dplyr::group_by(site, site_status) %>%
  dplyr::summarise(across(where(is.numeric), list(mean = ~mean(., na.rm = TRUE), 
                                                  sd = ~sd(., na.rm = TRUE),
                                                  se = ~sd(., na.rm = TRUE) / sqrt(sum(!is.na(.))),
                                                  n = ~sum(!is.na(.)))))

# Calculate confidence intervals using summarytools::conf_int()
#B_sem_df_filtered_summary <- B_sem_df_filtered_summary %>%
#  mutate_at(vars(-group_cols()), list(conf_int = ~conf_int(.)))

# Print the summary
B_sem_df_filtered_summary


#### Figure SXX - Metacommunity metrics  #### 
B_meta_sem_df_filtered<- B_meta_sem_df %>% 
  dplyr::filter(network_ID %in% networks_of_interest$network_ID)

island_status_list<-B_meta_sem_df_filtered %>% 
  dplyr::group_by(Island,site_status) %>% 
  dplyr::summarize()



## ANOVA for this
B_meta_sem_df_filtered$Island<-factor(B_meta_sem_df_filtered$Island, levels = c("Anacapa","SCI","SRI","SMI"))

# add replicates to x-axis
#meta_count_levels<-B_meta_sem_df_filtered %>% group_by(site_status,Island) %>% summarize(count = length(distance_inv))



# Stability (dispersion)
summary(aov(distance_inv~site_status*Island,B_meta_sem_df_filtered))
TukeyHSD(aov(distance_inv~site_status*Island,B_meta_sem_df_filtered))

meta_stab_lmer<-lmer(distance_inv~site_status*Island + (1|network_ID),B_meta_sem_df_filtered, REML = F)
anova(meta_stab_lmer)
ranova(meta_stab_lmer)

meta_stab_PW <- emmeans(meta_stab_lmer, list(pairwise ~ site_status * Island), adjust = "tukey")
meta_stab_cld <- cld(meta_stab_PW$`emmeans of site_status, Island`, Letters = letters)

facet_labels <- B_meta_sem_df_filtered %>%
  dplyr::group_by(Island, site_status) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::mutate(meta_xlabs = paste(site_status, "\n(N=", count, ")", sep = ""))

B_meta_sem_df_filtered_plotting <- B_meta_sem_df_filtered %>%
  left_join(facet_labels, by = c("Island", "site_status"))

meta_stab_cld<-merge(meta_stab_cld,facet_labels)


ANOVA_stab_meta<-ggplot(B_meta_sem_df_filtered_plotting, aes(x = factor(meta_xlabs, levels = unique(meta_xlabs)), y = distance_inv, color = site_status)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  scale_color_manual(values = c(MPA = "red", reference = "blue", all = "black")) +
  facet_grid(~Island, scales = "free_x") +
  theme_bw() +
  removeGrid() +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    #axis.title.x = element_blank(),
  #      axis.text.x = element_blank()) +
  #labs(y = "mean predator-prey\n interaction strength (cube-root)") +
  labs(x = "site status", y = "network stability",color = "Metacom. status") +
  #scale_x_discrete(labels = meta_xlabs) +
  geom_text(data = meta_stab_cld,  # DataFrame with x, y, and letters
            aes(x = meta_xlabs, y = -0.14, label = .group))


# Synchrony 
meta_sync_lmer<-lmer(synchrony~site_status*Island + (1|network_ID),B_meta_sem_df_filtered, REML = F)
anova(meta_sync_lmer)
ranova(meta_sync_lmer)

meta_sync_PW <- emmeans(meta_sync_lmer, list(pairwise ~ Island* site_status), adjust = "tukey")
meta_sync_cld <- cld(meta_sync_PW$`emmeans of Island, site_status`, Letters = letters)

meta_sync_PW_is <- emmeans(meta_sync_lmer, list(pairwise ~ site_status), adjust = "tukey")
pairs(meta_sync_PW_is)

summary(aov(synchrony~site_status*Island + (1|network_ID),B_meta_sem_df_filtered))
summary(aov(synchrony~site_status*Island,B_meta_sem_df_filtered))


ANOVA_sync_meta<-ggplot(B_meta_sem_df_filtered, aes(x = site_status, y = synchrony, color = site_status)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  scale_color_manual(values = c(MPA = "red", reference = "blue", all = "black")) +
  facet_grid(~Island) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  removeGrid() +
  #labs(y = "mean predator-prey\n interaction strength (cube-root)") +
  labs(y = "prey synchrony",color = "Metacom. status") +
  geom_text(data = meta_sync_cld,  # DataFrame with x, y, and letters
            aes(x = site_status, y = 1.1, label = .group))




# Predation - no difference in B for predator-prey between MPAs and refs; yes islands
# Anacapa shows more positive predation coefficients
meta_pred_lmer<-lmer(pred_B_3root~site_status*Island + (1|network_ID),B_meta_sem_df_filtered, REML = F)
anova(meta_pred_lmer)
ranova(meta_pred_lmer)

summary(aov(pred_B_3root~site_status*Island,B_meta_sem_df_filtered))
summary(aov(mean_cons_B~site_status*Island,B_meta_sem_df_filtered)) # untransformed is significant

meta_pred_PW <- emmeans(meta_pred_lmer, list(pairwise ~ Island * site_status), adjust = "tukey")
meta_pred_cld <- cld(meta_pred_PW$`emmeans of Island, site_status`, Letters = letters)


ANOVA_predation_meta<-ggplot(B_meta_sem_df_filtered, aes(x = site_status, y = pred_B_3root_inv, color = site_status)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  scale_color_manual(values = c(MPA = "red", reference = "blue", all = "black")) +
  facet_grid(~Island) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  removeGrid() +
  #labs(y = "mean predator-prey\n interaction strength (cube-root)") +
  labs(y = expression("predation " * beta), color = "Metacom. status") +
  geom_text(data = meta_pred_cld,  # DataFrame with x, y, and letters
            aes(x = site_status, y = 2.5, label = .group))



# Competition -  difference in B for competition between MPAs and refs; also varies among islands
# difference in competitive dynamics stronger at SCI and SRI
meta_prey_lmer<-lmer(prey_B_3root~site_status*Island + (1|network_ID),B_meta_sem_df_filtered, REML = F)
anova(meta_prey_lmer)
ranova(meta_prey_lmer)

summary(aov(prey_B_3root~site_status*Island,B_meta_sem_df_filtered))
#summary(aov(mean_comp_B~site_status*Island,B_meta_sem_df_filtered)) # untransformed is significant
#TukeyHSD(aov(prey_B_3root~site_status*Island,B_meta_sem_df_filtered))

meta_prey_PW <- emmeans(meta_prey_lmer, list(pairwise ~ Island * site_status), adjust = "tukey")
meta_prey_cld <- cld(meta_prey_PW$`emmeans of Island, site_status`, Letters = letters)

ANOVA_competition_meta<-ggplot(B_meta_sem_df_filtered, aes(x = site_status, y = prey_B_3root_inv, color = site_status)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c(MPA = "red", reference = "blue", all = "black")) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  facet_grid(~Island) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  removeGrid() +
  #labs(y = "mean interspecific competition \n interaction strength (cube-root)") +
  labs(y = expression("competiton " * beta), color = "Metacom. status") +
  geom_text(data = meta_prey_cld,  # DataFrame with x, y, and letters
            aes(x = site_status, y = 2, label = .group))



metacom_compar<-ANOVA_competition_meta + ANOVA_predation_meta + ANOVA_sync_meta + ANOVA_stab_meta + plot_layout(ncol = 1, guides = "collect") & plot_annotation(tag_levels = "a") & theme(legend.position = "top")

ggsave("./Manuscript_scripts/MS_figures/Fig_S2_metacom_metrics.pdf",
       plot = metacom_compar,
       width = 8,
       height = 9)

ggsave("./Manuscript_scripts/MS_figures/Figure_S2.png",
       plot = metacom_compar,
       width = 8,
       height = 9)



 
 # (5) MDS for pred-prey communities ####

 pred_prey_mds_df<- read.csv("./Manuscript_scripts/Exported_data/mds_points_community.csv")
 
 mds_net<-merge(pred_prey_mds_df,site_table) 
 
 mds_net %>% filter(year == 2015) %>%
   ggplot(aes(x=MDS1,y=MDS2)) + 
   #geom_errorbarh(aes(xmax = mean_MDS1 + se_MDS1, xmin = mean_MDS1 - se_MDS1), color = "grey") + # horizontal error bars (SE) 
   #geom_errorbar(aes(ymax = mean_MDS2 + se_MDS2, ymin = mean_MDS2 - se_MDS2), color = "grey") + # vertical error bars (SE) 
   geom_point(aes(fill = site_status, alpha = year),pch = 21,size=6) +  # actual points
   #geom_path(aes(color=site_status), show.legend = F) + # path that connects each point in time
   theme_bw() + # theme
   removeGrid() + # no grid please
   #geom_text(aes(label = year)) +
   theme(text = element_text(size = 14, color = "black"), # theme formatting 
         axis.text.x = element_text(color="black"), 
         axis.text.y = element_text(color="black")) +
   labs(fill = "MPA site status", x = "MDS1",y="MDS2") + # labels 
   scale_colour_manual(values=c("MPA" = "red","reference" = "blue"))+ # manual coloring
   scale_fill_manual(values=c("MPA" = "red","reference" = "blue"))+ # manual coloring
   scale_alpha_continuous(name = "year", range = c(0.1, 1)) +
   facet_wrap(~Island, ncol = 2)
 
 
 #### MDS  with averages and error bars
 ## For plotting means +/- SE on objects
 se_fn <- function(x) sd(x) / sqrt(length(x)) # Create own function
 
 summarized_mds <- mds_net %>% 
   group_by(year,site_status,Island) %>%
   dplyr::summarize(mean_MDS1 = mean(MDS1),
                    mean_MDS2 = mean(MDS2),
                    se_MDS1 = se_fn(MDS1),
                    se_MDS2 = se_fn(MDS2))
 
 
 # plot
 
 summarized_mds$Island<-factor(summarized_mds$Island, levels = c("Anacapa", "SCI", "SRI","SMI"))
 
 
 summarized_MDS<-ggplot(summarized_mds[order(summarized_mds$year),], aes(x=mean_MDS1,y=mean_MDS2)) + 
   geom_errorbarh(aes(xmax = mean_MDS1 + se_MDS1, xmin = mean_MDS1 - se_MDS1), color = "grey") + # horizontal error bars (SE) 
   geom_errorbar(aes(ymax = mean_MDS2 + se_MDS2, ymin = mean_MDS2 - se_MDS2), color = "grey") + # vertical error bars (SE) 
   geom_point(aes(fill = site_status, alpha = year),pch = 21,size=6) +  # actual points
   geom_path(aes(color=site_status), show.legend = F) + # path that connects each point in time
   theme_bw() + # theme
   removeGrid() + # no grid please
   theme(text = element_text(size = 14, color = "black"), # theme formatting 
         axis.text.x = element_text(color="black"), 
         axis.text.y = element_text(color="black")) +
   labs(fill = "MPA site status", x = "MDS1",y="MDS2") + # labels 
   scale_colour_manual(values=c("MPA" = "red","reference" = "blue"))+ # manual coloring
   scale_fill_manual(values=c("MPA" = "red","reference" = "blue"))+ # manual coloring
   scale_alpha_continuous(name = "year", range = c(0.1, 1)) +
   facet_wrap(~Island, ncol = 2) # facet for another treatment level
 
 
 
 
 # make each color a scale by time
 
 ggsave(
   #"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/Manuscript_scripts/MS_figures/community_MDS.pdf",
    "./Manuscript_scripts/MS_figures/Figure_2_MDS.pdf",
       plot = summarized_MDS,
        width = 11,
        height = 8)
 
 ggsave(
   #"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/CI_asynchrony_MPA/Manuscript_scripts/MS_figures/community_MDS.pdf",
   "./Manuscript_scripts/MS_figures/Figure_2.png",
   plot = summarized_MDS,
   width = 11,
   height = 8,
   dpi = 300)
 
 
 
 

 # From reviewer ####
 
 # calculate percent overlap of trophic networks
 
 head(trophic_network)
 
 TN_prey_count<-trophic_network %>% group_by(network_ID) %>%
   summarize(no_prey = n())
 
 overlap_data <- trophic_network %>%
   group_by(resource_ID) %>%
   summarise(
     unique_networks = n_distinct(network_ID),
     total_networks = n_distinct(trophic_network$network_ID),
     percent_overlap = (unique_networks / total_networks) * 100
   )
 
 
 # this is the mean ± SE percent overlap across trophic networks 
 overlap_data %>% summarize(mean = mean(percent_overlap),
                            se = sd(percent_overlap) / sqrt(n()))
 
 
 
 #### Package reporting #####
 # The following is a .bib download of all the packages used for analyses
 
 library(knitr)
 # get all packages for citation
 loaded_packages <- search()[-1]
 
 # Get session information
 session_info <- sessionInfo()
 
 # Extract loaded package names
 loaded_packages <- names(session_info$otherPkgs)
 
 write_bib(loaded_packages, file = "R-GUI-pkgs.bib")
 
 
 
 # END ####