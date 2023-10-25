# Script for SEM path diagrams

# Community --> (1) across sites and networks; (2) across sites for each network
# Metacom --> (3) across metacoms and networks; (4) across metacoms for each network

library(ggExtra)
library(tidyverse)
library(viridis)
#library(scale_brewer)
library(ggraph)
library(cowplot)
library(tidygraph)
library(igraph)
library(ggraph)


# Set coordinate positions
coord_positions <- data.frame(name = c("MPA site status","CV SST","prey richness",
                                       "competition","predation","prey synchrony",
                                       "network stability"),
                              y = c(15,15,13,
                                    10,10,8,
                                    5),
                              x = c(1.5,4,6,
                                    1.5,3,5.5,
                                    4))



## (1) Com - all ####

# rename 
SEM_all_segment_plot <- SEM_all_plotting %>% 
  mutate(from = fct_recode(rhs,
                           "CV SST" = "CV_SST", 
                           "Hab. PC1"= "PC1",
                           "MPA site status" = "MPA",
                           "prey richness" = "count_comp_presnt",
                           "competition" = "prey_B_3root_inv",
                           "predation" = "pred_B_3root_inv",
                           "prey synchrony" = "synchrony",
                           "network stability" = "distance_inv",
                           "consumer TP" = "consumer_PA",
                           "resource TP" = "resource_PA"),
         to = fct_recode(lhs,
                         "CV SST" = "CV_SST", 
                         "Hab. PC1"= "PC1",
                         "MPA site status" = "MPA",
                         "prey richness" = "count_comp_presnt",
                         "competition" = "prey_B_3root_inv",
                         "predation" = "pred_B_3root_inv",
                         "prey synchrony" = "synchrony",
                         "network stability" = "distance_inv",
                         "consumer TP" = "consumer_PA",
                         "resource TP" = "resource_PA"))


int_covariance<-SEM_all_segment_plot %>% 
  filter(from == "predation" & to == "competition") %>% 
  mutate(to = fct_recode(to, "predation" = "competition"),
         from = fct_recode(from, "competition" = "predation"))

SEM_all_segment_plot_ready<-rbind(SEM_all_segment_plot,int_covariance)

hs_graph <- as_tbl_graph(SEM_all_segment_plot_ready, directed = T)
vertex.attributes(hs_graph)

# make fig
arrow_cov = arrow(length=unit(3,"cm"),
      ends="both",
      type = "closed")

arrow_norm = arrow(length=unit(3,"cm"),
                    ends="last",
                    type = "closed")

sem_path_all<-
  ggraph(hs_graph, coord_positions) + # success!!!!
  #ggraph(hs_graph,layout = "sugiyama") + # success!!!
  theme_bw() +
  removeGrid() +
    #coord_fixed() +
  geom_edge_link(aes(edge_color = ifelse(!is.na(sig) & sig == "yes", est.std, NA), # changed from est
                     alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5),
                     #arrow = ifelse(op == "~~", arrow_cov, arrow_norm)
                     ),
                 arrow = arrow(length = unit(3, 'mm'),angle = 20,
                               ends = "last",
                               type = "open"),
                 start_cap = circle(12, 'mm'),
                 end_cap = square(25, 'mm'),
                 width = 3,
                 #check_overlap = T,
                 angle_calc = 'none',
                 label_colour = "white",
                 #label_dodge = unit(0, 'mm'),
                 label_size = 4) +
  #geom_edge_link(arrow=ar,aes(start_cap=circle(6, 'mm'),
  #                            end_cap=square(10, 'mm'))) +
  labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable") +
  geom_node_point(aes(fill = name), fill = "black",size = 30, pch =21) +

  #scale_fill_manual(limits = c('MPA site status','CV SST', 'Hab. PC1','resource TP','Prey rich.', "consumer TP",
  #                             'predation','competition','prey synchrony','network stability'),
  #                    values = fill_colors_com$fill) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.74,0.07),
        legend.direction = "horizontal",
        legend.key.width=unit(0.8, 'cm'),
        legend.background = element_blank()) +
  guides(alpha = "none", fill = "none") +
  scale_edge_linetype(guide = "none") +
  scale_edge_alpha(guide = 'none') +
  scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                              values = scales::rescale(c(-0.5, -0.15,0, 0.15,0.5)),
                              limits = c(-0.5, 0.5), oob = scales::squish) +
  geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "white") +
  lims(x = c(1,6.5),
       y = c(2,17))


ggsave("./Manuscript_scripts/MS_figures/Figure_4.pdf",
       plot = sem_path_all,
       width = 6,
       height = 4)

ggsave("./Manuscript_scripts/MS_figures/Figure_4.png",
       plot = sem_path_all,
       width = 6,
       height = 4)







# Metacommunity ####

# Set coordinate positions
coord_positions_metacom <- data.frame(name = as.factor(c("MPA metacom","Reference metacom","CV SST","prey richness",
                                                         "competition","predation","prey synchrony",
                                                         "network stability")),
                                      y = c(15,15,15,13,
                                            9,9,8,
                                            4),
                                      x = c(2,4,6,8,
                                            1.5,4,7,
                                            5))




SEM_meta_all_segment_plot <- SEM_meta_all_plotting %>% 
  mutate(from = forcats::fct_recode(rhs,
                           "CV SST" = "CV_SST", 
                           #"Hab. PC1"= "PC1",
                           "MPA metacom" = "MPA",
                           "Reference metacom" = "reference",
                           #"All sites metacom" = "all",
                           "prey richness" = "count_comp_presnt",
                           "competition" = "prey_B_3root_inv",
                           "predation" = "pred_B_3root_inv",
                           "prey synchrony" = "synchrony",
                           #"network stability" = "distance_inv",
                           #"consumer TP" = "consumer_PA",
                           #"resource TP" = "resource_PA"
                           ),
         to = forcats::fct_recode(lhs,
                         #"CV SST" = "CV_SST", 
                         #"Hab. PC1"= "PC1",
                         #"MPA metacom" = "MPA",
                         #"Reference metacom" = "reference",
                         #"All sites metacom" = "all",
                         "prey richness" = "count_comp_presnt",
                         "competition" = "prey_B_3root_inv",
                         "predation" = "pred_B_3root_inv",
                         "prey synchrony" = "synchrony",
                         "network stability" = "distance_inv",
                         #"consumer TP" = "consumer_PA",
                         #"resource TP" = "resource_PA"
                         ))


int_covariance_meta<-SEM_meta_all_segment_plot %>% 
  filter(from == "predation" & to == "competition") %>% 
  mutate(to = fct_recode(to, "predation" = "competition"),
         from = fct_recode(from, "competition" = "predation"))

SEM_meta_all_segment_plot_ready<-rbind(SEM_meta_all_segment_plot,int_covariance_meta)

# convert to graph object 
hs_meta_all_graph <- as_tbl_graph(SEM_meta_all_segment_plot_ready, directed = T)

name_order<-data.frame(name = vertex.attributes(hs_meta_all_graph))

coord_positions_metacom_orderd<-coord_positions_metacom[order(match(coord_positions_metacom[,1],name_order[,1])),]


sem_path_meta_all<-
  ggraph(hs_meta_all_graph,coord_positions_metacom_orderd) +
  theme_bw() +
  removeGrid() +
    geom_edge_link(aes(edge_color =  ifelse(!is.na(sig) & sig == "yes", est.std, NA), # changed from est
                       alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5)),
                   arrow = arrow(length = unit(3, 'mm'),angle = 20,
                                 ends = "last",
                                 type = "open"),
                   start_cap = circle(12, 'mm'),
                   end_cap = square(25, 'mm'),
                   width = 3,
                   angle_calc = 'none',
                   label_colour = "white",
                   label_size = 4) +
    labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable") +
    geom_node_point(aes(fill = name), fill = "black",size = 30, pch =21) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = c(0.78,0.08),
          legend.direction = "horizontal",
          legend.background = element_blank()) +
    guides(alpha = "none", fill = "none") +
    scale_edge_linetype(guide = "none") +
    scale_edge_alpha(guide = 'none') +
    scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                                values = scales::rescale(c(-0.6, -0.05,0, 0.05,0.6)),
                                limits = c(-0.8, 0.8), oob = scales::squish) +
    geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "white") +
  lims(x = c(.5,9),
       y = c(1,17))




ggsave("./Manuscript_scripts/MS_figures/Figure_5.pdf",
       plot = sem_path_meta_all,
       width = 6,
       height = 4)

ggsave("./Manuscript_scripts/MS_figures/Figure_5.png",
       plot = sem_path_meta_all,
       width = 6,
       height = 4)


## Metacom - all - supplement comparison ####

coord_positions_metacom_MvR <- data.frame(name = as.factor(c("MPA metacom","All sites metacom","CV SST","prey richness",
                                                         "competition","predation","prey synchrony",
                                                         "network stability")),
                                      y = c(15,15,15,13,
                                            9,9,8,
                                            4),
                                      x = c(2,4,6,8,
                                            1.5,4,7,
                                            5))




SEM_meta_all_segment_MvR_plot <- SEM_meta_all_plotting_MvR %>% 
  mutate(from = forcats::fct_recode(rhs,
                                    "CV SST" = "CV_SST", 
                                    #"Hab. PC1"= "PC1",
                                    "MPA metacom" = "MPA",
                                    "All sites metacom" = "all",
                                    "prey richness" = "count_comp_presnt",
                                    "competition" = "prey_B_3root_inv",
                                    "predation" = "pred_B_3root_inv",
                                    "prey synchrony" = "synchrony",
                                    #"network stability" = "distance_inv",
                                    #"consumer TP" = "consumer_PA",
                                    #"resource TP" = "resource_PA"
  ),
  to = forcats::fct_recode(lhs,
                           #"CV SST" = "CV_SST", 
                           #"Hab. PC1"= "PC1",
                           #"MPA metacom" = "MPA",
                           #"Reference metacom" = "reference",
                           #"All sites metacom" = "all",
                           "prey richness" = "count_comp_presnt",
                           "competition" = "prey_B_3root_inv",
                           "predation" = "pred_B_3root_inv",
                           "prey synchrony" = "synchrony",
                           "network stability" = "distance_inv",
                           #"consumer TP" = "consumer_PA",
                           #"resource TP" = "resource_PA"
  ))


int_covariance_MvR_meta<-SEM_meta_all_segment_MvR_plot %>% 
  filter(from == "predation" & to == "competition") %>% 
  mutate(to = fct_recode(to, "predation" = "competition"),
         from = fct_recode(from, "competition" = "predation"))

SEM_meta_all_segment_MvR_plot_ready<-rbind(SEM_meta_all_segment_MvR_plot,int_covariance_meta)

# convert to graph object 
hs_meta_all_MvR_graph <- as_tbl_graph(SEM_meta_all_segment_MvR_plot_ready, directed = T)

name_order_MvR<-data.frame(name = vertex.attributes(hs_meta_all_MvR_graph))

coord_positions_metacom_MvR_orderd<-coord_positions_metacom_MvR[order(match(coord_positions_metacom_MvR[,1],name_order_MvR[,1])),]


sem_path_meta_all_MvR<-
  ggraph(hs_meta_all_MvR_graph,coord_positions_metacom_MvR_orderd) +
  theme_bw() +
  removeGrid() +
  geom_edge_link(aes(edge_color =  ifelse(!is.na(sig) & sig == "yes", est.std, NA), 
                     alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5)),
                 arrow = arrow(length = unit(3, 'mm'),angle = 20,
                               ends = "last",
                               type = "open"),
                 start_cap = circle(12, 'mm'),
                 end_cap = square(25, 'mm'),
                 width = 3,
                 angle_calc = 'none',
                 label_colour = "white",
                 label_size = 4) +
  labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable") +
  geom_node_point(aes(fill = name), fill = "black",size = 30, pch =21) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.78,0.08),
        legend.direction = "horizontal",
        legend.background = element_blank()) +
  guides(alpha = "none", fill = "none") +
  scale_edge_linetype(guide = "none") +
  scale_edge_alpha(guide = 'none') +
  scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                              values = scales::rescale(c(-0.6, -0.05,0, 0.05,0.6)),
                              limits = c(-0.8, 0.8), oob = scales::squish) +
  geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "white") +
  lims(x = c(.5,9),
       y = c(1,17))




ggsave("./Manuscript_scripts/MS_figures/Fig_S3_MvR.pdf",
       plot = sem_path_meta_all_MvR,
       width = 6,
       height = 4)


ggsave("./Manuscript_scripts/MS_figures/Figure_S3_MvR.png",
       plot = sem_path_meta_all_MvR,
       width = 6,
       height = 4)





# In black for presentations ####

sem_path_all_black<-
  ggraph(hs_graph, coord_positions) + # success!!!!
  #ggraph(hs_graph,layout = "sugiyama") + # success!!!
  theme_bw() +
  removeGrid() +
  #coord_fixed() +
  geom_edge_link(aes(edge_color = ifelse(!is.na(sig) & sig == "yes", est.std, NA), # changed from est
                     alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5),
                     #arrow = ifelse(op == "~~", arrow_cov, arrow_norm)
  ),
  arrow = arrow(length = unit(3, 'mm'),angle = 20,
                ends = "last",
                type = "open"),
  start_cap = circle(12, 'mm'),
  end_cap = square(25, 'mm'),
  width = 3,
  #check_overlap = T,
  angle_calc = 'none',
  label_colour = "white",
  #label_dodge = unit(0, 'mm'),
  label_size = 4) +
  #geom_edge_link(arrow=ar,aes(start_cap=circle(6, 'mm'),
  #                            end_cap=square(10, 'mm'))) +
  labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable") +
  geom_node_point(aes(fill = name), fill = "white",size = 30, pch =21) +
  
  #scale_fill_manual(limits = c('MPA site status','CV SST', 'Hab. PC1','resource TP','Prey rich.', "consumer TP",
  #                             'predation','competition','prey synchrony','network stability'),
  #                    values = fill_colors_com$fill) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.74,0.07),
        legend.direction = "horizontal",
        legend.key.width=unit(0.8, 'cm'),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "black", color  =  "black"),  
        panel.border = element_rect(fill = NA, color = "black"),  
        plot.background = element_rect(color = "black", fill = "black"),
        panel.margin = unit(0.5, "lines"),
        strip.background = element_rect(fill = "grey30", color = "grey10"),
        text = element_text(color = "white")) +
  guides(alpha = "none", fill = "none") +
  scale_edge_linetype(guide = "none") +
  scale_edge_alpha(guide = 'none') +
  scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                              values = scales::rescale(c(-0.5, -0.15,0, 0.15,0.5)),
                              limits = c(-0.5, 0.5), oob = scales::squish) +
  geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "black") +
  lims(x = c(1,6.5),
       y = c(2,17))


ggsave("./Manuscript_scripts/MS_figures/Fig_3_black.pdf",
       plot = sem_path_all_black,
       width = 6,
       height = 4)




sem_path_meta_all_black<-
  #ggraph(hs_meta_all_graph, layout = "sugiyama") +
  ggraph(hs_meta_all_graph,coord_positions_metacom_orderd) +
  theme_bw() +
  removeGrid() +
  geom_edge_link(aes(edge_color =  ifelse(!is.na(sig) & sig == "yes", est.std, NA), # changed from est
                     #width =  ifelse(!is.na(sig) & sig == "yes", est.std, NA),
                     #linetype = ifelse(!is.na(sig) & sig == "yes", "dotted", 'solid'),
                     #label = ifelse(!is.na(sig) & sig == "yes", ifelse(est.std > 0, "+", "-"),""), # changed from est
                     alpha = ifelse(!is.na(sig) & sig == "yes", 1,0.5)),
                 arrow = arrow(length = unit(3, 'mm'),angle = 20,
                               ends = "last",
                               type = "open"),
                 start_cap = circle(12, 'mm'),
                 end_cap = square(25, 'mm'),
                 width = 3,
                 #check_overlap = T,
                 angle_calc = 'none',
                 label_colour = "white",
                 #label_dodge = unit(0, 'mm'),
                 label_size = 4) +
  #geom_edge_link(arrow=ar,aes(start_cap=circle(6, 'mm'),
  #                            end_cap=square(10, 'mm'))) +
  labs(edge_color = "Coef. estimate", linetype = "none", fill = "variable") +
  geom_node_point(aes(fill = name), fill = "white",size = 30, pch =21) +
  
  #scale_fill_manual(limits = c('MPA site status','CV SST', 'Hab. PC1','resource TP','Prey rich.', "consumer TP",
  #                             'predation','competition','prey synchrony','network stability'),
  #                    values = fill_colors_com$fill) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.78,0.08),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        panel.background = element_rect(fill = "black", color  =  "black"),  
        panel.border = element_rect(fill = NA, color = "black"),  
        plot.background = element_rect(color = "black", fill = "black"),
        panel.margin = unit(0.5, "lines"),
        strip.background = element_rect(fill = "grey30", color = "grey10"),
        text = element_text(color = "white")) +
  guides(alpha = "none", fill = "none") +
  scale_edge_linetype(guide = "none") +
  scale_edge_alpha(guide = 'none') +
  scale_edge_colour_gradientn(colors = c(low = "blue",mid = "white",high = "red"), 
                              values = scales::rescale(c(-0.6, -0.05,0, 0.05,0.6)),
                              limits = c(-0.8, 0.8), oob = scales::squish) +
  geom_node_text(aes(label = stringr::str_wrap(name, 8)), color = "black") +
  lims(x = c(.5,9),
       y = c(1,17))




ggsave("./Manuscript_scripts/MS_figures/Fig_5_black.pdf",
       plot = sem_path_meta_all_black,
       width = 6,
       height = 4)

# END #
