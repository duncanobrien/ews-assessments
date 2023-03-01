########################################################################################################################
# Figure 1 #
########################################################################################################################
require(tidyverse)
require(foreach)
require(patchwork)

source("Data/kinneret_plankton_data.R") #raw data not provided. See references in Data Availability statement and repository README
source("Data/kinneret_environmental_data.R")

source("Data/windermere_plankton_data.R")
source("Data/windermere_environmental_data.R")

rm(list = ls()[!ls() %in% c("plank_env.data.yr","phyto_env.windyrdata")])

source("Code/threshold_gam.R")

############ 
#Prepare plankton data
############ 
#for Lake Kinneret and Windermere, find total plankton density per year across trophic levels and create environmental PCA
state.kin.dat <- data.frame("lake" = "Kinneret",
                            "tot_density" =scale(rowSums(scale(log1p(plank_env.data.yr[5:46,c(2,4:77)])))),
                            "community" = prcomp(scale(plank_env.data.yr[5:46,c(2,4:77)]))$x[,1],
                            "env" = prcomp(scale(plank_env.data.yr[5:46,c("Temp","Nitrate","P_ort")]))$x[,1],
                            "env2" = prcomp(scale(plank_env.data.yr[5:46,c("Temp","Nitrate","P_ort")]))$x[,2],
                            "env2" = prcomp(scale(plank_env.data.yr[5:46,c("Temp","Nitrate","P_ort")]))$x[,2],
                            "phyto_density" = scale(log1p(rowSums(plank_env.data.yr[5:46,c(2,4:46)]))),
                            "zoo_density" =scale(log1p(rowSums(plank_env.data.yr[5:46,c(47:77)]))),
                            "date" = as.numeric(plank_env.data.yr$Date[5:46])) #unable to use first 4 years due to missing P_ort data

state.wind.dat <- data.frame("lake" = "Windermere",
                             "tot_density" =scale(rowSums(scale(log1p(phyto_env.windyrdata[,2:22])))),
                             "community" = prcomp(scale(phyto_env.windyrdata[,2:22]))$x[,1],
                             "env" = prcomp(scale(phyto_env.windyrdata[,c("TEMP","NO3N","TOTP")]))$x[,1],
                             "env2" = prcomp(scale(phyto_env.windyrdata[,c("TEMP","NO3N","TOTP")]))$x[,2],
                             "phyto_density" = scale(log1p(rowSums(phyto_env.windyrdata[,2:18]))),
                             "zoo_density" =scale(log1p(rowSums(phyto_env.windyrdata[,20:23]))),
                             "date" = as.numeric(phyto_env.windyrdata$Date)) 

############ 
#Fit temporal TGAMs
############ 
lake_temporal <- foreach::foreach(j = list(state.kas.dat,state.kin.dat,state.wind.dat), 
                                  .combine = "rbind") %do% {
                                    
                                    lapply(c("tot_density","phyto_density","zoo_density","community"), 
                                           FUN = function(i){
                                             
                                             cont_formula <-  as.formula(paste(i,"~ s(date, bs='tp', k=3)"))
                                             thresh_formula <- as.formula(paste(i,"~ 1"))
                                             dens_gam <- compare_gam(data = j,
                                                                     cont_formula = formula(cont_formula, method = "REML"),
                                                                     thresh_formula =  formula(thresh_formula, method = "REML"),
                                                                     thresh.var = "date",
                                                                     thresh.range = c(0.15,0.85), by = 1, k=3) #quartile range for threshold within 15-85%. Prevents bias towards thresholds at ends of time series
                                             
                                             best_gam <- predict_best_gam(object=dens_gam) %>%
                                               dplyr::rename("metric.val" = i) 
                                             
                                             return(best_gam)
                                             
                                           }) |> `names<-`(c("tot_density","phyto_density","zoo_density","community")) |>
                                      data.table::rbindlist(idcol = "metric",use.names = T) |>
                                      dplyr::mutate(threshold = tidyr::replace_na(as.character(threshold),"pre"),
                                                    start.date = dplyr::first(thresh.var),
                                                    last.date = dplyr::last(thresh.var),
                                                    lake = j$lake[1]
                                      )
                                  }

ss_time_plot_dat1 <- subset(lake_temporal, metric %in% c("phyto_density","zoo_density")) |>
  mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton"))) #only interested in phytoplankton and zooplankton densities

############ 
#Fit environmental TGAMs
############ 
lake_state_space <- foreach::foreach(j = list(state.kas.dat,state.wind.dat), 
                                     .combine = "rbind") %do% {
                                       
                                       lapply(c("tot_density","phyto_density","zoo_density","community"), 
                                              FUN = function(i){
                                                
                                                cont_formula <-  as.formula(paste(i,"~ s(env, bs='tp', k=3)"))
                                                thresh_formula <- as.formula(paste(i,"~ 1"))
                                                dens_gam <- compare_gam(data = j,
                                                                        cont_formula = formula(cont_formula, method = "REML"),
                                                                        thresh_formula =  formula(thresh_formula, method = "REML"),
                                                                        thresh.var = "date", expl.var = "env",
                                                                        thresh.range = c(0.15,0.85), by = 1, k=3)
                                                
                                                best_gam <- predict_best_gam(object=dens_gam) %>%
                                                  dplyr::rename("metric.val" = i) 
                                                
                                                return(best_gam)
                                                
                                              }) |> `names<-`(c("tot_density","phyto_density","zoo_density","community")) |>
                                         data.table::rbindlist(idcol = "metric",use.names = T) |>
                                         dplyr::mutate(threshold = tidyr::replace_na(as.character(threshold),"pre"),
                                                       start.date = dplyr::first(thresh.var),
                                                       last.date = dplyr::last(thresh.var),
                                                       lake = j$lake[1]
                                         )
                                     }
ss_plot_dat1 <- subset(lake_state_space, metric %in% c("phyto_density","zoo_density")) |>
  mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton"))) #only interested in phytoplankton and zooplankton densities

############ 
#Create figure
############ 
fig1_plot_dat <- rbind(ss_plot_dat1 |>
                         pivot_longer(env,names_to = "var_metric",values_to = "var.value") |>
                         mutate(var.value = var.value/2), #shrink for plotting purposes
                       ss_time_plot_dat1|>
                         pivot_longer(date,names_to = "var_metric",values_to = "var.value")) |>
  mutate(var_metric = factor(ifelse(var_metric == "env","Environment","Year"),levels = c("Year","Environment")))


fig1 <- ggplot(data = fig1_plot_dat, aes(x=var.value,y=metric.val)) + 
  geom_point(aes(x=var.value, y = metric.val))+
  geom_path(aes(x=var.value, y = metric.val)) +  
  xlab("Expanatory variable") + ylab("Scaled metric score")+ 
  ggh4x::facet_nested(lake + metric ~var_metric,scales = "free_x",
                      labeller = label_value,
                      switch = "x",
                      strip = ggh4x::strip_nested(size="constant",bleed=F, 
                                                  background_y = ggh4x::elem_list_rect(fill = c("#D6D6D6","#D6D6D6",rep("white",4))),
                                                  background_x = ggh4x::elem_list_rect(fill = c(rep("white",2))),
                                                  by_layer_y = F)) +
  geom_line(data = filter(fig1_plot_dat, threshold=="pre" & var_metric == "Environment"),aes(x=var.value,y=fit), col="blue",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="pre"& var_metric == "Environment"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "Environment"),aes(x=var.value,y=fit), col="red",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "Environment"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  geom_point(data = filter(fig1_plot_dat, thresh.var == start.date & var_metric == "Environment"),aes(x=var.value, y = metric.val,col="Start date"))+
  geom_point(data = filter(fig1_plot_dat, thresh.var == last.date & var_metric == "Environment"),aes(x=var.value, y = metric.val,col="End date"))+
  geom_point(data = filter(fig1_plot_dat,transition == "trans" & var_metric == "Environment"),aes(x=var.value, y = metric.val,col="Transition\ndates"))+
  ggrepel::geom_text_repel(data = filter(fig1_plot_dat,transition == "trans" & metric == "Phytoplankton" & var_metric == "Environment"), aes(x=var.value, y = metric.val,label=thresh.var),force =1,nudge_x=1.5,nudge_y=-0.25,segment.linetype=2,min.segment.length = 0.1)+
  ggrepel::geom_text_repel(data = filter(fig1_plot_dat,transition == "trans" & metric == "Zooplankton" & var_metric == "Environment"), aes(x=var.value, y = metric.val,label=thresh.var),force =1,nudge_x=-1.5,nudge_y=0.5,segment.linetype=2,min.segment.length = 0.1)+
  geom_line(data = filter(fig1_plot_dat, threshold=="pre" & var_metric == "Year"),aes(x=var.value,y=fit), col="blue",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="pre"& var_metric == "Year"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "Year"),aes(x=var.value,y=fit), col="red",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "Year"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  geom_point(data = filter(fig1_plot_dat, thresh.var == start.date & var_metric == "Year"),aes(x=var.value, y = metric.val,col="Start date"))+
  geom_point(data = filter(fig1_plot_dat, thresh.var == last.date & var_metric == "Year"),aes(x=var.value, y = metric.val,col="End date"))+
  geom_point(data = filter(fig1_plot_dat,transition == "trans" & var_metric == "Year"),aes(x=var.value, y = metric.val,col="Transition\ndates"))+
  ggrepel::geom_text_repel(data = filter(fig1_plot_dat,transition == "trans" & metric == "Phytoplankton" & var_metric == "Year"), aes(x=var.value, y = metric.val,label=thresh.var),force =1,nudge_x=1.5,nudge_y=-0.75,segment.linetype=2,min.segment.length = 0.1)+
  ggrepel::geom_text_repel(data = filter(fig1_plot_dat,transition == "trans" & metric == "Zooplankton" & var_metric == "Year"), aes(x=var.value, y = metric.val,label=thresh.var),force =1.5,nudge_x=-7.5,nudge_y=0.25,segment.linetype=2,min.segment.length = 0.1)+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates"), name = NULL, guide = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  theme(legend.text=element_text(size=12),
        strip.background.x = element_rect(linewidth = 0),
        strip.placement.y = "outside",
        legend.position = "top",
        panel.border = element_rect(linewidth = 1,fill = "transparent"))

ggplot2::ggsave("Figures/figure2.png",
                fig1,
                width = 6,height=6.5,dpi=200)

