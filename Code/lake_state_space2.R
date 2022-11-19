#############################################################################
# Preamble
#############################################################################
require(tidyverse)

source("/Users/ul20791/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_environmental_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Zurich/Data/zurich_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Zurich/Data/LZ_environmental_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Madison/Data/madison_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Madison/Data/madison_environmental_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Windermere/Data/windermere_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Windermere/Data/windermere_environmental_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Kasumigaura/Data/kasumigaura_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Kasumigaura/Data/Kasumigaura_environmental_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Zurich/Data/UZ_plankton_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Madison/Data/monona_plankton_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Leven/Data/leven_plankton_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Washington/Data/washington_plankton_data.R")

rm(list = ls()[!ls() %in% c("plank_env.data.mth","plank_env.data.yr","plank_env.kasmthdata","plank_env.kasyrdata",
              "plank_env.levemthdata","plank_env.leveyrdata","plank_env.LZmthdata","plank_env.LZyrdata",
              "plank_env.madmthdata","plank_env.madyrdata","plank_env.monmthdata","plank_env.monyrdata",
              "plank_env.UZmthdata","plank_env.UZyrdata","plank_env.washmthdata","plank_env.washyrdata",
              "phyto_env.windmthdata","phyto_env.windyrdata")]
)

source("/Users/ul20791/Desktop/Academia/PhD/Repositories/ews-assessments/ews-assessments/Code/threshold_gam.R")


#############################################################################
#Prepare Lake Data
#############################################################################

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

state.kas.dat <- data.frame("lake" = "Kasumigaura",
                            "tot_density" =scale(rowSums(scale(log1p(plank_env.kasyrdata[1:38,2:154])))),
                            "community" = prcomp(scale(plank_env.kasyrdata[1:38,2:154]))$x[,1],
                            "env" = prcomp(scale(plank_env.kasyrdata[1:38,c("wtemp","NO3N","totP")]))$x[,1],
                            "env2" = prcomp(scale(plank_env.kasyrdata[1:38,c("wtemp","NO3N","totP")]))$x[,2],
                            "phyto_density" = scale(log1p(rowSums(plank_env.kasyrdata[1:38,2:119]))),
                            "zoo_density" =scale(log1p(rowSums(plank_env.kasyrdata[1:38,120:154]))),
                            "date" = as.numeric(plank_env.kasyrdata$date[1:38])) #trim last outlier year

state.leve.dat <- data.frame("lake" = "Loch Leven",
                             "tot_density" =scale(rowSums(scale(log1p(plank_env.leveyrdata[,2:10])))),
                             "community" = prcomp(scale(plank_env.leveyrdata[,2:10]))$x[,1],
                             "env" = prcomp(scale(plank_env.leveyrdata[,c("Temp","TOTP")]))$x[,1],
                             "env2" = prcomp(scale(plank_env.leveyrdata[,c("Temp","TOTP")]))$x[,2],
                             "phyto_density" = scale(log1p(rowSums(plank_env.leveyrdata[,2:5]))),
                             "zoo_density" =scale(log1p(rowSums(plank_env.leveyrdata[,6:10]))),
                             "date" = as.numeric(plank_env.leveyrdata$Date)) 

state.LZ.dat <- data.frame("lake" = "Lower Zurich",
                           "tot_density" =scale(rowSums(scale(log1p(plank_env.LZyrdata[,5:206])))),
                            "community" = prcomp(scale(plank_env.LZyrdata[,5:206]))$x[,1],
                            "env" = prcomp(scale(plank_env.LZyrdata[,c("mean.t","mean.po4","NO3_N")]))$x[,1],
                           "env2" = prcomp(scale(plank_env.LZyrdata[,c("mean.t","mean.po4","NO3_N")]))$x[,2],
                           "phyto_density" = scale(log1p(rowSums(plank_env.LZyrdata[,5:179]))),
                            "zoo_density" =scale(log1p(rowSums(plank_env.LZyrdata[,180:206]))),
                            "date" = as.numeric(plank_env.LZyrdata$date)) 

state.UZ.dat <- data.frame("lake" = "Upper Zurich",
                           "tot_density" =scale(rowSums(scale(log1p(plank_env.UZyrdata[,2:95])))),
                           "community" = prcomp(scale(plank_env.UZyrdata[,2:95]))$x[,1],
                           "env" = prcomp(scale(plank_env.UZyrdata[,c("mean.t","mean.po4")]))$x[,1],
                           "env2" = prcomp(scale(plank_env.UZyrdata[,c("mean.t","mean.po4")]))$x[,2],
                           "phyto_density" = scale(log1p(rowSums(plank_env.UZyrdata[,2:69]))),
                           "zoo_density" =scale(log1p(rowSums(plank_env.UZyrdata[,70:95]))),
                           "date" = as.numeric(plank_env.UZyrdata$date)) 

state.mad.dat <- data.frame("lake" = "Mendota",
                            "tot_density" =scale(rowSums(scale(log1p(plank_env.madyrdata[,2:222])))),
                           "community" = prcomp(scale(plank_env.madyrdata[,2:222]))$x[,1],
                           "env" = prcomp(scale(plank_env.madyrdata[,c("wtemp","totP","NO3N")]))$x[,1],
                           "env2" = prcomp(scale(plank_env.madyrdata[,c("wtemp","totP","NO3N")]))$x[,2],
                           "phyto_density" = scale(log1p(rowSums(plank_env.madyrdata[,2:198]))),
                           "zoo_density" =scale(log1p(rowSums(plank_env.madyrdata[,199:222]))),
                           "date" = as.numeric(plank_env.madyrdata$date)) 

state.mon.dat <- data.frame("lake" = "Monona",
                            "tot_density" =scale(rowSums(scale(log1p(plank_env.monyrdata[1:19,2:202])))),
                            "community" = prcomp(scale(plank_env.monyrdata[1:19,2:202]))$x[,1],
                            "env" = prcomp(scale(plank_env.monyrdata[1:19,c("wtemp","totP","NO3N")]))$x[,1],
                            "env2" = prcomp(scale(plank_env.monyrdata[1:19,c("wtemp","totP","NO3N")]))$x[,2],
                            "phyto_density" = scale(log1p(rowSums(plank_env.monyrdata[1:19,2:180]))),
                            "zoo_density" =scale(log1p(rowSums(plank_env.monyrdata[1:19,181:202]))),
                            "date" = as.numeric(plank_env.monyrdata$date[1:19]))  #trim anomalous last data point

state.wash.dat <- data.frame("lake" = "Washington",
                             "tot_density" =scale(rowSums(scale(log1p(plank_env.washyrdata[,4:18])))),
                            "community" = prcomp(scale(plank_env.washyrdata[,4:18]))$x[,1],
                            "env" = prcomp(scale(plank_env.washyrdata[,c("temp","TP")]))$x[,1],
                            "env2" = prcomp(scale(plank_env.washyrdata[,c("temp","TP")]))$x[,1],
                            "phyto_density" = scale(log1p(rowSums(plank_env.washyrdata[,4:9]))),
                            "zoo_density" =scale(log1p(rowSums(plank_env.washyrdata[,10:18]))),
                            "date" = as.numeric(plank_env.washyrdata$date)) 

#############################################################################
#Perform State Space Fitting - Environmentals
#############################################################################
lake_state_space <- foreach::foreach(j = list(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                                                state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat), 
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

ss_plot_dat <- subset(lake_state_space, metric %in% c("phyto_density","zoo_density"))

p_env <- ggplot(data = ss_plot_dat, aes(x=env,y=metric.val)) + 
  geom_point(aes(x=env, y = metric.val))+
  geom_path(aes(x=env, y = metric.val)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.kin, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + 
  facet_grid(metric~lake,scales = "free_y") +
  geom_line(data = filter(ss_plot_dat, threshold=="pre"),aes(x=env,y=fit), col="blue",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_plot_dat, threshold=="pre"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(ss_plot_dat, threshold=="post"),aes(x=env,y=fit), col="red",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_plot_dat, threshold=="post"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  geom_point(data = filter(ss_plot_dat, thresh.var == start.date),aes(x=env, y = metric.val,col="Start date"))+
  geom_point(data = filter(ss_plot_dat, thresh.var == last.date),aes(x=env, y = metric.val,col="End date"))+
  ggrepel::geom_text_repel(data = filter(ss_plot_dat,transition == "trans"), aes(x=env, y = metric.val,label=thresh.var),force =1,nudge_x=3,nudge_y=-0.5,segment.linetype=2)+
  geom_point(data = filter(ss_plot_dat,transition == "trans"),aes(x=env, y = metric.val,col="Transition\ndates"))+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates","1993"), name = NULL)+  
  xlim(-6,6)


lake_state_space2 <- foreach::foreach(j = list(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                                              state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat), 
                                     .combine = "rbind") %do% {
                                       
                                       lapply(c("tot_density","phyto_density","zoo_density","community"), 
                                              FUN = function(i){
                                                
                                                cont_formula <-  as.formula(paste(i,"~ s(env2, bs='tp', k=10)"))
                                                thresh_formula <- as.formula(paste(i,"~ 1"))
                                                dens_gam <- compare_gam(data = j,
                                                                        cont_formula = formula(cont_formula, method = "REML"),
                                                                        thresh_formula =  formula(thresh_formula, method = "REML"),
                                                                        thresh.var = "date", 
                                                                        thresh.range = c(0.15,0.85),
                                                                        expl.var = "env2")
                                                
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

ss_plot_dat2 <- subset(lake_state_space2, metric %in% c("phyto_density","zoo_density"))

ggplot(data = ss_plot_dat2, aes(x=env2,y=metric.val)) + 
  geom_point(aes(x=env2, y = metric.val))+
  geom_path(aes(x=env2, y = metric.val)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.kin, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + 
  facet_grid(metric~lake,scales = "free_y") +
  geom_line(data = filter(ss_plot_dat2, threshold=="pre"),aes(x=env2,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_plot_dat2, threshold=="pre"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(ss_plot_dat2, threshold=="post"),aes(x=env2,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_plot_dat2, threshold=="post"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_point(data = filter(ss_plot_dat2, thresh.var == start.date),aes(x=env2, y = metric.val,col="Start date"))+
  geom_point(data = filter(ss_plot_dat2, thresh.var == last.date),aes(x=env2, y = metric.val,col="End date"))+
  ggrepel::geom_text_repel(data = filter(ss_plot_dat2,transition == "trans"), aes(x=env2, y = metric.val,label=thresh.var),force =1,nudge_x=3,nudge_y=-0.5,segment.linetype=2)+
  geom_point(data = filter(ss_plot_dat2,transition == "trans"),aes(x=env2, y = metric.val,col="Transition\ndates"))+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates","1993"), name = NULL)+  
  xlim(-6,6)


tot_plankton_ss <- rbind(lake_state_space |>  dplyr::mutate(driver = "PC1"),
                         lake_state_space2 |>  dplyr::rename(env = env2) |> dplyr::mutate(driver = "PC2")) |>
  dplyr::filter(metric == "tot_density")


ggplot(data = tot_plankton_ss, aes(x=env,y=metric.val)) + 
  geom_point(aes(x=env, y = metric.val))+
  geom_path(aes(x=env, y = metric.val)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.kin, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + 
  facet_grid(driver~lake,scales = "free_y") +
  geom_line(data = filter(tot_plankton_ss, threshold=="pre"),aes(x=env,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(tot_plankton_ss, threshold=="pre"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(tot_plankton_ss, threshold=="post"),aes(x=env,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(tot_plankton_ss, threshold=="post"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_point(data = filter(tot_plankton_ss, thresh.var == start.date),aes(x=env, y = metric.val,col="Start date"))+
  geom_point(data = filter(tot_plankton_ss, thresh.var == last.date),aes(x=env, y = metric.val,col="End date"))+
  ggrepel::geom_text_repel(data = filter(tot_plankton_ss,transition == "trans"), aes(x=env, y = metric.val,label=thresh.var),force =1,nudge_x=3,nudge_y=-0.5,segment.linetype=2)+
  geom_point(data = filter(tot_plankton_ss,transition == "trans"),aes(x=env, y = metric.val,col="Transition\ndates"))+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates","1993"), name = NULL)+  
  xlim(-6,6)

#############################################################################
#Perform State Space Fitting - Temporal
#############################################################################

lake_temporal <- foreach::foreach(j = list(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                                              state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat), 
                                     .combine = "rbind") %do% {
                                       
                                       lapply(c("tot_density","phyto_density","zoo_density","community"), 
                                              FUN = function(i){
                                                
                                                cont_formula <-  as.formula(paste(i,"~ s(date, bs='tp', k=3)"))
                                                thresh_formula <- as.formula(paste(i,"~ 1"))
                                                dens_gam <- compare_gam(data = j,
                                                                        cont_formula = formula(cont_formula, method = "REML"),
                                                                        thresh_formula =  formula(thresh_formula, method = "REML"),
                                                                        thresh.var = "date",
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

ss_time_plot_dat <- subset(lake_temporal, metric %in% c("phyto_density","zoo_density"))


p_temporal <- ggplot(data = ss_time_plot_dat, aes(x=date,y=metric.val)) + 
  geom_point(aes(x=date, y = metric.val))+
  geom_path(aes(x=date, y = metric.val)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.kin, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Year") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + 
  facet_grid(metric~lake,scales = "free_y") +
  geom_line(data = filter(ss_time_plot_dat, threshold=="pre"),aes(x=date,y=fit), col="blue",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_time_plot_dat, threshold=="pre"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(ss_time_plot_dat, threshold=="post"),aes(x=date,y=fit), col="red",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_time_plot_dat, threshold=="post"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(ss_time_plot_dat,transition == "trans"), aes(x=date, y = metric.val,label=thresh.var),force =1,nudge_x=3,nudge_y=-0.5,segment.linetype=2)+
  geom_point(data = filter(ss_time_plot_dat,date == start.date),aes(x=date, y = metric.val,col="Start date"))+
  geom_point(data = filter(ss_time_plot_dat,date == last.date),aes(x=date, y = metric.val,col="End date"))+
  geom_point(data = filter(ss_time_plot_dat,transition == "trans"),aes(x=date, y = metric.val,col="Transition\ndates"))+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates","1993"), name = NULL) +
  scale_x_continuous(breaks=seq(1970,2015,by=15),limits = c(1960,2020))

#############################################################################
#Compare State Spaces
#############################################################################
state_transition_dates <- subset(lake_state_space,transition == "trans") |>
  rename(state_date = date)
temporal_transition_dates <- subset(lake_temporal,transition == "trans") |>
  rename(temporal_date = date)

transition_dates <- left_join(state_transition_dates,temporal_transition_dates,by = lake) |>
  mutate(date_match = ifelse(state_date == temporal_date, TRUE, FALSE))

ggplot2::ggsave("/Users/ul20791/Downloads/lake_state_spaces.pdf",
       ggpubr::ggarrange(p_temporal + ggtitle("a) time series"),
                         p_env +  ggtitle("b) state space"),
                         common.legend = T,ncol=1,nrow=2,legend = "right")
       ,width = 18,height=10,dpi=144)

