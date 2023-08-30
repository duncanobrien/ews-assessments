#############################################################################
# Preamble
#############################################################################

require(tidyverse)
require(foreach)
require(patchwork)

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

source("Code/threshold_gam.R")

#############################################################################
#Prepare Lake Data
#############################################################################

state.kin.dat <- data.frame("lake" = "Kinneret",
                            "tot_density" =scale(rowSums(scale(log1p(plank_env.data.yr[5:46,c(2,4:77)])))),
                            "community" = prcomp(scale(plank_env.data.yr[5:46,c(2,4:77)]))$x[,1],
                            "env" = prcomp(scale(plank_env.data.yr[5:46,c("Temp","Nitrate","P_ort")]))$x[,1],
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
#Identify bi-modality
#############################################################################
require(LaplacesDemon)
require(diptest)
require(mousetrap)

lake_bimod <- foreach::foreach(j = list(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                                        state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat), 
                               .combine = "rbind") %do% {
                                 
                                 lapply(c("tot_density","phyto_density","zoo_density","community"),function(i){
                                   
                                   modes <- LaplacesDemon::Modes(j[,i])$modes
                                   
                                   return(data.frame("metric" = i,
                                                     "is_bimodal" =  LaplacesDemon::is.bimodal(j[,i]),
                                                     "mode1" = modes[1],
                                                     "mode2" = ifelse(length(modes)==1,NA,modes[2]),
                                                     "modality_coef" = mousetrap::bimodality_coefficient(j[,i]))) #Pfister et al., 2013
                                 }) |>
                                   data.table::rbindlist() |>
                                   dplyr::mutate(lake = j$lake[1])
                               }

bimod_plot_datS1 <- rbind(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                          state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat) |>
  dplyr::select(lake,phyto_density,zoo_density) |>
  tidyr::pivot_longer(-lake,names_to = "metric",values_to = "metric.value") |>
  dplyr::mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton")))


bimod_modes <- subset(lake_bimod, metric %in% c("phyto_density","zoo_density")) |>
  dplyr::mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton"))) |>
  dplyr::group_by(metric) |>
  dplyr::mutate(y = ifelse(metric == "Phytoplankton",c(0.525,0.5),c(0.725,0.7)),
                mode1 = round(mode1,digits = 2),
                mode2 = round(mode2,digits = 2),
                mode_sig = ifelse(modality_coef > 0.5, paste0(modes,"*"),paste0(modes)))
 
#############################################################################
#TGAM Figure 1
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

ss_time_plot_dat1 <- subset(lake_temporal, metric %in% c("phyto_density","zoo_density") & 
                              lake %in% c("Kinneret","Windermere")) |>
  mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton")))

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

ss_plot_dat1 <- subset(lake_state_space, metric %in% c("phyto_density","zoo_density") & 
                         lake %in% c("Kinneret","Windermere")) |>
  mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton")))

lake_bimod <- foreach::foreach(j = list(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                                        state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat), 
                               .combine = "rbind") %do% {
                                 
                                 lapply(c("tot_density","phyto_density","zoo_density","community"),function(i){
                                   
                                   modes <- LaplacesDemon::Modes(j[,i])$modes
                                   
                                   return(data.frame("metric" = i,
                                                     "is_bimodal" =  LaplacesDemon::is.bimodal(j[,i]),
                                                     "mode1" = modes[1],
                                                     "mode2" = ifelse(length(modes)==1,NA,modes[2]),
                                                     "modality_coef" = mousetrap::bimodality_coefficient(j[,i])))
                                   
                                 }) |>
                                   data.table::rbindlist() |>
                                   dplyr::mutate(lake = j$lake[1])
                               }

bimod_plot_datS1 <- rbind(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                          state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat) |>
  dplyr::select(lake,phyto_density,zoo_density) |>
  tidyr::pivot_longer(-lake,names_to = "metric",values_to = "metric.value") |>
  dplyr::mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton")))

bimod_modes <- subset(lake_bimod, metric %in% c("phyto_density","zoo_density")) |>
  dplyr::mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton"))) |>
  dplyr::group_by(metric) |>
  dplyr::mutate(y = ifelse(metric == "Phytoplankton",c(0.525,0.5),c(0.725,0.7)),
                mode1 = round(mode1,digits = 2),
                mode2 = round(mode2,digits = 2)) 

fig1_plot_dat <- rbind(ss_plot_dat1 |>
                         pivot_longer(env,names_to = "var_metric",values_to = "var.value") |>
                         mutate(var.value = var.value/2), #shrink for plotting purposes
                       ss_time_plot_dat1|>
                         pivot_longer(date,names_to = "var_metric",values_to = "var.value") ) |>
  mutate(var_metric = factor(ifelse(var_metric == "env","Environment",
                                    ifelse(var_metric == "date","Year","Bimodality")),
                             levels = c("Year","Environment","Bimodality")))

fig1_plot_dat_bimod <- subset(bimod_plot_datS1, 
                              lake %in% c("Kinneret","Windermere"))|>
  mutate(var_metric = "bimodal") |>
  mutate(var_metric = factor(ifelse(var_metric == "env","Environment",
                                    ifelse(var_metric == "date","Year","Bimodality")),
                             levels = c("Year","Environment","Bimodality"))) |>
  right_join(subset(bimod_modes,
                    lake %in% c("Kinneret","Windermere")),by = c("lake","metric"),
             relationship = "many-to-many") |>
  pivot_longer(mode1:mode2,values_to = "modes")|>
  dplyr::mutate(mode_sig = ifelse(modality_coef > 0.5, paste0(modes,"*"),paste0(modes)))
  


fig1 <- ggplot(data = fig1_plot_dat, aes(y=metric.val)) + 
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
  geom_line(data = filter(fig1_plot_dat, threshold=="pre" & var_metric == "Environment"),aes(x=var.value,y=fit), col="blue",linewidth=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="pre"& var_metric == "Environment"),aes(x=var.value,ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "Environment"),aes(x=var.value,y=fit), col="red",linewidth=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "Environment"),aes(x=var.value,ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  geom_point(data = filter(fig1_plot_dat, thresh.var == start.date & var_metric == "Environment"),aes(x=var.value, y = metric.val,col="Start date"))+
  geom_point(data = filter(fig1_plot_dat, thresh.var == last.date & var_metric == "Environment"),aes(x=var.value, y = metric.val,col="End date"))+
  geom_point(data = filter(fig1_plot_dat,transition == "trans" & var_metric == "Environment"),aes(x=var.value, y = metric.val,col="Transition\ndates"))+
  ggrepel::geom_text_repel(data = filter(fig1_plot_dat,transition == "trans" & metric == "Phytoplankton" & var_metric == "Environment"), aes(x=var.value, y = metric.val,label=thresh.var),force =1,nudge_x=1.5,nudge_y=-0.25,segment.linetype=2,min.segment.length = 0.1)+
  ggrepel::geom_text_repel(data = filter(fig1_plot_dat,transition == "trans" & metric == "Zooplankton" & var_metric == "Environment"), aes(x=var.value, y = metric.val,label=thresh.var),force =1,nudge_x=-1.5,nudge_y=0.5,segment.linetype=2,min.segment.length = 0.1)+
  geom_line(data = filter(fig1_plot_dat, threshold=="pre" & var_metric == "Year"),aes(x=var.value,y=fit), col="blue",linewidth=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="pre"& var_metric == "Year"),aes(x=var.value,ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "Year"),aes(x=var.value,y=fit), col="red",linewidth=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "Year"),aes(x=var.value,ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  geom_point(data = filter(fig1_plot_dat, thresh.var == start.date & var_metric == "Year"),aes(x=var.value, y = metric.val,col="Start date"))+
  geom_point(data = filter(fig1_plot_dat, thresh.var == last.date & var_metric == "Year"),aes(x=var.value, y = metric.val,col="End date"))+
  geom_point(data = filter(fig1_plot_dat,transition == "trans" & var_metric == "Year"),aes(x=var.value, y = metric.val,col="Transition\ndates"))+
  ggrepel::geom_text_repel(data = filter(fig1_plot_dat,transition == "trans" & metric == "Phytoplankton" & var_metric == "Year"), aes(x=var.value, y = metric.val,label=thresh.var),force =1,nudge_x=1.5,nudge_y=-0.75,segment.linetype=2,min.segment.length = 0.1)+
  ggrepel::geom_text_repel(data = filter(fig1_plot_dat,transition == "trans" & metric == "Zooplankton" & var_metric == "Year"), aes(x=var.value, y = metric.val,label=thresh.var),force =1.5,nudge_x=-7.5,nudge_y=0.25,segment.linetype=2,min.segment.length = 0.1)+
  geom_density(data = fig1_plot_dat_bimod,aes(y = metric.value)) +
  geom_hline(data = fig1_plot_dat_bimod,aes(yintercept = modes),linetype = "dashed",colour = "black")+
  geom_text(data = fig1_plot_dat_bimod,aes(y=modes+0.5,x=0.05,label = mode_sig),size = 5)+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates"), name = NULL, guide = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  theme(legend.text=element_text(size=12),
        strip.background.x = element_rect(linewidth = 0),
        strip.placement.y = "outside",
        legend.position = "top",
        panel.border = element_rect(linewidth = 1,fill = "transparent"))

ggplot2::ggsave("/Users/ul20791/Desktop/Academia/Papers/OBrien_et_al_PlankEWS/lake_state_spaces_fig.png",
                fig1,
                width = 6,height=6.5,dpi=200)


p_temporal1 + 
  p_env1 +  
  patchwork::plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A')

#############################################################################
#TGAM Figure Supp 1
#############################################################################

ss_time_plot_datS1 <- subset(lake_temporal, metric %in% c("phyto_density","zoo_density")) |>
  mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton")))

ss_plot_datS1 <- subset(lake_state_space, metric %in% c("phyto_density","zoo_density")) |>
  mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton")))

p_temporalS1 <- ggplot(data = ss_time_plot_datS1, aes(x=date,y=metric.val)) + 
  geom_point(aes(x=date, y = metric.val))+
  geom_path(aes(x=date, y = metric.val)) +  
  xlab("Year") + ylab("Scaled metric score")+ 
  theme_bw() + 
  facet_grid(metric~lake,scales = "free_y") +
  geom_line(data = filter(ss_time_plot_datS1, threshold=="pre"),aes(x=date,y=fit), col="blue",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_time_plot_datS1, threshold=="pre"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(ss_time_plot_datS1, threshold=="post"),aes(x=date,y=fit), col="red",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_time_plot_datS1, threshold=="post"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  geom_point(data = filter(ss_time_plot_datS1,date == start.date),aes(x=date, y = metric.val,col="Start date"))+
  geom_point(data = filter(ss_time_plot_datS1,date == last.date),aes(x=date, y = metric.val,col="End date"))+
  geom_point(data = filter(ss_time_plot_datS1,transition == "trans"),aes(x=date, y = metric.val,col="Transition\ndates"))+
  ggrepel::geom_text_repel(data = filter(ss_time_plot_datS1,transition == "trans" & !(lake == "Windermere" & metric == "Zooplankton")), aes(x=date, y = metric.val,label=thresh.var),force =1,nudge_x=3,nudge_y=-0.5,segment.linetype=2)+
  ggrepel::geom_text_repel(data = filter(ss_time_plot_datS1,transition == "trans"& lake == "Windermere" & metric == "Zooplankton"), aes(x=date, y = metric.val,label=thresh.var),force =1,nudge_x=-1.5,nudge_y=-0.25,segment.linetype=2)+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates","1993"), name = NULL, guide = guide_legend(override.aes = list(size = 5)))+  
  scale_x_continuous(breaks=seq(1970,2015,by=15),limits = c(1960,2020))+
  theme_classic() +
  theme(legend.text=element_text(size=12),
        strip.background = element_rect(fill = "#D6D6D6"),
        strip.placement.y = "outside",
        legend.position = "top",
        panel.border = element_rect(linewidth = 1,fill = "transparent"))

p_envS1 <- ggplot(data = ss_plot_datS1, aes(x=env,y=metric.val)) + 
  geom_point(aes(x=env, y = metric.val))+
  geom_path(aes(x=env, y = metric.val)) +  
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  theme_bw() + 
  facet_grid(metric~lake,scales = "free_y") +
  geom_line(data = filter(ss_plot_datS1, threshold=="pre"),aes(x=env,y=fit), col="blue",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_plot_datS1, threshold=="pre"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(ss_plot_datS1, threshold=="post"),aes(x=env,y=fit), col="red",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_plot_datS1, threshold=="post"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  geom_point(data = filter(ss_plot_datS1, thresh.var == start.date),aes(x=env, y = metric.val,col="Start date"))+
  geom_point(data = filter(ss_plot_datS1, thresh.var == last.date),aes(x=env, y = metric.val,col="End date"))+
  geom_point(data = filter(ss_plot_datS1,transition == "trans"),aes(x=env, y = metric.val,col="Transition\ndates"))+
  ggrepel::geom_text_repel(data = filter(ss_plot_datS1,transition == "trans" & metric == "Phytoplankton"), aes(x=env, y = metric.val,label=thresh.var),force =1,nudge_x=1.5,nudge_y=-0.25,segment.linetype=2)+
  ggrepel::geom_text_repel(data = filter(ss_plot_datS1,transition == "trans" & metric == "Zooplankton"), aes(x=env, y = metric.val,label=thresh.var),force =1,nudge_x=-1.5,segment.linetype=2)+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates","1993"), name = NULL, guide = guide_legend(override.aes = list(size = 5)))+  
  xlim(-6,6) +
  theme_classic() +
  theme(legend.text=element_text(size=12),
        strip.background = element_rect(fill = "#D6D6D6"),
        strip.placement.y = "outside",
        legend.position = "top",
        panel.border = element_rect(linewidth = 1,fill = "transparent"))

p_bimodal <- ggplot(right_join(bimod_plot_datS1,bimod_modes,relationship = "many-to-many") |>
                      group_by(metric,lake)) +
  geom_density(aes(y = metric.value)) +
  #geom_histogram(aes(y = metric.value),binwidth = 0.5) +
  geom_hline(aes(yintercept = modes),linetype = "dashed",colour = "black")+
  geom_text(aes(y=modes+0.25,x=0.2,label = mode_sig),size = 5)+
  facet_grid(metric~lake,scales = "free_x") +
  xlab("Density") + ylab("Scaled metric score")+ 
  theme_classic() +
  theme(legend.text=element_text(size=12),
        strip.background = element_rect(fill = "#D6D6D6"),
        strip.placement.y = "outside",
        legend.position = "top",
        panel.border = element_rect(linewidth = 1,fill = "transparent"))


ggplot2::ggsave("/Users/ul20791/Desktop/Academia/Papers/OBrien_et_al_PlankEWS/lake_state_spaces_supp_fig.png",
                ggpubr::ggarrange(p_temporalS1 + ggtitle("Time series"),
                                  p_envS1 +  ggtitle("State space"),
                                  p_bimodal +  ggtitle("Bimodality"),
                                  labels="AUTO",font.label = list(face="plain"),
                                  common.legend = T,ncol=1,nrow=3,legend = "top")
                ,width = 16,height=9,dpi=300)

#############################################################################
#Compare State Spaces
#############################################################################

state_transition_dates <- subset(lake_state_space,transition == "trans") |>
  rename(state_date = thresh.var)
temporal_transition_dates <- subset(lake_temporal,transition == "trans") |>
  rename(temporal_date = date)
bimod_dates <- lake_bimod |>
  dplyr::mutate(is_bimodal = ifelse(modality_coef > 0.5 & all(!is.na(c(mode1,mode2))),TRUE,FALSE))
  
transition_dates <-expand.grid(metric = unique(c(lake_state_space$metric,lake_temporal$metric)),
              lake = unique(c(lake_state_space$lake,lake_temporal$lake))) |>
  left_join(state_transition_dates |>
            select(metric,lake,state_date),by = c("lake","metric"),
            relationship = "many-to-many") |>
  left_join(temporal_transition_dates |>
            select(metric,lake,temporal_date),by = c("lake","metric"),
            relationship = "many-to-many") |>
  group_by(metric,lake) |>
  slice_head(n=1) |>
  mutate(date_match = ifelse(state_date == temporal_date, TRUE, FALSE),
         threshold_date = ifelse(isTRUE(date_match),state_date,NA)) |>
  left_join(bimod_dates,by = c("lake","metric"))

write.csv(transition_dates,"Data/transition_dates.csv")

ggplot2::ggsave("Figures/figure_S1.pdf",
       ggpubr::ggarrange(p_temporal + ggtitle("a) time series"),
                         p_env +  ggtitle("b) state space"),
                         common.legend = T,ncol=1,nrow=2,legend = "right")
       ,width = 18,height=10,dpi=144)
