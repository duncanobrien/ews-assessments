########################################################################################################################
# Figure S1 #
########################################################################################################################

require(tidyverse)
require(foreach)
require(patchwork)

source("Data/kinneret_plankton_data.R") #raw data not provided. See references in Data Availability statement and repository README
source("Data/kinneret_environmental_data.R")

source("Data/zurich_plankton_data.R")
source("Data/LZ_environmental_data.R")

source("Data/madison_plankton_data.R")
source("Data/madison_environmental_data.R")

source("Data/windermere_plankton_data.R")
source("Data/windermere_environmental_data.R")

source("Data/kasumigaura_plankton_data.R")
source("Data/Kasumigaura_environmental_data.R")

source("Data/UZ_plankton_data.R")

source("Data/monona_plankton_data.R")

source("Data/leven_plankton_data.R")

source("Data/washington_plankton_data.R")

rm(list = ls()[!ls() %in% c("plank_env.data.mth","plank_env.data.yr","plank_env.kasmthdata","plank_env.kasyrdata",
                            "plank_env.levemthdata","plank_env.leveyrdata","plank_env.LZmthdata","plank_env.LZyrdata",
                            "plank_env.madmthdata","plank_env.madyrdata","plank_env.monmthdata","plank_env.monyrdata",
                            "plank_env.UZmthdata","plank_env.UZyrdata","plank_env.washmthdata","plank_env.washyrdata",
                            "phyto_env.windmthdata","phyto_env.windyrdata")]
)

source("Code/threshold_gam.R")

############ 
#Prepare plankton data
############ 

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

############ 
#Fit temporal TGAMs
############ 

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


############ 
#Fit environmental TGAMs
############ 

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

############ 
#Create figureS1
############ 

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


ggplot2::ggsave("Figures/figure_S1.png",
                ggpubr::ggarrange(p_temporalS1 + ggtitle("Time series"),
                                  p_envS1 +  ggtitle("State space"),
                                  labels="AUTO",font.label = list(face="plain"),
                                  common.legend = T,ncol=1,nrow=2,legend = "top")
                ,width = 16,height=9,dpi=300)
