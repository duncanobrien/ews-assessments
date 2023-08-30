########################################################################################################################
# Figure 4 #
########################################################################################################################

require(tidyverse)
require(foreach)
require(ggh4x)
require(LaplacesDemon)
require(mousetrap)
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

#for Lake Kinneret, Washington and Windermere, find total plankton density per year across trophic levels and create environmental PCA
state.kin.dat <- data.frame("lake" = "Kinneret",
                            "tot_density" =scale(rowSums(scale(log1p(plank_env.data.yr[5:46,c(2,4:77)])))),
                            "env" = prcomp(scale(plank_env.data.yr[5:46,c("Temp","Nitrate","P_ort")]))$x[,1],
                            "phyto_density" = scale(log1p(rowSums(plank_env.data.yr[5:46,c(2,4:46)]))),
                            "zoo_density" =scale(log1p(rowSums(plank_env.data.yr[5:46,c(47:77)]))),
                            "date" = as.numeric(plank_env.data.yr$Date[5:46])) #unable to use first 4 years due to missing P_ort data

state.wind.dat <- data.frame("lake" = "Windermere",
                             "tot_density" =scale(rowSums(scale(log1p(phyto_env.windyrdata[,2:22])))),
                             "env" = prcomp(scale(phyto_env.windyrdata[,c("TEMP","NO3N","TOTP")]))$x[,1],
                             "phyto_density" = scale(log1p(rowSums(phyto_env.windyrdata[,2:18]))),
                             "zoo_density" =scale(log1p(rowSums(phyto_env.windyrdata[,20:23]))),
                             "date" = as.numeric(phyto_env.windyrdata$Date)) 

state.wash.dat <- data.frame("lake" = "Washington",
                             "tot_density" =scale(rowSums(scale(log1p(plank_env.washyrdata[,4:18])))),
                             "env" = prcomp(scale(plank_env.washyrdata[,c("temp","TP")]))$x[,1],
                             "phyto_density" = scale(log1p(rowSums(plank_env.washyrdata[,4:9]))),
                             "zoo_density" =scale(log1p(rowSums(plank_env.washyrdata[,10:18]))),
                             "date" = as.numeric(plank_env.washyrdata$date)) 

############ 
#Fit temporal TGAMs
############ 

lake_temporal <- foreach::foreach(j = list(state.kin.dat,state.wind.dat,state.wash.dat), 
                                  .combine = "rbind") %do% {
                                    
                                    lapply(c("tot_density","phyto_density","zoo_density"), 
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
                                             
                                           }) |> `names<-`(c("tot_density","phyto_density","zoo_density")) |>
                                      data.table::rbindlist(idcol = "metric",use.names = T) |>
                                      dplyr::mutate(threshold = tidyr::replace_na(as.character(threshold),"pre"),
                                                    start.date = dplyr::first(thresh.var),
                                                    last.date = dplyr::last(thresh.var),
                                                    lake = j$lake[1]
                                      )
                                  }

ss_time_plot_dat1 <- subset(lake_temporal, metric %in% c("phyto_density","zoo_density")) |>
  dplyr::mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton"))) |> #only interested in phytoplankton and zooplankton densities
  dplyr::filter((metric == "Phytoplankton" & lake == "Kinneret")| 
           (metric == "Phytoplankton" & lake == "Windermere") |
           (metric == "Zooplankton" & lake == "Washington"))

############ 
#Fit environmental TGAMs
############ 

lake_state_space <- foreach::foreach(j = list(state.kin.dat,state.wind.dat,state.wash.dat), 
                                     .combine = "rbind") %do% {
                                       
                                       lapply(c("tot_density","phyto_density","zoo_density"), 
                                              FUN = function(i){
                                                
                                                cont_formula <-  as.formula(paste(i,"~ s(env, bs='tp', k=3)"))
                                                thresh_formula <- as.formula(paste(i,"~ 1"))
                                                dens_gam <- compare_gam(data = j,
                                                                        cont_formula = formula(cont_formula, method = "REML"),
                                                                        thresh_formula =  formula(thresh_formula, method = "REML"),
                                                                        thresh.var = "date", expl.var = "env",
                                                                        thresh.range = c(0.2,0.75), by = 1, k=3)
                                                
                                                best_gam <- predict_best_gam(object=dens_gam) %>%
                                                  dplyr::rename("metric.val" = i) 
                                                
                                                return(best_gam)
                                                
                                              }) |> `names<-`(c("tot_density","phyto_density","zoo_density")) |>
                                         data.table::rbindlist(idcol = "metric",use.names = T) |>
                                         dplyr::mutate(threshold = tidyr::replace_na(as.character(threshold),"pre"),
                                                       start.date = dplyr::first(thresh.var),
                                                       last.date = dplyr::last(thresh.var),
                                                       lake = j$lake[1]
                                         )
                                     }
ss_plot_dat1 <- subset(lake_state_space, metric %in% c("phyto_density","zoo_density")) |>
  dplyr::mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton"))) |> #only interested in phytoplankton and zooplankton densities
  dplyr::filter((metric == "Phytoplankton" & lake == "Kinneret")| 
           (metric == "Phytoplankton" & lake == "Windermere") |
           (metric == "Zooplankton" & lake == "Washington"))


############ 
#Identify bi-modality
############ 

lake_bimod <- foreach::foreach(j = list(state.kin.dat,state.wind.dat,state.wash.dat), 
                               .combine = "rbind") %do% {
                                 
                                 lapply(c("tot_density","phyto_density","zoo_density"),function(i){
                                   
                                   
                                   return(data.frame("metric" = i,
                                                     "is_bimodal" =  LaplacesDemon::is.bimodal(j[,i]),
                                                     "modes" = LaplacesDemon::Modes(j[,i])$modes,
                                                     "modality_coef" = mousetrap::bimodality_coefficient(j[,i]))) #Pfister et al., 2013
                                 }) |>
                                   data.table::rbindlist() |>
                                   dplyr::mutate(lake = j$lake[1])
                               }

bimod_plot_datS1 <- rbind(state.kin.dat,state.wind.dat,state.wash.dat) |>
  dplyr::select(lake,phyto_density,zoo_density) |>
  tidyr::pivot_longer(-lake,names_to = "metric",values_to = "metric.value") |>
  dplyr::mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton"))) |>
  dplyr::filter((metric == "Phytoplankton" & lake == "Kinneret")| 
           (metric == "Phytoplankton" & lake == "Windermere") |
           (metric == "Zooplankton" & lake == "Washington"))


bimod_modes <- subset(lake_bimod, metric %in% c("phyto_density","zoo_density")) |>
  dplyr::mutate(metric = factor(metric, labels = c("Phytoplankton","Zooplankton"))) |>
  dplyr::filter((metric == "Phytoplankton" & lake == "Kinneret")| 
                  (metric == "Phytoplankton" & lake == "Windermere") |
                  (metric == "Zooplankton" & lake == "Washington")) |>
  dplyr::group_by(metric) |>
  dplyr::mutate(y = ifelse(metric == "Phytoplankton",c(0.525,0.5),c(0.725,0.7)),
                modes = round(modes,digits = 2),
                mode_sig = ifelse(modality_coef > 0.5, paste0(modes,"*"),paste0(modes)))

############ 
#Create figure
############ 

fig1_plot_dat <- rbind(ss_plot_dat1 |>
                         pivot_longer(env,names_to = "var_metric",values_to = "var.value") |>
                         mutate(var.value = var.value/2), #shrink for plotting purposes
                       ss_time_plot_dat1|>
                         pivot_longer(date,names_to = "var_metric",values_to = "var.value")) |>
  mutate(var_metric = factor(ifelse(var_metric == "env","Environment","Year"),levels = c("Year","Environment"), labels = c("Time series","State space")),
         title = factor(paste(lake,metric,sep = ": "),levels = c("Washington: Zooplankton","Kinneret: Phytoplankton","Windermere: Phytoplankton")),
         scenario = case_when(
           lake == "Washington" ~ "Abrupt shift",
           lake == "Kinneret" ~ "Critical transition",
           lake == "Windermere" ~ "No transition"
         ))


fig2a <- ggplot(data = subset(fig1_plot_dat,var_metric == "Time series"), aes(x=var.value,y=metric.val)) + 
  geom_point(aes(x=var.value, y = metric.val))+
  geom_path(aes(x=var.value, y = metric.val)) +  
  xlab("Year") + ylab("Scaled metric score") +
  ggh4x::facet_nested(title + scenario ~ var_metric, scales = "free_x") +
  geom_line(data = filter(fig1_plot_dat, threshold=="pre" & var_metric == "Time series"),aes(x=var.value,y=fit), col="blue",linewidth=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="pre"& var_metric == "Time series"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "Time series"),aes(x=var.value,y=fit), col="red",linewidth=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "Time series"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  geom_point(data = filter(fig1_plot_dat, thresh.var == start.date & var_metric == "Time series"),aes(x=var.value, y = metric.val,col="Start date"))+
  geom_point(data = filter(fig1_plot_dat, thresh.var == last.date & var_metric == "Time series"),aes(x=var.value, y = metric.val,col="End date"))+
  geom_point(data = filter(fig1_plot_dat,transition == "trans" & var_metric == "Time series"),aes(x=var.value, y = metric.val,col="Transition\ndates"))+
  ggrepel::geom_text_repel(data = filter(fig1_plot_dat,transition == "trans" & var_metric == "Time series"), aes(x=var.value, y = metric.val,label=thresh.var),force =1.5,nudge_x=5,nudge_y=-0.75,segment.linetype=2,min.segment.length = 0.1)+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates"), name = NULL, guide = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  scale_x_continuous(breaks = seq(1970,2010,20)) +
  theme(legend.text=element_text(size=14),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size=14),
        axis.title = element_text(size=14),
        legend.position = "top",
        panel.border = element_rect(linewidth = 1,fill = "transparent"))


fig2b <- ggplot(data = subset(fig1_plot_dat,var_metric == "State space"), aes(x=var.value,y=metric.val)) + 
  geom_point(aes(x=var.value, y = metric.val))+
  geom_path(aes(x=var.value, y = metric.val)) +  
  xlab("Environment score") + ylab("Scaled metric score") +
  ggh4x::facet_nested(title + scenario ~ var_metric, scales = "free_x") +
  geom_line(data = filter(fig1_plot_dat, threshold=="pre" & var_metric == "State space"),aes(x=var.value,y=fit), col="blue",linewidth=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="pre"& var_metric == "State space"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "State space"),aes(x=var.value,y=fit), col="red",linewidth=0.8, linetype = "solid")+
  geom_ribbon(data = filter(fig1_plot_dat, threshold=="post"& var_metric == "State space"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  geom_point(data = filter(fig1_plot_dat, thresh.var == start.date & var_metric == "State space"),aes(x=var.value, y = metric.val,col="Start date"))+
  geom_point(data = filter(fig1_plot_dat, thresh.var == last.date & var_metric == "State space"),aes(x=var.value, y = metric.val,col="End date"))+
  geom_point(data = filter(fig1_plot_dat,transition == "trans" & var_metric == "State space"),aes(x=var.value, y = metric.val,col="Transition\ndates"))+
  ggrepel::geom_text_repel(data = filter(fig1_plot_dat,transition == "trans" & var_metric == "State space"), aes(x=var.value, y = metric.val,label=thresh.var),force =1,nudge_x=0.55,nudge_y=-0.35,segment.linetype=2,min.segment.length = 0.1)+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates"), name = NULL, guide = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  theme(legend.text=element_text(size=14),
        strip.text.y = element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size=14),
        axis.title.x= element_text(size=14),
        legend.position = "top",
        panel.border = element_rect(linewidth = 1,fill = "transparent"))

fig1_plot_dat2 <- dplyr::right_join(bimod_plot_datS1,bimod_modes,
                                    relationship = "many-to-many") |>
  dplyr::group_by(metric,lake) |>
  dplyr::mutate(var_metric = "Bimodality",
                title = factor(paste(lake,metric,sep = ": "),levels = c("Washington: Zooplankton","Kinneret: Phytoplankton","Windermere: Phytoplankton"),
                               labels = c("Wash: Zoo","Kin: Phyto","Wind: Phyto")),
                scenario = case_when(
                  lake == "Washington" ~ "Abrupt shift",
                  lake == "Kinneret" ~ "Critical transition",
                  lake == "Windermere" ~ "No transition"
                  ))
                
fig2c <- ggplot(fig1_plot_dat2) +
  geom_density(aes(y = metric.value)) +
  geom_hline(aes(yintercept = modes),linetype = "dashed",colour = "black")+
  geom_text(aes(y=modes+0.35,x=0.2,label = mode_sig),size = 5)+
  ggh4x::facet_nested(title + scenario ~ var_metric, scales = "free_x",
                      labeller = label_value,
                      strip = ggh4x::strip_nested(size="constant",bleed=F, 
                                                  background_y = ggh4x::elem_list_rect(fill = c(rep("#D6D6D6",3),rep("white",3))),
                                                  background_x = ggh4x::elem_list_rect(fill = c(rep("white",3))),
                                                  by_layer_y = F)) +
  xlab("Density") + ylab("Scaled metric score")+ 
  coord_cartesian(xlim = c(0.05,0.6))+
  theme_classic() +
  theme(legend.text=element_text(size=14),
        #strip.background.x= element_rect(linewidth = 0),
        strip.placement.y = "outside",
        axis.title.y=element_blank(),
        strip.text = element_text(size=14),
        axis.title.x= element_text(size=14),
        panel.border = element_rect(linewidth = 1,fill = "transparent"))

ggplot2::ggsave("Figures/figure_4a.pdf",
                ggpubr::ggarrange(fig2a,fig2b,fig2c,
                                  ncol = 3,common.legend = T),
                width = 10,height=6,dpi=300)

############ 
#Create figure 4b
############

ggplot2::ggsave("Figures/figure_4b.pdf",
                ggplot(state.kas.dat,aes(x=date,y=phyto_density)) +
                  geom_point(col="#000080") +
                  geom_vline(xintercept = 2010.5, linetype = "longdash",color = "grey40")+
                  xlab("Year") + ylab("Daphnid\ndensity (ind/ml)") +
                  theme_classic() +
                  theme(panel.background = element_rect(fill='transparent'), 
                         plot.background = element_rect(fill='transparent', color=NA), 
                         legend.background = element_rect(fill='transparent'),
                         legend.box.background = element_rect(fill='transparent') ),
                width = 3.5,height=2,dpi=300,bg='transparent')

ggplot2::ggsave("Figures/figure_4c.pdf",
                ggplot(state.kas.dat,aes(x=date,y=zoo_density)) +
                  geom_point(col="#552200") +
                  geom_vline(xintercept = 2010.5, linetype = "longdash",color = "grey40")+
                  xlab("Year") + ylab("Cyclopoid\ndensity (ind/ml)") +
                  theme_classic() +
                  theme(panel.background = element_rect(fill='transparent'), 
                        plot.background = element_rect(fill='transparent', color=NA), 
                        legend.background = element_rect(fill='transparent'),
                        legend.box.background = element_rect(fill='transparent') ),
                width = 3.5,height=2,dpi=300,bg='transparent')

