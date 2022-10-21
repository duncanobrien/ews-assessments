require(Matrix)
require(tidyverse)
require(foreach)

################################################################################################################
## Kinneret State Space ##
################################################################################################################

source("/Users/ul20791/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_environmental_data.R")
#phyo, zoo + combined raw concentrations at week, month and year groupings

state.kin.dat <- data.frame("density" =scale(log1p(rowSums(plank_env.data.yr[5:46,c(2,4:77)]))),
                            "community" = prcomp(scale(plank_env.data.yr[5:46,c(2,4:77)]))$x[,1],
                            "PC1" = prcomp(scale(plank_env.data.yr[5:46,c("Temp","Nitrate","P_ort","Water.Level")]))$x[,1],
                            #"PC1" = prcomp(scale(plank_env.data.yr[5:46,c("Temp","Nitrate","P_ort")]))$x[,1],
                            "date" = as.numeric(plank_env.data.yr$Date[5:46])) #unable to use first 4 years due to missing P_ort data
mid.kin<-median(state.kin.dat$date)

kin_state_gam_yr <-foreach(i = c("density","community")) %do%{
  res <- list() #result out
  best <- list() #assess best fitting model
  thresh.values<-c()
  
  sub.dat <- data.frame(date=state.kin.dat$date,PC = state.kin.dat$PC1,metric=state.kin.dat[,paste(i)])
  
  gamm.cont<- mgcv::gam(metric~ s(PC, bs="tp", k=10) , data = sub.dat,
                        family = gaussian(), method = "REML")
  cont.summary <- summary(gamm.cont)
  cont.out <- data.frame(rbind(cont.summary$p.table,cont.summary$s.table),
                         method = rep("continuous",length.out = 2),
                         metric = rep(paste(i),length.out = 2))
  
  res$cont$mod <- gamm.cont
  res$cont$summary <- cont.out
  res$cont$AIC <- gamm.cont$gcv.ubre
  res$cont$method <- "continuous"
  res$cont$threshold <- NA
  best$cont <-gamm.cont$gcv.ubre
  if(res$cont$summary$Pr...t..[2] < 0.1 && res$cont$summary$Estimate[2] >1.00){
    
    thresh.values<-round(seq(stats::quantile(state.kin.dat$date[],0.2),
                             stats::quantile(state.kin.dat$date[],0.8),by=1)) #round to nearest whole number
    aic.thresh<-thresh.values*NA
    
    for(jj in 1:length(thresh.values)){
      sub.dat$tmp.thresh <- factor(ifelse(sub.dat$date<thresh.values[jj],'pre','post'), levels = c("pre","post"))
      gam.test<- mgcv::gam(metric~ tmp.thresh + s(PC, bs="tp", k=5,by=tmp.thresh)-1, data = sub.dat, 
                           family = gaussian(), method = "REML")
      aic.thresh[jj] <- gam.test$gcv.ubre
    }
    sub.dat$threshold<-factor(ifelse(sub.dat$date<thresh.values[order(aic.thresh)][1],'pre','post'),levels = c("pre","post"))
    
    gamm.thresh<- mgcv::gam(metric~ threshold + s(PC, bs="tp", k=5,by=as.factor(threshold))-1, 
                            data = sub.dat, 
                            family = gaussian(), method = "REML")  
    thresh.summary <- summary(gamm.thresh)
    thresh.out <- data.frame(rbind(thresh.summary$p.table,thresh.summary$s.table),
                             method = rep("threshold",length.out = 4),
                             FD = rep(paste(i),length.out = 4))
    
    res$thresh$mod <- gamm.thresh
    res$thresh$summary <- thresh.out
    res$thresh$AIC <- gamm.thresh$gcv.ubre
    res$thresh$method <- "threshold"
    res$thresh$threshold <- thresh.values[base::order(aic.thresh)][1]
    res$thresh$threshold.gcv <-cbind("thresh.date" = thresh.values[base::order(aic.thresh)],"gcv" = aic.thresh[base::order(aic.thresh)])
    best$thresh <-gamm.thresh$gcv.ubre
  }else{res$thresh$mod <- "NS"
  res$thresh$AIC <- Inf
  }
  return(res)
}|>
  `names<-`(c("density","community"))

kin.state.gam.yr.pred <- as.data.frame(foreach(i = c("density","community"),.combine="cbind",.packages=c('tidyverse',"magrittr")) %do%{
  if(abs(kin_state_gam_yr[[i]]$cont$AIC-kin_state_gam_yr[[i]]$thresh$AIC) <=6){
    best <- 1
  }else{
    best <- which.min(c(kin_state_gam_yr[[i]]$cont$AIC,kin_state_gam_yr[[i]]$thresh$AIC))
  }
  best.mod <- kin_state_gam_yr[[i]][[best]] #subset best model (continuous vs threshold based upon AIC)
  sub.dat <- data.frame(date=state.kin.dat$date,PC = state.kin.dat$PC1,metric=state.kin.dat[,i]) %>%
    mutate(threshold = ifelse(date>= best.mod$threshold,'post','pre')) %>%
    mutate(transition = ifelse(threshold == "pre" & lead(threshold) == "post",'trans',
                               ifelse(threshold == "post" & lag(threshold) == "pre",'trans',      
                                      'no.trans')))
  pred <- predict(best.mod$mod, newdata = sub.dat, type = 'link', se.fit = TRUE,
                  exclude = c("s(month)"))
  
  tmp <- cbind(sub.dat,pred) %>%
    set_colnames(c(paste(i,"date",sep = "_"),paste(i,"PC",sep = "_"),paste(i,"value",sep = "_"),
                   paste(i,"threshold",sep = "_"),paste(i,"transition",sep = "_"),paste(i,"fit",sep = "_"),
                   paste(i,"CI",sep = "_")))
  return(tmp)
}) %>% #drop accum column caused by .init specification
  pivot_longer(everything(), 
               names_to = c("metric",".value"), 
               names_pattern = "(.*)_(.*)") %>%
  mutate(threshold = replace_na(as.character(threshold),"pre"),
         lake = "Kinneret") #replace NAs in threshold variable that appear for continuous GAMs

kin.map <- ggplot(data =kin.state.gam.yr.pred,aes(x=PC,y=value)) + 
  geom_point(aes(x=PC, y = value))+
  geom_path(aes(x=PC, y = value)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.kin, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + facet_wrap(~metric,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  geom_line(data = filter(kin.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(kin.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(kin.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(kin.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(kin.state.gam.yr.pred,transition == "trans"), aes(x=PC, y = value,label=date),force =2,nudge_x=2,segment.linetype=2)+
  xlim(-4,4)

################################################################################################################
## Mendota State Space ##
################################################################################################################

source("/Users/ul20791/Desktop/Academia/PhD/Data/Madison/Data/madison_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Madison/Data/madison_environmental_data.R")

state.mad.dat <- data.frame("density" =scale(log1p(rowSums(plank_env.madyrdata[,2:222]))),
                            "community" = prcomp(scale(plank_env.madyrdata[,c(2:222)]))$x[,1],
                            "PC1" = prcomp(scale(plank_env.madyrdata[,c("wtemp","totP","NO3N","water.lvl")]))$x[,1],
                            #"PC1" = prcomp(scale(plank_env.madyrdata[,c("wtemp","totP","NO3N")]))$x[,1],
                            "date" = as.numeric(plank_env.madyrdata$date)) #unable to use first 4 years due to missing P_ort data
mid.mad<-median(state.mad.dat$date)

mad_state_gam_yr <-foreach(i = c("density","community")) %do%{
  res <- list() #result out
  best <- list() #assess best fitting model
  thresh.values<-c()
  
  sub.dat <- data.frame(date=state.mad.dat$date,PC = state.mad.dat$PC1,metric=state.mad.dat[,paste(i)])
  
  gamm.cont<- mgcv::gam(metric~ s(PC, bs="tp", k=10) , data = sub.dat,
                        family = gaussian(), method = "REML")
  cont.summary <- summary(gamm.cont)
  cont.out <- data.frame(rbind(cont.summary$p.table,cont.summary$s.table),
                         method = rep("continuous",length.out = 2),
                         metric = rep(paste(i),length.out = 2))
  
  res$cont$mod <- gamm.cont
  res$cont$summary <- cont.out
  res$cont$AIC <- gamm.cont$gcv.ubre
  res$cont$method <- "continuous"
  res$cont$threshold <- NA
  best$cont <-gamm.cont$gcv.ubre
  if(res$cont$summary$Pr...t..[2] < 0.1 && res$cont$summary$Estimate[2] >1.00){
    
    thresh.values<-round(seq(stats::quantile(state.mad.dat$date[],0.2),
                             stats::quantile(state.mad.dat$date[],0.8),by=1)) #round to nearest whole number
    aic.thresh<-thresh.values*NA
    
    for(jj in 1:length(thresh.values)){
      sub.dat$tmp.thresh <- factor(ifelse(sub.dat$date<thresh.values[jj],'pre','post'), levels = c("pre","post"))
      gam.test<- mgcv::gam(metric~ tmp.thresh + s(PC, bs="tp", k=5,by=tmp.thresh)-1, data = sub.dat, 
                           family = gaussian(), method = "REML")
      aic.thresh[jj] <- gam.test$gcv.ubre
    }
    sub.dat$threshold<-factor(ifelse(sub.dat$date<thresh.values[order(aic.thresh)][1],'pre','post'),levels = c("pre","post"))
    
    gamm.thresh<- mgcv::gam(metric~ threshold + s(PC, bs="tp", k=5,by=as.factor(threshold))-1, 
                            data = sub.dat, 
                            family = gaussian(), method = "REML")  
    thresh.summary <- summary(gamm.thresh)
    thresh.out <- data.frame(rbind(thresh.summary$p.table,thresh.summary$s.table),
                             method = rep("threshold",length.out = 4),
                             FD = rep(paste(i),length.out = 4))
    
    res$thresh$mod <- gamm.thresh
    res$thresh$summary <- thresh.out
    res$thresh$AIC <- gamm.thresh$gcv.ubre
    res$thresh$method <- "threshold"
    res$thresh$threshold <- thresh.values[base::order(aic.thresh)][1]
    res$thresh$threshold.gcv <-cbind("thresh.date" = thresh.values[base::order(aic.thresh)],"gcv" = aic.thresh[base::order(aic.thresh)])
    best$thresh <-gamm.thresh$gcv.ubre
  }else{res$thresh$mod <- "NS"
  res$thresh$AIC <- Inf
  }
  return(res)
} |>
  `names<-`(c("density","community"))

mad.state.gam.yr.pred <- as.data.frame(foreach(i = c("density","community"),.combine="cbind",.packages=c('tidyverse',"magrittr")) %do%{
  if(abs(mad_state_gam_yr[[i]]$cont$AIC-mad_state_gam_yr[[i]]$thresh$AIC) <=6){
    best <- 1
  }else{
    best <- which.min(c(mad_state_gam_yr[[i]]$cont$AIC,mad_state_gam_yr[[i]]$thresh$AIC))
  }
  best.mod <- mad_state_gam_yr[[i]][[best]] #subset best model (continuous vs threshold based upon AIC)
  sub.dat <- data.frame(date=state.mad.dat$date,PC = state.mad.dat$PC1,metric=state.mad.dat[,i]) %>%
    mutate(threshold = ifelse(date>= best.mod$threshold,'post','pre')) %>%
    mutate(transition = ifelse(threshold == "pre" & lead(threshold) == "post",'trans',
                               ifelse(threshold == "post" & lag(threshold) == "pre",'trans',      
                                      'no.trans')))
  pred <- predict(best.mod$mod, newdata = sub.dat, type = 'link', se.fit = TRUE,
                  exclude = c("s(month)"))
  
  tmp <- cbind(sub.dat,pred) %>%
    set_colnames(c(paste(i,"date",sep = "_"),paste(i,"PC",sep = "_"),paste(i,"value",sep = "_"),
                   paste(i,"threshold",sep = "_"),paste(i,"transition",sep = "_"),paste(i,"fit",sep = "_"),
                   paste(i,"CI",sep = "_")))
  return(tmp)
}) %>% #drop accum column caused by .init specification
  pivot_longer(everything(), 
               names_to = c("metric",".value"), 
               names_pattern = "(.*)_(.*)") %>%
  mutate(threshold = replace_na(as.character(threshold),"pre"),
         lake = "Mendota") #replace NAs in threshold variable that appear for continuous GAMs

mad.map <- ggplot(data =mad.state.gam.yr.pred,aes(x=PC,y=value)) + 
  geom_point(aes(x=PC, y = value))+
  #geom_path(aes(x=PC, y = value)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.mad, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + facet_wrap(~metric,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  geom_line(data = filter(mad.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(mad.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(mad.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(mad.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(mad.state.gam.yr.pred,transition == "trans"), aes(x=PC, y = value,label=date),force =2,nudge_x=2,segment.linetype=2)+
  xlim(-4,4)

################################################################################################################
## Kasumigaura State Space ##
################################################################################################################

source("/Users/ul20791/Desktop/Academia/PhD/Data/Kasumigaura/Data/kasumigaura_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Kasumigaura/Data/Kasumigaura_environmental_data.R")

state.kas.dat <- data.frame("density" = scale(log1p(rowSums(plank_env.kasyrdata[,2:154]))),
                            "community" = prcomp(scale(plank_env.kasyrdata[,2:154]))$x[,1],
                            "PC1" = prcomp(scale(plank_env.kasyrdata[,c("wtemp","totP","NO3N","water.lvl")]))$x[,1],
                            #"PC1" = prcomp(scale(plank_env.kasyrdata[,c("wtemp","totP","NO3N")]))$x[,1],
                            "date" = as.numeric(plank_env.kasyrdata$date)) #unable to use first 4 years due to missing P_ort data
mid.kas<-median(state.kas.dat$date)

kas_state_gam_yr <-foreach(i = c("density","community")) %do%{
  res <- list() #result out
  best <- list() #assess best fitting model
  thresh.values<-c()
  
  sub.dat <- data.frame(date=state.kas.dat$date,PC = state.kas.dat$PC1,metric=state.kas.dat[,paste(i)])
  
  gamm.cont<- mgcv::gam(metric~ s(PC, bs="tp", k=10) , data = sub.dat,
                        family = gaussian(), method = "REML")
  cont.summary <- summary(gamm.cont)
  cont.out <- data.frame(rbind(cont.summary$p.table,cont.summary$s.table),
                         method = rep("continuous",length.out = 2),
                         metric = rep(paste(i),length.out = 2))
  
  res$cont$mod <- gamm.cont
  res$cont$summary <- cont.out
  res$cont$AIC <- gamm.cont$gcv.ubre
  res$cont$method <- "continuous"
  res$cont$threshold <- NA
  best$cont <-gamm.cont$gcv.ubre
  if(res$cont$summary$Pr...t..[2] < 0.1 && res$cont$summary$Estimate[2] >1.00){
    
    thresh.values<-round(seq(stats::quantile(state.kas.dat$date[],0.2),
                             stats::quantile(state.kas.dat$date[],0.8),by=1)) #round to nearest whole number
    aic.thresh<-thresh.values*NA
    
    for(jj in 1:length(thresh.values)){
      sub.dat$tmp.thresh <- factor(ifelse(sub.dat$date<thresh.values[jj],'pre','post'), levels = c("pre","post"))
      gam.test<- mgcv::gam(metric~ tmp.thresh + s(PC, bs="tp", k=5,by=tmp.thresh)-1, data = sub.dat, 
                           family = gaussian(), method = "REML")
      aic.thresh[jj] <- gam.test$gcv.ubre
    }
    sub.dat$threshold<-factor(ifelse(sub.dat$date<thresh.values[order(aic.thresh)][1],'pre','post'),levels = c("pre","post"))
    
    gamm.thresh<- mgcv::gam(metric~ threshold + s(PC, bs="tp", k=5,by=as.factor(threshold))-1, 
                            data = sub.dat, 
                            family = gaussian(), method = "REML")  
    thresh.summary <- summary(gamm.thresh)
    thresh.out <- data.frame(rbind(thresh.summary$p.table,thresh.summary$s.table),
                             method = rep("threshold",length.out = 4),
                             FD = rep(paste(i),length.out = 4))
    
    res$thresh$mod <- gamm.thresh
    res$thresh$summary <- thresh.out
    res$thresh$AIC <- gamm.thresh$gcv.ubre
    res$thresh$method <- "threshold"
    res$thresh$threshold <- thresh.values[base::order(aic.thresh)][1]
    res$thresh$threshold.gcv <-cbind("thresh.date" = thresh.values[base::order(aic.thresh)],"gcv" = aic.thresh[base::order(aic.thresh)])
    best$thresh <-gamm.thresh$gcv.ubre
  }else{res$thresh$mod <- "NS"
  res$thresh$AIC <- Inf
  }
  return(res)
}|>
  `names<-`(c("density","community"))

kas.state.gam.yr.pred <- as.data.frame(foreach(i = c("density","community"),.combine="cbind",.packages=c('tidyverse',"magrittr")) %do%{
  if(abs(kas_state_gam_yr[[i]]$cont$AIC-kas_state_gam_yr[[i]]$thresh$AIC) <=6){
    best <- 1
  }else{
    best <- which.min(c(kas_state_gam_yr[[i]]$cont$AIC,kas_state_gam_yr[[i]]$thresh$AIC))
  }
  best.mod <- kas_state_gam_yr[[i]][[best]] #subset best model (continuous vs threshold based upon AIC)
  sub.dat <- data.frame(date=state.kas.dat$date,PC = state.kas.dat$PC1,metric=state.kas.dat[,i]) %>%
    mutate(threshold = ifelse(date>= best.mod$threshold,'post','pre')) %>%
    mutate(transition = ifelse(threshold == "pre" & lead(threshold) == "post",'trans',
                               ifelse(threshold == "post" & lag(threshold) == "pre",'trans',      
                                      'no.trans')))
  pred <- predict(best.mod$mod, newdata = sub.dat, type = 'link', se.fit = TRUE,
                  exclude = c("s(month)"))
  
  tmp <- cbind(sub.dat,pred) %>%
    set_colnames(c(paste(i,"date",sep = "_"),paste(i,"PC",sep = "_"),paste(i,"value",sep = "_"),
                   paste(i,"threshold",sep = "_"),paste(i,"transition",sep = "_"),paste(i,"fit",sep = "_"),
                   paste(i,"CI",sep = "_")))
  return(tmp)
}) %>% #drop accum column caused by .init specification
  pivot_longer(everything(), 
               names_to = c("metric",".value"), 
               names_pattern = "(.*)_(.*)") %>%
  mutate(threshold = replace_na(as.character(threshold),"pre"),
         lake = "Kasumigaura") #replace NAs in threshold variable that appear for continuous GAMs

kas.map <- ggplot(data =kas.state.gam.yr.pred,aes(x=PC,y=value)) + 
  geom_point(aes(x=PC, y = value))+
  #geom_path(aes(x=PC, y = value)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.kas, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  theme_bw() + facet_wrap(~metric,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  geom_line(data = filter(kas.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(kas.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(kas.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(kas.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(kas.state.gam.yr.pred,transition == "trans"), aes(x=PC, y = value,label=date),force =2,nudge_x=2,segment.linetype=2)+
  xlim(-4,4)

################################################################################################################
## Lower Zurich State Space ##
################################################################################################################

source("/Users/ul20791/Desktop/Academia/PhD/Data/Zurich/Data/zurich_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Zurich/Data/LZ_environmental_data.R")

state.LZ.dat <- data.frame("density" = scale(log1p(rowSums(plank_env.LZyrdata[,5:206]))),
                            "community" = prcomp(scale(plank_env.LZyrdata[,5:206]))$x[,1],
                            "PC1" = prcomp(scale(plank_env.LZyrdata[,c("mean.t","mean.po4","NO3_N")]))$x[,1],
                            "date" = as.numeric(plank_env.LZyrdata$date)) #unable to use first 4 years due to missing P_ort data
mid.LZ<-median(state.LZ.dat$date)

LZ_state_gam_yr <-foreach(i = c("density","community")) %do%{
  res <- list() #result out
  best <- list() #assess best fitting model
  thresh.values<-c()
  
  sub.dat <- data.frame(date=state.LZ.dat$date,PC = state.LZ.dat$PC1,metric=state.LZ.dat[,paste(i)])
  
  gamm.cont<- mgcv::gam(metric~ s(PC, bs="tp", k=10) , data = sub.dat,
                        family = gaussian(), method = "REML")
  cont.summary <- summary(gamm.cont)
  cont.out <- data.frame(rbind(cont.summary$p.table,cont.summary$s.table),
                         method = rep("continuous",length.out = 2),
                         metric = rep(paste(i),length.out = 2))
  
  res$cont$mod <- gamm.cont
  res$cont$summary <- cont.out
  res$cont$AIC <- gamm.cont$gcv.ubre
  res$cont$method <- "continuous"
  res$cont$threshold <- NA
  best$cont <-gamm.cont$gcv.ubre
  if(res$cont$summary$Pr...t..[2] < 0.1 && res$cont$summary$Estimate[2] >1.00){
    
    thresh.values<-round(seq(stats::quantile(state.LZ.dat$date[],0.2),
                             stats::quantile(state.LZ.dat$date[],0.8),by=1)) #round to nearest whole number
    aic.thresh<-thresh.values*NA
    
    for(jj in 1:length(thresh.values)){
      sub.dat$tmp.thresh <- factor(ifelse(sub.dat$date<thresh.values[jj],'pre','post'), levels = c("pre","post"))
      gam.test<- mgcv::gam(metric~ tmp.thresh + s(PC, bs="tp", k=5,by=tmp.thresh)-1, data = sub.dat, 
                           family = gaussian(), method = "REML")
      aic.thresh[jj] <- gam.test$gcv.ubre
    }
    sub.dat$threshold<-factor(ifelse(sub.dat$date<thresh.values[order(aic.thresh)][1],'pre','post'),levels = c("pre","post"))
    
    gamm.thresh<- mgcv::gam(metric~ threshold + s(PC, bs="tp", k=5,by=as.factor(threshold))-1, 
                            data = sub.dat, 
                            family = gaussian(), method = "REML")  
    thresh.summary <- summary(gamm.thresh)
    thresh.out <- data.frame(rbind(thresh.summary$p.table,thresh.summary$s.table),
                             method = rep("threshold",length.out = 4),
                             FD = rep(paste(i),length.out = 4))
    
    res$thresh$mod <- gamm.thresh
    res$thresh$summary <- thresh.out
    res$thresh$AIC <- gamm.thresh$gcv.ubre
    res$thresh$method <- "threshold"
    res$thresh$threshold <- thresh.values[base::order(aic.thresh)][1]
    res$thresh$threshold.gcv <-cbind("thresh.date" = thresh.values[base::order(aic.thresh)],"gcv" = aic.thresh[base::order(aic.thresh)])
    best$thresh <-gamm.thresh$gcv.ubre
  }else{res$thresh$mod <- "NS"
  res$thresh$AIC <- Inf
  }
  return(res)
}|>
  `names<-`(c("density","community"))

LZ.state.gam.yr.pred <- as.data.frame(foreach(i = c("density","community"),.combine="cbind",.packages=c('tidyverse',"magrittr")) %do%{
  if(abs(LZ_state_gam_yr[[i]]$cont$AIC-LZ_state_gam_yr[[i]]$thresh$AIC) <=6){
    best <- 1
  }else{
    best <- which.min(c(LZ_state_gam_yr[[i]]$cont$AIC,LZ_state_gam_yr[[i]]$thresh$AIC))
  }
  best.mod <- LZ_state_gam_yr[[i]][[best]] #subset best model (continuous vs threshold based upon AIC)
  sub.dat <- data.frame(date=state.LZ.dat$date,PC = state.LZ.dat$PC1,metric=state.LZ.dat[,i]) %>%
    mutate(threshold = ifelse(date>= best.mod$threshold,'post','pre')) %>%
    mutate(transition = ifelse(threshold == "pre" & lead(threshold) == "post",'trans',
                               ifelse(threshold == "post" & lag(threshold) == "pre",'trans',      
                                      'no.trans')))
  pred <- predict(best.mod$mod, newdata = sub.dat, type = 'link', se.fit = TRUE,
                  exclude = c("s(month)"))
  
  tmp <- cbind(sub.dat,pred) %>%
    set_colnames(c(paste(i,"date",sep = "_"),paste(i,"PC",sep = "_"),paste(i,"value",sep = "_"),
                   paste(i,"threshold",sep = "_"),paste(i,"transition",sep = "_"),paste(i,"fit",sep = "_"),
                   paste(i,"CI",sep = "_")))
  return(tmp)
}) %>% #drop accum column caused by .init specification
  pivot_longer(everything(), 
               names_to = c("metric",".value"), 
               names_pattern = "(.*)_(.*)") %>%
  mutate(threshold = replace_na(as.character(threshold),"pre"),
         lake = "Lower Zurich") #replace NAs in threshold variable that appear for continuous GAMs

LZ.map <- ggplot(data =LZ.state.gam.yr.pred,aes(x=PC,y=value)) + 
  geom_point(aes(x=PC, y = value))+
  #geom_path(aes(x=PC, y = value)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.LZ, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + facet_wrap(~metric,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  geom_line(data = filter(LZ.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(LZ.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(LZ.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(LZ.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(LZ.state.gam.yr.pred,transition == "trans"), aes(x=PC, y = value,label=date),force =2,nudge_x=2,segment.linetype=2)+
  xlim(-4,4)

################################################################################################################
## Upper Zurich State Space ##
################################################################################################################

source("/Users/ul20791/Desktop/Academia/PhD/Data/Zurich/Data/UZ_plankton_data.R")

state.UZ.dat <- data.frame("density" = scale(log1p(rowSums(plank_env.UZyrdata[,5:95]))),
                           "community" = prcomp(scale(plank_env.UZyrdata[,5:95]))$x[,1],
                           "PC1" = prcomp(scale(plank_env.UZyrdata[,c("mean.t","mean.po4")]))$x[,1],
                           "date" = as.numeric(plank_env.UZyrdata$date)) #unable to use first 4 years due to missing P_ort data
mid.UZ<-median(state.UZ.dat$date)

UZ_state_gam_yr <-foreach(i = c("density","community")) %do%{
  res <- list() #result out
  best <- list() #assess best fitting model
  thresh.values<-c()
  
  sub.dat <- data.frame(date=state.UZ.dat$date,PC = state.UZ.dat$PC1,metric=state.UZ.dat[,paste(i)])
  
  gamm.cont<- mgcv::gam(metric~ s(PC, bs="tp", k=10) , data = sub.dat,
                        family = gaussian(), method = "REML")
  cont.summary <- summary(gamm.cont)
  cont.out <- data.frame(rbind(cont.summary$p.table,cont.summary$s.table),
                         method = rep("continuous",length.out = 2),
                         metric = rep(paste(i),length.out = 2))
  
  res$cont$mod <- gamm.cont
  res$cont$summary <- cont.out
  res$cont$AIC <- gamm.cont$gcv.ubre
  res$cont$method <- "continuous"
  res$cont$threshold <- NA
  best$cont <-gamm.cont$gcv.ubre
  if(res$cont$summary$Pr...t..[2] < 0.1 && res$cont$summary$Estimate[2] >1.00){
    
    thresh.values<-round(seq(stats::quantile(state.UZ.dat$date[],0.2),
                             stats::quantile(state.UZ.dat$date[],0.8),by=1)) #round to nearest whole number
    aic.thresh<-thresh.values*NA
    
    for(jj in 1:length(thresh.values)){
      sub.dat$tmp.thresh <- factor(ifelse(sub.dat$date<thresh.values[jj],'pre','post'), levels = c("pre","post"))
      gam.test<- mgcv::gam(metric~ tmp.thresh + s(PC, bs="tp", k=5,by=tmp.thresh)-1, data = sub.dat, 
                           family = gaussian(), method = "REML")
      aic.thresh[jj] <- gam.test$gcv.ubre
    }
    sub.dat$threshold<-factor(ifelse(sub.dat$date<thresh.values[order(aic.thresh)][1],'pre','post'),levels = c("pre","post"))
    
    gamm.thresh<- mgcv::gam(metric~ threshold + s(PC, bs="tp", k=5,by=as.factor(threshold))-1, 
                            data = sub.dat, 
                            family = gaussian(), method = "REML")  
    thresh.summary <- summary(gamm.thresh)
    thresh.out <- data.frame(rbind(thresh.summary$p.table,thresh.summary$s.table),
                             method = rep("threshold",length.out = 4),
                             FD = rep(paste(i),length.out = 4))
    
    res$thresh$mod <- gamm.thresh
    res$thresh$summary <- thresh.out
    res$thresh$AIC <- gamm.thresh$gcv.ubre
    res$thresh$method <- "threshold"
    res$thresh$threshold <- thresh.values[base::order(aic.thresh)][1]
    res$thresh$threshold.gcv <-cbind("thresh.date" = thresh.values[base::order(aic.thresh)],"gcv" = aic.thresh[base::order(aic.thresh)])
    best$thresh <-gamm.thresh$gcv.ubre
  }else{res$thresh$mod <- "NS"
  res$thresh$AIC <- Inf
  }
  return(res)
}|>
  `names<-`(c("density","community"))

UZ.state.gam.yr.pred <- as.data.frame(foreach(i = c("density","community"),.combine="cbind",.packages=c('tidyverse',"magrittr")) %do%{
  if(abs(UZ_state_gam_yr[[i]]$cont$AIC-UZ_state_gam_yr[[i]]$thresh$AIC) <=6){
    best <- 1
  }else{
    best <- which.min(c(UZ_state_gam_yr[[i]]$cont$AIC,UZ_state_gam_yr[[i]]$thresh$AIC))
  }
  best.mod <- UZ_state_gam_yr[[i]][[best]] #subset best model (continuous vs threshold based upon AIC)
  sub.dat <- data.frame(date=state.UZ.dat$date,PC = state.UZ.dat$PC1,metric=state.UZ.dat[,i]) %>%
    mutate(threshold = ifelse(date>= best.mod$threshold,'post','pre')) %>%
    mutate(transition = ifelse(threshold == "pre" & lead(threshold) == "post",'trans',
                               ifelse(threshold == "post" & lag(threshold) == "pre",'trans',      
                                      'no.trans')))
  pred <- predict(best.mod$mod, newdata = sub.dat, type = 'link', se.fit = TRUE,
                  exclude = c("s(month)"))
  
  tmp <- cbind(sub.dat,pred) %>%
    set_colnames(c(paste(i,"date",sep = "_"),paste(i,"PC",sep = "_"),paste(i,"value",sep = "_"),
                   paste(i,"threshold",sep = "_"),paste(i,"transition",sep = "_"),paste(i,"fit",sep = "_"),
                   paste(i,"CI",sep = "_")))
  return(tmp)
}) %>% #drop accum column caused by .init specification
  pivot_longer(everything(), 
               names_to = c("metric",".value"), 
               names_pattern = "(.*)_(.*)") %>%
  mutate(threshold = replace_na(as.character(threshold),"pre"),
         lake = "Upper Zurich") #replace NAs in threshold variable that appear for continuous GAMs

UZ.map <- ggplot(data =UZ.state.gam.yr.pred,aes(x=PC,y=value)) + 
  geom_point(aes(x=PC, y = value))+
  #geom_path(aes(x=PC, y = value)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.UZ, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + facet_wrap(~metric,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  geom_line(data = filter(UZ.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(UZ.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(UZ.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(UZ.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(UZ.state.gam.yr.pred,transition == "trans"), aes(x=PC, y = value,label=date),force =2,nudge_x=2,segment.linetype=2)+
  xlim(-4,4)

################################################################################################################
## Windermere State Space ##
################################################################################################################

source("/Users/ul20791/Desktop/Academia/PhD/Data/Windermere/Data/windermere_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Windermere/Data/windermere_environmental_data.R")

state.wind.dat <- data.frame("density" =scale(log1p(rowSums(phyto_env.windyrdata[,c(2:18,20:23)]))),
                           "community" = prcomp(scale(phyto_env.windyrdata[,c(2:18,20:23)]))$x[,1],
                           "PC1" = prcomp(scale(phyto_env.windyrdata[,c("TEMP","TOTP","NO3N")]))$x[,1],
                           "date" = as.numeric(phyto_env.windyrdata$Date)) #unable to use first 4 years due to missing P_ort data
mid.wind<-median(state.wind.dat$date)

wind_state_gam_yr <-foreach(i = c("density","community")) %do%{
  res <- list() #result out
  best <- list() #assess best fitting model
  thresh.values<-c()
  
  sub.dat <- data.frame(date=state.wind.dat$date,PC = state.wind.dat$PC1,metric=state.wind.dat[,paste(i)])
  
  gamm.cont<- mgcv::gam(metric~ s(PC, bs="tp", k=10) , data = sub.dat,
                        family = gaussian(), method = "REML")
  cont.summary <- summary(gamm.cont)
  cont.out <- data.frame(rbind(cont.summary$p.table,cont.summary$s.table),
                         method = rep("continuous",length.out = 2),
                         metric = rep(paste(i),length.out = 2))
  
  res$cont$mod <- gamm.cont
  res$cont$summary <- cont.out
  res$cont$AIC <- gamm.cont$gcv.ubre
  res$cont$method <- "continuous"
  res$cont$threshold <- NA
  best$cont <-gamm.cont$gcv.ubre
  if(res$cont$summary$Pr...t..[2] < 0.1 && res$cont$summary$Estimate[2] >1.00){
    
    thresh.values<-round(seq(stats::quantile(state.wind.dat$date[],0.2),
                             stats::quantile(state.wind.dat$date[],0.8),by=1)) #round to nearest whole number
    aic.thresh<-thresh.values*NA
    
    for(jj in 1:length(thresh.values)){
      sub.dat$tmp.thresh <- factor(ifelse(sub.dat$date<thresh.values[jj],'pre','post'), levels = c("pre","post"))
      gam.test<- mgcv::gam(metric~ tmp.thresh + s(PC, bs="tp", k=5,by=tmp.thresh)-1, data = sub.dat, 
                           family = gaussian(), method = "REML")
      aic.thresh[jj] <- gam.test$gcv.ubre
    }
    sub.dat$threshold<-factor(ifelse(sub.dat$date<thresh.values[order(aic.thresh)][1],'pre','post'),levels = c("pre","post"))
    
    gamm.thresh<- mgcv::gam(metric~ threshold + s(PC, bs="tp", k=5,by=as.factor(threshold))-1, 
                            data = sub.dat, 
                            family = gaussian(), method = "REML")  
    thresh.summary <- summary(gamm.thresh)
    thresh.out <- data.frame(rbind(thresh.summary$p.table,thresh.summary$s.table),
                             method = rep("threshold",length.out = 4),
                             FD = rep(paste(i),length.out = 4))
    
    res$thresh$mod <- gamm.thresh
    res$thresh$summary <- thresh.out
    res$thresh$AIC <- gamm.thresh$gcv.ubre
    res$thresh$method <- "threshold"
    res$thresh$threshold <- thresh.values[base::order(aic.thresh)][1]
    res$thresh$threshold.gcv <-cbind("thresh.date" = thresh.values[base::order(aic.thresh)],"gcv" = aic.thresh[base::order(aic.thresh)])
    best$thresh <-gamm.thresh$gcv.ubre
  }else{res$thresh$mod <- "NS"
  res$thresh$AIC <- Inf
  }
  return(res)
}|>
  `names<-`(c("density","community"))

wind.state.gam.yr.pred <- as.data.frame(foreach(i = c("density","community"),.combine="cbind",.packages=c('tidyverse',"magrittr")) %do%{
  if(abs(wind_state_gam_yr[[i]]$cont$AIC-wind_state_gam_yr[[i]]$thresh$AIC) <=6){
    best <- 1
  }else{
    best <- which.min(c(wind_state_gam_yr[[i]]$cont$AIC,wind_state_gam_yr[[i]]$thresh$AIC))
  }
  best.mod <- wind_state_gam_yr[[i]][[best]] #subset best model (continuous vs threshold based upon AIC)
  sub.dat <- data.frame(date=state.wind.dat$date,PC = state.wind.dat$PC1,metric=state.wind.dat[,i]) %>%
    mutate(threshold = ifelse(date>= best.mod$threshold,'post','pre')) %>%
    mutate(transition = ifelse(threshold == "pre" & lead(threshold) == "post",'trans',
                               ifelse(threshold == "post" & lag(threshold) == "pre",'trans',      
                                      'no.trans')))
  pred <- predict(best.mod$mod, newdata = sub.dat, type = 'link', se.fit = TRUE,
                  exclude = c("s(month)"))
  
  tmp <- cbind(sub.dat,pred) %>%
    set_colnames(c(paste(i,"date",sep = "_"),paste(i,"PC",sep = "_"),paste(i,"value",sep = "_"),
                   paste(i,"threshold",sep = "_"),paste(i,"transition",sep = "_"),paste(i,"fit",sep = "_"),
                   paste(i,"CI",sep = "_")))
  return(tmp)
}) %>% #drop accum column caused by .init specification
  pivot_longer(everything(), 
               names_to = c("metric",".value"), 
               names_pattern = "(.*)_(.*)") %>%
  mutate(threshold = replace_na(as.character(threshold),"pre"),
         lake = "Windermere") #replace NAs in threshold variable that appear for continuous GAMs

wind.map <- ggplot(data =wind.state.gam.yr.pred,aes(x=PC,y=value)) + 
  geom_point(aes(x=PC, y = value))+
  #geom_path(aes(x=PC, y = value)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.wind, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + facet_wrap(~metric,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  geom_line(data = filter(wind.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(wind.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(wind.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(wind.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(wind.state.gam.yr.pred,transition == "trans"), aes(x=PC, y = value,label=date),force =2,nudge_x=2,segment.linetype=2)+
  xlim(-4,4)

################################################################################################################
## Monona State Space ##
################################################################################################################
source("/Users/ul20791/Desktop/Academia/PhD/Data/Madison/Data/monona_plankton_data.R")

state.mon.dat <- data.frame("density" =scale(log1p(rowSums(plank_env.monyrdata[,2:202]))),
                            "community" = prcomp(scale(plank_env.monyrdata[,c(2:202)]))$x[,1],
                            "PC1" = prcomp(scale(plank_env.monyrdata[,c("wtemp","totP","NO3N","water.lvl")]))$x[,1],
                            #"PC1" = prcomp(scale(plank_env.monyrdata[,c("wtemp","totP","NO3N")]))$x[,1],
                            "date" = as.numeric(plank_env.monyrdata$date)) #unable to use first 4 years due to missing P_ort data
mid.mon<-median(state.mon.dat$date)

mon_state_gam_yr <-foreach(i = c("density","community")) %do%{
  res <- list() #result out
  best <- list() #assess best fitting model
  thresh.values<-c()
  
  sub.dat <- data.frame(date=state.mon.dat$date,PC = state.mon.dat$PC1,metric=state.mon.dat[,paste(i)])
  
  gamm.cont<- mgcv::gam(metric~ s(PC, bs="tp", k=10) , data = sub.dat,
                        family = gaussian(), method = "REML")
  cont.summary <- summary(gamm.cont)
  cont.out <- data.frame(rbind(cont.summary$p.table,cont.summary$s.table),
                         method = rep("continuous",length.out = 2),
                         metric = rep(paste(i),length.out = 2))
  
  res$cont$mod <- gamm.cont
  res$cont$summary <- cont.out
  res$cont$AIC <- gamm.cont$gcv.ubre
  res$cont$method <- "continuous"
  res$cont$threshold <- NA
  best$cont <-gamm.cont$gcv.ubre
  if(res$cont$summary$Pr...t..[2] < 0.1 && res$cont$summary$Estimate[2] >1.00){
    
    thresh.values<-round(seq(stats::quantile(state.mon.dat$date[],0.2),
                             stats::quantile(state.mon.dat$date[],0.8),by=1)) #round to nearest whole number
    aic.thresh<-thresh.values*NA
    
    for(jj in 1:length(thresh.values)){
      sub.dat$tmp.thresh <- factor(ifelse(sub.dat$date<thresh.values[jj],'pre','post'), levels = c("pre","post"))
      gam.test<- mgcv::gam(metric~ tmp.thresh + s(PC, bs="tp", k=5,by=tmp.thresh)-1, data = sub.dat, 
                           family = gaussian(), method = "REML")
      aic.thresh[jj] <- gam.test$gcv.ubre
    }
    sub.dat$threshold<-factor(ifelse(sub.dat$date<thresh.values[order(aic.thresh)][1],'pre','post'),levels = c("pre","post"))
    
    gamm.thresh<- mgcv::gam(metric~ threshold + s(PC, bs="tp", k=5,by=as.factor(threshold))-1, 
                            data = sub.dat, 
                            family = gaussian(), method = "REML")  
    thresh.summary <- summary(gamm.thresh)
    thresh.out <- data.frame(rbind(thresh.summary$p.table,thresh.summary$s.table),
                             method = rep("threshold",length.out = 4),
                             FD = rep(paste(i),length.out = 4))
    
    res$thresh$mod <- gamm.thresh
    res$thresh$summary <- thresh.out
    res$thresh$AIC <- gamm.thresh$gcv.ubre
    res$thresh$method <- "threshold"
    res$thresh$threshold <- thresh.values[base::order(aic.thresh)][1]
    res$thresh$threshold.gcv <-cbind("thresh.date" = thresh.values[base::order(aic.thresh)],"gcv" = aic.thresh[base::order(aic.thresh)])
    best$thresh <-gamm.thresh$gcv.ubre
  }else{res$thresh$mod <- "NS"
  res$thresh$AIC <- Inf
  }
  return(res)
} |>
  `names<-`(c("density","community"))

mon.state.gam.yr.pred <- as.data.frame(foreach(i = c("density","community"),.combine="cbind",.packages=c('tidyverse',"magrittr")) %do%{
  if(abs(mon_state_gam_yr[[i]]$cont$AIC-mon_state_gam_yr[[i]]$thresh$AIC) <=6){
    best <- 1
  }else{
    best <- which.min(c(mon_state_gam_yr[[i]]$cont$AIC,mon_state_gam_yr[[i]]$thresh$AIC))
  }
  best.mod <- mon_state_gam_yr[[i]][[best]] #subset best model (continuous vs threshold based upon AIC)
  sub.dat <- data.frame(date=state.mon.dat$date,PC = state.mon.dat$PC1,metric=state.mon.dat[,i]) %>%
    mutate(threshold = ifelse(date>= best.mod$threshold,'post','pre')) %>%
    mutate(transition = ifelse(threshold == "pre" & lead(threshold) == "post",'trans',
                               ifelse(threshold == "post" & lag(threshold) == "pre",'trans',      
                                      'no.trans')))
  pred <- predict(best.mod$mod, newdata = sub.dat, type = 'link', se.fit = TRUE,
                  exclude = c("s(month)"))
  
  tmp <- cbind(sub.dat,pred) %>%
    set_colnames(c(paste(i,"date",sep = "_"),paste(i,"PC",sep = "_"),paste(i,"value",sep = "_"),
                   paste(i,"threshold",sep = "_"),paste(i,"transition",sep = "_"),paste(i,"fit",sep = "_"),
                   paste(i,"CI",sep = "_")))
  return(tmp)
}) %>% #drop accum column caused by .init specification
  pivot_longer(everything(), 
               names_to = c("metric",".value"), 
               names_pattern = "(.*)_(.*)") %>%
  mutate(threshold = replace_na(as.character(threshold),"pre"),
         lake = "Monona") #replace NAs in threshold variable that appear for continuous GAMs

mon.map <- ggplot(data =mon.state.gam.yr.pred,aes(x=PC,y=value)) + 
  geom_point(aes(x=PC, y = value))+
  #geom_path(aes(x=PC, y = value)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.mon, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + facet_wrap(~metric,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  geom_line(data = filter(mon.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(mon.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(mon.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(mon.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(mon.state.gam.yr.pred,transition == "trans"), aes(x=PC, y = value,label=date),force =2,nudge_x=2,segment.linetype=2)+
  xlim(-4,4)

################################################################################################################
## Leven State Space ##
################################################################################################################
source("/Users/ul20791/Desktop/Academia/PhD/Data/Leven/Data/leven_plankton_data.R")

state.leve.dat <- data.frame("density" =scale(log1p(rowSums(plank_env.leveyrdata[,2:10]))),
                            "community" = prcomp(scale(plank_env.leveyrdata[,c(2:10)]))$x[,1],
                            "PC1" = prcomp(scale(plank_env.leveyrdata[,c("Temp","TOTP")]))$x[,1],
                            #"PC1" = prcomp(scale(plank_env.leveyrdata[,c("wtemp","totP","NO3N")]))$x[,1],
                            "date" = as.numeric(plank_env.leveyrdata$Date)) #unable to use first 4 years due to missing P_ort data
mid.leve<-median(state.leve.dat$date)

leve_state_gam_yr <-foreach(i = c("density","community")) %do%{
  res <- list() #result out
  best <- list() #assess best fitting model
  thresh.values<-c()
  
  sub.dat <- data.frame(date=state.leve.dat$date,PC = state.leve.dat$PC1,metric=state.leve.dat[,paste(i)])
  
  gamm.cont<- mgcv::gam(metric~ s(PC, bs="tp", k=10) , data = sub.dat,
                        family = gaussian(), method = "REML")
  cont.summary <- summary(gamm.cont)
  cont.out <- data.frame(rbind(cont.summary$p.table,cont.summary$s.table),
                         method = rep("continuous",length.out = 2),
                         metric = rep(paste(i),length.out = 2))
  
  res$cont$mod <- gamm.cont
  res$cont$summary <- cont.out
  res$cont$AIC <- gamm.cont$gcv.ubre
  res$cont$method <- "continuous"
  res$cont$threshold <- NA
  best$cont <-gamm.cont$gcv.ubre
  if(res$cont$summary$Pr...t..[2] < 0.1 && res$cont$summary$Estimate[2] >1.00){
    
    thresh.values<-round(seq(stats::quantile(state.leve.dat$date[],0.2),
                             stats::quantile(state.leve.dat$date[],0.8),by=1)) #round to nearest whole number
    aic.thresh<-thresh.values*NA
    
    for(jj in 1:length(thresh.values)){
      sub.dat$tmp.thresh <- factor(ifelse(sub.dat$date<thresh.values[jj],'pre','post'), levels = c("pre","post"))
      gam.test<- mgcv::gam(metric~ tmp.thresh + s(PC, bs="tp", k=5,by=tmp.thresh)-1, data = sub.dat, 
                           family = gaussian(), method = "REML")
      aic.thresh[jj] <- gam.test$gcv.ubre
    }
    sub.dat$threshold<-factor(ifelse(sub.dat$date<thresh.values[order(aic.thresh)][1],'pre','post'),levels = c("pre","post"))
    
    gamm.thresh<- mgcv::gam(metric~ threshold + s(PC, bs="tp", k=5,by=as.factor(threshold))-1, 
                            data = sub.dat, 
                            family = gaussian(), method = "REML")  
    thresh.summary <- summary(gamm.thresh)
    thresh.out <- data.frame(rbind(thresh.summary$p.table,thresh.summary$s.table),
                             method = rep("threshold",length.out = 4),
                             FD = rep(paste(i),length.out = 4))
    
    res$thresh$mod <- gamm.thresh
    res$thresh$summary <- thresh.out
    res$thresh$AIC <- gamm.thresh$gcv.ubre
    res$thresh$method <- "threshold"
    res$thresh$threshold <- thresh.values[base::order(aic.thresh)][1]
    res$thresh$threshold.gcv <-cbind("thresh.date" = thresh.values[base::order(aic.thresh)],"gcv" = aic.thresh[base::order(aic.thresh)])
    best$thresh <-gamm.thresh$gcv.ubre
  }else{res$thresh$mod <- "NS"
  res$thresh$AIC <- Inf
  }
  return(res)
} |>
  `names<-`(c("density","community"))

leve.state.gam.yr.pred <- as.data.frame(foreach(i = c("density","community"),.combine="cbind",.packages=c('tidyverse',"magrittr")) %do%{
  if(abs(leve_state_gam_yr[[i]]$cont$AIC-leve_state_gam_yr[[i]]$thresh$AIC) <=6){
    best <- 1
  }else{
    best <- which.min(c(leve_state_gam_yr[[i]]$cont$AIC,leve_state_gam_yr[[i]]$thresh$AIC))
  }
  best.mod <- leve_state_gam_yr[[i]][[best]] #subset best model (continuous vs threshold based upon AIC)
  sub.dat <- data.frame(date=state.leve.dat$date,PC = state.leve.dat$PC1,metric=state.leve.dat[,i]) %>%
    mutate(threshold = ifelse(date>= best.mod$threshold,'post','pre')) %>%
    mutate(transition = ifelse(threshold == "pre" & lead(threshold) == "post",'trans',
                               ifelse(threshold == "post" & lag(threshold) == "pre",'trans',      
                                      'no.trans')))
  pred <- predict(best.mod$mod, newdata = sub.dat, type = 'link', se.fit = TRUE,
                  exclude = c("s(leveth)"))
  
  tmp <- cbind(sub.dat,pred) %>%
    set_colnames(c(paste(i,"date",sep = "_"),paste(i,"PC",sep = "_"),paste(i,"value",sep = "_"),
                   paste(i,"threshold",sep = "_"),paste(i,"transition",sep = "_"),paste(i,"fit",sep = "_"),
                   paste(i,"CI",sep = "_")))
  return(tmp)
}) %>% #drop accum column caused by .init specification
  pivot_longer(everything(), 
               names_to = c("metric",".value"), 
               names_pattern = "(.*)_(.*)") %>%
  mutate(threshold = replace_na(as.character(threshold),"pre"),
         lake = "Loch Leven") #replace NAs in threshold variable that appear for continuous GAMs

leve.map <- ggplot(data =leve.state.gam.yr.pred,aes(x=PC,y=value)) + 
  geom_point(aes(x=PC, y = value))+
  #geom_path(aes(x=PC, y = value)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.leve, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + facet_wrap(~metric,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  geom_line(data = filter(leve.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(leve.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(leve.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(leve.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(leve.state.gam.yr.pred,transition == "trans"), aes(x=PC, y = value,label=date),force =2,nudge_x=2,segment.linetype=2)+
  xlim(-4,4)

################################################################################################################
## Washington State Space ##
################################################################################################################

source("/Users/ul20791/Desktop/Academia/PhD/Data/Washington/Data/washington_plankton_data.R")
#phyo, zoo + combined raw concentrations at week, month and year groupings

state.wash.dat <- data.frame("density" =scale(log1p(rowSums(plank_env.washyrdata[,4:18]))),
                            "community" = prcomp(scale(plank_env.washyrdata[,4:18]))$x[,1],
                            "PC1" = prcomp(scale(plank_env.washyrdata[,c("temp","TP")]))$x[,1],
                            "date" = as.numeric(plank_env.washyrdata$date)) #unable to use first 4 years due to missing P_ort data
mid.wash<-median(state.wash.dat$date)

wash_state_gam_yr <-foreach(i = c("density","community")) %do%{
  res <- list() #result out
  best <- list() #assess best fitting model
  thresh.values<-c()
  
  sub.dat <- data.frame(date=state.wash.dat$date,PC = state.wash.dat$PC1,metric=state.wash.dat[,paste(i)])
  
  gamm.cont<- mgcv::gam(metric~ s(PC, bs="tp", k=10) , data = sub.dat,
                        family = gaussian(), method = "REML")
  cont.summary <- summary(gamm.cont)
  cont.out <- data.frame(rbind(cont.summary$p.table,cont.summary$s.table),
                         method = rep("continuous",length.out = 2),
                         metric = rep(paste(i),length.out = 2))
  
  res$cont$mod <- gamm.cont
  res$cont$summary <- cont.out
  res$cont$AIC <- gamm.cont$gcv.ubre
  res$cont$method <- "continuous"
  res$cont$threshold <- NA
  best$cont <-gamm.cont$gcv.ubre
  if(res$cont$summary$Pr...t..[2] < 0.1 && res$cont$summary$Estimate[2] >1.00){
    
    thresh.values<-round(seq(stats::quantile(state.wash.dat$date[],0.2),
                             stats::quantile(state.wash.dat$date[],0.8),by=1)) #round to nearest whole number
    aic.thresh<-thresh.values*NA
    
    for(jj in 1:length(thresh.values)){
      sub.dat$tmp.thresh <- factor(ifelse(sub.dat$date<thresh.values[jj],'pre','post'), levels = c("pre","post"))
      gam.test<- mgcv::gam(metric~ tmp.thresh + s(PC, bs="tp", k=5,by=tmp.thresh)-1, data = sub.dat, 
                           family = gaussian(), method = "REML")
      aic.thresh[jj] <- gam.test$gcv.ubre
    }
    sub.dat$threshold<-factor(ifelse(sub.dat$date<thresh.values[order(aic.thresh)][1],'pre','post'),levels = c("pre","post"))
    
    gamm.thresh<- mgcv::gam(metric~ threshold + s(PC, bs="tp", k=5,by=as.factor(threshold))-1, 
                            data = sub.dat, 
                            family = gaussian(), method = "REML")  
    thresh.summary <- summary(gamm.thresh)
    thresh.out <- data.frame(rbind(thresh.summary$p.table,thresh.summary$s.table),
                             method = rep("threshold",length.out = 4),
                             FD = rep(paste(i),length.out = 4))
    
    res$thresh$mod <- gamm.thresh
    res$thresh$summary <- thresh.out
    res$thresh$AIC <- gamm.thresh$gcv.ubre
    res$thresh$method <- "threshold"
    res$thresh$threshold <- thresh.values[base::order(aic.thresh)][1]
    res$thresh$threshold.gcv <-cbind("thresh.date" = thresh.values[base::order(aic.thresh)],"gcv" = aic.thresh[base::order(aic.thresh)])
    best$thresh <-gamm.thresh$gcv.ubre
  }else{res$thresh$mod <- "NS"
  res$thresh$AIC <- Inf
  }
  return(res)
}|>
  `names<-`(c("density","community"))

wash.state.gam.yr.pred <- as.data.frame(foreach(i = c("density","community"),.combine="cbind",.packages=c('tidyverse',"magrittr")) %do%{
  if(abs(wash_state_gam_yr[[i]]$cont$AIC-wash_state_gam_yr[[i]]$thresh$AIC) <=6){
    best <- 1
  }else{
    best <- which.min(c(wash_state_gam_yr[[i]]$cont$AIC,wash_state_gam_yr[[i]]$thresh$AIC))
  }
  best.mod <- wash_state_gam_yr[[i]][[best]] #subset best model (continuous vs threshold based upon AIC)
  sub.dat <- data.frame(date=state.wash.dat$date,PC = state.wash.dat$PC1,metric=state.wash.dat[,i]) %>%
    mutate(threshold = ifelse(date>= best.mod$threshold,'post','pre')) %>%
    mutate(transition = ifelse(threshold == "pre" & lead(threshold) == "post",'trans',
                               ifelse(threshold == "post" & lag(threshold) == "pre",'trans',      
                                      'no.trans')))
  pred <- predict(best.mod$mod, newdata = sub.dat, type = 'link', se.fit = TRUE,
                  exclude = c("s(month)"))
  
  tmp <- cbind(sub.dat,pred) %>%
    set_colnames(c(paste(i,"date",sep = "_"),paste(i,"PC",sep = "_"),paste(i,"value",sep = "_"),
                   paste(i,"threshold",sep = "_"),paste(i,"transition",sep = "_"),paste(i,"fit",sep = "_"),
                   paste(i,"CI",sep = "_")))
  return(tmp)
}) %>% #drop accum column caused by .init specification
  pivot_longer(everything(), 
               names_to = c("metric",".value"), 
               names_pattern = "(.*)_(.*)") %>%
  mutate(threshold = replace_na(as.character(threshold),"pre"),
         lake = "Washington") #replace NAs in threshold variable that appear for continuous GAMs

wash.map <- ggplot(data =wash.state.gam.yr.pred,aes(x=PC,y=value)) + 
  geom_point(aes(x=PC, y = value))+
  geom_path(aes(x=PC, y = value)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.wash, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + facet_wrap(~metric,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  geom_line(data = filter(wash.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(wash.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(wash.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(wash.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(wash.state.gam.yr.pred,transition == "trans"), aes(x=PC, y = value,label=date),force =2,nudge_x=2,segment.linetype=2)+
  xlim(-4,4)

################################################################################################################
## Combine Lakes + Plot ##
################################################################################################################
all.lakes.state.gam.yr.pred <- rbind(kin.state.gam.yr.pred,kas.state.gam.yr.pred,mad.state.gam.yr.pred,LZ.state.gam.yr.pred,wind.state.gam.yr.pred,UZ.state.gam.yr.pred,mon.state.gam.yr.pred,leve.state.gam.yr.pred,wash.state.gam.yr.pred)

write.csv(all.lakes.state.gam.yr.pred,file =  "Data/all.lakes.state.gam.yr.pred.csv")
all.lakes.state.gam.yr.pred <- read.csv(file =  "Data/all.lakes.state.gam.yr.pred.csv")[,-1]

all.lakes.state.gam.yr.pred <- all.lakes.state.gam.yr.pred |>
  group_by(metric,lake)|>
  mutate(start.date = first(date),
         last.date = last(date))
  
ggsave(ggplot(data =all.lakes.state.gam.yr.pred,aes(x=PC,y=value)) + 
  geom_point(aes(x=PC, y = value))+
  #geom_path(aes(x=PC, y = value)) +  
  #scale_color_gradient2(aes(x=PC, y = value,col = date), midpoint=mid.wind, low="blue", mid="grey",
  #                     high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + 
  facet_grid(metric~lake,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  geom_line(data = filter(all.lakes.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(all.lakes.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_line(data = filter(all.lakes.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(all.lakes.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  geom_point(data = filter(all.lakes.state.gam.yr.pred, date == start.date),aes(x=PC, y = value,col="Start date"))+
  geom_point(data = filter(all.lakes.state.gam.yr.pred, date == last.date),aes(x=PC, y = value,col="End date"))+
  ggrepel::geom_text_repel(data = filter(all.lakes.state.gam.yr.pred,transition == "trans"), aes(x=PC, y = value,label=date),force =1,nudge_x=3,nudge_y=-0.5,segment.linetype=2)+
  geom_point(data = filter(all.lakes.state.gam.yr.pred,transition == "trans"),aes(x=PC, y = value,col="Transition\ndates"))+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates"), name = NULL)+
  xlim(-6,6),
  filename = "Results/lake_results/state_space_fig.pdf", width = 10,height = 6, dpi= 144)

ts_state_plot_dat <- all.lakes.state.gam.yr.pred |>
  group_by(lake,metric) |>
  mutate(scaled_date = scales::rescale(date, to=c(0,10)))

ggsave(ggplot(data =ts_state_plot_dat,aes(x=PC,y=value,col=scaled_date)) + 
  geom_point(aes(x=PC, y = value))+
  geom_path(aes(x=PC, y = value)) +  
  scale_color_gradient2(aes(x=PC, y = value,col = as.numeric(scaled_date)), midpoint=5, low="blue", mid="grey",
                       high="red",guide = "colourbar")+
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  #guides(colourbar=guide_legend(title="New Legend Title"))+
  theme_bw() + 
  facet_grid(metric~lake,scales = "free_y",labeller = labeller(metric = c(community = "Community",density = "Density"))) +
  #geom_line(data = filter(all.lakes.state.gam.yr.pred, threshold=="pre"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  #geom_ribbon(data = filter(all.lakes.state.gam.yr.pred, threshold=="pre"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  #geom_line(data = filter(all.lakes.state.gam.yr.pred, threshold=="post"),aes(x=PC,y=fit), col="black",size=0.8, linetype = "solid")+
  #geom_ribbon(data = filter(all.lakes.state.gam.yr.pred, threshold=="post"),aes(ymin = fit - (2 * CI),ymax = fit + (2 * CI)  ), fill = "grey", col="grey",alpha = 0.2)+
  #geom_point(data = filter(all.lakes.state.gam.yr.pred, date == start.date),aes(x=PC, y = value,col="Start date"))+
  #geom_point(data = filter(all.lakes.state.gam.yr.pred, date == last.date),aes(x=PC, y = value,col="End date"))+
  ggrepel::geom_text_repel(data = filter(ts_state_plot_dat,date %in% paste(c(start.date,last.date))),
                           aes(x=PC, y = value,label=date),force =1,nudge_x=3,nudge_y=-0.5,segment.linetype=2,col="black")+
  #geom_point(data = filter(all.lakes.state.gam.yr.pred,transition == "trans"),aes(x=PC, y = value,col="Transition\ndates"))+
  #scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates"), name = NULL)+
  xlim(-6,6),
  filename = "Results/lake_results/lotsa_lakes_ts_state_space_fig.pdf", width = 12,height = 8, dpi= 144)

  