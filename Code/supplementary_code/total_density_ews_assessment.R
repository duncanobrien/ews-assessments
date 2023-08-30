##########################################################################################
## Preamble ##
##########################################################################################

require(dplyr)
require(EWSmethods)
require(foreach)
require(zoo)
require(stringr)

load("Data/wrangled_genus_plank_data.Rdata") #data wrangled to genus level and trimmed prior to TGAM estimated transition dates  

lapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
            wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
            kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
            wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat),
       function(x){
         
         if(any(stringr::str_length(x[,1])>4)){ #if monthly data, drop timeseries that are 0 for longer than a year
           x <- x |>
             select(-c(pca1,pca2))|>
             select(where(function(x){all(rle(x)$lengths[rle(x)$values == 0] < 12)}))
         }else if(all(stringr::str_length(x[,1])==4)){
           x <- x |>
             select(-c(pca1,pca2))|>
             select(where(function(x){all(rle(x)$lengths[rle(x)$values == 0] < 2)})) #if yearly data, drop timeseries that are 0 for longer than a year
         }
         return(x)
       }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat")) |>
  list2env(env = .GlobalEnv) #replace global env objects with wrangled

rm(list = ls()[!ls() %in% c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
                            "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat")])

################################################################################################################
## Univariate EWS Assessment Expanding ##
################################################################################################################

plank_ls <- list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
     wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
     kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
     wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))

exp_uni_tot <- pbapply::pblapply(seq_along(plank_ls),
                                     FUN = function(x){
                                       
                                       sub_dat <- data.frame("time" = plank_ls[[x]][,1]) |>
                                         dplyr::mutate(tot.density =  rowSums(plank_ls[[x]][,-1])) |>
                                         dplyr::mutate(phyto.density = switch(names(plank_ls)[x],
                                                                                "kin_yr_dat" = rowSums(plank_ls[[x]][,2:14]),
                                                                                "kin_mth_dat" = rowSums(plank_ls[[x]][,2:14]),
                                                                                "kas_yr_dat" = rowSums(plank_ls[[x]][,2:25]),
                                                                                "kas_mth_dat" = rowSums(plank_ls[[x]][,2:25]),
                                                                                "LZ_yr_dat" = rowSums(plank_ls[[x]][,2:45]),
                                                                                "LZ_mth_dat" = rowSums(plank_ls[[x]][,2:36]),
                                                                                "UZ_yr_dat" = rowSums(plank_ls[[x]][,2:59]),
                                                                                "UZ_mth_dat" = rowSums(plank_ls[[x]][,2:58]),
                                                                                "mad_yr_dat" = rowSums(plank_ls[[x]][,2:32]),
                                                                                "mad_mth_dat" = rowSums(plank_ls[[x]][,2:17]),
                                                                                "wind_yr_dat" = rowSums(plank_ls[[x]][,2:18]),
                                                                                "wind_mth_dat" = rowSums(plank_ls[[x]][,2:12]),
                                                                                "wash_yr_dat" = rowSums(plank_ls[[x]][,2:7]),
                                                                                "wash_mth_dat" = rowSums(plank_ls[[x]][,2:7]),
                                                                                "leve_yr_dat" = rowSums(plank_ls[[x]][,2:4]),
                                                                                "leve_mth_dat" = rowSums(plank_ls[[x]][,2:4]),
                                                                                "mon_yr_dat" = rowSums(plank_ls[[x]][,2:45]),
                                                                                "mon_mth_dat" = rowSums(plank_ls[[x]][,2:29])
                                         ),
                                         zoo.density = switch(names(plank_ls)[x],
                                                                "kin_yr_dat" = rowSums(plank_ls[[x]][,15:32]),
                                                                "kin_mth_dat" = rowSums(plank_ls[[x]][,15:22]),
                                                                "kas_yr_dat" = rowSums(plank_ls[[x]][,26:34]),
                                                                "kas_mth_dat" = rowSums(plank_ls[[x]][,26:34]),
                                                                "LZ_yr_dat" = rowSums(plank_ls[[x]][,46:51]),
                                                                "LZ_mth_dat" = rowSums(plank_ls[[x]][,37:41]),
                                                                "UZ_yr_dat" = rowSums(plank_ls[[x]][,60:65]),
                                                                "UZ_mth_dat" = rowSums(plank_ls[[x]][,59:63]),
                                                                "mad_yr_dat" = rowSums(plank_ls[[x]][,33:38]),
                                                                "mad_mth_dat" = rowSums(plank_ls[[x]][,18:22]),
                                                                "wind_yr_dat" = rowSums(plank_ls[[x]][,20:23]),
                                                                "wind_mth_dat" = rowSums(plank_ls[[x]][,14:17]),
                                                                "wash_yr_dat" = rowSums(plank_ls[[x]][,8:15]),
                                                                "wash_mth_dat" = rowSums(plank_ls[[x]][,8:13]),
                                                                "leve_yr_dat" = rowSums(plank_ls[[x]][,5:7]),
                                                                "leve_mth_dat" = rowSums(plank_ls[[x]][,5:6]),
                                                                "mon_yr_dat" = rowSums(plank_ls[[x]][,46:53]),
                                                                "mon_mth_dat" = rowSums(plank_ls[[x]][,30:36])
                                         )
                                         )

                                       out <- lapply(c("none","linear","gaussian","loess"),function(j){ #detrend 
                                         if(j != "none"){
                                           sub_dat <- EWSmethods::detrend_ts(sub_dat,method = j,span=0.5)   
                                           sub_dat[,-1] <- sapply(colnames(sub_dat)[-1],FUN = function(i){
                                             sub_dat[,i] <- sub_dat[,i] + abs(min(sub_dat[,i])) #set 0 as minimum
                                           })
                                         }
                                         
                                         out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                           
                                           print(paste(names(plank_ls)[x],j,k,sep="_"))
                                           if(k != "none"){
                                             
                                             if(any(stringr::str_length(sub_dat[,1])>4)){ #deseason if monthly data
                                               sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat[,1]))
                                               sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = k,order = "ymd")
                                               sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
                                               
                                               if(k == "decompose"){
                                                 sub_dat <- na.omit(sub_dat) #remove NAs at head and tail of df
                                               }
                                             }
                                           }
                                           
                                           time_vec <- colnames(sub_dat)[1]
                                           
                                           foreach(i=colnames(sub_dat[,-1]), .combine='rbind',.verbose = F) %do%{
                                             
                                             out <- EWSmethods::uniEWS(sub_dat[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",
                                                                       burn_in = ceiling(dim(sub_dat)[1]*0.5),
                                                                       threshold = 1)$EWS
                                             
                                             data.frame(cbind(data_source = paste(i),
                                                              out))
                                           }
                                           
                                          
                                           #return(data.frame(rbind(out.det$raw)))
                                           
                                         } ) |>
                                           `names<-`(c("none","average","decompose","stl")) |>
                                           data.table::rbindlist(idcol="deseason_meth")
                                         
                                         return(out_inner)
                                       }) |>
                                         `names<-`(c("none","linear","gaussian","loess")) |>
                                         data.table::rbindlist(idcol="detrend_meth")
                                       
                                       return(out)
                                     }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "expanding",
    troph_level = NA) |>
  select(-name)


################################################################################################################
## Permuted Univariate EWS Assessment Rolling ##
################################################################################################################

plank_ls <- list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
                 wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
                 kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
                 wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))

perm_roll_uni_tot <- lapply(seq_along(plank_ls),
                              FUN = function(x){
                                
                                sub_dat <- data.frame("time" = plank_ls[[x]][,1]) |>
                                  dplyr::mutate(tot.density =  rowSums(plank_ls[[x]][,-1])) |>
                                  dplyr::mutate(phyto.density = switch(names(plank_ls)[x],
                                                                       "kin_yr_dat" = rowSums(plank_ls[[x]][,2:14]),
                                                                       "kin_mth_dat" = rowSums(plank_ls[[x]][,2:14]),
                                                                       "kas_yr_dat" = rowSums(plank_ls[[x]][,2:25]),
                                                                       "kas_mth_dat" = rowSums(plank_ls[[x]][,2:25]),
                                                                       "LZ_yr_dat" = rowSums(plank_ls[[x]][,2:45]),
                                                                       "LZ_mth_dat" = rowSums(plank_ls[[x]][,2:36]),
                                                                       "UZ_yr_dat" = rowSums(plank_ls[[x]][,2:59]),
                                                                       "UZ_mth_dat" = rowSums(plank_ls[[x]][,2:58]),
                                                                       "mad_yr_dat" = rowSums(plank_ls[[x]][,2:32]),
                                                                       "mad_mth_dat" = rowSums(plank_ls[[x]][,2:17]),
                                                                       "wind_yr_dat" = rowSums(plank_ls[[x]][,2:18]),
                                                                       "wind_mth_dat" = rowSums(plank_ls[[x]][,2:12]),
                                                                       "wash_yr_dat" = rowSums(plank_ls[[x]][,2:7]),
                                                                       "wash_mth_dat" = rowSums(plank_ls[[x]][,2:7]),
                                                                       "leve_yr_dat" = rowSums(plank_ls[[x]][,2:4]),
                                                                       "leve_mth_dat" = rowSums(plank_ls[[x]][,2:4]),
                                                                       "mon_yr_dat" = rowSums(plank_ls[[x]][,2:45]),
                                                                       "mon_mth_dat" = rowSums(plank_ls[[x]][,2:29])
                                  ),
                                  zoo.density = switch(names(plank_ls)[x],
                                                       "kin_yr_dat" = rowSums(plank_ls[[x]][,15:32]),
                                                       "kin_mth_dat" = rowSums(plank_ls[[x]][,15:22]),
                                                       "kas_yr_dat" = rowSums(plank_ls[[x]][,26:34]),
                                                       "kas_mth_dat" = rowSums(plank_ls[[x]][,26:34]),
                                                       "LZ_yr_dat" = rowSums(plank_ls[[x]][,46:51]),
                                                       "LZ_mth_dat" = rowSums(plank_ls[[x]][,37:41]),
                                                       "UZ_yr_dat" = rowSums(plank_ls[[x]][,60:65]),
                                                       "UZ_mth_dat" = rowSums(plank_ls[[x]][,59:63]),
                                                       "mad_yr_dat" = rowSums(plank_ls[[x]][,33:38]),
                                                       "mad_mth_dat" = rowSums(plank_ls[[x]][,18:22]),
                                                       "wind_yr_dat" = rowSums(plank_ls[[x]][,20:23]),
                                                       "wind_mth_dat" = rowSums(plank_ls[[x]][,14:17]),
                                                       "wash_yr_dat" = rowSums(plank_ls[[x]][,8:15]),
                                                       "wash_mth_dat" = rowSums(plank_ls[[x]][,8:13]),
                                                       "leve_yr_dat" = rowSums(plank_ls[[x]][,5:7]),
                                                       "leve_mth_dat" = rowSums(plank_ls[[x]][,5:6]),
                                                       "mon_yr_dat" = rowSums(plank_ls[[x]][,46:53]),
                                                       "mon_mth_dat" = rowSums(plank_ls[[x]][,30:36])
                                  )
                                  )
                                
                                out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                  
                                  if(j != "none"){
                                    sub_dat <- EWSmethods::detrend_ts(sub_dat,method = j,span = 0.5) 
                                    sub_dat[,-1] <- sapply(colnames(sub_dat)[-1],FUN = function(i){
                                      sub_dat[,i] <- sub_dat[,i] + abs(min(sub_dat[,i])) #set 0 as minimum
                                    })
                                  }
                                  
                                  out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                    
                                    if(k != "none"){
                                      
                                      if(any(stringr::str_length(sub_dat[,1])>4)){
                                        sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat[,1]))
                                        sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = k,order = "ymd")
                                        sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
                                        
                                        if(k == "decompose"){
                                          sub_dat <- na.omit(sub_dat) #remove NAs at head and tail of df
                                        }
                                      }
                                    }
                                    
                                    time_vec <- colnames(sub_dat)[1]
                                    print(paste(names(plank_ls)[x],j,k,sep="_"))
                                    
                                    foreach(i=colnames(sub_dat[,-1]), .combine='rbind',.verbose = F) %do%{
                                      
                                      # out <-  perm_rollEWS(data = x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),
                                      #                      winsize = 50, perm.meth = "replacement", variate = "uni",
                                      #                      iter = 500,mc.cores=5)   
                                      
                                      out <-  EWSmethods::perm_rollEWS(data = sub_dat[,c(time_vec,i)],metrics = c("ar1","SD","skew"),
                                                                       winsize = 50, perm.meth = "sample", variate = "uni",
                                                                       iter = 500)$EWS   
                                      
                                      data.frame(cbind(data_source = paste(i),
                                                       out$cor))
                                    }
                                  }) |>
                                    `names<-`(c("none","average","decompose","stl")) |>
                                    data.table::rbindlist(idcol="deseason_meth")
                                  
                                  return(out_inner)
                                }) |>
                                  `names<-`(c("none","linear","gaussian","loess")) |>
                                  data.table::rbindlist(idcol="detrend_meth")
                                return(out)
                                
                              }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "rolling",
    troph_level = NA) |>
  select(-name)

################################################################################################################
## EWSNet ##
################################################################################################################

EWSmethods::ewsnet_init("EWSNET_env", auto=T)
EWSmethods::ewsnet_reset(auto=T)

plank_ls <- list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
                 wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
                 kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
                 wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))

ewsnet_tot <- pbapply::pblapply(seq_along(plank_ls),
                                  FUN = function(x){
                                    
                                    sub_dat <- data.frame("time" = plank_ls[[x]][,1]) |>
                                      dplyr::mutate(tot.density =  rowSums(plank_ls[[x]][,-1])) |>
                                      dplyr::mutate(phyto.density = switch(names(plank_ls)[x],
                                                                           "kin_yr_dat" = rowSums(plank_ls[[x]][,2:14]),
                                                                           "kin_mth_dat" = rowSums(plank_ls[[x]][,2:14]),
                                                                           "kas_yr_dat" = rowSums(plank_ls[[x]][,2:25]),
                                                                           "kas_mth_dat" = rowSums(plank_ls[[x]][,2:25]),
                                                                           "LZ_yr_dat" = rowSums(plank_ls[[x]][,2:45]),
                                                                           "LZ_mth_dat" = rowSums(plank_ls[[x]][,2:36]),
                                                                           "UZ_yr_dat" = rowSums(plank_ls[[x]][,2:59]),
                                                                           "UZ_mth_dat" = rowSums(plank_ls[[x]][,2:58]),
                                                                           "mad_yr_dat" = rowSums(plank_ls[[x]][,2:32]),
                                                                           "mad_mth_dat" = rowSums(plank_ls[[x]][,2:17]),
                                                                           "wind_yr_dat" = rowSums(plank_ls[[x]][,2:18]),
                                                                           "wind_mth_dat" = rowSums(plank_ls[[x]][,2:12]),
                                                                           "wash_yr_dat" = rowSums(plank_ls[[x]][,2:7]),
                                                                           "wash_mth_dat" = rowSums(plank_ls[[x]][,2:7]),
                                                                           "leve_yr_dat" = rowSums(plank_ls[[x]][,2:4]),
                                                                           "leve_mth_dat" = rowSums(plank_ls[[x]][,2:4]),
                                                                           "mon_yr_dat" = rowSums(plank_ls[[x]][,2:45]),
                                                                           "mon_mth_dat" = rowSums(plank_ls[[x]][,2:29])
                                      ),
                                      zoo.density = switch(names(plank_ls)[x],
                                                           "kin_yr_dat" = rowSums(plank_ls[[x]][,15:32]),
                                                           "kin_mth_dat" = rowSums(plank_ls[[x]][,15:22]),
                                                           "kas_yr_dat" = rowSums(plank_ls[[x]][,26:34]),
                                                           "kas_mth_dat" = rowSums(plank_ls[[x]][,26:34]),
                                                           "LZ_yr_dat" = rowSums(plank_ls[[x]][,46:51]),
                                                           "LZ_mth_dat" = rowSums(plank_ls[[x]][,37:41]),
                                                           "UZ_yr_dat" = rowSums(plank_ls[[x]][,60:65]),
                                                           "UZ_mth_dat" = rowSums(plank_ls[[x]][,59:63]),
                                                           "mad_yr_dat" = rowSums(plank_ls[[x]][,33:38]),
                                                           "mad_mth_dat" = rowSums(plank_ls[[x]][,18:22]),
                                                           "wind_yr_dat" = rowSums(plank_ls[[x]][,20:23]),
                                                           "wind_mth_dat" = rowSums(plank_ls[[x]][,14:17]),
                                                           "wash_yr_dat" = rowSums(plank_ls[[x]][,8:15]),
                                                           "wash_mth_dat" = rowSums(plank_ls[[x]][,8:13]),
                                                           "leve_yr_dat" = rowSums(plank_ls[[x]][,5:7]),
                                                           "leve_mth_dat" = rowSums(plank_ls[[x]][,5:6]),
                                                           "mon_yr_dat" = rowSums(plank_ls[[x]][,46:53]),
                                                           "mon_mth_dat" = rowSums(plank_ls[[x]][,30:36])
                                      )
                                      )
                                    
                                    out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                      
                                      if(j != "none"){
                                        sub_dat <- detrend_ts(sub_dat,method = j,span = 0.5) 
                                        sub_dat[,-1] <- sapply(colnames(sub_dat)[-1],FUN = function(i){
                                          sub_dat[,i] <- sub_dat[,i] + abs(min(sub_dat[,i])) #set 0 as minimum
                                        })
                                      }
                                      
                                      out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                        
                                        if(k != "none"){
                                          if(any(stringr::str_length(sub_dat[,1])>4)){
                                            sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat[,1]))
                                            sub_dat <- deseason_ts(sub_dat,increment = "month", method = k,order = "ymd")
                                            sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
                                            
                                            if(k == "decompose"){
                                              sub_dat <- na.omit(sub_dat) #remove NAs at head and tail of df
                                            }
                                          }
                                        }
                                        
                                        time_vec <- colnames(sub_dat)[1]
                                        
                                        foreach(i=colnames(sub_dat[,-1]), .combine='rbind',.verbose = F) %do%{
                                          
                                          print(paste(i,j,k,collapse = "_"))
                                          
                                          ews.tmp.scaled <- EWSmethods::ewsnet_predict(c(sub_dat[,i]),scaling = T, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                          ews.tmp.unscaled <- EWSmethods::ewsnet_predict(c(sub_dat[,i]),scaling = F, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                          
                                          data.frame(data_source = paste(i),
                                                     model_ensemble = "25",
                                                     predictionScaled = ews.tmp.scaled$pred,
                                                     no_transScaled = ews.tmp.scaled$no_trans_prob,
                                                     smth_transScaled =  ews.tmp.scaled$smooth_trans_prob,
                                                     crt_transScaled = ews.tmp.scaled$critical_trans_prob,
                                                     predictionUnscaled = ews.tmp.unscaled$pred,
                                                     no_transUnscaled = ews.tmp.unscaled$no_trans_prob,
                                                     smth_transUnscaled =  ews.tmp.unscaled$smooth_trans_prob,
                                                     crt_transUnscaled = ews.tmp.unscaled$critical_trans_prob)
                                          
                                        } 
                                      }) |>
                                        `names<-`(c("none","average","decompose","stl")) |>
                                        data.table::rbindlist(idcol="deseason_meth")
                                      
                                      return(out_inner)
                                    }) |>
                                      `names<-`(c("none","linear","gaussian","loess")) |>
                                      data.table::rbindlist(idcol="detrend_meth")
                                    return(out)
                                    
                                  }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "ML",
    troph_level = NA) |>
  select(-name)


################################################################################################################
## Combine ##
################################################################################################################

save(ewsnet_tot,exp_uni_tot,perm_roll_uni_tot, file = "Results/ews_tot_raw_data.RData")

################################################################################################################
## Fit models ##
################################################################################################################

require(Matrix)
require(tidybayes)
require(ggplot2)
require(data.table)
require(tidyverse)
require(brms)

load("Results/ews_tot_raw_data.RData")
load("Results/ews_tot_raw_datab.RData")

source("Code/extract_ews_pred_fn.R")

ewsnet_tot <- ewsnet_tot |> 
  mutate(troph_level = ifelse(grepl("phyto",data_source),"phytoplankton",
                              ifelse(grepl("zoo",data_source),"zooplankton",NA)))

exp_uni_tot <- exp_uni_tot |> 
  mutate(troph_level = ifelse(grepl("phyto",data_source),"phytoplankton",
                              ifelse(grepl("zoo",data_source),"zooplankton",NA)))

perm_roll_uni_tot <- perm_roll_uni_tot |> 
  mutate(troph_level = ifelse(grepl("phyto",data_source),"phytoplankton",
                              ifelse(grepl("zoo",data_source),"zooplankton",NA)))

lake_outcome_troph <- subset(ewsnet_tot, select = c(lake,res,troph_level)) |>
  dplyr::mutate(outcome = case_when(
    grepl("Kinneret", lake) & troph_level %in% c("phytoplankton",NA) ~ "trans",
    grepl("Kasumigaura", lake) & troph_level == c("zooplankton",NA) ~ "trans",
    grepl("Washington", lake) & troph_level == c("phytoplankton",NA) ~ "trans",
    grepl("Monona", lake) & troph_level == c("zooplankton",NA) ~ "trans",
    TRUE ~ "no.trans"
  )) |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_"))|>
  dplyr::select(reference,outcome)|>
  dplyr::distinct(reference,.keep_all = T)


suc_ewsnet_max_tot <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    extract_ews_pred(ews.data = subset(ewsnet_tot,detrend_meth == x & deseason_meth == j)|>
                       select(-c(detrend_meth,deseason_meth)),
                     sensitivity = NULL,
                     outcome = lake_outcome_troph,
                     method = "ML") |>
      mutate(detrend_meth = x,deseason_meth = j)
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()

suc_exp_uni_tot <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    out.mth <- extract_ews_pred(ews.data = subset(exp_uni_tot,detrend_meth == x & deseason_meth == j & res == "Monthly")|>
                                  select(-c(detrend_meth,deseason_meth)),
                                sensitivity = 2,
                                outcome = lake_outcome_troph,
                                method = "expanding") |>
      mutate(detrend_meth = x,deseason_meth = j)
    
    out.yr <- extract_ews_pred(ews.data = subset(exp_uni_tot,detrend_meth == x & deseason_meth == j & res == "Yearly")|>
                                 select(-c(detrend_meth,deseason_meth)),
                               sensitivity = 1,
                               outcome = lake_outcome_troph,
                               method = "expanding") |>
      mutate(detrend_meth = x,deseason_meth = j)
    
    return(rbind(out.mth,out.yr))
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()

suc_roll_perm_uni_tot <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    extract_ews_pred(ews.data = subset(perm_roll_uni_tot,detrend_meth == x & deseason_meth == j)|>
                       select(-c(detrend_meth,deseason_meth)),
                     sensitivity = 0.05,
                     outcome = lake_outcome_troph,
                     method = "rolling",
                     surrogate = TRUE) |>
      mutate(detrend_meth = x,deseason_meth = j)
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()


all_ews_data <- as.data.table(suc_ewsnet_max_tot) %>%
  .[,.SD[1], by = c("data_source","troph_level","scaling","lake","res","detrend_meth","deseason_meth")] %>%
  .[,c("data_source","troph_level","scaling","lake","res","method", "computation","prediction","detrend_meth","deseason_meth")] %>%
  setnames(old = "scaling",new= "metric.code") %>%
  rbind(as.data.table(suc_exp_uni_tot)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth")]) %>%
  rbind(as.data.table(suc_roll_perm_uni_tot)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth")]) %>%
  .[,prediction := ifelse(prediction %in% c("match","prior"),
                          1,0)] %>%
  .[,variate := ifelse(grepl("multivariate",method),"multivariate","univariate")]%>%
  .[,method_code := paste(variate,computation,sep = "_")] %>%
  .[,method_code := factor(method_code, levels = c("univariate_rolling","univariate_expanding","univariate_ML"))] %>%
  .[,computation := factor(computation, levels = c("rolling","expanding","ML"))] %>%
  .[,res := factor(res, levels = c("Yearly","Monthly"))]

ts_balance <- copy(all_ews_data) %>%
  .[,reference := paste(lake,res,troph_level,sep="_")] %>% #trials in terms of assessed time series
  merge(.,as.data.table(lake_outcome_troph),by="reference") %>%
  .[,data_ref := paste(data_source,lake,sep ="_")] %>%
  .[,.(data_ref = length(unique(data_ref))),by="outcome"] 

computation_data <- all_ews_data %>%
  .[,reference := paste(lake,res,troph_level,sep="_")] %>% #trials in terms of assessed time series
  merge(.,as.data.table(lake_outcome_troph),by="reference") %>%
  .[,offset := length(unique(data_source)), by =c("method_code","lake","res","troph_level","detrend_meth","deseason_meth") ] %>% #trials in terms of assessed time series
  .[,offset2 := length(unique(data_source))*length(unique(metric.code)), by =c("method_code","lake","res","troph_level","detrend_meth","deseason_meth")] %>% #trials in terms of assessed time series AND metrics 
  .[,.(total_success = sum(prediction),
       offset = unique(offset2),
       variate = unique(variate),
       reference = paste(lake,res,troph_level,sep="_"),
       outcome = unique(outcome),
       weights = unique(ifelse(outcome == "trans",(ts_balance$data_ref[2]),(ts_balance$data_ref[1])))
  ),
  by = c("lake", "res", "method_code","troph_level","detrend_meth","deseason_meth")] %>%
  .[,detrend_meth := factor(detrend_meth, levels = c("none","linear","loess","gaussian"))] %>%
  .[,deseason_meth := factor(deseason_meth, levels = c("none","average","decompose","stl"))] %>%
  as.data.frame()


overall_table <- lapply(unique(paste(all_ews_data$detrend_meth,all_ews_data$deseason_meth,
                                     sep = "_")), function(x){
                                       
                                       lapply(unique(all_ews_data$res),function(k){        
                                         gtsummary::tbl_strata(subset(all_ews_data,
                                                                      detrend_meth == gsub("_.*","",x) & 
                                                                        deseason_meth ==  gsub(".*_","",x) &
                                                                        res == k) |>
                                                                 dplyr::select(method_code,prediction) |>
                                                                 dplyr::mutate(prediction = paste(prediction)),
                                                               method_code,.tbl_fun = gtsummary::tbl_summary) |>
                                           as.data.frame() %>%
                                           rbind(colnames(.),.) |>
                                           t() |>
                                           as.data.frame() |>
                                           dplyr::slice(-1) |>
                                           `colnames<-`(c("Number of Trials","drop","Number of failures","Number of successes")) |>
                                           dplyr::select(-drop) |>
                                           dplyr::mutate(`Method code` = levels(all_ews_data$method_code),
                                                         Resolution = k, .before = `Number of Trials`) |>
                                           dplyr::mutate(`Detrend method` = gsub("_.*","",x), 
                                                         `Deseason method` = gsub(".*_","",x),.before = `Number of Trials`) |>
                                           dplyr::mutate(`Number of Trials` = gsub("\\D+","",`Number of Trials`))
                                       }) |> data.table::rbindlist()
                                     }) |> data.table::rbindlist() |>
  dplyr::arrange(Resolution,`Method code`)

metric_table <- lapply(unique(paste(all_ews_data$detrend_meth,all_ews_data$deseason_meth,
                                    sep = "_")), function(x){
                                      
                                      lapply(paste(unique(all_ews_data$method_code),rep(unique(all_ews_data$res),5),sep = ":"),function(k){
                                        
                                        sub_dat <- subset(all_ews_data,
                                                          detrend_meth == gsub("_.*","",x) & 
                                                            deseason_meth ==  gsub(".*_","",x) &
                                                            method_code == gsub(":.*","",k) &
                                                            res == gsub(".*:","",k)) |>
                                          dplyr::select(metric.code,prediction) |>
                                          dplyr::mutate(prediction = paste(prediction))
                                        
                                        gtsummary::tbl_strata(sub_dat,
                                                              metric.code,.tbl_fun =tbl_summary) |>
                                          as.data.frame() %>%
                                          rbind(colnames(.),.) |>
                                          t() |>
                                          as.data.frame() |>
                                          dplyr::slice(-1) |>
                                          `colnames<-`(c("Number of Trials","drop","Number of failures","Number of successes")) |>
                                          dplyr::select(-drop) |>
                                          dplyr::mutate(Indicator = unique(sub_dat$metric.code), .before = `Number of Trials`) |>
                                          dplyr::mutate(`Method code` = gsub(":.*","",k),
                                                        Resolution =  gsub(".*:","",k), 
                                                        .before = `Number of Trials`) |>
                                          dplyr::mutate(`Detrend method` = gsub("_.*","",x), 
                                                        `Deseason method` = gsub(".*_","",x),.before = `Number of Trials`) |>
                                          dplyr::mutate(`Number of Trials` = gsub("\\D+","",`Number of Trials`))
                                      }) |> data.table::rbindlist() |> dplyr::distinct()
                                    }
) |> data.table::rbindlist() |>
  dplyr::arrange(Resolution,Indicator)


metric_data <- all_ews_data %>%
  .[,offset := length(unique(data_source)), by =c("metric.code","lake","res","detrend_meth","deseason_meth","method_code") ] %>% #trials in terms of assessed time series
  .[,.(total_success = sum(prediction),
       offset = unique(offset),
       #ts_length = unique(ts_length),
       variate = unique(variate),
       reference = paste(lake,res,troph_level,sep="_")),by = c("lake", "res", "metric.code","detrend_meth","deseason_meth","method_code")] %>%
  .[,indicator := paste(metric.code,sub(".*_", "",method_code),sep = "_")] %>%
  merge(.,as.data.table(lake_outcome_troph),by="reference") %>%
  .[,detrend_meth := factor(detrend_meth, levels = c("none","linear","loess","gaussian"))] %>%
  .[,deseason_meth := factor(deseason_meth, levels = c("none","average","decompose","stl"))] %>%
  as.data.frame()


its <- 10000
thn <- 0.0005*its
wrmup <- 0.1*its

ews_mod_detrend_mth_uni_roll <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                          data = computation_data[computation_data$res == "Monthly" & grepl("univariate_rolling",computation_data$method_code)  & computation_data$lake %in% c("Kinneret", "Kasumigaura", "Monona", "Washington"),] |>
                                            mutate(combo_code = relevel(factor(paste(method_code,detrend_meth,deseason_meth,sep = "_")), ref = "univariate_rolling_none_none"))
                                          ,
                                          iter = 10000,
                                          thin =thn,
                                          warmup = 2000,
                                          prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                                   prior(normal(1, 1),class = sd)), #1.2 sd as corresponds to 1% and 99% success rate
                                          family = binomial(link = "logit"), 
                                          chains = 4,
                                          backend = "cmdstanr",
                                          control = list(adapt_delta =0.99, max_treedepth = 20,step_size = 0.01),
                                          seed = 12345, cores = 4,sample_prior = TRUE)
#linear_stl
ews_mod_detrend_mth_uni_exp <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                         data = computation_data[computation_data$res == "Monthly" & grepl("univariate_expanding",computation_data$method_code)  & computation_data$lake %in% c("Kinneret", "Kasumigaura", "Monona", "Washington"),] |>
                                           mutate(combo_code = relevel(factor(paste(method_code,detrend_meth,deseason_meth,sep = "_")), ref = "univariate_expanding_none_none"))
                                         ,
                                         iter = 10000,
                                         thin =thn,
                                         warmup = 2000,
                                         prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                                  prior(normal(1, 1),class = sd)), #1.2 sd as corresponds to 1% and 99% success rate
                                         family = binomial(link = "logit"), 
                                         chains = 4,
                                         backend = "cmdstanr",
                                         control = list(adapt_delta = .99, max_treedepth = 20,step_size = 0.01),
                                         seed = 12345, cores = 4,sample_prior = TRUE)
#gaussian-decompose
ews_mod_detrend_mth_ml <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                    data = computation_data[computation_data$res == "Monthly" & grepl("ML",computation_data$method_code)  & computation_data$lake %in% c("Kinneret", "Kasumigaura", "Monona", "Washington"),] |>
                                      mutate(combo_code = relevel(factor(paste(method_code,detrend_meth,deseason_meth,sep = "_")), ref = "univariate_ML_none_none"))
                                    ,
                                    iter = 10000,
                                    thin =thn,
                                    warmup = 2000,
                                    prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                             prior(normal(1, 1),class = sd)), #1.2 sd as corresponds to 1% and 99% success rate
                                    family = binomial(link = "logit"), 
                                    chains = 4,
                                    backend = "cmdstanr",
                                    control = list(adapt_delta = .99, max_treedepth = 20,step_size = 0.01),
                                    seed = 12345, cores = 4,sample_prior = TRUE)
#linear-none

ews_mod_tot_ind_mth_true <- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 + (1|ID1|lake) ), 
                                  data =  rbind(subset(metric_data, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling" & outcome == "trans"),
                                                subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "decompose" & method_code == "univariate_expanding" & outcome == "trans"),
                                                subset(metric_data, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "none" & method_code == "univariate_ML"& outcome == "trans")),
                                  iter = its,
                                  thin = thn,
                                  warmup = wrmup,
                                  prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                           prior(normal(1, 1),class = sd)),
                                  family = binomial(link = "logit"), 
                                  chains = 4,
                                  backend = "cmdstanr",
                                  control = list(adapt_delta = .975, max_treedepth = 20),
                                  seed = 12345, cores = 4,sample_prior = TRUE)

ews_mod_tot_ind_yr_true<- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 +(1|ID1|lake) ), 
                                data = rbind(subset(metric_data, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling" & outcome == "trans"),
                                            subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "decompose" & method_code == "univariate_expanding" & outcome == "trans"),
                                            subset(metric_data, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "none" & method_code == "univariate_ML"& outcome == "trans")),
                                iter = its,
                                thin = thn,
                                warmup = wrmup,
                                prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                         prior(normal(1, 1),class = sd)),
                                family = binomial(link = "logit"), 
                                chains = 4,
                                backend = "cmdstanr",
                                control = list(adapt_delta = .975, max_treedepth = 20),
                                seed = 12345, cores = 4,sample_prior = TRUE)

dat_tot_ind_true_trials <-  ews_mod_tot_ind_mth_true |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind(ews_mod_tot_ind_yr_true |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_indicator", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                 res = "Yearly")) |>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = forcats::fct_reorder(.variable,.value,.fun = mean,.desc = FALSE))  #reorder variable in to increasing ability

dat_tot_ind_true_halfeye <- ews_mod_tot_ind_mth_true |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind( ews_mod_tot_ind_yr_true |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_indicator", "", .variable)) |> 
           mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                  res = "Yearly"))|>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = factor(.variable, levels = levels(dat_tot_ind_true_trials$.variable))) |> #reorder variable in to increasing average ability based on dat_metric_true_trials
  mutate(.variable = sub("_.*", "",.variable))

ews_mod_tot_ind_mth_false <- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 + (1|ID1|lake)), 
                                   data = rbind(subset(metric_data, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling" & outcome == "no.trans"),
                                                subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "decompose" & method_code == "univariate_expanding" & outcome == "no.trans"),
                                                subset(metric_data, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "none" & method_code == "univariate_ML"& outcome == "no.trans")),
                                   iter = its,
                                   thin = thn,
                                   warmup = wrmup,
                                   prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                            prior(normal(1, 1),class = sd)),
                                   family = binomial(link = "logit"), 
                                   chains = 4,
                                   backend = "cmdstanr",
                                   control = list(adapt_delta = .99, max_treedepth = 20),
                                   seed = 12345, cores = 4,sample_prior = TRUE)

ews_mod_tot_ind_yr_false<- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 + (1|ID1|lake)), 
                                 data =  rbind(subset(metric_data, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling" & outcome == "no.trans"),
                                               subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "decompose" & method_code == "univariate_expanding" & outcome == "no.trans"),
                                               subset(metric_data, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "none" & method_code == "univariate_ML"& outcome == "no.trans")),
                                 iter = its,
                                 thin = thn,
                                 warmup = wrmup,
                                 prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                          prior(normal(1, 1),class = sd)),
                                 family = binomial(link = "logit"), 
                                 chains = 4,
                                 backend = "cmdstanr",
                                 control = list(adapt_delta = .99, max_treedepth = 20,stepsize = 0.01),
                                 seed = 12345, cores = 4,sample_prior = TRUE)

dat_tot_ind_false_trials <-  ews_mod_tot_ind_mth_false |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind(ews_mod_tot_ind_yr_false |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_indicator", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                 res = "Yearly")) |>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = factor(.variable,levels = levels(dat_tot_ind_true_trials$.variable)))  #reorder variable in to increasing ability based upon true postive results

dat_tot_ind_false_halfeye <- ews_mod_tot_ind_mth_false |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind( ews_mod_tot_ind_yr_false |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_indicator", "", .variable)) |> 
           mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                  res = "Yearly"))|>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = factor(.variable, levels = levels(dat_tot_ind_true_trials$.variable))) |> #reorder variable in to increasing average ability based on dat_metric_true_trials
  mutate(.variable = sub("_.*", "",.variable))

require(patchwork)

ggsave(
ggplot(data = dat_tot_ind_true_trials |>   
         mutate(.variable = sub("_.*", "",.variable),
                method_code = ifelse(method_code == "ML","EWSNet",method_code),
                reference = paste(variate,method_code,sep = "_")),
       aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= subset(dat_tot_ind_true_halfeye |>
                                      mutate( method_code = ifelse(method_code == "ML","EWSNet",method_code),
                                              reference = paste(variate,method_code,sep = "_")),
                                    method_code == "EWSNet"), alpha=0.5, aes(fill=reference),normalize = "panels") +
  tidybayes::stat_slab(data= subset(dat_tot_ind_true_halfeye |>
                                      mutate( method_code = ifelse(method_code == "ML","EWSNet",method_code),
                                              reference = paste(variate,method_code,sep = "_")
                                      ), variate != "EWSNet"),
                       aes(fill=reference,group=method_code),alpha=0.65,position = position_dodge(width=0.65),normalize = "panels") +
  scale_fill_manual(values=c("#529928", "#403F14","#bfbd3d"),
                    labels = c("EWSNet", "Univariate\nexpanding","Univariate\nrolling"))+
  tidybayes::geom_pointinterval(data = ~subset(.,method_code %in% c("expanding","rolling")),aes(xmin = .lower, xmax = .upper,group = method_code),colour = "black",position = position_dodge(width=0.65),interval_size_range = c(0.4, 1.2)) +
  tidybayes::geom_pointinterval(data = ~subset(.,method_code == "EWSNet"),aes(xmin = .lower, xmax = .upper),colour = "black",interval_size_range = c(0.4, 1.2)) +
  labs(x="True positive prediction probability", y = "Early warning signal indicator",
       fill = "EWS method",colour = "Computation") +
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-3,3))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'lines')) +
ggplot(data = dat_tot_ind_false_trials |>   
         mutate(.variable = sub("_.*", "",.variable),
                method_code = ifelse(method_code == "ML","EWSNet",method_code),
                reference = paste(variate,method_code,sep = "_")),
       aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= subset(dat_tot_ind_false_halfeye |>
                                      mutate( .variable = sub("_.*", "",.variable),
                                              method_code = ifelse(method_code == "ML","EWSNet",method_code),
                                              reference = paste(variate,method_code,sep = "_")),
                                    method_code == "EWSNet"), alpha=0.5, aes(fill=reference),normalize = "panels") +
  tidybayes::stat_slab(data= subset(dat_tot_ind_false_halfeye |>
                                      mutate( .variable = sub("_.*", "",.variable),
                                              method_code = ifelse(method_code == "ML","EWSNet",method_code),
                                              reference = paste(variate,method_code,sep = "_")
                                      ), variate != "EWSNet"),
                       aes(fill=reference,group=method_code),alpha=0.65,position = position_dodge(width=0.65),normalize = "panels") +
  
  scale_fill_manual(values=c("#529928", "#403F14","#bfbd3d"),
                    labels = c("EWSNet", "Univariate\nexpanding","Univariate\nrolling"))+
  tidybayes::geom_pointinterval(data = ~subset(.,method_code %in% c("expanding","rolling")),aes(xmin = .lower, xmax = .upper,group = method_code),colour = "black",position = position_dodge(width=0.65),interval_size_range = c(0.4, 1.2)) +
  tidybayes::geom_pointinterval(data = ~subset(.,method_code == "EWSNet"),aes(xmin = .lower, xmax = .upper),colour = "black",interval_size_range = c(0.4, 1.2)) +
  labs(x="True negative prediction probability", y = "Early warning signal indicator",
       fill = "EWS method",colour = "Computation") +
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-3,3))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  theme_bw()+
  theme(axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'lines')) +
  patchwork::plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A'),
filename = "/Users/ul20791/Desktop/Academia/PhD/Repositories/ews-assessments/tot_ews.pdf",width = 9,height = 6)


##########################################################################################
# Fit ROCs
##########################################################################################
all_roc_data <- as.data.table(suc_ewsnet_max_tot) %>%
  .[,.SD[1], by = c("data_source","troph_level","scaling","lake","res","detrend_meth","deseason_meth")] %>%
  .[,c("data_source","troph_level","scaling","lake","res","method", "computation","prediction","detrend_meth","deseason_meth","prob")] %>%
  setnames(old = "scaling",new= "metric.code") %>%
  setnames(old = "prob",new= "metric_score") %>%
  rbind(as.data.table(suc_exp_uni_tot)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth","prop_run")] %>% setnames(old = "prop_run",new= "metric_score")) %>%
  rbind(as.data.table(suc_roll_perm_uni_tot)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth","obs")] %>% setnames(old = "obs",new= "metric_score")) %>%
  .[,prediction := ifelse(prediction %in% c("match","prior"),
                          1,0)] %>%
  .[,variate := ifelse(grepl("multivariate",method),"multivariate","univariate")]%>%
  .[,method_code := paste(variate,computation,sep = "_")] %>%
  .[,method_code := factor(method_code, levels = c("univariate_rolling","univariate_expanding","univariate_ML","multivariate_rolling","multivariate_expanding"))] %>%
  .[,computation := factor(computation, levels = c("rolling","expanding","ML"))] %>%
  .[,variate := factor(variate, levels = c("univariate","multivariate"))] %>%
  .[,res := factor(res, levels = c("Yearly","Monthly"))]

roc_data <- all_roc_data %>%
  .[,reference := paste(lake,res,troph_level,sep="_")] %>% #trials in terms of assessed time series
  merge(.,as.data.table(lake_outcome_troph),by="reference") %>%
  .[,detrend_meth := factor(detrend_meth, levels = c("none","linear","loess","gaussian"))] %>%
  .[,deseason_meth := factor(deseason_meth, levels = c("none","average","decompose","stl"))] %>%
  as.data.frame()

data_grid <- expand.grid("method_code" = unique(roc_data$method_code),
                         "res" = unique(roc_data$res),
                         "detrend_meth" = unique(roc_data$detrend_meth),
                         "deseason_meth" = unique(roc_data$deseason_meth))

ews_roc <- lapply(1:NROW(data_grid),function(x){
  
  dat <- subset(roc_data,
                method_code == data_grid$method_code[x] &
                  detrend_meth == data_grid$detrend_meth[x] &
                  deseason_meth == data_grid$deseason_meth[x] &
                  res == data_grid$res[x]) |>
    dplyr::mutate(metric.code = as.character(metric.code)) |>
    tidyr::drop_na(metric_score)
  
  metrics <- c("all",unique(dat$metric.code))
  out.dat <- lapply(metrics,function(j){
    
    if(j != "all"){
      sub_dat <- subset(dat,metric.code == j) |>
        dplyr::mutate(prediction = case_when(prediction == 1 & outcome == "trans" ~ "trans",
                                             prediction != 1 & outcome == "trans" ~ "no.trans",
                                             prediction == 1 & outcome == "no.trans" ~ "no.trans",
                                             prediction != 1 & outcome == "no.trans" ~ "trans")) |>
        dplyr::mutate(dplyr::across(c(prediction,outcome), ~ifelse(.x == "trans",1,0))) |>         #convert from correct prediction flag to actual prediction label
        dplyr::mutate(dplyr::across(c(prediction,outcome),~ factor(.x,levels = c(0,1))))
    }else{
      sub_dat <- dat |>
        dplyr::mutate(metric.code = "all") |>
        dplyr::mutate(prediction = case_when(prediction == 1 & outcome == "trans" ~ "trans",
                                             prediction != 1 & outcome == "trans" ~ "no.trans",
                                             prediction == 1 & outcome == "no.trans" ~ "no.trans",
                                             prediction != 1 & outcome == "no.trans" ~ "trans"))|>
        dplyr::mutate(dplyr::across(c(prediction,outcome), ~ifelse(.x == "trans",1,0))) |>         #convert from correct prediction flag to actual prediction label
        dplyr::mutate(dplyr::across(c(prediction,outcome),~ factor(.x,levels = c(0,1))))
    }
    
    pred <- tryCatch(ROCR::prediction(sub_dat$metric_score, as.numeric(sub_dat$outcome)),
                     error = function(err){NA})
    # pred <- tryCatch(ROCR::prediction(as.numeric(sub_dat$prediction), as.numeric(sub_dat$outcome)),
    #                  error = function(err){NA})
    if(!inherits(pred,"prediction")){
      out <- data.frame("fpr" = NA,"tpr" = NA,
                        "auc" = NA, "pr_auc" = NA,
                        "f1_binary" = NA,
                        "f1_cont" = NA) |>
        cbind(sub_dat |> dplyr::select(method,computation,variate,method_code,detrend_meth,deseason_meth,res, metric.code) |>
                dplyr::distinct()) 
    }else{
      auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
      pr_auc <- ROCR::performance(pred,measure = "aucpr")@y.values[[1]]
      
      per <- ROCR::performance(pred,"tpr","fpr")
      f1 <- na.omit(tail(ROCR::performance(pred,"f")@y.values[[1]],n=1))
      
      c_mat <-caret::confusionMatrix(data = sub_dat$prediction,
                                     reference = sub_dat$outcome)
      
      out <- data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],
                        "auc" = auc, "pr_auc" = pr_auc,
                        "f1_binary" = c_mat$byClass[["F1"]],
                        "f1_cont" = f1,
                        "bal_acc" = c_mat$byClass[["Balanced Accuracy"]]) |>
        cbind(sub_dat |> dplyr::select(method,computation,variate,method_code,detrend_meth,deseason_meth,res, metric.code) |>
                dplyr::distinct()) 
    }
    return(out)
    
  }) |>
    data.table::rbindlist()
  return(out.dat)
}) |>
  data.table::rbindlist()

plot_data <- rbind(subset(ews_roc, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling"),
                   subset(ews_roc, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "decompose" & method_code == "univariate_expanding"),
                   subset(ews_roc, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "none" & method_code == "univariate_ML")) |>
  rbind(subset(ews_roc, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling"),
        subset(ews_roc, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "decompose" & method_code == "univariate_expanding"),
        subset(ews_roc, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "none" & method_code == "univariate_ML")) |>
  dplyr::mutate(variate = ifelse(grepl("ML",method_code),"EWSNet",ifelse(grepl("multivariate",method_code),"Multivariate","Univariate")),
                res = factor(res,levels = c("Monthly","Yearly"))) |>
  dplyr::select(metric.code,method_code,f1_cont,f1_binary,auc,pr_auc,bal_acc,res,variate) |>
  dplyr::distinct()


ggplot(subset(plot_data,metric.code != "all"),
       aes(x=bal_acc,y=metric.code))  +
  geom_vline(xintercept = 0.5,linetype = "dashed", colour="grey50") + 
  geom_point(data = subset(plot_data |>
                             dplyr::mutate(method_code = as.character(method_code),
                                           method_code = ifelse(grepl("ML",method_code),"EWSNet",method_code),
                                           reference = paste(variate,method_code,sep = "_")),
                           method_code == "EWSNet" & metric.code != "all"), 
             aes(fill = reference),size = 4, shape = 21,alpha = 0.8) +
  geom_point(data= subset(plot_data |>
                            mutate(method_code = as.character(method_code),
                                   method_code = ifelse(grepl("ML",method_code),"EWSNet",method_code),
                                   reference = paste(variate,method_code,sep = "_")
                            ), method_code != "EWSNet" & metric.code != "all"),
             aes(fill = reference),size = 4, shape = 21,alpha = 0.8) +
  scale_color_manual(values=c("#529928","#403F14","#bfbd3d"),
                     labels = c("EWSNet","Multivariate\nexpanding","Multivariate\nrolling","Univariate\nexpanding","Univariate\nrolling"),
                     guide = "none")+
  scale_fill_manual(values=c("#529928","#403F14","#bfbd3d"),
                    labels = c("EWSNet","Univariate\nexpanding","Univariate\nrolling"))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  labs(x="F1 statistic", y = "Early warning signal indicator",
       fill = "EWS method",colour = "Computation")+
  scale_x_continuous(breaks = c(0.2,0.5,0.8),limits = c(0,1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())



        