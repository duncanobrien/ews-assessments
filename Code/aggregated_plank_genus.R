require(EWSmethods)
require(Matrix)
require(tidyverse)
require(foreach) 
require(pbapply)

source("Code/perm_rollEWS_fn.R") 

EWSmethods::ewsnet_init("EWSNET_env", auto=T)
reticulate::py_config() #check for 'forced by use_python

EWSmethods::ewsnet_reset(auto=T)

################################################################################################################
## Load Data ##
################################################################################################################

source("/Users/ul20791/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_environmental_data.R")
#phyo, zoo + combined raw concentrations at week, month and year groupings

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

################################################################################################################
## Wangle Data ##

# save(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
#      kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat,
#     file =  "Data/wrangled_genus_plank_data.Rdata" )

#prewrangled load("Data/wrangled_genus_plank_data.Rdata")
################################################################################################################

kin_yr_dat <- plank_env.data.yr |> #drop environmentals
  dplyr::filter(as.numeric(Date) <= as.numeric(Date[25])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(Date:`Anuraeopsis fissa`) |> #drop environmentals
  pivot_longer(-Date,names_to = "species",values_to = "density") |>
  #mutate(troph_level = ifelse(species %in% colnames(plank_env.data.yr)[47:77],"Zooplankton","Phytoplankton")) |>
  mutate(genus = trimws(gsub("[0-9]|,","",species),which = "left",whitespace = "-"), #remove numbers and leading "-" from species name
         genus = gsub( " .*$", "", genus )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(Date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(Anabaena:Thermocyclops,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Anabaena:Thermocyclops) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(Date)) |>
  as.data.frame()

kin_mth_dat <- plank_env.data.mth |> 
  dplyr::filter(as.numeric(Date) <= as.numeric(Date[288])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(Date:`Anuraeopsis fissa`) |>
  tidyr::pivot_longer(-Date,names_to = "species",values_to = "density") |>
  mutate(genus = trimws(gsub("[0-9]|,","",species),which = "left",whitespace = "-"), #remove numbers and leading "-" from species name
         genus = gsub( " .*$", "", genus )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(Date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(Anabaena:Thermocyclops,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Anabaena:Thermocyclops) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(Date = base::as.Date(zoo::as.Date(Date)))|>
  #deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(zoo::as.yearmon(Date))) |>
  as.data.frame()

kas_yr_dat <- plank_env.kasyrdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[16])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date:Thermocyclops_taihokuensis) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(Acanthoceras:Thermocyclops,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Acanthoceras:Thermocyclops) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(date)) |>
  as.data.frame()

kas_mth_dat <- plank_env.kasmthdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[185])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(c(date,Acanthoceras_zachariasii:Thermocyclops_taihokuensis)) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(Acanthoceras:Thermocyclops,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Acanthoceras:Thermocyclops) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(date = base::as.Date(zoo::as.Date(date)))|>
  #deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  as.data.frame()

LZ_yr_dat <- plank_env.LZyrdata |>
  dplyr::filter(as.numeric(date) <= as.numeric(date[length(date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date,Microcystis_sp:Polyphemus_pediculus) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(Actinastrum:Willea,~log1p(.x))) |>
  ungroup() |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(Actinastrum:Willea,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Actinastrum:Willea) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(date)) |>
  as.data.frame()

LZ_mth_dat <- plank_env.LZmthdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[length(date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date,Microcystis_sp:Polyphemus_pediculus) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(Actinastrum:Willea,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Actinastrum:Willea) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(date = base::as.Date(zoo::as.Date(date))) |>
  #deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  as.data.frame()

mad_yr_dat <- plank_env.madyrdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[15]))|>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date:TROPOCYCLOPS_PRASINUS_MEXICANUS) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species ),
         genus = gsub( " .*$", "", genus )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(ACANTHOCYCLOPS:Uroglena,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(ACANTHOCYCLOPS:Uroglena) %>%
                                   dplyr::select_if(colSums(.) != 0))) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2))|>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(date)) |>
  as.data.frame()

mad_mth_dat <- plank_env.madmthdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[168]))|>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date,`Achnanthes minutissima`:TROPOCYCLOPS_PRASINUS_MEXICANUS) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species ),
         genus = gsub( " .*$", "", genus )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(ACANTHOCYCLOPS:Uroglena,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(ACANTHOCYCLOPS:Uroglena) %>%
                                   dplyr::select_if(colSums(.) != 0))) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(date = base::as.Date(zoo::as.Date(date))) |>
  #deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  as.data.frame()
#dplyr::select(-Aphanothece.clathrata) #randomly dropped species due to bug in svd

wind_yr_dat <- phyto_env.windyrdata |>
  dplyr::filter(as.numeric(Date) <= as.numeric(Date[length(Date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(Asterionella:TotCyclopoids,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Asterionella:TotCyclopoids) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  select(c(Date:TotCyclopoids,pca1,pca2)) |>
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(Date)) |>
  as.data.frame()

wind_mth_dat <- phyto_env.windmthdata |>
  dplyr::filter(as.numeric(Date) <= as.numeric(Date[length(Date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(Asterionella:TotCyclopoids,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Asterionella:TotCyclopoids) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  select(c(Date,Asterionella:TotCyclopoids,pca1,pca2))|>
  mutate(Date = base::as.Date(zoo::as.Date(Date))) |>
  #deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(zoo::as.yearmon(Date))) |>
  as.data.frame()

UZ_yr_dat <- plank_env.UZyrdata |>
  dplyr::filter(as.numeric(date) <= as.numeric(date[length(date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date,Ankistrodesmus_sp.:Polyphemus_pediculus) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(Achnanthes:Volvox,~log1p(.x))) |>
  ungroup() |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(Achnanthes:Volvox,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Achnanthes:Volvox) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(date)) |>
  as.data.frame()

UZ_mth_dat <- plank_env.UZmthdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[length(date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date,Ankistrodesmus_sp.:Polyphemus_pediculus) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(Achnanthes:Volvox,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Achnanthes:Volvox) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(date = base::as.Date(zoo::as.Date(date))) |>
  #deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  as.data.frame()

mon_yr_dat <- plank_env.monyrdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[length(date)*0.85]))|>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date:TROPOCYCLOPS_PRASINUS_MEXICANUS) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species ),
         genus = gsub( " .*$", "", genus )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(ACANTHOCYCLOPS:Woronichinia,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(ACANTHOCYCLOPS:Woronichinia) %>%
                                   dplyr::select_if(colSums(.) != 0))) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2))|>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(date)) |>
  as.data.frame()

mon_mth_dat <- plank_env.monmthdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[length(date)*0.85]))|>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date,`Achnanthes minutissima`:TROPOCYCLOPS_PRASINUS_MEXICANUS) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species ),
         genus = gsub( " .*$", "", genus )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(ACANTHOCYCLOPS:Woronichinia,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(ACANTHOCYCLOPS:Woronichinia) %>%
                                   dplyr::select_if(colSums(.) != 0))) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(date = base::as.Date(zoo::as.Date(date))) |>
  #deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  as.data.frame()

leve_yr_dat <- plank_env.leveyrdata |> 
  dplyr::filter(as.numeric(Date) <= as.numeric(Date[length(Date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(Date:Daphnia) |> #drop environmentals
  pivot_longer(-Date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species ),
         genus = gsub( " .*$", "", genus )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(Date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(Asterionella:Unicellular,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Asterionella:Unicellular) %>%
                                   dplyr::select_if(colSums(.) != 0))) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2))|>
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(Date)) |>
  as.data.frame()

leve_mth_dat <- plank_env.levemthdata |> 
  dplyr::filter(as.numeric(Date) <= as.numeric(Date[length(Date)*0.85]))|>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(Date,`Asterionella formosa`:`Eudiaptomus nauplii`) |> #drop environmentals
  pivot_longer(-Date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species ),
         genus = gsub( " .*$", "", genus )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(Date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  mutate(across(Asterionella:Unicellular,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Asterionella:Unicellular) %>%
                                   dplyr::select_if(colSums(.) != 0))) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(Date = base::as.Date(zoo::as.Date(Date))) |>
  #deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(zoo::as.yearmon(Date))) |>
  as.data.frame()

wash_yr_dat <- plank_env.washyrdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[length(date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(c(date,cryptomonas:nc.rotifers)) |> #drop environmentals
  mutate(across(cryptomonas:nc.rotifers,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(cryptomonas:nc.rotifers) %>%
                                   dplyr::select_if(colSums(.) != 0))) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2))|>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(date)) |>
  as.data.frame()

wash_mth_dat <- plank_env.washmthdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[length(date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(c(date,cryptomonas:nc.rotifers)) |> #drop environmentals
  mutate(across(cryptomonas:nc.rotifers,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(cryptomonas:nc.rotifers) %>%
                                   dplyr::select_if(colSums(.) != 0))) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2))|>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(date)) |>
  as.data.frame()

################################################################################################################
## EWSNet Assessment ##
################################################################################################################

genus_lake_ewsnet <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
                                            wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
                                            kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
                                            wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat),
                                 FUN = function(x){
                                   
                                   x = base::subset(x,select = -c(pca1,pca2))

                                   x <- EWSmethods::detrend_ts(x,method = "loess",span = 0.75)
                                   
                                   if(any(str_length(x[,1])>4)){
                                     x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                     x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                   }
                                   
                                   foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                     
                                     ews.tmp.scaled <- ewsnet_predict(c(x[,i]),scaling = T, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                     ews.tmp.unscaled <- ewsnet_predict(c(x[,i]),scaling = F, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                     
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
    method = "EWSNet",
    computation = "ML") |>
  select(-name)

write.csv(genus_lake_ewsnet,file = "Results/lake_results/genus_lake_ewsnet.csv")


genus_lake_phyto_ewsnet <- pbapply::pblapply(list(kin_yr_dat[,c(1:28)],kas_yr_dat[,c(1:63)],LZ_yr_dat[,c(1:119)],mad_yr_dat[,c(1:87)],
                                                  wind_yr_dat[,c(1:18)],wash_yr_dat[,c(1:7)],leve_yr_dat[,c(1:5)],UZ_yr_dat[,c(1:61)],mon_yr_dat[,c(1:90)],
                                                  kin_mth_dat[,c(1:28)],kas_mth_dat[,c(1:63)],LZ_mth_dat[,c(1:119)],mad_mth_dat[,c(1:86)],
                                                  wind_mth_dat[,c(1:18)],wash_mth_dat[,c(1:7)],leve_mth_dat[,c(1:5)],UZ_mth_dat[,c(1:61)],mon_mth_dat[,c(1:87)]),
                                       FUN = function(x){
                                         
                                         x <- EWSmethods::detrend_ts(x,method = "loess",span = 0.75)
                                         
                                         if(any(str_length(x[,1])>4)){
                                           x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                           x <- EWSmethods::deseason_ts(data=x,increment = "month", method = "stl",order = "ymd")
                                         }
                                         
                                         foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                           
                                           ews.tmp.scaled <- ewsnet_predict(c(x[,i]),scaling = T, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                           ews.tmp.unscaled <- ewsnet_predict(c(x[,i]),scaling = F, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                           
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
    method = "EWSNet",
    computation = "ML") |>
  mutate(troph_level = "phytoplankton")|>
  select(-name)

write.csv(genus_lake_phyto_ewsnet,file = "Results/lake_results/genus_lake_phyto_ewsnet.csv")

genus_lake_zoo_ewsnet <- pbapply::pblapply(list(kin_yr_dat[,c(1,29:50)],kas_yr_dat[,c(1,64:89)],LZ_yr_dat[,c(1,120:137)],mad_yr_dat[,c(1,88:105)],
                                                wind_yr_dat[,c(1,19:23)],wash_yr_dat[,c(1,8:16)],leve_yr_dat[,c(1,6:8)],UZ_yr_dat[,c(1,62:79)],mon_yr_dat[,c(1,91:107)],
                                                kin_mth_dat[,c(1,29:50)],kas_mth_dat[,c(1,64:89)],LZ_mth_dat[,c(1,120:137)],mad_mth_dat[,c(1,87:103)],
                                                wind_mth_dat[,c(1,19:23)],wash_mth_dat[,c(1,8:16)],leve_mth_dat[,c(1,6:7)],UZ_mth_dat[,c(1,62:79)],mon_mth_dat[,c(1,89:105)]),
                                             FUN = function(x){
                                               
                                               x <- EWSmethods::detrend_ts(x,method = "loess",span = 0.75)
                                               
                                               if(any(str_length(x[,1])>4)){
                                                 x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                                 x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                               }
                                               
                                               foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                                 
                                                 ews.tmp.scaled <- ewsnet_predict(c(x[,i]),scaling = T, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                                 ews.tmp.unscaled <- ewsnet_predict(c(x[,i]),scaling = F, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                                 
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
    method = "EWSNet",
    computation = "ML") |>
  mutate(troph_level = "zoooplankton")|>
  select(-name)

write.csv(genus_lake_zoo_ewsnet,file = "Results/lake_results/genus_lake_zoo_ewsnet.csv")


################################################################################################################
## Multivariate EWS Assessment Rolling ##
################################################################################################################

genus_lake_roll_multi_ews <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
                                                    wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
                                                    kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
                                                    wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat),
                                         FUN = function(x){

                                           sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                             mean(x[,i] == min(x[,i])) <= 0.49
                                           })]
                                           
                                           sub_dat_det <- EWSmethods::detrend_ts(subset(sub_dat, select=-c(pca1,pca2)),method = "loess",span = 0.75,degree=1) 
                                           
                                           # if(any(str_length(sub_dat_det[,1])>4)){
                                           #   sub_dat_det[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat_det[,1]))
                                           #   sub_dat_det <- EWSmethods::deseason_ts(sub_dat_det,increment = "month", method = "stl",order = "ymd")
                                           # }
                                           
                                           out.det <- EWSmethods::multiEWS(sub_dat_det,
                                                                       metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                       method = "rolling",ggplotIt = F,
                                                                       winsize = 50)
                                           
                                           out.undet <- EWSmethods::multiEWS(subset(sub_dat, select=-c(pca1,pca2)),
                                                                            metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                        "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                            method = "rolling",ggplotIt = F,
                                                                            winsize = 50)
                                           
                                           return(data.frame(cbind(out.det$cor,out.undet$cor)))
                                           
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
    method = "multivariate EWS",
    computation = "rolling") |>
  select(-name)

write.csv(genus_lake_roll_multi_ews,file = "Results/lake_results/genus_lake_roll_multi_ews.csv")


genus_lake_zoo_roll_multi_ews <- pbapply::pblapply(list(kin_yr_dat[,c(1,29:50)],kas_yr_dat[,c(1,64:89)],LZ_yr_dat[,c(1,120:137)],mad_yr_dat[,c(1,88:105)],
                                                          wind_yr_dat[,c(1,19:23)],wash_yr_dat[,c(1,8:16)],leve_yr_dat[,c(1,6:8)],UZ_yr_dat[,c(1,62:79)],mon_yr_dat[,c(1,91:107)],
                                                    kin_mth_dat[,c(1,29:50)],kas_mth_dat[,c(1,64:89)],LZ_mth_dat[,c(1,120:137)],mad_mth_dat[,c(1,87:103)],
                                                    wind_mth_dat[,c(1,19:23)],wash_mth_dat[,c(1,8:16)],leve_mth_dat[,c(1,6:7)],UZ_mth_dat[,c(1,62:79)],mon_mth_dat[,c(1,89:105)]),
                                               
                                                    FUN = function(x){
                                                 

                                                 sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                   mean(x[,i] == min(x[,i])) <= 0.49
                                                 })]
                                                 
                                                 
                                                 sub_dat_det <- EWSmethods::detrend_ts(sub_dat,method = "loess",span = 0.75,degree=1) 
                                                 
                                                 # if(any(str_length(sub_dat_det[,1])>4)){
                                                 #   sub_dat_det[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat_det[,1]))
                                                 #   sub_dat_det <- EWSmethods::deseason_ts(sub_dat_det,increment = "month", method = "stl",order = "ymd")
                                                 # }
                                                 
                                                 out.det <- EWSmethods::multiEWS(sub_dat_det,
                                                                                 metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                                 method = "rolling",ggplotIt = F,
                                                                                 winsize = 50)
                                                 
                                                 out.undet <- EWSmethods::multiEWS(sub_dat,
                                                                                   metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                               "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                                   method = "rolling",ggplotIt = F,
                                                                                   winsize = 50)
                                                 
                                                 return(data.frame(cbind(out.det$cor,out.undet$cor)))
                                                 
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
    method = "multivariate EWS",
    computation = "rolling") |>
  mutate(troph_level = "zooplankton")|>
  select(-name)

write.csv(genus_lake_zoo_roll_multi_ews,file = "Results/lake_results/genus_lake_zoo_roll_multi_ews.csv")

genus_lake_phyto_roll_multi_ews <- pbapply::pblapply(list(kin_yr_dat[,c(1:28)],kas_yr_dat[,c(1:63)],LZ_yr_dat[,c(1:119)],mad_yr_dat[,c(1:87)],
                                                        wind_yr_dat[,c(1:18)],wash_yr_dat[,c(1:7)],leve_yr_dat[,c(1:5)],UZ_yr_dat[,c(1:61)],mon_yr_dat[,c(1:90)],
                                                        kin_mth_dat[,c(1:28)],kas_mth_dat[,c(1:63)],LZ_mth_dat[,c(1:119)],mad_mth_dat[,c(1:86)],
                                                        wind_mth_dat[,c(1:18)],wash_mth_dat[,c(1:7)],leve_mth_dat[,c(1:5)],UZ_mth_dat[,c(1:61)],mon_mth_dat[,c(1:87)]),
                                                   
                                                   FUN = function(x){
                                                     
                                                     
                                                     sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                       mean(x[,i] == min(x[,i])) <= 0.47 #0.47 required because of Chroococcus
                                                     })]
                                                     
                                                     sub_dat_det <- EWSmethods::detrend_ts(sub_dat,method = "loess",span = 0.75,degree=1) 
                                                     
                                                     # if(any(str_length(sub_dat_det[,1])>4)){
                                                     #   sub_dat_det[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat_det[,1]))
                                                     #   sub_dat_det <- EWSmethods::deseason_ts(sub_dat_det,increment = "month", method = "stl",order = "ymd")
                                                     # }
                                                     
                                                     out.det <- EWSmethods::multiEWS(sub_dat_det,
                                                                                     metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                                     method = "rolling",ggplotIt = F,
                                                                                     winsize = 50)
                                                     
                                                     out.undet <- EWSmethods::multiEWS(data = sub_dat,
                                                                                       metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                                   "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                                       method = "rolling",ggplotIt = F,
                                                                                       winsize = 50)
                                                     
                                                     return(data.frame(cbind(out.det$cor,out.undet$cor)))
                                                     
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
    method = "multivariate EWS",
    computation = "rolling") |>
  mutate(troph_level = "phytoplankton")|>
  select(-name)

write.csv(genus_lake_phyto_roll_multi_ews,file = "Results/lake_results/genus_lake_phyto_roll_multi_ews.csv")

################################################################################################################
## Multivariate EWS Assessment Permuted Rolling ##
################################################################################################################
cl <- parallel::makeCluster(parallel::detectCores()-2)
doParallel::registerDoParallel(cl)

genus_lake_roll_multi_ews_perm <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
                                                         wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
                                                         kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
                                                         wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat),
                                                  FUN = function(x){
                                                    
                                                      sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                        mean(x[,i] == min(x[,i])) <= 0.49
                                                      })]

                                                      sub_dat_det <- EWSmethods::detrend_ts(subset(sub_dat, select=-c(pca1,pca2)),method = "loess",span = 0.75,degree=1) 
                                                      
                                                      # if(any(str_length(sub_dat_det[,1])>4)){
                                                      #   sub_dat_det[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat_det[,1]))
                                                      #   sub_dat_det <- EWSmethods::deseason_ts(sub_dat_det,increment = "month", method = "stl",order = "ymd")
                                                      # }
                                                      
                                                      out.det <- perm_rollEWS(data=sub_dat_det,perm.meth = "replacement", 
                                                                                      metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                                      winsize = 50,
                                                                               variate = "multi",
                                                                               iter = 500)
                                                      
                                                      out.undet <- perm_rollEWS(data = subset(sub_dat, select=-c(pca1,pca2)),perm.meth = "replacement", 
                                                                                        metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                                    "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                                winsize = 50,
                                                                                variate = "multi",
                                                                                iter = 500)
                                                      
                                                      return(data.frame(cbind(out.det$cor,out.undet$cor)))
                                                  }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "multivariate EWS",
    computation = "rolling") |>
  select(-name)

parallel::stopCluster(cl)
write.csv(genus_lake_roll_multi_ews_perm,file = "Results/lake_results/genus_lake_roll_multi_ews_perm.csv")


cl <- parallel::makeCluster(parallel::detectCores()-2)
doParallel::registerDoParallel(cl)
genus_lake_zoo_roll_multi_ews_perm <- pbapply::pblapply(list(kin_yr_dat[,c(1,29:50)],kas_yr_dat[,c(1,64:89)],LZ_yr_dat[,c(1,120:137)],mad_yr_dat[,c(1,88:105)],
                                                             wind_yr_dat[,c(1,19:23)],wash_yr_dat[,c(1,8:16)],leve_yr_dat[,c(1,6:8)],UZ_yr_dat[,c(1,62:79)],mon_yr_dat[,c(1,91:107)],
                                                             kin_mth_dat[,c(1,29:50)],kas_mth_dat[,c(1,64:89)],LZ_mth_dat[,c(1,120:137)],mad_mth_dat[,c(1,87:103)],
                                                             wind_mth_dat[,c(1,19:23)],wash_mth_dat[,c(1,8:16)],leve_mth_dat[,c(1,6:7)],UZ_mth_dat[,c(1,62:79)],mon_mth_dat[,c(1,89:105)]),
                                                          FUN = function(x){
                                                            
                                                            sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                              mean(x[,i] == min(x[,i])) <= 0.49
                                                            })]
                                                            
                                                            sub_dat_det <- EWSmethods::detrend_ts(sub_dat,method = "loess",span = 0.75,degree=1) 
                                                            
                                                            # if(any(str_length(sub_dat_det[,1])>4)){
                                                            #   sub_dat_det[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat_det[,1]))
                                                            #   sub_dat_det <- EWSmethods::deseason_ts(sub_dat_det,increment = "month", method = "stl",order = "ymd")
                                                            # }
                                                            
                                                            out.det <- perm_rollEWS(data=sub_dat_det,perm.meth = "replacement", 
                                                                                    metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                                    winsize = 50,
                                                                                    variate = "multi",
                                                                                    iter = 500)
                                                            
                                                            out.undet <- perm_rollEWS(sub_dat,perm.meth = "replacement", 
                                                                                      metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                                  "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                                      winsize = 50,
                                                                                      variate = "multi",
                                                                                      iter = 500)
                                                            
                                                            return(data.frame(cbind(out.det$cor,out.undet$cor)))
                                                            
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
    method = "multivariate EWS",
    computation = "rolling") |>
  mutate(troph_level = "zooplankton")|>
  select(-name)

parallel::stopCluster(cl)
write.csv(genus_lake_zoo_roll_multi_ews_perm,file = "Results/lake_results/genus_lake_zoo_roll_multi_ews_perm.csv")

cl <- parallel::makeCluster(parallel::detectCores()-2)
doParallel::registerDoParallel(cl)
genus_lake_phyto_roll_multi_ews_perm <- pbapply::pblapply(list(kin_yr_dat[,c(1:28)],kas_yr_dat[,c(1:63)],LZ_yr_dat[,c(1:119)],mad_yr_dat[,c(1:87)],
                                                         wind_yr_dat[,c(1:18)],wash_yr_dat[,c(1:7)],leve_yr_dat[,c(1:5)],UZ_yr_dat[,c(1:61)],mon_yr_dat[,c(1:90)],
                                                         kin_mth_dat[,c(1:28)],kas_mth_dat[,c(1:63)],LZ_mth_dat[,c(1:119)],mad_mth_dat[,c(1:86)],
                                                         wind_mth_dat[,c(1:18)],wash_mth_dat[,c(1:7)],leve_mth_dat[,c(1:5)],UZ_mth_dat[,c(1:61)],mon_mth_dat[,c(1:87)]),
                                                    FUN = function(x){
                                                      
                                                      sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                        mean(x[,i] == min(x[,i])) <= 0.49
                                                      })]
                                                      
                                                      sub_dat_det <- EWSmethods::detrend_ts(sub_dat,method = "loess",span = 0.75,degree=1) 
                                                      
                                                      # if(any(str_length(sub_dat_det[,1])>4)){
                                                      #   sub_dat_det[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat_det[,1]))
                                                      #   sub_dat_det <- EWSmethods::deseason_ts(sub_dat_det,increment = "month", method = "stl",order = "ymd")
                                                      # }
                                                      
                                                      out.det <- perm_rollEWS(data=sub_dat_det,perm.meth = "replacement", 
                                                                              metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                              winsize = 50,
                                                                              variate = "multi",
                                                                              iter = 500)
                                                      
                                                      out.undet <- perm_rollEWS(sub_dat,perm.meth = "replacement", 
                                                                                metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                            "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                                winsize = 50,
                                                                                variate = "multi",
                                                                                iter = 500)
                                                      
                                                      return(data.frame(cbind(out.det$cor,out.undet$cor)))
                                                      
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
    method = "multivariate EWS",
    computation = "rolling") |>
  mutate(troph_level = "phytoplankton")|>
  select(-name)

parallel::stopCluster(cl)
write.csv(genus_lake_phyto_roll_multi_ews_perm,file = "Results/lake_results/genus_lake_phyto_roll_multi_ews_perm.csv")



################################################################################################################
## Multivariate EWS Assessment Expanding ##
################################################################################################################

genus_lake_exp_multi_ews <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
                                                   wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
                                                   kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
                                                   wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat),
                                        FUN = function(x){

                                          sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                            mean(x[,i] == min(x[,i])) <= 0.49
                                          })]

                                          sub_dat_det <- EWSmethods::detrend_ts(subset(sub_dat, select=-c(pca1,pca2)),method = "loess",span = 0.75,degree=1) 
                                          
                                          # if(any(str_length(sub_dat_det[,1])>4)){
                                          #   sub_dat_det[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat_det[,1]))
                                          #   sub_dat_det <- EWSmethods::deseason_ts(sub_dat_det,increment = "month", method = "stl",order = "ymd")
                                          # }
                                          out.det <- EWSmethods::multiEWS(sub_dat_det,
                                                                          metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                          method = "expanding",ggplotIt = F,
                                                                          burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                          
                                          out.undet <- EWSmethods::multiEWS(subset(sub_dat, select=-c(pca1,pca2)),
                                                                            metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                        "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                            method = "expanding",ggplotIt = F,
                                                                            burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                          
                                          return(data.frame(cbind(out.det$raw,out.undet$raw)))
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
    method = "multivariate EWS",
    computation = "expanding") |>
  select(-name)

write.csv(genus_lake_exp_multi_ews,file = "Results/lake_results/genus_lake_exp_multi_ews.csv")

genus_lake_phyto_exp_multi_ews <- pbapply::pblapply(list(kin_yr_dat[,c(1:28)],kas_yr_dat[,c(1:63)],LZ_yr_dat[,c(1:119)],mad_yr_dat[,c(1:87)],
                                                         wind_yr_dat[,c(1:18)],wash_yr_dat[,c(1:7)],leve_yr_dat[,c(1:5)],UZ_yr_dat[,c(1:61)],mon_yr_dat[,c(1:90)],
                                                         kin_mth_dat[,c(1:28)],kas_mth_dat[,c(1:63)],LZ_mth_dat[,c(1:119)],mad_mth_dat[,c(1:86)],
                                                         wind_mth_dat[,c(1:18)],wash_mth_dat[,c(1:7)],leve_mth_dat[,c(1:5)],UZ_mth_dat[,c(1:61)],mon_mth_dat[,c(1:87)]),
                                                  FUN = function(x){
                                                    
                                                    sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                      mean(x[,i] == min(x[,i])) <= 0.49
                                                    })]
                                                    
                                                    sub_dat_det <- EWSmethods::detrend_ts(sub_dat,method = "loess",span = 0.75,degree=1) 
                                                    
                                                    # if(any(str_length(sub_dat_det[,1])>4)){
                                                    #   sub_dat_det[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat_det[,1]))
                                                    #   sub_dat_det <- EWSmethods::deseason_ts(sub_dat_det,increment = "month", method = "stl",order = "ymd")
                                                    # }
                                                    
                                                    out.det <- EWSmethods::multiEWS(sub_dat_det,
                                                                                    metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                                    method = "expanding",ggplotIt = F,
                                                                                    burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                                    
                                                    out.undet <- EWSmethods::multiEWS(sub_dat,
                                                                                      metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                                  "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                                      method = "expanding",ggplotIt = F,
                                                                                      burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                                    
                                                    return(data.frame(cbind(out.det$raw,out.undet$raw)))
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
    method = "multivariate EWS",
    computation = "expanding") |>
  mutate(troph_level = "phytoplankton")|>
  select(-name)

write.csv(genus_lake_phyto_exp_multi_ews,file = "Results/lake_results/genus_lake_phyto_exp_multi_ews.csv")

genus_lake_zoo_exp_multi_ews <- pbapply::pblapply(list(kin_yr_dat[,c(1,29:50)],kas_yr_dat[,c(1,64:89)],LZ_yr_dat[,c(1,120:137)],mad_yr_dat[,c(1,88:105)],
                                                   wind_yr_dat[,c(1,19:23)],wash_yr_dat[,c(1,8:16)],leve_yr_dat[,c(1,6:8)],UZ_yr_dat[,c(1,62:79)],mon_yr_dat[,c(1,91:107)],
                                                   kin_mth_dat[,c(1,29:50)],kas_mth_dat[,c(1,64:89)],LZ_mth_dat[,c(1,120:137)],mad_mth_dat[,c(1,87:103)],
                                                   wind_mth_dat[,c(1,19:23)],wash_mth_dat[,c(1,8:16)],leve_mth_dat[,c(1,6:7)],UZ_mth_dat[,c(1,62:79)],mon_mth_dat[,c(1,89:105)]),
                                              FUN = function(x){
                                                
                                                sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                  mean(x[,i] == min(x[,i])) <= 0.49
                                                })]
                                                
                                                sub_dat_det <- EWSmethods::detrend_ts(sub_dat,method = "loess",span = 0.75,degree=1) 
                                                # if(any(str_length(sub_dat_det[,1])>4)){
                                                #   sub_dat_det[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat_det[,1]))
                                                #   sub_dat_det <- EWSmethods::deseason_ts(sub_dat_det,increment = "month", method = "stl",order = "ymd")
                                                # }
                                                out.det <- EWSmethods::multiEWS(sub_dat_det,
                                                                                metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                                method = "expanding",ggplotIt = F,
                                                                                burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                                
                                                out.undet <- EWSmethods::multiEWS(sub_dat,
                                                                                  metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                              "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                                  method = "expanding",ggplotIt = F,
                                                                                  burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                                
                                                return(data.frame(cbind(out.det$raw,out.undet$raw)))
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
    method = "multivariate EWS",
    computation = "expanding") |>
  mutate(troph_level = "zooplankton")|>
  select(-name)

write.csv(genus_lake_zoo_exp_multi_ews,file = "Results/lake_results/genus_lake_zoo_exp_multi_ews.csv")

################################################################################################################
## Univariate EWS Assessment Rolling ##
################################################################################################################

genus_lake_roll_uni_ews <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
                                                  wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
                                                  kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
                                                  wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat),
                                       FUN = function(x){
                                         
                                         x <- EWSmethods::detrend_ts(subset(x, select=-c(pca1,pca2)),method = "loess",span = 0.75)                                           
                                         
                                         # if(any(str_length(x[,1])>4)){
                                         #   x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                         #   x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                         # }
                                         
                                         time_vec <- colnames(x)[1]
                                         
                                         foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                           
                                           out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                     winsize = 50)
                                           
                                           data.frame(cbind(data_source = paste(i),
                                                            out$cor))
                                           
                                         } 
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
    computation = "rolling") |>
  select(-name)

write.csv(genus_lake_roll_uni_ews,file = "Results/lake_results/genus_lake_roll_uni_ews.csv")


genus_lake_phyto_roll_uni_ews <- pbapply::pblapply(list(kin_yr_dat[,c(1:28)],kas_yr_dat[,c(1:63)],LZ_yr_dat[,c(1:119)],mad_yr_dat[,c(1:87)],
                                                        wind_yr_dat[,c(1:18)],wash_yr_dat[,c(1:7)],leve_yr_dat[,c(1:5)],UZ_yr_dat[,c(1:61)],mon_yr_dat[,c(1:90)],
                                                        kin_mth_dat[,c(1:28)],kas_mth_dat[,c(1:63)],LZ_mth_dat[,c(1:119)],mad_mth_dat[,c(1:86)],
                                                        wind_mth_dat[,c(1:18)],wash_mth_dat[,c(1:7)],leve_mth_dat[,c(1:5)],UZ_mth_dat[,c(1:61)],mon_mth_dat[,c(1:87)]),
                                             FUN = function(x){
                                               
                                               x <- EWSmethods::detrend_ts(x,method = "loess",span = 0.75)                                           
                                               # if(any(str_length(x[,1])>4)){
                                               #   x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                               #   x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                               # }
                                               time_vec <- colnames(x)[1]
                                               
                                               foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                                 
                                                 out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                           winsize = 50)
                                                 
                                                 data.frame(cbind(data_source = paste(i),
                                                                  out$cor))
                                               } 
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
    troph_level = "phytoplankton") |>
  select(-name)

write.csv(genus_lake_phyto_roll_uni_ews,file = "Results/lake_results/genus_lake_phyto_roll_uni_ews.csv")

genus_lake_zoo_roll_uni_ews <- pbapply::pblapply(list(kin_yr_dat[,c(1,29:50)],kas_yr_dat[,c(1,64:89)],LZ_yr_dat[,c(1,120:137)],mad_yr_dat[,c(1,88:105)],
                                                      wind_yr_dat[,c(1,19:23)],wash_yr_dat[,c(1,8:16)],leve_yr_dat[,c(1,6:8)],UZ_yr_dat[,c(1,62:79)],mon_yr_dat[,c(1,91:107)],
                                                      kin_mth_dat[,c(1,29:50)],kas_mth_dat[,c(1,64:89)],LZ_mth_dat[,c(1,120:137)],mad_mth_dat[,c(1,87:103)],
                                                      wind_mth_dat[,c(1,19:23)],wash_mth_dat[,c(1,8:16)],leve_mth_dat[,c(1,6:7)],UZ_mth_dat[,c(1,62:79)],mon_mth_dat[,c(1,89:105)]),
                                                   FUN = function(x){
                                                     
                                                     x <- EWSmethods::detrend_ts(x,method = "loess",span = 0.75)                                           
                                                     # if(any(str_length(x[,1])>4)){
                                                     #   x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                                     #   x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                                     # }
                                                     time_vec <- colnames(x)[1]
                                                     
                                                     foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                                       
                                                       out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                                 winsize = 50)
                                                       
                                                       data.frame(cbind(data_source = paste(i),
                                                                        out$cor))
                                                     } 
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
    troph_level = "zooplankton") |>
  select(-name)

write.csv(genus_lake_zoo_roll_uni_ews,file = "Results/lake_results/genus_lake_zoo_roll_uni_ews.csv")


################################################################################################################
## Univariate EWS Assessment Permuted Rolling ##
################################################################################################################

cl <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl)
genus_lake_roll_uni_ews_perm <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
                                                       wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
                                                       kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
                                                       wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat),
                                            FUN = function(x){
                                              
                                              x <- EWSmethods::detrend_ts(subset(x, select=-c(pca1,pca2)),method = "loess",span = 0.75)                                           
                                              # if(any(str_length(x[,1])>4)){
                                              #   x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                              #   x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                              # }
                                              time_vec <- colnames(x)[1]
                                              
                                              foreach(i=colnames(x[,-1]), .combine='rbind',.export = "perm_rollEWS",.verbose = F) %dopar%{

                                                out <- perm_rollEWS(data = x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),
                                                                    winsize = 50, perm.meth = "replacement", variate = "uni",
                                                                    iter = 500)                    
                                                
                                                data.frame(cbind(data_source = paste(i), out$cor))
                                              } 
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
    computation = "rolling") |>
  select(-name)

parallel::stopCluster(cl)
write.csv(genus_lake_roll_uni_ews_perm,file = "Results/lake_results/genus_lake_roll_uni_ews_perm.csv")

cl <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl)
genus_lake_phyto_roll_uni_ews_perm <- pbapply::pblapply(list(kin_yr_dat[,c(1:28)],kas_yr_dat[,c(1:63)],LZ_yr_dat[,c(1:119)],mad_yr_dat[,c(1:87)],
                                                             wind_yr_dat[,c(1:18)],wash_yr_dat[,c(1:7)],leve_yr_dat[,c(1:5)],UZ_yr_dat[,c(1:61)],mon_yr_dat[,c(1:90)],
                                                             kin_mth_dat[,c(1:28)],kas_mth_dat[,c(1:63)],LZ_mth_dat[,c(1:119)],mad_mth_dat[,c(1:86)],
                                                             wind_mth_dat[,c(1:18)],wash_mth_dat[,c(1:7)],leve_mth_dat[,c(1:5)],UZ_mth_dat[,c(1:61)],mon_mth_dat[,c(1:87)]),
                                                 FUN = function(x){
                                                   
                                                   x <- EWSmethods::detrend_ts(x,method = "loess",span = 0.75)                                           
                                                   # if(any(str_length(x[,1])>4)){
                                                   #   x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                                   #   x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                                   # }
                                                   time_vec <- colnames(x)[1]
                                                   
                                                   foreach(i=colnames(x[,-1]), .combine='rbind',.export = "perm_rollEWS",.verbose = F) %dopar%{
                                                     
                                                     out <- perm_rollEWS(data = x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),
                                                                         winsize = 50, perm.meth = "replacement", variate = "uni",
                                                                         iter = 500)                    
                                                     
                                                     data.frame(cbind(data_source = paste(i), out$cor))
                                                   } 
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
    troph_level = "phytoplankton") |>
  select(-name)

parallel::stopCluster(cl)
write.csv(genus_lake_phyto_roll_uni_ews_perm,file = "Results/lake_results/genus_lake_phyto_roll_uni_ews_perm.csv")

cl <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl)
genus_lake_zoo_roll_uni_ews_perm <- pbapply::pblapply(list(kin_yr_dat[,c(1,29:50)],kas_yr_dat[,c(1,64:89)],LZ_yr_dat[,c(1,120:137)],mad_yr_dat[,c(1,88:105)],
                                                           wind_yr_dat[,c(1,19:23)],wash_yr_dat[,c(1,8:16)],leve_yr_dat[,c(1,6:8)],UZ_yr_dat[,c(1,62:79)],mon_yr_dat[,c(1,91:107)],
                                                           kin_mth_dat[,c(1,29:50)],kas_mth_dat[,c(1,64:89)],LZ_mth_dat[,c(1,120:137)],mad_mth_dat[,c(1,87:103)],
                                                           wind_mth_dat[,c(1,19:23)],wash_mth_dat[,c(1,8:16)],leve_mth_dat[,c(1,6:7)],UZ_mth_dat[,c(1,62:79)],mon_mth_dat[,c(1,89:105)]),
                                                        FUN = function(x){
                                                          
                                                          x <- EWSmethods::detrend_ts(x,method = "loess",span = 0.75)                                           
                                                          # if(any(str_length(x[,1])>4)){
                                                          #   x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                                          #   x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                                          # }
                                                          time_vec <- colnames(x)[1]
                                                          
                                                          foreach(i=colnames(x[,-1]), .combine='rbind',.export = "perm_rollEWS",.verbose = F) %dopar%{
                                                            
                                                            out <- perm_rollEWS(data = x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),
                                                                                winsize = 50, perm.meth = "replacement", variate = "uni",
                                                                                iter = 500)                    
                                                            
                                                            data.frame(cbind(data_source = paste(i), out$cor))
                                                          } 
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
    troph_level = "phytoplankton") |>
  select(-name)

parallel::stopCluster(cl)
write.csv(genus_lake_zoo_roll_uni_ews_perm,file = "Results/lake_results/genus_lake_zoo_roll_uni_ews_perm.csv")

################################################################################################################
## Univariate EWS Assessment Expanding ##
################################################################################################################

genus_lake_exp_uni_ews <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
                                                 wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
                                                 kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
                                                 wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat),
                                      FUN = function(x){
                                        
                                        x <- EWSmethods::detrend_ts(subset(x, select=-c(pca1,pca2)),method = "loess",span = 0.75)                                           
                                        # if(any(str_length(x[,1])>4)){
                                        #   x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                        #   x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                        # }
                                        time_vec <- colnames(x)[1]
                                        
                                        foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                          
                                          out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                    burn_in = ceiling(dim(x)[1]*0.50))
                                          
                                          data.frame(cbind(data_source = paste(i),
                                                           out))
                                        } 
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
    computation = "expanding") |>
  select(-name)

write.csv(genus_lake_exp_uni_ews,file = gzfile("Results/lake_results/genus_lake_exp_uni_ews.csv.gz"))

genus_lake_phyto_exp_uni_ews <- pbapply::pblapply(list(kin_yr_dat[,c(1:28)],kas_yr_dat[,c(1:63)],LZ_yr_dat[,c(1:119)],mad_yr_dat[,c(1:87)],
                                                 wind_yr_dat[,c(1:18)],wash_yr_dat[,c(1:7)],leve_yr_dat[,c(1:5)],UZ_yr_dat[,c(1:61)],mon_yr_dat[,c(1:90)],
                                                 kin_mth_dat[,c(1:28)],kas_mth_dat[,c(1:63)],LZ_mth_dat[,c(1:119)],mad_mth_dat[,c(1:86)],
                                                 wind_mth_dat[,c(1:18)],wash_mth_dat[,c(1:7)],leve_mth_dat[,c(1:5)],UZ_mth_dat[,c(1:61)],mon_mth_dat[,c(1:87)]),
                                            FUN = function(x){
                                              
                                              x <- EWSmethods::detrend_ts(x,method = "loess",span = 0.75)                                           
                                              # if(any(str_length(x[,1])>4)){
                                              #   x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                              #   x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                              # }
                                              time_vec <- colnames(x)[1]
                                              
                                              foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                                
                                                out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                          burn_in = ceiling(dim(x)[1]*0.50))
                                                
                                                data.frame(cbind(data_source = paste(i),
                                                                 out))
                                              } 
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
    troph_level = "phytoplankton") |>
  select(-name)

write.csv(genus_lake_phyto_exp_uni_ews,file = gzfile("Results/lake_results/genus_lake_phyto_exp_uni_ews.csv.gz"))

genus_lake_zoo_exp_uni_ews <- pbapply::pblapply(list(kin_yr_dat[,c(1,29:50)],kas_yr_dat[,c(1,64:89)],LZ_yr_dat[,c(1,120:137)],mad_yr_dat[,c(1,88:105)],
                                                     wind_yr_dat[,c(1,19:23)],wash_yr_dat[,c(1,8:16)],leve_yr_dat[,c(1,6:8)],UZ_yr_dat[,c(1,62:79)],mon_yr_dat[,c(1,91:107)],
                                                     kin_mth_dat[,c(1,29:50)],kas_mth_dat[,c(1,64:89)],LZ_mth_dat[,c(1,120:137)],mad_mth_dat[,c(1,87:103)],
                                                     wind_mth_dat[,c(1,19:23)],wash_mth_dat[,c(1,8:16)],leve_mth_dat[,c(1,6:7)],UZ_mth_dat[,c(1,62:79)],mon_mth_dat[,c(1,89:105)]),
                                                  FUN = function(x){
                                                    
                                                    x <- EWSmethods::detrend_ts(x,method = "loess",span = 0.75)                                           
                                                    # if(any(str_length(x[,1])>4)){
                                                    #   x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                                    #   x <- EWSmethods::deseason_ts(x,increment = "month", method = "stl",order = "ymd")
                                                    # }
                                                    time_vec <- colnames(x)[1]
                                                    
                                                    foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                                      
                                                      out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                                burn_in = ceiling(dim(x)[1]*0.50))
                                                      
                                                      data.frame(cbind(data_source = paste(i),
                                                                       out))
                                                    } 
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
    troph_level = "zooplankton") |>
  select(-name)

write.csv(genus_lake_zoo_exp_uni_ews,file = gzfile("Results/lake_results/genus_lake_zoo_exp_uni_ews.csv.gz"))

