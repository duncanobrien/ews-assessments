##########################################################################################
## Preamble ##
##########################################################################################

require(Matrix)
require(tidyverse)

transition_dates <- read.csv("Data/transition_dates.csv")[,-1] |>
  filter(metric %in% c("phyto_density","zoo_density"))

################################################################################################################
## Load Data ##
################################################################################################################

source("Data/kinneret_plankton_data.R") #raw data not provided. See references in Data Availability statement and repository README
source("Data/kinneret_environmental_data.R")
#phyo, zoo + combined raw concentrations at week, month and year groupings

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

################################################################################################################
## Wangle Data ##
################################################################################################################

kin_yr_dat <- plank_env.data.yr |> #drop environmentals
  dplyr::filter(as.numeric(Date) <= as.numeric(na.omit(transition_dates$threshold_date[transition_dates$lake == "Kinneret"]))) |>
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
  #mutate(across(Anabaena:Thermocyclops,~log1p(.x))) |>
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
  dplyr::filter(as.numeric(Date) <= as.numeric(na.omit(transition_dates$threshold_date[transition_dates$lake == "Kinneret"]))) |>
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
  #mutate(across(Anabaena:Thermocyclops,~log1p(.x))) |>
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
  dplyr::filter(as.numeric(date) <= as.numeric(na.omit(transition_dates$threshold_date[transition_dates$lake == "Kasumigaura"]))) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date:Thermocyclops_taihokuensis) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  #mutate(across(Acanthoceras:Thermocyclops,~log1p(.x))) |>
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
  dplyr::filter(as.numeric(date) <= as.numeric(na.omit(transition_dates$threshold_date[transition_dates$lake == "Kasumigaura"]))) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(c(date,Acanthoceras_zachariasii:Thermocyclops_taihokuensis)) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
  #mutate(across(Acanthoceras:Thermocyclops,~log1p(.x))) |>
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
  ungroup() |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  #mutate(across(Actinastrum:Willea,~log1p(.x))) |>
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
  #mutate(across(Actinastrum:Willea,~log1p(.x))) |>
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
  #mutate(across(ACANTHOCYCLOPS:Uroglena,~log1p(.x))) |>
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
  #mutate(across(ACANTHOCYCLOPS:Uroglena,~log1p(.x))) |>
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
  #mutate(across(Asterionella:TotCyclopoids,~log1p(.x))) |>
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
  #mutate(across(Asterionella:TotCyclopoids,~log1p(.x))) |>
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
  #mutate(across(Achnanthes:Volvox,~log1p(.x))) |>
  ungroup() |>
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
  #mutate(across(Achnanthes:Volvox,~log1p(.x))) |>
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
  dplyr::filter(as.numeric(date) <= as.numeric(na.omit(transition_dates$threshold_date[transition_dates$lake == "Monona"])))|>
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
  #mutate(across(ACANTHOCYCLOPS:Woronichinia,~log1p(.x))) |>
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
  dplyr::filter(as.numeric(date) <= as.numeric(na.omit(transition_dates$threshold_date[transition_dates$lake == "Monona"])))|>
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
  #mutate(across(ACANTHOCYCLOPS:Woronichinia,~log1p(.x))) |>
  ungroup() |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Achnanthes:TROPOCYCLOPS) %>%
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
  #mutate(across(Asterionella:Unicellular,~log1p(.x))) |>
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
 # mutate(across(Asterionella:Unicellular,~log1p(.x))) |>
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
  dplyr::filter(as.numeric(date) <= as.numeric(na.omit(transition_dates$threshold_date[transition_dates$lake == "Washington"]))) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(c(date,cryptomonas:nc.rotifers)) |> #drop environmentals
  #mutate(across(cryptomonas:nc.rotifers,~log1p(.x))) |>
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
  dplyr::filter(as.numeric(date) <= as.numeric(na.omit(transition_dates$threshold_date[transition_dates$lake == "Washington"]))) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(c(date,cryptomonas:nc.rotifers)) |> #drop environmentals
  #mutate(across(cryptomonas:nc.rotifers,~log1p(.x))) |>
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
## Save out ##
################################################################################################################

save(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
     kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat,
    file =  "Data/wrangled_genus_plank_data.Rdata" )
