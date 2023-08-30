require(dplyr)
require(EWSmethods)
require(foreach)
require(zoo)
require(stringr)

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

################################################################################################################
## Wangle Data ##
################################################################################################################

kin_mth_dat <- plank_env.data.mth |> 
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
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(zoo::as.yearmon(Date))) |>
  mutate(phyto_density = log1p(rowSums(across(Anabaena:Rhodomonas))),
         zoo_density = log1p(rowSums(across(`copepotide-i`:Anuraeopsis)))) |>
  as.data.frame()

kas_mth_dat <- plank_env.kasmthdata |> 
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(c(date,Acanthoceras_zachariasii:Thermocyclops_taihokuensis)) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
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
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  mutate(phyto_density = log1p(rowSums(across(Acanthoceras:Westella))),
         zoo_density = log1p(rowSums(across(Alona:Thermocyclops)))) |>
  as.data.frame()

LZ_mth_dat <- plank_env.LZmthdata |> 
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
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Actinastrum:Willea) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  mutate(date = base::as.Date(zoo::as.Date(date))) |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  mutate(phyto_density = log1p(rowSums(across(Microcystis:Heteronema))),
         zoo_density = log1p(rowSums(across(Alona:Polyphemus)))) |>
  as.data.frame()

mad_mth_dat <- plank_env.madmthdata |> 
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
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  mutate(phyto_density = log1p(rowSums(across(Achnanthes:Woronichinia))),
         zoo_density = log1p(rowSums(across(ACANTHOCYCLOPS:TROPOCYCLOPS)))) |>
  as.data.frame()

wind_mth_dat <- phyto_env.windmthdata |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
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
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(zoo::as.yearmon(Date))) |>
  mutate(phyto_density = log1p(rowSums(across(Asterionella:Monoraphidium))),
         zoo_density = log1p(rowSums(across(Bosmina:TotCyclopoids)))) |>
  as.data.frame()

UZ_mth_dat <- plank_env.UZmthdata |> 
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(date,Ankistrodesmus_sp.:Polyphemus_pediculus) |> #drop environmentals
  pivot_longer(-date,names_to = "species",values_to = "density") |>
  mutate(genus = gsub( "_.*$", "", species )) |> #then shrink to just genus
  mutate(key = forcats::fct_inorder(genus))|> #set order 
  group_by(date,key) |>
  summarise(density = sum(density)) |>
  rename("genus" = key) |>
  pivot_wider(names_from = "genus",values_from  = "density") |>
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
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  mutate(phyto_density = log1p(rowSums(across(Ankistrodesmus:Gomphosphaeria))),
         zoo_density = log1p(rowSums(across(Alona:Polyphemus)))) |>
  as.data.frame()

mon_mth_dat <- plank_env.monmthdata |> 
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
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  mutate(phyto_density = log1p(rowSums(across(Achnanthes:Woronichinia))),
         zoo_density = log1p(rowSums(across(ACANTHOCYCLOPS:TROPOCYCLOPS)))) |>
  as.data.frame()

leve_mth_dat <- plank_env.levemthdata |> 
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
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(zoo::as.yearmon(Date))) |>
  mutate(phyto_density = log1p(rowSums(across(Asterionella:Aulacoseira))),
         zoo_density = log1p(rowSums(across(Cyclops:Eudiaptomus)))) |>
  as.data.frame()

wash_mth_dat <- plank_env.washmthdata |> 
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  dplyr::select(c(date,cryptomonas:nc.rotifers)) |> #drop environmentals
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
  mutate(phyto_density = log1p(rowSums(across(cryptomonas:other.algae))),
         zoo_density = log1p(rowSums(across(conochilus:nc.rotifers)))) |>
  as.data.frame()

#############################################################################
#Prepare Lake Data
#############################################################################

state.kin.dat <- kin_mth_dat |>
  dplyr::select(-c(pca1,pca2)) |>
  dplyr::rename(date = Date) |>
  dplyr::filter(date >=1975) |>
  dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(date))) |>
  dplyr::mutate(dplyr::across(-date,~c(scale(log1p(.x)))))|>
  EWSmethods::deseason_ts(increment = "month",method = "average",order = "ymd") |>
  dplyr::mutate(date = as.numeric(zoo::as.yearmon(date))) |>
  dplyr::left_join(plank_env.data.mth |> dplyr::rename(date = Date) |>
                     dplyr::mutate(date = as.numeric(date)) |>
                     dplyr::filter(date >=1975) |>
                     tidyr::nest(data = everything())|>
                     dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                                      dplyr::select(c("Temp","Nitrate","P_ort")) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~zoo::na.approx(.x))) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~c(scale(.x)))))) |>
                     
                     dplyr::mutate(env =purrr::map(tmp, ~prcomp(.x)$x[,1])) |>
                     dplyr::select(-tmp)|>
                     tidyr::unnest(c(data,env)) |>
                     dplyr::select(date,env),
                   by="date") |>
  tidyr::pivot_longer(-c(date,env),names_to = "spp",values_to = "density") |>
  dplyr::mutate(lake = "Kinneret",
                troph_level = ifelse(spp %in% c("phyto_density","zoo_density"),NA,ifelse(
                  spp %in% colnames(kin_mth_dat)[2:30],
                  "Phytoplankton","Zooplankton")))

state.wind.dat <- wind_mth_dat |>
  dplyr::select(-c(pca1,pca2)) |>
  dplyr::rename(date = Date) |>
  dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(date))) |>
  dplyr::mutate(dplyr::across(-date,~c(scale(log1p(.x)))))|>
  EWSmethods::deseason_ts(increment = "month",method = "average",order = "ymd") |>
  dplyr::mutate(date = as.numeric(zoo::as.yearmon(date))) |>
  dplyr::left_join(phyto_env.windmthdata |> dplyr::rename(date = Date) |>
                     dplyr::mutate(date = as.numeric(date)) |>
                     tidyr::nest(data = everything())|>
                     dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                                      dplyr::select(c("TEMP","NO3N","TOTP")) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~zoo::na.approx(.x))) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~c(scale(.x)))))) |>
                     
                     dplyr::mutate(env =purrr::map(tmp, ~prcomp(.x)$x[,1])) |>
                     dplyr::select(-tmp)|>
                     tidyr::unnest(c(data,env)) |>
                     dplyr::select(date,env),
                   by="date") |>
  tidyr::pivot_longer(-c(date,env),names_to = "spp",values_to = "density") |>
  dplyr::mutate(lake = "Windermere",
                troph_level = ifelse(spp %in% c("phyto_density","zoo_density"),NA,ifelse(
                  spp %in% colnames(wind_mth_dat)[2:18],
                  "Phytoplankton","Zooplankton")))

state.kas.dat <- kas_mth_dat |>
  dplyr::select(-c(pca1,pca2)) |>
  dplyr::filter(date <2019) |>
  dplyr::mutate(dplyr::across(-date,~c(scale(log1p(.x)))))|>
  dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(date))) |>
  EWSmethods::deseason_ts(increment = "month",method = "average",order = "ymd") |>
  dplyr::mutate(date = as.numeric(zoo::as.yearmon(date))) |>
  dplyr::left_join(plank_env.kasmthdata |> 
                     dplyr::filter(date <2019) |>
                     dplyr::mutate(date = as.numeric(date)) |>
                     tidyr::nest(data = everything())|>
                     dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                                      dplyr::select(c("wtemp","NO3N","totP")) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~zoo::na.approx(.x))) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~c(scale(.x)))))) |>
                     
                     dplyr::mutate(env =purrr::map(tmp, ~prcomp(.x)$x[,1])) |>
                     dplyr::select(-tmp)|>
                     tidyr::unnest(c(data,env)) |>
                     dplyr::select(date,env),
                   by="date") |>
  tidyr::pivot_longer(-c(date,env),names_to = "spp",values_to = "density") |>
  dplyr::mutate(lake = "Kasumigaura",
                troph_level = ifelse(spp %in% c("phyto_density","zoo_density"),NA,ifelse(
                  spp %in% colnames(kas_mth_dat)[2:99],
                  "Phytoplankton","Zooplankton")))

state.leve.dat <- leve_mth_dat |>
  dplyr::select(-c(pca1,pca2)) |>
  dplyr::rename(date = Date) |>
  dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(date))) |>
  dplyr::mutate(dplyr::across(-date,~c(scale(log1p(.x)))))|>
  EWSmethods::deseason_ts(increment = "month",method = "average",order = "ymd") |>
  dplyr::mutate(date = as.numeric(zoo::as.yearmon(date))) |>
  dplyr::left_join(plank_env.levemthdata |> dplyr::rename(date = Date) |>
                     dplyr::mutate(date = as.numeric(date)) |>
                     tidyr::nest(data = everything())|>
                     dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                                      dplyr::select(c("Temp","TOTP")) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~zoo::na.approx(.x))) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~c(scale(.x)))))) |>
                     
                     dplyr::mutate(env =purrr::map(tmp, ~prcomp(.x)$x[,1])) |>
                     dplyr::select(-tmp)|>
                     tidyr::unnest(c(data,env)) |>
                     dplyr::select(date,env),
                   by="date") |>
  tidyr::pivot_longer(-c(date,env),names_to = "spp",values_to = "density") |>
  dplyr::mutate(lake = "Loch Leven",
                troph_level = ifelse(spp %in% c("phyto_density","zoo_density"),NA,ifelse(
                  spp %in% colnames(leve_mth_dat)[2:5],
                  "Phytoplankton","Zooplankton")))

state.LZ.dat <- LZ_mth_dat |>
  dplyr::select(-c(pca1,pca2)) |>
  dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(date))) |>
  dplyr::mutate(dplyr::across(-date,~c(scale(log1p(.x)))))|>
  EWSmethods::deseason_ts(increment = "month",method = "average",order = "ymd") |>
  dplyr::mutate(date = as.numeric(zoo::as.yearmon(date))) |>
  dplyr::left_join(plank_env.LZmthdata |> 
                     dplyr::mutate(date = as.numeric(date)) |>
                     tidyr::nest(data = everything())|>
                     dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                                      dplyr::select(c("mean.t","mean.po4","NO3_N")) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~zoo::na.approx(.x))) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~c(scale(.x)))))) |>
                     
                     dplyr::mutate(env =purrr::map(tmp, ~prcomp(.x)$x[,1])) |>
                     dplyr::select(-tmp)|>
                     tidyr::unnest(c(data,env)) |>
                     dplyr::select(date,env),
                   by="date") |>
  dplyr::mutate(env = zoo::na.approx(env)) |> #for singular missing years
  tidyr::pivot_longer(-c(date,env),names_to = "spp",values_to = "density") |>
  dplyr::mutate(lake = "Lower Zurich",
                troph_level = ifelse(spp %in% c("phyto_density","zoo_density"),NA,ifelse(
                  spp %in% colnames(LZ_mth_dat)[2:122],
                  "Phytoplankton","Zooplankton")))

state.UZ.dat <- UZ_mth_dat |>
  dplyr::select(-c(pca1,pca2)) |>
  dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(date))) |>
  dplyr::mutate(dplyr::across(-date,~c(scale(log1p(.x)))))|>
  EWSmethods::deseason_ts(increment = "month",method = "average",order = "ymd") |>
  dplyr::mutate(date = as.numeric(zoo::as.yearmon(date))) |>
  dplyr::left_join(plank_env.UZmthdata |> 
                     dplyr::mutate(date = as.numeric(date)) |>
                     tidyr::nest(data = everything())|>
                     dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                                      dplyr::select(c("mean.t","mean.po4")) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~zoo::na.approx(.x))) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~c(scale(.x)))))) |>
                     
                     dplyr::mutate(env =purrr::map(tmp, ~prcomp(.x)$x[,1])) |>
                     dplyr::select(-tmp)|>
                     tidyr::unnest(c(data,env)) |>
                     dplyr::select(date,env),
                   by="date") |>
  dplyr::mutate(env = zoo::na.approx(env)) |> #for singular missing years
  tidyr::pivot_longer(-c(date,env),names_to = "spp",values_to = "density") |>
  dplyr::mutate(lake = "Upper Zurich",
                troph_level = ifelse(spp %in% c("phyto_density","zoo_density"),NA,ifelse(
                  spp %in% colnames(UZ_mth_dat)[2:61],
                  "Phytoplankton","Zooplankton")))

state.mad.dat <- mad_mth_dat |>
  dplyr::select(-c(pca1,pca2)) |>
  dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(date))) |>
  dplyr::mutate(dplyr::across(-date,~c(scale(log1p(.x)))))|>
  EWSmethods::deseason_ts(increment = "month",method = "average",order = "ymd") |>
  dplyr::mutate(date = as.numeric(zoo::as.yearmon(date))) |>
  dplyr::left_join(plank_env.madmthdata |> 
                     dplyr::mutate(date = as.numeric(date)) |>
                     tidyr::nest(data = everything())|>
                     dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                                      dplyr::select(c("wtemp","totP","NO3N")) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~zoo::na.approx(.x))) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~c(scale(.x)))))) |>
                     
                     dplyr::mutate(env =purrr::map(tmp, ~prcomp(.x)$x[,1])) |>
                     dplyr::select(-tmp)|>
                     tidyr::unnest(c(data,env)) |>
                     dplyr::select(date,env),
                   by="date") |>
  tidyr::pivot_longer(-c(date,env),names_to = "spp",values_to = "density") |>
  dplyr::mutate(lake = "Mendota",
                troph_level = ifelse(spp %in% c("phyto_density","zoo_density"),NA,ifelse(
                  spp %in% colnames(mad_mth_dat)[2:98],
                  "Phytoplankton","Zooplankton")))

state.mon.dat <- mon_mth_dat |>
  dplyr::select(-c(pca1,pca2)) |>
  dplyr::filter(date <2018) |>
  dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(date))) |>
  dplyr::mutate(dplyr::across(-date,~c(scale(log1p(.x)))))|>
  EWSmethods::deseason_ts(increment = "month",method = "average",order = "ymd") |>
  dplyr::mutate(date = as.numeric(zoo::as.yearmon(date))) |>
  dplyr::left_join(plank_env.monmthdata |>  dplyr::filter(date <2018) |>
                     dplyr::mutate(date = as.numeric(date)) |>
                     tidyr::nest(data = everything())|>
                     dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                                      dplyr::select(c("wtemp","totP","NO3N")) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~zoo::na.approx(.x))) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~c(scale(.x)))))) |>
                     
                     dplyr::mutate(env =purrr::map(tmp, ~prcomp(.x)$x[,1])) |>
                     dplyr::select(-tmp)|>
                     tidyr::unnest(c(data,env)) |>
                     dplyr::select(date,env),
                   by="date") |>
  tidyr::pivot_longer(-c(date,env),names_to = "spp",values_to = "density") |>
  dplyr::mutate(lake = "Monona",
                troph_level = ifelse(spp %in% c("phyto_density","zoo_density"),NA,ifelse(
                  spp %in% colnames(mon_mth_dat)[2:91],
                  "Phytoplankton","Zooplankton")))

state.wash.dat <- wash_mth_dat |>
  dplyr::select(-c(pca1,pca2)) |>
  dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(date))) |>
  dplyr::mutate(dplyr::across(-date,~c(scale(log1p(.x)))))|>
  EWSmethods::deseason_ts(increment = "month",method = "average",order = "ymd") |>
  dplyr::mutate(date = as.numeric(zoo::as.yearmon(date))) |>
  dplyr::left_join(plank_env.washmthdata |>  
                     dplyr::mutate(date = as.numeric(date)) |>
                     tidyr::nest(data = everything())|>
                     dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                                      dplyr::select(c("temp","TP")) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~zoo::na.approx(.x))) |>
                                                      dplyr::mutate(dplyr::across(dplyr::everything(),~c(scale(.x)))))) |>
                     
                     dplyr::mutate(env =purrr::map(tmp, ~prcomp(.x)$x[,1])) |>
                     dplyr::select(-tmp)|>
                     tidyr::unnest(c(data,env)) |>
                     dplyr::select(date,env),
                   by="date") |>
  tidyr::pivot_longer(-c(date,env),names_to = "spp",values_to = "density") |>
  dplyr::mutate(lake = "Washington",
                troph_level = ifelse(spp %in% c("phyto_density","zoo_density"),NA,ifelse(
                  spp %in% colnames(wash_mth_dat)[2:7],
                  "Phytoplankton","Zooplankton")))

#############################################################################
#Perform State Space Fitting - Environmentals
#############################################################################

lake_state_space_mth <- foreach::foreach(j = list(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                                               state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat), 
                                      .combine = "rbind") %do% {
                                        
                                        lapply(c("phyto_density","zoo_density"), 
                                               FUN = function(i){
                                                
                                                 cont_formula <-  as.formula(paste("density ~ s(env, bs='tp', k=6)"))
                                                 thresh_formula <- as.formula(paste("density ~ 1"))
                                                 sub_dat <- subset(j,spp == i) |>
                                                   tidyr::drop_na(density)
                                                 
                                                 dens_gam <- compare_gam(data = sub_dat,
                                                                         cont_formula = formula(cont_formula, method = "REML"),
                                                                         thresh_formula =  formula(thresh_formula, method = "REML"),
                                                                         thresh.var = "date", expl.var = "env",
                                                                         thresh.range = c(0.25,0.75), by = 1/12, k=3)
                                                 
                                                 best_gam <- predict_best_gam(object=dens_gam) |> 
                                                   dplyr::select(density,threshold,thresh.var,env,thresh.var,
                                                                 transition,fit,ci) |>
                                                   dplyr::mutate("spp" = i,
                                                                 "troph_level" = sub_dat$troph_level[1]) 
                                                 
                                                 
                                                 return(best_gam)
                                                 
                                               })  |> 
                                          data.table::rbindlist() |>
                                          dplyr::group_by(spp) |>
                                          dplyr::mutate(threshold = tidyr::replace_na(as.character(threshold),"pre"),
                                                        start.date = dplyr::first(thresh.var),
                                                        last.date = dplyr::last(thresh.var),
                                                        lake = j$lake[1] 
                                          )
                                      }

ss_plot_dat <- subset(lake_state_space_mth, spp %in% c("phyto_density","zoo_density")) |>
  dplyr::mutate(spp = factor(spp, levels = c("phyto_density","zoo_density"),labels = c("Phytoplankton","Zooplankton")))


ggplot(data = ss_plot_dat, aes(x=env,y=density)) + 
  geom_point(aes(x=env, y = density))+
  geom_path(aes(x=env, y = density)) +  
  xlab("Environmental stressor") + ylab("Scaled metric score")+ 
  theme_bw() + 
  facet_grid(spp~lake,scales = "free_y") +
  geom_line(data = filter(ss_plot_dat, threshold=="pre"),aes(x=env,y=fit), col="blue",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_plot_dat, threshold=="pre"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(ss_plot_dat, threshold=="post"),aes(x=env,y=fit), col="red",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_plot_dat, threshold=="post"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  geom_point(data = filter(ss_plot_dat, thresh.var == start.date),aes(x=env, y = density,col="Start date"))+
  geom_point(data = filter(ss_plot_dat, thresh.var == last.date),aes(x=env, y = density,col="End date"))+
  ggrepel::geom_text_repel(data = filter(ss_plot_dat,transition == "trans"), aes(x=env, y = density,label=thresh.var),force =1,nudge_x=3,nudge_y=-0.5,segment.linetype=2)+
  geom_point(data = filter(ss_plot_dat,transition == "trans"),aes(x=env, y = density,col="Transition\ndates"))+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates","1993"), name = NULL)+  
  xlim(-6,6)

#############################################################################
#Perform State Space Fitting - Temporal
#############################################################################

lake_temporal_mth <- foreach::foreach(j = list(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                                           state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat), 
                                  .combine = "rbind") %do% {
                                    
                                lapply(c("phyto_density","zoo_density"), 
                                           FUN = function(i){
                                             
                                             cont_formula <-  as.formula(paste("density ~ s(date, bs='tp', k=6)"))
                                             thresh_formula <- as.formula(paste("density ~ 1"))
                                             sub_dat <- subset(j,spp == i) |>
                                               tidyr::drop_na(density)
                                             
                                             dens_gam <- compare_gam(data = sub_dat,
                                                                     cont_formula = formula(cont_formula, method = "REML"),
                                                                     thresh_formula =  formula(thresh_formula, method = "REML"),
                                                                     thresh.var = "date",
                                                                     thresh.range = c(0.25,0.75), by = 1/12, k=3)
                                             
                                             best_gam <- predict_best_gam(object=dens_gam) |> 
                                               dplyr::select(density,threshold,date,thresh.var,
                                                             transition,fit,ci) |>
                                               dplyr::mutate("spp" = i,
                                                             "troph_level" = sub_dat$troph_level[1]) 
                                               
                                             
                                             return(best_gam)
                                             
                                           })  |> 
                                     data.table::rbindlist() |>
                                     dplyr::group_by(spp) |>
                                      dplyr::mutate(threshold = tidyr::replace_na(as.character(threshold),"pre"),
                                                    start.date = dplyr::first(thresh.var),
                                                    last.date = dplyr::last(thresh.var),
                                                    lake = j$lake[1] 
                                      )
                                  }

ss_time_plot_dat <- subset(lake_temporal_mth, spp %in% c("phyto_density","zoo_density")) |>
  dplyr::mutate(spp = factor(spp, levels = c("phyto_density","zoo_density"),labels = c("Phytoplankton","Zooplankton")))


ggplot(data = ss_time_plot_dat, aes(x=date,y=density)) + 
  geom_point(aes(x=date, y = density))+
  geom_path(aes(x=date, y = density)) +  
  xlab("Year") + ylab("Scaled density score")+ 
  theme_bw() + 
  facet_grid(spp~lake,scales = "free_y") +
  geom_line(data = filter(ss_time_plot_dat, threshold=="pre"),aes(x=date,y=fit), col="blue",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_time_plot_dat, threshold=="pre"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#A1B4FE", col="#A1B4FE",alpha = 0.2)+
  geom_line(data = filter(ss_time_plot_dat, threshold=="post"),aes(x=date,y=fit), col="red",size=0.8, linetype = "solid")+
  geom_ribbon(data = filter(ss_time_plot_dat, threshold=="post"),aes(ymin = fit - (1.96 * ci),ymax = fit + (1.96 * ci)  ), fill = "#FFA6B9", col="#FFA6B9",alpha = 0.2)+
  ggrepel::geom_text_repel(data = filter(ss_time_plot_dat,transition == "trans"), aes(x=date, y = density,label=thresh.var),force =1,nudge_x=3,nudge_y=-0.5,segment.linetype=2)+
  geom_point(data = filter(ss_time_plot_dat,date == start.date),aes(x=date, y = density,col="Start date"))+
  geom_point(data = filter(ss_time_plot_dat,date == last.date),aes(x=date, y = density,col="End date"))+
  geom_point(data = filter(ss_time_plot_dat,transition == "trans"),aes(x=date, y = density,col="Transition\ndates"))+
  scale_colour_manual(values = c("blue","red","#FFE823"),breaks = c("Start date", "End date", "Transition\ndates","1993"), name = NULL) +
  scale_x_continuous(breaks=seq(1970,2015,by=15),limits = c(1960,2020))

#############################################################################
#Identify bi-modality
#############################################################################
require(LaplacesDemon)
require(diptest)
require(mousetrap)

lake_bimod_mth <- foreach::foreach(j = list(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                                        state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat), 
                               .combine = "rbind") %do% {
                                 
                                 lapply(c("phyto_density","zoo_density"),function(i){
                                   
                                   sub_dat <- subset(j,spp == i) |>
                                     tidyr::drop_na(density)

                                   modes <- LaplacesDemon::Modes(sub_dat[["density"]])$modes
                                   
                                   return(data.frame("spp" = i,
                                                     "troph_level" = sub_dat$troph_level[1],
                                                     "is_bimodal" =  LaplacesDemon::is.bimodal(sub_dat[["density"]]),
                                                     "mode1" = modes[1],
                                                     "mode2" = ifelse(length(modes)==1,NA,modes[2]),
                                                     "modality_coef" = mousetrap::bimodality_coefficient(sub_dat[["density"]])))

                                 }) |>
                                   data.table::rbindlist() |>
                                   dplyr::mutate(lake = j$lake[1])
                               }

bimod_plot_datS1 <- rbind(state.kas.dat,state.kin.dat,state.leve.dat,state.LZ.dat,
                          state.mad.dat,state.mon.dat,state.UZ.dat,state.wash.dat,state.wind.dat) |>
  dplyr::select(lake,spp,density) |>
  dplyr::filter(spp %in% c("phyto_density","zoo_density")) |>
  dplyr::mutate(spp = factor(spp, levels = c("phyto_density","zoo_density"),labels = c("Phytoplankton","Zooplankton")))

bimod_modes_mth <- lake_bimod_mth |>
  dplyr::mutate(spp = factor(spp, levels = c("phyto_density","zoo_density"),labels = c("Phytoplankton","Zooplankton"))) |>
  dplyr::group_by(spp) |>
  dplyr::mutate(mode1 = round(mode1,digits = 2),
                mode2 = round(mode2,digits = 2)
               )

ggplot(bimod_plot_datS1) +
  geom_density(aes(y = density)) +
  geom_hline(data = bimod_modes_mth,
             aes(yintercept = mode1),linetype = "dashed",colour = "black")+
  geom_hline(data = bimod_modes_mth,
             aes(yintercept = mode2),linetype = "dashed",colour = "black")+
  facet_grid(spp~lake,scales = "free_y") +
  #ylim(c(-1,1)) +
  theme_classic()


#############################################################################
#Compare State Spaces
#############################################################################

state_transition_mth_dates <- subset(lake_state_space_mth,transition == "trans") |>
  rename(state_date = thresh.var)
temporal_transition_mth_dates <- subset(lake_temporal_mth,transition == "trans") |>
  rename(temporal_date = date)
bimod_mth_dates <- lake_bimod_mth |>
  ungroup() |>
  dplyr::mutate(is_bimodal = ifelse(modality_coef > 0.5 & all(!is.na(c(mode1,mode2))),TRUE,FALSE))

transition_mth_dates <-expand.grid(spp = unique(c(lake_state_space_mth$spp,lake_temporal_mth$spp)),
                               lake = unique(c(lake_state_space_mth$lake,lake_temporal_mth$lake))) |>
  left_join(state_transition_mth_dates |>
              select(spp,lake,state_date),by = c("lake","spp"),
            relationship = "many-to-many") |>
  left_join(temporal_transition_mth_dates |>
              select(spp,lake,temporal_date),by = c("lake","spp"),
            relationship = "many-to-many") |>
  group_by(spp,lake) |>
  slice_head(n=1) |>
  mutate(date_match = ifelse(state_date == temporal_date, TRUE, FALSE),
         threshold_date = ifelse(isTRUE(date_match),state_date,NA)) |>
  left_join(bimod_mth_dates,by = c("lake","spp"))

write.csv(transition_mth_dates,"Results/supplementary_info/monthly_transition_dates.csv")

