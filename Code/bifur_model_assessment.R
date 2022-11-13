########################################################################################
## Preamble ##
########################################################################################
bifur.mod.scripts <- list.files(path ="Code/bifur_models", pattern="*.R",full.names = T) 
purrr::walk(bifur.mod.scripts, source) # source silently

load("Data/wrangled_genus_plank_data.Rdata")
load("Data/whole_wrangled_genus_plank_data.Rdata")

kin_density <- data.frame(seq_along(as.numeric(kin_mth_dat$Date[kin_mth_dat$Date >1980 & kin_mth_dat$Date < 1994])),rowSums(kin_mth_dat[kin_mth_dat$Date >1980 & kin_mth_dat$Date < 1994,2:50]))
kin_density_yr <- data.frame(seq_along(as.numeric(kin_yr_dat$Date[kin_yr_dat$Date <1994])),rowSums(kin_yr_dat[kin_yr_dat$Date<1994,2:50]))

wind_density <- data.frame(seq_along(wind_mth_dat$Date[wind_mth_dat$Date > 1984]),rowSums(wind_mth_dat[wind_mth_dat$Date > 1984,2:23]))
wind_density_yr <- data.frame(seq_along(wind_yr_dat$Date[wind_yr_dat$Date > 1987]),rowSums(wind_yr_dat[wind_yr_dat$Date > 1987,2:23]))

mad_density <- data.frame(seq_along(mad_mth_dat$date),rowSums(mad_mth_dat[,2:103]))
mad_density_yr <- data.frame(seq_along(mad_yr_dat$date),rowSums(mad_yr_dat[,2:103]))

kas_density <- data.frame(seq_along(kas_mth_dat$date[kas_yr_dat$date<1998]),rowSums(kas_mth_dat[kas_yr_dat$date<1998,2:89]))
kas_density_yr <- data.frame(seq_along(kas_yr_dat$date[kas_yr_dat$date<1998]),rowSums(kas_yr_dat[kas_yr_dat$date<1998,2:89]))

LZ_density <- data.frame(seq_along(LZ_mth_dat$date),rowSums(LZ_mth_dat[,2:137]))
LZ_density_yr <- data.frame(seq_along(LZ_yr_dat$date),rowSums(LZ_yr_dat[,2:137]))

UZ_density <- data.frame(seq_along(UZ_mth_dat$date[UZ_yr_dat$date<1993]),rowSums(UZ_mth_dat[UZ_yr_dat$date<1993,2:79]))
UZ_density_yr <- data.frame(seq_along(UZ_yr_dat$date[UZ_yr_dat$date<1993]),rowSums(UZ_yr_dat[UZ_yr_dat$date<1993,2:79]))

leve_density <- data.frame(seq_along(leve_mth_dat$Date),rowSums(leve_mth_dat[,2:8]))
leve_density_yr <- data.frame(seq_along(leve_yr_dat$Date),rowSums(leve_yr_dat[,2:8]))

mon_density <- data.frame(seq_along(mon_mth_dat$date[mon_yr_dat$date<2014]),rowSums(mon_mth_dat[mon_yr_dat$date<2014,2:107]))
mon_density_yr <- data.frame(seq_along(mon_yr_dat$date[mon_yr_dat$date<2014]),rowSums(mon_yr_dat[mon_yr_dat$date<2014,2:107]))

wash_density <- data.frame(seq_along(wash_mth_dat$date[wash_yr_dat$date<1976]),rowSums(wash_mth_dat[wash_yr_dat$date<1976,2:16]))
wash_density_yr <- data.frame(seq_along(wash_yr_dat$date[wash_yr_dat$date<1976]),rowSums(wash_yr_dat[wash_yr_dat$date<1976,2:16]))

####################################################################################################
## Observed Fits
####################################################################################################
kinOU <- fit_bif_mod(kin_density,model = "OU", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
kinSN <- fit_bif_mod(kin_density,model = "LSN", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
kinTC <- fit_bif_mod(kin_density,model = "LTC", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
kinsuperPF <- fit_bif_mod(kin_density,model = "superLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
kinsubPF <- fit_bif_mod(kin_density,model = "subLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))

kin.observedSN <- -2*(kinOU$loglik - kinSN$loglik)
kin.observedTC<- -2*(kinOU$loglik - kinTC$loglik)
kin.observedsuperPF <- -2*(kinOU$loglik - kinsuperPF$loglik)
kin.observedsubPF <- -2*(kinOU$loglik - kinsubPF$loglik)
knitr::kable(cbind(kin.observedSN,kin.observedTC,kin.observedsuperPF,kin.observedsubPF))

windOU <- fit_bif_mod(wind_density,model = "OU", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
windSN <- fit_bif_mod(wind_density,model = "LSN", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
windTC <- fit_bif_mod(wind_density,model = "LTC", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
windsuperPF <- fit_bif_mod(wind_density,model = "superLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
windsubPF <- fit_bif_mod(wind_density,model = "subLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))

wind.observedSN <- -2*(windOU$loglik - windSN$loglik)
wind.observedTC<- -2*(kinOU$loglik - windTC$loglik)
wind.observedsuperPF <- -2*(windOU$loglik - windsuperPF$loglik)
wind.observedsubPF <- -2*(windOU$loglik - windsubPF$loglik)
knitr::kable(cbind(wind.observedSN,wind.observedTC,wind.observedsuperPF,wind.observedsubPF))

#, method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0)
madOU <- fit_bif_mod(mad_density,model = "OU", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
madSN <- fit_bif_mod(mad_density,model = "LSN", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
madTC <- fit_bif_mod(mad_density,model = "LTC", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
madsuperPF <- fit_bif_mod(mad_density,model = "superLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
madsubPF <- fit_bif_mod(mad_density,model = "subLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))

mad.observedSN <- -2*(madOU$loglik - madSN$loglik)
mad.observedTC<- -2*(madOU$loglik - madTC$loglik)
mad.observedsuperPF <- -2*(madOU$loglik - madsuperPF$loglik)
mad.observedsubPF <- -2*(madOU$loglik - madsubPF$loglik)
knitr::kable(cbind(mad.observedSN,mad.observedTC,mad.observedsuperPF,mad.observedsubPF))

kasOU <- fit_bif_mod(kas_density,model = "OU", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
kasSN <- fit_bif_mod(kas_density,model = "LSN", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
kasTC <- fit_bif_mod(kas_density,model = "LTC", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
kassuperPF <- fit_bif_mod(kas_density,model = "superLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
kassubPF <- fit_bif_mod(kas_density,model = "subLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))

kas.observedSN <- -2*(kasOU$loglik - kasSN$loglik)
kas.observedTC<- -2*(kasOU$loglik - kasTC$loglik)
kas.observedsuperPF <- -2*(kasOU$loglik - kassuperPF$loglik)
kas.observedsubPF <- -2*(kasOU$loglik - kassubPF$loglik)
knitr::kable(cbind(kas.observedSN,kas.observedTC,kas.observedsuperPF,kas.observedsubPF))

LZOU <- fit_bif_mod(LZ_density,model = "OU", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
LZSN <- fit_bif_mod(LZ_density,model = "LSN", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
LZTC <- fit_bif_mod(LZ_density,model = "LTC", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
LZsuperPF <- fit_bif_mod(LZ_density,model = "superLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
LZsubPF <- fit_bif_mod(LZ_density,model = "subLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))

LZ.observedSN <- -2*(LZOU$loglik - LZSN$loglik)
LZ.observedTC<- -2*(LZOU$loglik - LZTC$loglik)
LZ.observedsuperPF <- -2*(LZOU$loglik - LZsuperPF$loglik)
LZ.observedsubPF <- -2*(LZOU$loglik - LZsubPF$loglik)
knitr::kable(cbind(LZ.observedSN,LZ.observedTC,LZ.observedsuperPF,LZ.observedsubPF))

UZOU <- fit_bif_mod(UZ_density,model = "OU", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
UZSN <- fit_bif_mod(UZ_density,model = "LSN", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
UZTC <- fit_bif_mod(UZ_density,model = "LTC", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
UZsuperPF <- fit_bif_mod(UZ_density,model = "superLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
UZsubPF <- fit_bif_mod(UZ_density,model = "subLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))

UZ.observedSN <- -2*(UZOU$loglik - UZSN$loglik)
UZ.observedTC<- -2*(UZOU$loglik - UZTC$loglik)
UZ.observedsuperPF <- -2*(UZOU$loglik - UZsuperPF$loglik)
UZ.observedsubPF <- -2*(UZOU$loglik - UZsubPF$loglik)
knitr::kable(cbind(UZ.observedSN,UZ.observedTC,UZ.observedsuperPF,UZ.observedsubPF))

washOU <- fit_bif_mod(wash_density,model = "OU", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
washSN <- fit_bif_mod(wash_density,model = "LSN", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
washTC <- fit_bif_mod(wash_density,model = "LTC", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
washsuperPF <- fit_bif_mod(wash_density,model = "superLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
washsubPF <- fit_bif_mod(wash_density,model = "subLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))

wash.observedSN <- -2*(washOU$loglik - washSN$loglik)
wash.observedTC<- -2*(washOU$loglik - washTC$loglik)
wash.observedsuperPF <- -2*(washOU$loglik - washsuperPF$loglik)
wash.observedsubPF <- -2*(washOU$loglik - washsubPF$loglik)
knitr::kable(cbind(wash.observedSN,wash.observedTC,wash.observedsuperPF,wash.observedsubPF))

leveOU <- fit_bif_mod(mon_density,model = "OU", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
leveSN <- fit_bif_mod(mon_density,model = "LSN", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
leveTC <- fit_bif_mod(mon_density,model = "LTC", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
levesuperPF <- fit_bif_mod(mon_density,model = "superLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
levesubPF <- fit_bif_mod(mon_density,model = "subLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))

leve.observedSN <- -2*(leveOU$loglik - leveSN$loglik)
leve.observedTC<- -2*(leveOU$loglik - leveTC$loglik)
leve.observedsuperPF <- -2*(leveOU$loglik - levesuperPF$loglik)
leve.observedsubPF <- -2*(leveOU$loglik - levesubPF$loglik)
knitr::kable(cbind(leve.observedSN,leve.observedTC,leve.observedsuperPF,leve.observedsubPF))

monOU <- fit_bif_mod(mon_density,model = "OU", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
monSN <- fit_bif_mod(mon_density,model = "LSN", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
monTC <- fit_bif_mod(mon_density,model = "LTC", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
monsuperPF <- fit_bif_mod(mon_density,model = "superLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))
monsubPF <- fit_bif_mod(mon_density,model = "subLPF", method ="L-BFGS-B", lower = c(-Inf,0,-Inf,0))

mon.observedSN <- -2*(monOU$loglik - monSN$loglik)
mon.observedTC<- -2*(monOU$loglik - monTC$loglik)
mon.observedsuperPF <- -2*(monOU$loglik - monsuperPF$loglik)
mon.observedsubPF <- -2*(monOU$loglik - monsubPF$loglik)
knitr::kable(cbind(mon.observedSN,mon.observedTC,mon.observedsuperPF,mon.observedsubPF))

####################################################################################################
## Simulate Replicates
####################################################################################################
kin.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(kinOU, kinSN)},mc.cores=9)
kin.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(kinOU, kinTC)},mc.cores=9)
kin.repssuperPF <- pbmcapply::pbmclapply(1:500, function(i){compare(kinOU, kinsuperPF)},mc.cores=9)

kin.lrSN <- lik_ratios(kin.repsSN) |> dplyr::mutate(model = "SN")
kin.lrTC <- lik_ratios(kin.repsTC) |> dplyr::mutate(model = "TC")
kin.lrPF <- lik_ratios(kin.repssuperPF) |> dplyr::mutate(model = "PF")

kin.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                           "obs" = c(-2*(kinOU$loglik - kinSN$loglik),-2*(kinOU$loglik - kinTC$loglik),-2*(kinOU$loglik - kinsuperPF$loglik)))

kin_lr <- dplyr::left_join(rbind(kin.lrSN,kin.lrTC,kin.lrPF ),kin.observed.dat,by="model") |>
  dplyr::mutate(lake = "Kinneret")

ggplot(kin_lr) + 
  geom_density(aes(value, fill=simulation), alpha=0.8,size=0.3) + 
  geom_vline(aes(xintercept=obs),linetype = "dashed",col="grey20") + 
  scale_fill_manual(values = c(c("#A1B4FE","#FFE7A1")))+
  facet_wrap(~model,scales = "free_x") +theme_bw() +
  xlim(-50,50)


wind.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(windOU, windSN)},mc.cores=9)
wind.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(windOU, windTC)},mc.cores=9)
wind.repssuperPF <- pbmcapply::pbmclapply(1:500, function(i){compare(windOU, windsuperPF)},mc.cores=9)

wind.lrSN <- lik_ratios(wind.repsSN) |> dplyr::mutate(model = "SN")
wind.lrTC <- lik_ratios(wind.repsTC) |> dplyr::mutate(model = "TC")
wind.lrPF <- lik_ratios(wind.repssuperPF) |> dplyr::mutate(model = "PF")

wind.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                               "obs" = c(-2*(windOU$loglik - windSN$loglik),-2*(windOU$loglik - windTC$loglik),-2*(windOU$loglik - windsuperPF$loglik)))

wind_lr <- dplyr::left_join(rbind(wind.lrSN,wind.lrTC,wind.lrPF),wind.observed.dat,by="model") |>
  dplyr::mutate(lake = "Windermere")


kin.observed.dat <- data.frame("model" = c("SN","TC"),
                               "obs" = c(-2*(kinOU$loglik - kinSN$loglik),-2*(kinOU$loglik - kinTC$loglik)))

kin_lr <- dplyr::left_join(rbind(kin.lrSN,kin.lrTC),kin.observed.dat,by="model") |>
  dplyr::mutate(lake = "Kinneret")

wind.observed.dat <- data.frame("model" = c("SN","TC"),
                                "obs" = c(-2*(windOU$loglik - windSN$loglik),-2*(windOU$loglik - windTC$loglik)))

wind_lr <- dplyr::left_join(rbind(wind.lrSN,wind.lrTC),wind.observed.dat,by="model") |>
  dplyr::mutate(lake = "Windermere")


ggplot(rbind(kin_lr,wind_lr)) + 
  geom_density(aes(value, fill=simulation), alpha=0.8,size=0.3) + 
  geom_vline(aes(xintercept=obs),linetype = "dashed",col="grey20") + 
  scale_fill_manual(values = c(c("#A1B4FE","#FFE7A1")))+
  facet_grid(lake~model,scales = "free_x") +
  xlim(-80,80)+
  theme_bw()


mad.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(madOU, madSN)},mc.cores=9)
mad.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(madOU, madTC)},mc.cores=9)
mad.repssuperPF <- pbmcapply::pbmclapply(1:500, function(i){compare(madOU, madsuperPF)},mc.cores=9)

mad.lrSN <- lik_ratios(mad.repsSN) |> dplyr::mutate(model = "SN")
mad.lrTC <- lik_ratios(mad.repsTC) |> dplyr::mutate(model = "TC")
mad.lrPF <- lik_ratios(mad.repssuperPF) |> dplyr::mutate(model = "PF")

mad.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                               "obs" = c(-2*(madOU$loglik - madSN$loglik),-2*(madOU$loglik - madTC$loglik),-2*(madOU$loglik - madsuperPF$loglik)))

mad_lr <- dplyr::left_join(rbind(mad.lrSN,mad.lrTC,
                                 mad.lrPF),mad.observed.dat,by="model") |>
  dplyr::mutate(lake = "Mendota")


kas.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(kasOU, kasSN)},mc.cores=9)
kas.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(kasOU, kasTC)},mc.cores=9)
kas.repssuperPF <- pbmcapply::pbmclapply(1:500, function(i){compare(kasOU, kassuperPF)},mc.cores=9)

kas.lrSN <- lik_ratios(kas.repsSN) |> dplyr::mutate(model = "SN")
kas.lrTC <- lik_ratios(kas.repsTC) |> dplyr::mutate(model = "TC")
kas.lrPF <- lik_ratios(kas.repssuperPF) |> dplyr::mutate(model = "PF")

kas.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                               "obs" = c(-2*(kasOU$loglik - kasSN$loglik),-2*(kasOU$loglik - kasTC$loglik),-2*(kasOU$loglik - kassuperPF$loglik)))

kas_lr <- dplyr::left_join(rbind(kas.lrSN,kas.lrTC,
                                 kas.lrPF),kas.observed.dat,by="model") |>
  dplyr::mutate(lake = "Kasumigaura")

LZ.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(LZOU, LZSN)},mc.cores=9)
LZ.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(LZOU, LZTC)},mc.cores=9)
LZ.repssuperPF <- pbmcapply::pbmclapply(1:500, function(i){compare(LZOU, LZsuperPF)},mc.cores=9)

LZ.lrSN <- lik_ratios(LZ.repsSN) |> dplyr::mutate(model = "SN")
LZ.lrTC <- lik_ratios(LZ.repsTC) |> dplyr::mutate(model = "TC")
LZ.lrPF <- lik_ratios(LZ.repssuperPF) |> dplyr::mutate(model = "PF")

LZ.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                               "obs" = c(-2*(LZOU$loglik - LZSN$loglik),-2*(LZOU$loglik - LZTC$loglik),-2*(LZOU$loglik - LZsuperPF$loglik)))

LZ_lr <- dplyr::left_join(rbind(LZ.lrSN,LZ.lrTC,
                                LZ.lrPF),
                          LZ.observed.dat,by="model") |>
  dplyr::mutate(lake = "Lower Zurich")

UZ.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(UZOU, UZSN)},mc.cores=9)
UZ.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(UZOU, UZTC)},mc.cores=9)
UZ.repssuperPF <- pbmcapply::pbmclapply(1:500, function(i){compare(UZOU, UZsuperPF)},mc.cores=9)

UZ.lrSN <- lik_ratios(UZ.repsSN) |> dplyr::mutate(model = "SN")
UZ.lrTC <- lik_ratios(UZ.repsTC) |> dplyr::mutate(model = "TC")
UZ.lrPF <- lik_ratios(UZ.repssuperPF) |> dplyr::mutate(model = "PF")

UZ.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                              "obs" = c(-2*(UZOU$loglik - UZSN$loglik),-2*(UZOU$loglik - UZTC$loglik),-2*(UZOU$loglik - UZsuperPF$loglik)))

UZ_lr <- dplyr::left_join(rbind(UZ.lrSN,UZ.lrTC,
                                UZ.lrPF),
                          UZ.observed.dat,by="model") |>
  dplyr::mutate(lake = "Upper Zurich")

wash.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(washOU, washSN)},mc.cores=9)
wash.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(washOU, washTC)},mc.cores=9)
wash.repssuperPF <- pbmcapply::pbmclapply(1:500, function(i){compare(washOU, washsuperPF)},mc.cores=9)

wash.lrSN <- lik_ratios(wash.repsSN) |> dplyr::mutate(model = "SN")
wash.lrTC <- lik_ratios(wash.repsTC) |> dplyr::mutate(model = "TC")
wash.lrPF <- lik_ratios(wash.repssuperPF) |> dplyr::mutate(model = "PF")

wash.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                              "obs" = c(-2*(washOU$loglik - washSN$loglik),-2*(washOU$loglik - washTC$loglik),-2*(washOU$loglik - washsuperPF$loglik)))

wash_lr <- dplyr::left_join(rbind(wash.lrSN,wash.lrTC,
                                wash.lrPF),
                          wash.observed.dat,by="model") |>
  dplyr::mutate(lake = "Washington")

leve.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(leveOU, leveSN)},mc.cores=9)
leve.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(leveOU, leveTC)},mc.cores=9)
leve.repssuperPF <- pbmcapply::pbmclapply(1:500, function(i){compare(leveOU, levesuperPF)},mc.cores=9)

leve.lrSN <- lik_ratios(leve.repsSN) |> dplyr::mutate(model = "SN")
leve.lrTC <- lik_ratios(leve.repsTC) |> dplyr::mutate(model = "TC")
leve.lrPF <- lik_ratios(leve.repssuperPF) |> dplyr::mutate(model = "PF")

leve.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                                "obs" = c(-2*(leveOU$loglik - leveSN$loglik),-2*(leveOU$loglik - leveTC$loglik),-2*(leveOU$loglik - levesuperPF$loglik)))

leve_lr <- dplyr::left_join(rbind(leve.lrSN,leve.lrTC,
                                  leve.lrPF),
                            leve.observed.dat,by="model") |>
  dplyr::mutate(lake = "Loch Leven")

mon.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(monOU, monSN)},mc.cores=9)
mon.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(monOU, monTC)},mc.cores=9)
mon.repssuperPF <- pbmcapply::pbmclapply(1:500, function(i){compare(monOU, monsuperPF)},mc.cores=9)

mon.lrSN <- lik_ratios(mon.repsSN) |> dplyr::mutate(model = "SN")
mon.lrTC <- lik_ratios(mon.repsTC) |> dplyr::mutate(model = "TC")
mon.lrPF <- lik_ratios(mon.repssuperPF) |> dplyr::mutate(model = "PF")

mon.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                                "obs" = c(-2*(monOU$loglik - monSN$loglik),-2*(monOU$loglik - monTC$loglik),-2*(monOU$loglik - monsuperPF$loglik)))

mon_lr <- dplyr::left_join(rbind(mon.lrSN,mon.lrTC,
                                  mon.lrPF),
                            mon.observed.dat,by="model") |>
  dplyr::mutate(lake = "Monona")

x <- wind_lr
obs = wind.observed.dat

x <- rbind(kin.lrSN,kin.lrTC,
           kin.lrPF)
obs = kin.observed.dat

find_quartile(x=rbind(kin.lrSN,kin.lrTC,kin.lrPF),quartiles = c(0.25,0.75),obs=kin.observed.dat)
find_quartile(x=rbind(kas.lrSN,kas.lrTC,kas.lrPF),quartiles = c(0.25,0.75),obs=kas.observed.dat)
find_quartile(x = rbind(LZ.lrSN,LZ.lrTC,LZ.lrPF),quartiles = c(0.25,0.75),obs = LZ.observed.dat)
find_quartile(x=rbind(mad.lrSN,mad.lrTC,mad.lrPF),quartiles = c(0.25,0.75),obs = mad.observed.dat)
find_quartile(x = rbind(wind.lrSN,wind.lrTC,wind.lrPF),quartiles = c(0.25,0.75),obs = wind.observed.dat)

ggplot(rbind(kin_lr,wind_lr,mad_lr,LZ_lr,leve_lr,wash_lr,kas_lr,UZ_lr,mon_lr)) + 
  geom_density(aes(value, fill=simulation), alpha=0.8,size=0.3) + 
  geom_vline(aes(xintercept=obs),linetype = "dashed",col="grey20") + 
  xlim(-30,30)+
  scale_fill_manual(values = c(c("#A1B4FE","#FFE7A1")))+
  facet_grid(lake~model,scales = "free_x") +theme_bw()


ggplot(rbind(kin_lr,wind_lr,mad_lr,LZ_lr,wash_lr,kas_lr,UZ_lr) ) + 
  geom_density(aes(value, fill=simulation), alpha=0.8,size=0.3) + 
  geom_vline(aes(xintercept=obs),linetype = "dashed",col="grey20") + 
  xlim(-30,30)+
  scale_fill_manual(values = c(c("#A1B4FE","#FFE7A1")))+
  facet_wrap(~model,scales = "free_x") +theme_bw()

tt <- lapply(unique(kin_lr$model),function(x){
  
  roc_data(dat = subset(kin_lr,model == x)) %>%
    dplyr::mutate(model = x)
  
}) |> data.table::rbindlist()

ggplot(tt) + geom_line(aes(False.positives, True.positives)) + 
  facet_wrap(~model) +theme_bw()


