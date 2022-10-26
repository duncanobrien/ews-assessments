########################################################################################
## Preamble ##
########################################################################################
bifur.mod.scripts <- list.files(path ="Code/bifur_models", pattern="*.R",full.names = T) 
purrr::walk(bifur.mod.scripts, source) # source silently

load("Data/wrangled_genus_plank_data.Rdata")

########################################################################################
## Perform Assessments ##
########################################################################################
kin_density <- data.frame(seq_along(as.numeric(kin_mth_dat$Date[kin_mth_dat$Date >1980])),rowSums(kin_mth_dat[kin_mth_dat$Date >1980,2:50]))

wind_density <- data.frame(seq_along(wind_mth_dat$Date[wind_mth_dat$Date > 1984]),rowSums(wind_mth_dat[wind_mth_dat$Date > 1984,2:23]))

mad_density <- data.frame(seq_along(mad_mth_dat$date),rowSums(mad_mth_dat[,2:103]))

kas_density <- data.frame(seq_along(kas_mth_dat$date),rowSums(kas_mth_dat[,2:89]))

LZ_density <- data.frame(seq_along(LZ_mth_dat$date),rowSums(LZ_mth_dat[,2:137]))

kinOU <- fit_bif_mod(kin_density,model = "OU")
kinSN <- fit_bif_mod(kin_density,model = "LSN")
kinTC <- fit_bif_mod(kin_density,model = "LTC")
kinPF <- fit_bif_mod(kin_density,model = "LPF")

kin.observedSN <- -2*(kinOU$loglik - kinSN$loglik)
kin.observedTC<- -2*(kinOU$loglik - kinTC$loglik)
kin.observedPF <- -2*(kinOU$loglik - kinPF$loglik)
knitr::kable(cbind(kin.observedSN,kin.observedTC,kin.observedPF))

windOU <- fit_bif_mod(wind_density,model = "OU")
windSN <- fit_bif_mod(wind_density,model = "LSN")
windTC <- fit_bif_mod(wind_density,model = "LTC")
windPF <- fit_bif_mod(wind_density,model = "LPF")

wind.observedSN <- -2*(windOU$loglik - windSN$loglik)
wind.observedTC<- -2*(kinOU$loglik - windTC$loglik)
wind.observedPF <- -2*(windOU$loglik - windPF$loglik)
knitr::kable(cbind(wind.observedSN,wind.observedTC,wind.observedPF))

madOU <- fit_bif_mod(mad_density,model = "OU")
madSN <- fit_bif_mod(mad_density,model = "LSN")
madTC <- fit_bif_mod(mad_density,model = "LTC")
madPF <- fit_bif_mod(mad_density,model = "LPF")

mad.observedSN <- -2*(madOU$loglik - madSN$loglik)
mad.observedTC<- -2*(madOU$loglik - madTC$loglik)
mad.observedPF <- -2*(madOU$loglik - madPF$loglik)

kasOU <- fit_bif_mod(kas_density,model = "OU")
kasSN <- fit_bif_mod(kas_density,model = "LSN")
kasTC <- fit_bif_mod(kas_density,model = "LTC")
kasPF <- fit_bif_mod(kas_density,model = "LPF")

kas.observedSN <- -2*(kasOU$loglik - kasSN$loglik)
kas.observedTC<- -2*(kasOU$loglik - kasTC$loglik)
kas.observedPF <- -2*(kasOU$loglik - kasPF$loglik)

LZOU <- fit_bif_mod(LZ_density,model = "OU")
LZSN <- fit_bif_mod(LZ_density,model = "LSN")
LZTC <- fit_bif_mod(LZ_density,model = "LTC")
LZPF <- fit_bif_mod(LZ_density,model = "LPF")

LZ.observedSN <- -2*(LZOU$loglik - LZSN$loglik)
LZ.observedTC<- -2*(LZOU$loglik - LZTC$loglik)
LZ.observedPF <- -2*(LZOU$loglik - LZPF$loglik)

kin.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(kinOU, kinSN)},mc.cores=9)
kin.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(kinOU, kinTC)},mc.cores=9)
kin.repsPF <- pbmcapply::pbmclapply(1:500, function(i){compare(kinOU, kinPF)},mc.cores=9)

kin.lrSN <- lik_ratios(kin.repsSN) |> dplyr::mutate(model = "SN")
kin.lrTC <- lik_ratios(kin.repsTC) |> dplyr::mutate(model = "TC")
kin.lrPF <- lik_ratios(kin.repsPF) |> dplyr::mutate(model = "PF")

kin.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                           "obs" = c(-2*(kinOU$loglik - kinSN$loglik),-2*(kinOU$loglik - kinTC$loglik),-2*(kinOU$loglik - kinPF$loglik)))

kin_lr <- dplyr::left_join(rbind(kin.lrSN,kin.lrTC,kin.lrPF),observed.dat,by="model") |>
  dplyr::mutate(lake = "Kinneret")

wind.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(windOU, windSN)},mc.cores=9)
wind.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(windOU, windTC)},mc.cores=9)
wind.repsPF <- pbmcapply::pbmclapply(1:500, function(i){compare(windOU, windPF)},mc.cores=9)

wind.lrSN <- lik_ratios(wind.repsSN) |> dplyr::mutate(model = "SN")
wind.lrTC <- lik_ratios(wind.repsTC) |> dplyr::mutate(model = "TC")
wind.lrPF <- lik_ratios(wind.repsPF) |> dplyr::mutate(model = "PF")

wind.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                               "obs" = c(-2*(windOU$loglik - windSN$loglik),-2*(windOU$loglik - windTC$loglik),-2*(windOU$loglik - windPF$loglik)))

wind_lr <- dplyr::left_join(rbind(wind.lrSN,wind.lrTC,wind.lrPF),observed.dat,by="model") |>
  dplyr::mutate(lake = "Windermere")

mad.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(madOU, madSN)},mc.cores=9)
mad.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(madOU, madTC)},mc.cores=9)
mad.repsPF <- pbmcapply::pbmclapply(1:500, function(i){compare(madOU, madPF)},mc.cores=9)

mad.lrSN <- lik_ratios(mad.repsSN) |> dplyr::mutate(model = "SN")
mad.lrTC <- lik_ratios(mad.repsTC) |> dplyr::mutate(model = "TC")
mad.lrPF <- lik_ratios(mad.repsPF) |> dplyr::mutate(model = "PF")

mad.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                               "obs" = c(-2*(madOU$loglik - madSN$loglik),-2*(madOU$loglik - madTC$loglik),-2*(madOU$loglik - madPF$loglik)))

mad_lr <- dplyr::left_join(rbind(mad.lrSN,mad.lrTC,mad.lrPF),observed.dat,by="model") |>
  dplyr::mutate(lake = "Mendota")


kas.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(kasOU, kasSN)},mc.cores=9)
kas.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(kasOU, kasTC)},mc.cores=9)
kas.repsPF <- pbmcapply::pbmclapply(1:500, function(i){compare(kasOU, kasPF)},mc.cores=9)

kas.lrSN <- lik_ratios(kas.repsSN) |> dplyr::mutate(model = "SN")
kas.lrTC <- lik_ratios(kas.repsTC) |> dplyr::mutate(model = "TC")
kas.lrPF <- lik_ratios(kas.repsPF) |> dplyr::mutate(model = "PF")

kas.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                               "obs" = c(-2*(kasOU$loglik - kasSN$loglik),-2*(kasOU$loglik - kasTC$loglik),-2*(kasOU$loglik - kasPF$loglik)))

kas_lr <- dplyr::left_join(rbind(kas.lrSN,kas.lrTC,kas.lrPF),observed.dat,by="model") |>
  dplyr::mutate(lake = "Kasumigaura")

LZ.repsSN <- pbmcapply::pbmclapply(1:500, function(i){compare(LZOU, LZSN)},mc.cores=9)
LZ.repsTC <- pbmcapply::pbmclapply(1:500, function(i){compare(LZOU, LZTC)},mc.cores=9)
LZ.repsPF <- pbmcapply::pbmclapply(1:500, function(i){compare(LZOU, LZPF)},mc.cores=9)

LZ.lrSN <- lik_ratios(LZ.repsSN) |> dplyr::mutate(model = "SN")
LZ.lrTC <- lik_ratios(LZ.repsTC) |> dplyr::mutate(model = "TC")
LZ.lrPF <- lik_ratios(LZ.repsPF) |> dplyr::mutate(model = "PF")

LZ.observed.dat <- data.frame("model" = c("SN","TC","PF"),
                               "obs" = c(-2*(LZOU$loglik - LZSN$loglik),-2*(LZOU$loglik - LZTC$loglik),-2*(LZOU$loglik - LZPF$loglik)))

LZ_lr <- dplyr::left_join(rbind(LZ.lrSN,LZ.lrTC,LZ.lrPF),observed.dat,by="model") |>
  dplyr::mutate(lake = "Lower Zurich")





ggplot(rbind(kin_lr,wind_lr,mad_lr,kas_lr,LZ_lr)) + 
  geom_density(aes(value, fill=simulation), alpha=0.8,size=0.3) + 
  geom_vline(aes(xintercept=obs),linetype = "dashed",col="grey20") + 
  scale_fill_manual(values = c(c("#A1B4FE","#FFE7A1")))+
  facet_grid(lake~model,scales = "free_x") +theme_bw()

