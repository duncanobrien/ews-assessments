find_quartile <- function(x,obs,quartiles = c(0.025,0.975),zero_tol = 0.5, selection.meth = "min.log.lik"){
  
  require(magrittr)
  selection.meth <- match.arg(selection.meth,choices = c("max.median.diff","min.log.lik","min.IQR"))
  
  if(length(quartiles) != 2){
    stop("quartiles must contain a lower and upper quartile")
  }
  if(quartiles[1] <! quartiles[1]){
    stop("the first entry in quartiles must be smaller than the second")
  }
  
  out_df <- x %>%
    dplyr::group_by(model)%>%
    dplyr::mutate(mw_p = wilcox.test(value ~ simulation,alternative = "two.sided")$p.value) %>%
    nest(data = -c(model)) %>%
    dplyr::mutate(data = purrr::map(data,
                                    ~.x %>%
                                      dplyr::group_by(simulation) %>%
                                      dplyr::filter(value < quantile(value,0.99,type = 8) &
                                                      value > quantile(value,0.01,type = 8)) %>%
                                      dplyr::group_by(rep) %>%
                                      dplyr::filter(abs(value) < abs(2e+3))
                                    ))%>%
    unnest(data) %>%
    dplyr::group_by(model,simulation)  %>%
    dplyr::group_modify(function(x,...){
      dens <-ecdf(x[["value"]]) #estimate Empirical Cumulative Distribution Function
      percent_zero <- length(x[["value"]][x[["value"]]==0])/length(x[["value"]]) #count number of zeroes in vale
      return(data.frame(x,"median" = median(x[["value"]]),"LQ" = quantile(dens,c(quartiles[1]),type = 8) %>% `names<-`(NULL),"UQ" = quantile(dens,c(quartiles[2]),type = 8)%>% `names<-`(NULL), "percent_zero" = percent_zero))
    }) %>% #return quantiles and percentage of zeroes
    ungroup() %>%
    dplyr::group_by(rep,model)%>%
    dplyr::filter(simulation[1] == "null" & simulation[2] == "test") %>% #filter to reps where have both null and test
    arrange(model,rep,simulation)%>%
    dplyr::group_by(model) %>%
    dplyr::select(-c(value,rep)) %>%
    dplyr::slice_head(n=2) %>% #only keep head
    dplyr::left_join(obs , by = "model") %>% #join observed value
    dplyr::mutate(reject_model = ifelse(obs < LQ | obs > UQ, TRUE,FALSE)) %>% #identify whether obs is within CI
    dplyr::group_by(model) %>%
    dplyr::mutate(sim_choice = dplyr::case_when(
      mw_p[1] > 0.025 ~ "null",
      all(reject_model == TRUE) ~ "null",
      #all(reject_model  == FALSE) ~ model[1],
      all(reject_model  == FALSE) ~ "null",
      reject_model[2]  == FALSE & reject_model[1]  == TRUE ~ model[1],
      reject_model[2]  == FALSE & percent_zero[2] > zero_tol ~ "null",
      TRUE ~ "null")) #initial classification of model rejection WITHIN each null test

  outcome <- subset(out_df,sim_choice != "null") #filter to just non-rejected models

  if(dim(outcome)[1] == 0){ #if none are test models, return null
    outcome <- "null"
    
  }else{ #else classify base upon whether only one non-rejected model is returned vs, if multiple, return model with smallest CI
    outcome <-  dplyr::case_when(
      length(unique(outcome$sim_choice)) == 1 ~ unique(outcome$sim_choice),
      length(unique(outcome$sim_choice)) >= 1 & selection.meth == "max.median.diff" ~ 
        outcome %>%  dplyr::group_by(model) %>% 
        dplyr::summarise(median_diff = diff(median)) %>%
        dplyr::filter(median_diff == max(median_diff)) %>%
        dplyr::pull(model),
      length(unique(outcome$sim_choice)) >= 1  & selection.meth == "min.IQR" ~ 
        subset(outcome,simulation == "test")$sim_choice[which.min(subset(outcome,simulation == "test")$UQ - 
                                                                    subset(outcome,simulation == "test")$LQ)],
      length(unique(outcome$sim_choice)) >= 1 & selection.meth == "min.log.lik" ~ 
        subset(outcome,simulation == "test")$sim_choice[which.max(subset(outcome,simulation == "test")$obs)])
  }

  

    # outcome <- dplyr::case_when(isTRUE(all(out_df$reject_model)) ~ "null",
    #                           isFALSE(subset(out_df,model == "SN" & simulation == "test")$reject_model) & 
    #                             all((subset(out_df,model != "SN") %>% 
    #                                      dplyr::group_by(model) %>%
    #                                    mutate(group_truth = any(reject_model)))$group_truth)  ~ "SN",
    #                           isFALSE(subset(out_df,model == "TC" & simulation == "test")$reject_model) & 
    #                             all((subset(out_df,model != "TC") %>% 
    #                                    dplyr::group_by(model) %>%
    #                                    mutate(group_truth = any(reject_model)))$group_truth)  ~ "TC",
    #                           isFALSE(subset(out_df,model == "PF" & simulation == "test")$reject_model) & 
    #                             all((subset(out_df,model != "PF") %>% 
    #                                    dplyr::group_by(model) %>%
    #                                    mutate(group_truth = any(reject_model)))$group_truth) ~ "PF",
    #                           TRUE ~ "null")
  
    return(list("x" = as.data.frame(out_df), "outcome" = outcome[1]))
}
