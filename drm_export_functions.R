
retrieve_results <- function(model, assays){
  #get coefs and calculate GI50
  coefs <- setNames(
    model$coefficients,
    c("hill", "min_value", "max_value", "IC50")
  )
  
  GI50 <- with(
    as.list(coefs),
    exp(
      log(IC50) + (1 / hill) * log(max_value / (max_value - 2 * min_value))
    )
  )
  #Collect IC50 and IC90 with estimated standard errors
  SE<- data.frame(ED(model, c(50, 90), display = FALSE))
  assays$IC50<- as.numeric(SE[1,1])
  assays$IC50_SE<- as.numeric(SE[1,2])
  assays$hill_slope <- as.numeric(coefs[1])
  assays$min <- as.numeric(coefs[2])
  assays$max <- as.numeric(coefs[3])
  assays$GI50<- GI50
  
  return(assays)
}

export_results <- function(batch_id, results_all,curves_all,
                           data_all,plots_all, grouped_all){
  
  
  colnames(results_all) <- c("plate_id","position_id","format","replicates","drug",           
                             "cell","starting_uM","dilution_factor","vol_ul","batch","IC50","IC50_SE",
                             "hill_slope","min","max","GI50") 
  
  colnames(curves_all) <- c( 'dose','pr', 'pmin', 'pmax', 'cell', 'drug',
                             'starting_uM', 'plate_id','pos')
  colnames(data_all) <- c( 'dose', 'variable', 'value', 'cell', 'drug', 
                           'starting_uM', 'plate_id','pos')
  
  #Save results, data and predicted curves in csv results file
  write.csv(curves_all,paste0(batch_id,"_curves.csv"), row.names = FALSE)
  write.csv(data_all,paste0(batch_id,"_assays.csv"), row.names = FALSE)
  write.csv(results_all,paste0(batch_id,"_results.csv"), row.names = FALSE)
  
  #all individual plots into pdf 
  plots <- marrangeGrob(plots_all, nrow = 3, ncol = 3)           
  ind_name<- paste0('all_plots_', batch_id,'.pdf')
  ggsave(ind_name, plots, width = 11, height = 8.5, units = "in")
  
  #generate and save pdf of all grouped plots
  group_name<- paste0('all_plots_grouped_', batch_id,'.pdf')
  plots <- marrangeGrob(grouped_all, nrow = 2, ncol = 2)
  ggsave(group_name, plots, width = 11, height = 8.5, units = "in")
  
}  
  

