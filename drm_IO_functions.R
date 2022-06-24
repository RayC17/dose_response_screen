
import_plate <- function(file.name, assays){
  message <- NULL
  file_id <-NULL
  data <- NULL
  tryCatch(
    data <- as.data.frame(read.csv(file.name, sep = ",",header = FALSE,
                                    skip = 10, stringsAsFactors = FALSE)),
    error=function(e) NULL)
  
  if(is.null(data)||nrow(data)==0){
      message <- 'plate readings are empty!'
  }
  tag <- as.data.frame(read.csv(file.name, sep = ",",header = FALSE,
                                    nrows = 3, stringsAsFactors = FALSE))[3,1]

  if(nchar(tag) < 5){
    message <- 'ID1 empty, no plate ID!'
    
  }else{
    file_id <- gsub(" ","", strsplit(tag, " ")[[1]][2], fixed = TRUE)
    
    if(!any(file_id == assays$plate_id)  ){
      message <-'plate id does not match meta file'}
  }

  
  plate <- list(data = data, file_id = file_id, message = message)
  return(plate)
}


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
  
  if (as.numeric(SE[1,1]) > as.numeric(coefs[3])){
    assays$IC50 <-  as.numeric(coefs[3]) 
  }else {assays$IC50<- as.numeric(SE[1,1])}
  
 
  assays$IC50_SE<- as.numeric(SE[1,2])
  assays$hill_slope <- as.numeric(coefs[1])
  assays$min <- as.numeric(coefs[2])
  assays$max <- as.numeric(coefs[3])
  assays$GI50<- GI50
  
 
  
  return(assays)
}

export_results <- function(batch_id, results_all,curves_all,
                           data_all,plots_all, grouped_all){
  
  
  colnames(results_all) <- c("plate_id","position_id","format","replicates",'index',"compound",           
                             "cell","starting_uM","dilution_factor","batch","file_name","IC50","IC50_SE",
                             "hill_slope","min","max","GI50") 
  
  colnames(curves_all) <- c( 'dose','pr', 'pmin', 'pmax', 'cell', 'compound',
                             'starting_uM', 'plate_id','pos','index')
  colnames(data_all) <- c( 'dose', 'variable', 'value', 'cell', 'compound', 
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
  

