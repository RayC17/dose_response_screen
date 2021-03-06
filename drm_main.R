
#------- import the packages ---------------------------------------------------

use_these_packages <- c("dplyr", "ggplot2", "data.table", "drc", 
                        "gridExtra", "plyr", "scales")
new_packages <- use_these_packages[!(use_these_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  install.packages(new_packages)
  sapply(use_these_packages, require, character.only = TRUE)
} else {
  sapply(use_these_packages, require, character.only = TRUE)
}

#-------- source the tool functions---------------------------------------------
#utility scripts with jobs broken down into functions
use_these_utilities <- c("~/dose_response_screen/drm_data_structures.R",
                         "~/dose_response_screen/drm_plotting_functions.R",
                         "~/dose_response_screen/drm_IO_functions.R")
sapply(use_these_utilities, source)

#--------- set job parameters --------------------------------------------------
#identifer used as prefix for results output
batch_id <- 'screen-03-22'
#csv format containing assay data, see project folder /example_meta.csv
meta <- "screen_03-22-meta.csv"  
#set you're working directory to a folder containing raw 'TRno' prefixed csv
# from the plate reader and the meta file
# sessions > set working directory > choose directory

#-------- Initialise main and set job running ------------------------------------------------------

main <- function(meta, batch_id){
  cat('    ---------------------------------------------------------------------\n')
  cat(paste0('    | Drug Assay Batch: ', batch_id, '\n'))
  cat('    ---------------------------------------------------------------------\n')
  
  file.names <- dir(pattern ="TRno|automated")
  assays <- as.data.frame(read.csv(meta, sep = ",",header = TRUE, stringsAsFactors = FALSE))
  
  all_Plots <- list()
  results_all <- data.frame(matrix(ncol = 17, nrow = 0))
  curves_all <- data.frame(matrix(ncol = 8, nrow = 0))
  data_all <- data.frame(matrix(ncol = 8, nrow = 0))
  c<-1
  #---------loop through files ---------------------------------------------------
  for (i in 1:length(file.names)){
    plate <- import_plate(file.names[i], assays) 
  
    if (!is.null(plate$message)){
      cat(paste0('    | ', file.names[i], ': ',plate$message, '. Skipping plate\n'))
    
      cat('    ---------------------------------------------------------------------\n')
      i<-i+1
      next
    }
    cat(paste0('    | ',file.names[i], ': ID1 = ',plate$file_id,'.\n'))
    plate_assays <-  assays %>% filter(plate_id  == plate$file_id)
    plate_assays$file_name <- file.names[i]
   
    #--------loop through plate---------------------------------------------------
    for (j in 1:nrow(plate_assays)){
      cat(paste0('    | ',plate_assays$position_id[j],
                 ' is ', plate_assays$compound[j], ' treated ', plate_assays$cell[j],'.\n'))
      
      # get plate layout info for each assay
      locate <- get_locations(plate_assays[j,], plate$data)
      
      # retrieve data from plate for assay
      data <- get_assay_data(plate$data, locate, plate_assays[j,])
      
    
      
      #fitting function
      control <- drmc(errorm = TRUE,noMessage = TRUE, warnVal = -1)

      possibleError <- tryCatch(
        LL.4 <- drm(data = data,value~dose,fct=LL.4(fixed=c(NA,NA, NA, NA)),
                    na.action = na.omit, control = control),
        error=function(e) e,
        silent = TRUE
      )
   
      if(inherits(possibleError, "error")) {
        cat('    | Convergence failed. The model was not fit!\n')
        plate_assays$IC50[j] <- 'Convergence failed. The model was not fit! '
        results_all[nrow(results_all) + 1,] <- c(plate_assays[j,],NA,NA,NA,NA,NA)
        all_Plots[[c]] <-failed_plots(data , plate_assays[j,])
        c<-c+1
        next
      }
      cat(paste0('    | model fit succesffuly.\n'))
      #get model coefficients and combine in one table
      model_results <- retrieve_results(LL.4 , plate_assays[j,] )
     
      results_all[nrow(results_all) + 1,] <- model_results 
      
      #generates predited curves from model and plots each assay
      plot_fit <- individual_plots(data,LL.4 , model_results)
     
      all_Plots[[c]] <- plot_fit$plot
      c <- c+1
      
      #combine all predicted curves together 
      curves_all <- rbind(curves_all,plot_fit$curve)
      #combine all input data together
      data_all<- rbind(data_all, plot_fit$data)
    
    }#-----end loop for plate/file -----------------------------------------------
    cat('    ---------------------------------------------------------------------\n')
  }#-------end of loop for files, final exports-----------------------------------
  
  #takes all data and generate plots grouped by compound
  all_grouped_plots <-  grouped_plots(data_all, curves_all)
  #exports all results

  export_results(batch_id,results_all,curves_all,
                 data_all, all_Plots,all_grouped_plots)
}
#sapply(use_these_utilities, source)

main(meta, batch_id)          





