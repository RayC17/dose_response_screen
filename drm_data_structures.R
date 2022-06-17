

get_locations <- function(assay, plate){
  #here are the four available formats
  
  triplicate_384 <- data.frame(position_name = c('pos_1','pos_2','pos_3','pos_4','pos_5',
                                                 'pos_6','pos_7','pos_8','pos_9','pos_10'),
                               set_rows_t = c(2,5,8,11,14,2,5,8,11,14),
                               set_rows_b = c(4,7,10,13,16,4,7,10,13,16),
                               set_cols_l = c(2,2,2,2,2,14,14,14,14,14),
                               set_cols_r = c(11,11,11,11,11,23,23,23,23,23))
  
  duplicate_384 <- data.frame(position_name = c('pos_1','pos_2','pos_3','pos_4','pos_5',
                                                'pos_6','pos_7','pos_8','pos_9','pos_10',
                                                'pos_11','pos_12','pos_13','pos_14'),
                              set_rows_t = c(2,4,6,8,10,12,14,2,4,6,8,10,12,14),
                              set_rows_b = c(3,5,7,9,11,13,15,3,5,7,9,11,13,15),
                              set_cols_l = c(2,2,2,2,2,2,2,14,14,14,14,14,14,14),
                              set_cols_r = c(11,11,11,11,11,11,11,23,23,23,23,23,23,23))
  
  triplicate_96 <- data.frame(position_name = c('pos_1','pos_2'),
                              set_rows_t = c(2,5),
                              set_rows_b = c(4,7),
                              set_cols_l = c(1,1),
                              set_cols_r = c(10,10))
  
  duplicate_96 <- data.frame(position_name = c('pos_1','pos_2','pos_3'),
                             set_rows_t = c(2,4,6),
                             set_rows_b = c(3,5,7),
                             set_cols_l = c(1,1,1),
                             set_cols_r = c(10,10,10))
  
  locate<-NULL

  if (assay$format == '384w'){

    if (assay$replicates == 3){
      locate <- filter(triplicate_384, position_name == assay$position_id)
    }
    if (assay$replicates == 2){
      locate <- filter(duplicate_384, position_name == assay$position_id)
    }
    locate$blank <-  mean(as.numeric(as.character(plate[locate$set_rows_t:locate$set_rows_b,13])))
    locate$no_drug <- mean(as.numeric(as.character(plate[locate$set_rows_t:locate$set_rows_b,12])))
    
  }
  if (assay$format == '96w'){

    if (assay$replicates == 3){
      locate <- filter(triplicate_96, position_name == assay$position_id)
    }
    if (assay$replicates == 2){
      locate <- filter(duplicate_96, position_name == assay$position_id)
    }
    locate$blank <-  mean(as.numeric(as.character(plate[locate$set_rows_t:locate$set_rows_b,12])))
    locate$no_drug <- mean(as.numeric(as.character(plate[locate$set_rows_t:locate$set_rows_b,11])))
  }
  
  return(locate)
}


get_assay_data <- function(data, location, assays_info){
  #retrieve assay specific data
  data_set <- as.data.frame(t(data[location$set_rows_t:location$set_rows_b,
                                   location$set_cols_l:location$set_cols_r]))
  
  #add no drug and blank data
  data_set[nrow(data_set) + 1,] = c( rep(location$no_drug,assays_info$replicates))
  data_set[] <- lapply(data_set, function(x) as.numeric(as.character(x)))
  data_set_blank <- data_set %>% mutate_all(.funs = list(~. - location$blank))
  
  #Create percentage dataframe
 
  data_pct <- data.frame(sapply(names(data_set_blank), 
      function(x) { data_set_blank[paste0(x, "_pct")] <<- ((data_set_blank[x] /                                                     data_set_blank[nrow(data_set_blank),x])*100) }))
  
  #create dose series and add to data
  dilutionF <- as.numeric(as.character(assays_info$dilution_factor))
  startDose <- as.numeric(as.character(assays_info$starting_uM))
  dose <- c()
  
  for (q in 1:nrow(data_pct )){
    dose[q] <- startDose
    startDose <- startDose/dilutionF
    
  }

  data_pct  <- cbind(dose = dose, data_pct )
  
  #rename columns
  if (ncol(data_pct ) == 3){
    colnames(data_pct ) <- c("dose", "R1", "R2")
  } 
  else { 
    colnames(data_pct ) <- c("dose", "R1", "R2", "R3")
  }
  
  #convert to long format
  data_pct  <- reshape2::melt(data_pct , id.vars = "dose")

  return(data_pct)

}






