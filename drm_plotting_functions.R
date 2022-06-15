
individual_plots <- function(data, model, assay){

  ####Use model to generate fitted curve for plotting
  MaxD <- max(data$dose)
  MinD <- min(data$dose)
  MaxR <- max(data$value)
  MinR <- min(data$value)
  #generate predicted curves 
  fit <- expand.grid(dose=exp(seq(log(MinD), log(MaxD), length=100))) 
  #use model to make predictions for new dose scale
  pv <- predict(model, newdata=fit, interval="confidence") 
  #add results to dataframe
  fit$pr <- pv[,1]
  fit$pmin <- pv[,2]
  fit$pmax <- pv[,3]
  fit<- data.frame(fit)
  
  extra_var <- data.frame(cell = rep(assay$cell, length(data)),
                          drug = rep(assay$drug, length(data)),
                          starting_uM = rep(assay$starting_uM, length(data)),
                          plate_id = rep(assay$plate_id, length(data)),
                          position_id = rep(assay$position_id, length(data)))
  print(length(data))
  print(data)
  data_exp <- cbind(data, extra_var) 
  print(data)
  extra_var <- data.frame(cell = rep(assay$cell, length(fit)),
                          drug = rep(assay$drug, length(fit)),
                          starting_uM = rep(assay$starting_uM, length(fit)),
                          plate_id = rep(assay$plate_id, length(fit)),
                          position_id = rep(assay$position_id, length(fit)))

  fit <- cbind(fit, extra_var) 
  
  fit_2 <-fit
  title = paste(assay$drug, ' treated ', assay$cell)
  subtitle = paste(assay$format, ', ',assay$plate_id,', ',assay$position_id)
  IC50 <- as.numeric(as.character(assay$IC50))
  conc <- 'uM'

  MinD <- min(data$dose)-1000
  if (assay$IC50  > MaxD | MinR >45){IC50<-paste0(">",MaxD)}
  else if (IC50 < 0.1){
    data$dose <- data$dose*1000
    fit_2$dose <- fit_2$dose*1000
    IC50 <- IC50 * 1000
    conc <- 'nM'
 
    MinD <- min(data$dose)
    MaxD <- max(data$dose)
  }
 
  IC50 <- as.numeric(as.character(IC50))
  IC50 <- round(IC50, 2)
  label <- paste0("IC50\n",as.character(IC50), conc)
  print(data)
  #Generate plot
  p <- ggplot(data, aes(x = dose, y = value)) +
    geom_point(size = 0.5) +
    geom_ribbon(data=fit_2, aes(x=dose, y=pr, ymin=pmin, ymax=pmax), alpha=0.2) +
    geom_line(data=fit_2, aes(x=dose, y=pr))+
    geom_vline(aes(xintercept=IC50), size = 0.4, linetype="dotted", color="blue")+
    annotate("text", x = MaxD, y = 120, label = paste0(label), colour = 'blue', size=2)+
    labs(title = title, subtitle = subtitle, x = paste0('log(Drug[',conc,'])'), y = "Growth (%)") +
    scale_x_log10(limits=c(MinD,MaxD))+
    ylim(-5,120)+ theme_classic(base_size = 8)
  
  #print(data_exp)
  plot_fit = list("plot" = p, "curve" = fit, "data" = data_exp)
  return(plot_fit)
}


failed_plots <- function(dataP, info){
  #when fitting fails still include a plot of the data
  MaxD <- max(dataP[1])
  MinD <- min(dataP[1])
  MaxR <- max(dataP[3])
  
  p <- ggplot(dataP, aes(x = dose, y = value)) +
    geom_point(size = 0.5)+
    labs(x = 'log(Drug[uM])', y = "Growth (%)")+
    ggtitle(info[1, 1])+ #plot.title = element_text(hjust = 0.5)+
    scale_x_log10(limits=c(MinD-1000,MaxD))+
    ylim(-5,120)+ theme_classic(base_size = 8)
  
  return(p)
  
}


grouped_plots <- function(all_assays,all_curves){
  
  all_grouped_plots <- list()
  drugs <- unique(all_assays$drug)
  h<-1

  for (i in drugs){
    fit <- subset(all_curves, drug ==i)
    
    MaxD <- max(fit$dose)
    MinD <- min(fit$dose)
    

    fit$cell<-as.factor(fit$cell)
    fit$plate_id <- as.factor(fit$plate_id)
    fit$position_id <- as.factor(fit$position_id)
    
    p <- ggplot(fit, aes(x = dose, y = value,  colour =plate_id,  group=interaction(plate_id, position_id))) +
      geom_line(data=fit, aes(x=dose, y=pr, linetype = cell), size = 0.4)+
      #geom_ribbon(data=fit, aes(x=conc, y=pr, ymin=pmin, ymax=pmax), alpha=0.1, linetype = 0) +
      labs(title = i, x = 'log(Drug[uM])', y = "Growth (%)")+
      #scale_linetype(guide = "none") +
      scale_x_log10(labels = scales::comma, limits=c(MinD-1000,MaxD))+
      ylim(-5,110)+ 
      theme_classic(base_size = 8)+
      theme(plot.margin = unit(c(1,0.5,1,0.5), "cm"))
  
 
      
    
    all_grouped_plots[[h]] <- p
    h <- h+1

  }
  return(all_grouped_plots)
}
