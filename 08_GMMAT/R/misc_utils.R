time.since <- function(start_time, units = "minutes"){
  if(units == "seconds"){
    divide_by <- 1
  }
  if(units == "minutes"){
    divide_by <- 60
  }
  if(units == "hours"){
    divide_by <- 3600
  }
  elapsed_time <- ((proc.time() - start_time)[3])/divide_by
  elapsed_time_to_print <- prettyNum(elapsed_time, digits = 2)
  time_out <- paste(elapsed_time_to_print, units, sep=" ")
  time_out
}

find_FDR <- function(p_vals, alpha = 0.05){
  p_vals <- as.data.table(p_vals)
  p_vals <- as.vector(p_vals[order(p_vals[,1])])
  p_adjusted <- p.adjust(unlist(p_vals[,1]),
                         method = "BH", n = nrow(p_vals))
  FDR_p_val_threshold <- as.numeric(p_vals[length(subset(p_adjusted,
                                                         p_adjusted < alpha))])
  FDR_p_val_threshold
}
