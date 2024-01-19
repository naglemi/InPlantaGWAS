getwd()
setwd("/mnt/data/NSF_GWAS/notebooks/InPlantaGWAS/09_Data_mining/")

library(data.table)
library(qqman)
library(ggplot2)
library(foreach)
library(ggpubr)
library(InterMineR)
library(tools)
library(dplyr)

data <- fread("/mnt/data/NSF_GWAS/notebooks/InPlantaGWAS/13_GMMAT_ART/5-OUT_Table_results_by_method.csv", na.strings="")

outdir <- "Results_GMMAT_ART/"

for(batch in batches_to_exclude){
  data[, batch] <- NA
}

titles <- fread("batch_titles_a2.csv", header = FALSE)

colnames(titles) <- c("tag", "title")

im <- initInterMine(mine=listMines()["PhytoMine"])

cl <- parallel::makeForkCluster(24)
doParallel::registerDoParallel(cl)

for(i in 1:nrow(data)){
  print(i)
  this_trait <- data[i, ]

  p_g_list <- plot_stack(this_trait,
                         outdir = outdir)

  this_trait <- na.omit(unlist(this_trait))
  this_raw_trait_name <- this_trait['raw_trait']
  this_raw_trait_path <- this_trait['raw_trait_path']
  these_results_this_trait <- this_trait[3:length(this_trait)]

  plot_list <- p_g_list[[1]]

  p_stack <- ggarrange(plotlist = plot_list, ncol = 1, heights = 25)

  gene_dataframe_list <- p_g_list[[2]]
  
  # if(debug_counter == 1){
  #   browser()
  # }

  gene_dataframe_combined <- dplyr::bind_rows(gene_dataframe_list)
  
  output_dir <- paste0(outdir, "MethodsStacked/")

  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  p_stack <- annotate_figure(p_stack,
                             top = text_grob(paste0("Trait: ",
                                                    this_raw_trait_name)))

  output_file_prefix <- paste0(
    output_dir,
    this_raw_trait_name)

  this_plot_out_path <- paste0(output_file_prefix,
                               ".png")

  this_gene_table_out_path <- paste0(output_file_prefix,
                                     ".csv")

  ggsave(this_plot_out_path, plot = p_stack,
         width = 18, height = 27)

  fwrite(gene_dataframe_combined, file = this_gene_table_out_path)
  debug_counter <- debug_counter + 1
}

