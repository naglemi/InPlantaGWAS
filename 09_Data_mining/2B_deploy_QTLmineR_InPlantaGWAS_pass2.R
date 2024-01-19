getwd()
setwd("/mnt/data/NSF_GWAS/notebooks/InPlantaGWAS/11_Data_mining/")

library(data.table)
library(qqman)
library(ggplot2)
library(foreach)
library(ggpubr)
library(InterMineR)
library(tools)
library(dplyr)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("InterMineR")

data <- fread("1_Table_results_by_method_a2.csv", na.strings="")
data <- data[, 1:29] # only keep GEMMA and GMMAT for second pass

batches_to_exclude <- c("GEMMA_boxcox_outa1",
                        "GMMAT_threshold-0.198412874212136.binary_logitlink-batch3_maf05andmaf01_mixperhaps",
                        "GMMAT_nothreshold_duplicates_binarized.binary_logitlink-batch3_maf05andmaf01_mixperhaps",
                        #"GMMAT_binarized_logitlink-batch2_maf01_geno10",
                        "GMMAT_binarized_logitlink-batch3_maf05andmaf01_mixperhaps")

outdir <- "Results_pass2/"

for(batch in batches_to_exclude){
  data[, batch] <- NA
}

titles <- fread("batch_titles_a2.csv", header = FALSE)
colnames(titles) <- c("tag", "title")

im <- initInterMine(mine=listMines()["PhytoMine"])

cl <- parallel::makeForkCluster(24)
doParallel::registerDoParallel(cl)

#for(i in 9:22){
foreach(i = 1:nrow(data)) %dopar% {
#for(i in 17:nrow(data)){ #37:nrow(data)){
#for(i in 2:3){
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

