load_covariates <- function(path){
    Covar <- fread(path)
    Covar <- Covar[,2:ncol(Covar)]
    Covar <- as.matrix(Covar)
    Covar
}

filter_phenotypes <- function(phenotype, covariates){
    # If any genotypes have NA for any covariate, they shold also have NA for the phenotype being tested.
    indices_missing_value <- which(is.na(covariates[,1]))
    phenotype[indices_missing_value, 3] <- NA
    phenotype
}

load_kinship <- function(kinship_path, labels){
    kinship <- as.matrix(fread(kinship_path))
    colnames(kinship) <- rownames(kinship) <- as.character(labels)
    
    if(!matrixcalc::is.positive.semi.definite(kinship)){
        stop("Error! This kinship matrix is not positive semi-definite.")
    }
    
    kinship
}

score_test <- function(phenotype, kinship, Covar, gds_path, outfile, ncores){
    start_time <- proc.time()
    model0 <- NULL
    
    if (!is.null(Covar)){
        model0 <- glmmkin(V3 ~ as.matrix(Covar), 
                          data = phenotype, kins = kinship, id = "V1",
                          family = binomial(link = logit),
                          verbose = TRUE#,
                          #maxiter = 1000000
                         )
    }
    
    if (is.null(Covar)){
        model0 <- glmmkin(V3 ~ 1, 
                          data = phenotype, kins = kinship, id = "V1",
                          family = binomial(link = logit),
                          verbose = TRUE#,
                          #maxiter = 1000000
                         )
    }
    
    message("Null model complete")
    
    glmm.score(model0, infile = gds_path, outfile = outfile, ncores=ncores)
    
    time_complete <- time.since(start_time)
    
    message(paste("Finished running score tests using",
               ncores, "cores in",
               time_complete,
                 "and wrote results to",
                 outfile))
}

filter_SNPs <- function(x,
                        #missing_threshold = 0.1,
                        contig = TRUE,
                        n_SNPs = 100){

    if(contig == TRUE){ # Filter to SNPs on assembled chromosomes only
        x$CHR <- suppressWarnings(as.numeric(as.character(x$CHR)))
        x <- na.omit(x)
    }

    #x <- x[which(x$MISSRATE < missing_threshold), ]
    
    x <- x[order(x, x$PVAL), ]
    
    if(is.infinite(n_SNPs)){
        n_SNPs <- nrow(x)
    }
    
    #x <- x[which(x$PVAL < pthresh), ]
    
    x <- x[1:min(n_SNPs, nrow(x)), ]

    x <- x[which(x$PVAL <= 1e-4), ]

    print("Head of filtered SNPs: ")
    print(x)
    
    x
}

wald_test <- function(phenotype, Covar, kinship, gds_path, snps, outfile){
    
    start_time <- proc.time()
    
    # We can only pass a whole matrix as a covariate (with clean code that
    #   generalizes to any # of covariates) if we make the covariate matrix
    #   a global variable.


    #print(as.list(.GlobalEnv))
    if (!is.null(Covar)){
        assign("covar", Covar, envir = .GlobalEnv)
        wald_out <- glmm.wald(fixed = V3 ~ covar,
                              kins = kinship,
                              id = "V1",
                              family = binomial(link = "logit"),
                              #infile = "1323_cohort_maf0.05_defaultmissingrates.SeqArrayGDS",
                              infile = gds_path,
                              data = phenotype,
                              #out.fn = outfile,
                              snps = snps)
    }
    if (is.null(Covar)){
        wald_out <- glmm.wald(fixed = V3 ~ 1,
                              kins = kinship,
                              id = "V1",
                              family = binomial(link = "logit"),
                              #infile = "1323_cohort_maf0.05_defaultmissingrates.SeqArrayGDS",
                              infile = gds_path,
                              data = phenotype,
                              #out.fn = outfile,
                              snps = snps)
    }

    time_complete <- time.since(start_time)

    fwrite(wald_out, outfile)

    message(paste("Finished running Wald tests for",
                  length(snps), "SNPs in",
                   time_complete,
                     "and wrote results to",
                     outfile))
    wald_out
}

score_Wald_scatter <- function(wald_out, score_results, phenotype_path,
                               outfile){
    results_merged <- merge(wald_out, score_results, all.x = TRUE, all.y = FALSE,
                        by = c("SNP", "CHR", "POS", "REF", "ALT", "N"))

    nonconverged <- results_merged[which(results_merged$converged_Wald == FALSE), ]

    print(paste("N SNPs for which model did not converge:", nrow(nonconverged)))
    print(nonconverged$SNP)

    converged <- results_merged[which(results_merged$converged_Wald == TRUE), ]

    ggplot(converged, aes(x = -log(p_score_test, base = 10), y = -log(p_Wald, base = 10))) + 
    geom_point() + 
    geom_text_repel(label = stringr::str_split_fixed(converged$SNP, "_", 2)[,1], nudge_y = 0.1) +
    geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
    xlab("p-value from score test (-log scale)") +
    ylab("p-value from Wald test (-log scale)") +
    ggtitle(paste("Comparison of p-values from score test and Wald tests\n",
               "Phenotype: ",
                  tools::file_path_sans_ext(
                      tools::file_path_sans_ext(
                          basename(phenotype_path)))))
    
    ggsave(outfile)
}


GMMAT_manhattan <- function(wald_out, score_results, phenotype_path, outdir, alpha = 0.05, contig = TRUE){
    
    wd <- getwd()
    
    if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    
    single_chr_plot_dir <- paste0(tools::file_path_sans_ext(
                      tools::file_path_sans_ext(
                          basename(phenotype_path))), "_SingleChrPlots")
    
    if(!dir.exists(single_chr_plot_dir)) dir.create(single_chr_plot_dir, recursive = TRUE)
    
    setwd(outdir)
    wald_out <- wald_out[which(wald_out$converged_Wald == TRUE), ]
    results_merged <- merge(wald_out, score_results, all.x = TRUE, all.y = TRUE,
                        by = c("SNP", "CHR", "POS", "REF", "ALT", "N"))
    
    if(contig == TRUE){ # Filter to SNPs on assembled chromosomes only
        results_merged$CHR <- suppressWarnings(as.numeric(as.character(results_merged$CHR)))
        results_merged <- results_merged[which(!is.na(results_merged$CHR)), ]
        #results_merged <- na.omit(results_merged)
    }

    # Do it this way because we want to keep those with NA 
    #   Many have NA since we did score test only and no Wald test if they are not in
    #   top N SNPs
    
    # Keep contiguous chromosomes only

    results_merged.4 <- results_merged[,c(1, 2, 3, 10, 16)]

    FDR <- find_FDR(as.data.table(results_merged.4$p_score_test), alpha)
    BF <- alpha / nrow(results_merged.4)

    # To speed up plotting when we have millions of SNPs
    results_merged.4.subset <- results_merged.4[which(results_merged.4$p_score_test < 0.01), ]

    print("Making whole genome Manhattan plot")
    if(is.na(FDR)){
        CMplot(results_merged.4.subset, plot.type="m",multracks=TRUE,threshold=c(BF),threshold.lty=c(1), 
            threshold.lwd=c(1), threshold.col=c("black"), amplify=FALSE,bin.size=1e6,
            chr.den.col=c("darkgreen", "yellow", "red"), signal.col=NULL,
            signal.cex=1, file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,
            highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4, highlight.col=NULL,
              memo = tools::file_path_sans_ext(
                      tools::file_path_sans_ext(
                          basename(phenotype_path))))
    }

    if(!is.na(FDR)){
        CMplot(results_merged.4.subset, plot.type="m",multracks=TRUE,threshold=c(BF,FDR),threshold.lty=c(1,2), 
            threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=FALSE,bin.size=1e6,
            chr.den.col=c("darkgreen", "yellow", "red"), signal.col=NULL,
            signal.cex=1, file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,
            highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4, highlight.col=NULL,
              memo = tools::file_path_sans_ext(
                      tools::file_path_sans_ext(
                          basename(phenotype_path))))
    }
    
    chromosomes <- levels(factor(results_merged.4.subset$CHR))
    
    setwd(paste0("../../", single_chr_plot_dir))

    for(chromosome in chromosomes){

        this_chr_subset <- results_merged.4.subset[which(results_merged.4.subset$CHR == chromosome), ]
        
        if(nrow(this_chr_subset) > 1){
            print(paste("Making plot for chromosome", chromosome))
            
            if(is.na(FDR)){
            CMplot(this_chr_subset, plot.type="m",multracks=TRUE,threshold=c(BF),threshold.lty=c(1), 
                   threshold.lwd=c(1), threshold.col=c("black"), amplify=FALSE,bin.size=1e6,
                   chr.den.col=c("darkgreen", "yellow", "red"), signal.col=NULL,
                   signal.cex=1, file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,
                   highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4, highlight.col=NULL,
                   memo = paste0(tools::file_path_sans_ext(
                       tools::file_path_sans_ext(
                              basename(phenotype_path))), "Chr", chromosome))
            }
            
            if(!is.na(FDR)){
                CMplot(this_chr_subset, plot.type="m",multracks=TRUE,threshold=c(BF,FDR),threshold.lty=c(1,2), 
                       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=FALSE,bin.size=1e6,
                       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=NULL,
                       signal.cex=1, file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,
                       highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4, highlight.col=NULL,
                       memo = paste0(tools::file_path_sans_ext(
                          tools::file_path_sans_ext(
                              basename(phenotype_path))), "Chr", chromosome))
            }
        }
        

        
    }
    setwd(wd)
}
 
GMMAT_workflow <- function(phenotype_path, covariate_path, kinship_path, gds_path, ncores, outdir, alpha = 0.05, n_SNPs = 100){
    
  Covar <- load_covariates(covariate_path)
  phenotype <- fread(phenotype_path)
  phenotype <- filter_phenotypes(phenotype, covariates = Covar)
  kinship <- load_kinship(kinship_path, labels = phenotype$V1)
  
  output_prefix <- paste0(outdir,
                          file_path_sans_ext(file_path_sans_ext(basename(phenotype_path))),
                          "_logitlink")
  
  score_outfile <- paste0(output_prefix,
                          ".score",
                          ".glmm")
  
  Wald_outfile <- paste0(output_prefix,
                         ".Wald",
                         ".glmm")
  
  scorevswald_outfile <- paste0(output_prefix,
                                ".scorevsWald",
                                ".png")

# Perform score test on all SNPs ------------------------------------------

    message("Starting score test")
    
    if(!file.exists(score_outfile)){
        score_test(phenotype, kinship,
                   Covar, gds_path,
                   score_outfile, ncores)
    }
  


# Filter to top SNPs from score results -----------------------------------

  message("Filtering SNPs")
  
  score_results <- fread(score_outfile)
    
  print("Head of score results: ")
  print(head(score_results))
  score_results_filtered <- filter_SNPs(x = score_results,
                                        #missing_threshold = 0.1,
                                        contig = TRUE,
                                        n_SNPs = n_SNPs)
  snps <- score_results_filtered$SNP
    
  message(paste0("Got ", length(snps), " snps"))
  

# Perform Wald test on top SNPs from score test ---------------------------

  message("Starting Wald test")
    
  Sys.setenv(MKL_NUM_THREADS = 1)
  
  phenotype <- na.omit(phenotype)
    
  # Do we need to reload Covar here to avoid error of
  #  object `Covar` not found? Very strange...
  Covar <- load_covariates(covariate_path)
  Covar <- na.omit(Covar)
  Covar <- Covar[which(!is.na(phenotype$V3)), ]
  
  wald_out <- wald_test(phenotype = phenotype,
                        Covar = Covar,
                        kinship = kinship,
                        gds_path = gds_path,
                        snps = snps,
                        outfile = Wald_outfile)

# Make plots --------------------------------------------------------------

  
#   colnames(score_results)[4:ncol(score_results)] <- 
#     c("REF", "ALT", "N", "MISSRATE", "AF_score_test", "SCORE", "VAR_score_test", "p_score_test")
  
#   colnames(wald_out)[7:11] <- c("AF_Wald", "Beta_Wald", "SE_Wald", "p_Wald", "converged_Wald")
  
  
  
#   score_Wald_scatter(wald_out = wald_out,
#                      score_results = score_results,
#                      phenotype_path = phenotype_path,
#                      outfile = scorevswald_outfile)
  
#   GMMAT_manhattan(wald_out, score_results, phenotype_path = phenotype_path, outdir = outdir)
  
}


GMMAT_workflow.newbeta <- function(phenotype_path, covariate_path, kinship_path, gds_path, ncores, outdir,
                           alpha = 0.05, n_SNPs = 100, pthresh = 1){
    
  Covar <- load_covariates(covariate_path)
  phenotype <- fread(phenotype_path)
  phenotype <- filter_phenotypes(phenotype, covariates = Covar)
  kinship <- load_kinship(kinship_path, labels = phenotype$V1)
  
  output_prefix <- paste0(outdir,
                          file_path_sans_ext(file_path_sans_ext(basename(phenotype_path))),
                          "_logitlink")
  
  score_outfile <- paste0(output_prefix,
                          ".score",
                          ".glmm")
  
  Wald_outfile <- paste0(output_prefix,
                         ".Wald",
                         ".glmm")
  
  scorevswald_outfile <- paste0(output_prefix,
                                ".scorevsWald",
                                ".png")

# Perform score test on all SNPs ------------------------------------------

    if(!file.exists(score_outfile)){
        score_test(phenotype, kinship,
                   Covar, gds_path,
                   score_outfile, ncores)
    }
  


# Filter to top SNPs from score results -----------------------------------

  
  score_results <- fread(score_outfile)
  score_results_filtered <- filter_SNPs.newbeta(x = score_results,
                                        missing_threshold = 0.1,
                                        contig = TRUE,
                                        n_SNPs = n_SNPs,
                                        pthresh = pthresh)
  snps <- score_results_filtered$SNP
  

# Perform Wald test on top SNPs from score test ---------------------------


  Sys.setenv(MKL_NUM_THREADS = 1)
  
  phenotype <- na.omit(phenotype)
    
  # Do we need to reload Covar here to avoid error of
  #  object `Covar` not found? Very strange...
  Covar <- load_covariates(covariate_path)
  Covar <- na.omit(Covar)
  Covar <- Covar[which(!is.na(phenotype$V3)), ]
  
  wald_out <- wald_test(phenotype = phenotype,
                        Covar = Covar,
                        kinship = kinship,
                        gds_path = gds_path,
                        snps = snps,
                        outfile = Wald_outfile)

# Make plots --------------------------------------------------------------

  
#   colnames(score_results)[4:ncol(score_results)] <- 
#     c("REF", "ALT", "N", "MISSRATE", "AF_score_test", "SCORE", "VAR_score_test", "p_score_test")
  
#   colnames(wald_out)[7:11] <- c("AF_Wald", "Beta_Wald", "SE_Wald", "p_Wald", "converged_Wald")
  
  
  
#   score_Wald_scatter(wald_out = wald_out,
#                      score_results = score_results,
#                      phenotype_path = phenotype_path,
#                      outfile = scorevswald_outfile)
  
#   GMMAT_manhattan(wald_out, score_results, phenotype_path = phenotype_path, outdir = outdir)
  
}

GMMAT_workflow_modified <- function(phenotype_path, covariate_path, kinship_path, gds_path, ncores,
                           alpha = 0.05, n_SNPs = 100){
    
  Covar <- load_covariates(covariate_path)
  phenotype <- fread(phenotype_path)
  phenotype <- filter_phenotypes(phenotype, covariates = Covar)
  kinship <- load_kinship(kinship_path, labels = phenotype$V1)
  
  output_prefix <- paste0("Results/batch2_maf01_geno10/",
                          file_path_sans_ext(file_path_sans_ext(basename(phenotype_path))),
                          "_logitlink")
  
  score_outfile <- paste0(output_prefix,
                          ".score",
                          ".glmm")
  
  Wald_outfile <- paste0(output_prefix,
                         ".Wald",
                         ".glmm")
  
  scorevswald_outfile <- paste0(output_prefix,
                                ".scorevsWald",
                                ".png")

# Perform score test on all SNPs ------------------------------------------

  # Was done earlier. Now we are just doing wald test for top 1k SNPs.
  #score_test(phenotype, kinship, Covar, gds_path, score_outfile, ncores)

# Filter to top SNPs from score results -----------------------------------

  
  score_results <- fread(score_outfile)
  score_results_filtered <- filter_SNPs(x = score_results,
                                        #missing_threshold = 0.1,
                                        contig = TRUE,
                                        n_SNPs = n_SNPs,
                                        pthresh = pthresh)
  snps <- score_results_filtered$SNP
    
  print(paste0("Running Wald test for ", nrow(score_results_filtered), " SNPs below ", pthresh))
  

# Perform Wald test on top SNPs from score test ---------------------------


  Sys.setenv(MKL_NUM_THREADS = 1)
  
  phenotype <- na.omit(phenotype)
    
  # Do we need to reload Covar here to avoid error of
  #  object `Covar` not found? Very strange...
  Covar <- load_covariates(covariate_path)
  Covar <- na.omit(Covar)
  Covar <- Covar[which(!is.na(phenotype$V3)), ]
  
  wald_out <- wald_test(phenotype = phenotype,
                        Covar = Covar,
                        kinship = kinship,
                        gds_path = gds_path,
                        snps = snps,
                        outfile = Wald_outfile)

# Make plots --------------------------------------------------------------

  
#   colnames(score_results)[4:ncol(score_results)] <- 
#     c("REF", "ALT", "N", "MISSRATE", "AF_score_test", "SCORE", "VAR_score_test", "p_score_test")
  
#   colnames(wald_out)[7:11] <- c("AF_Wald", "Beta_Wald", "SE_Wald", "p_Wald", "converged_Wald")
  
  
  
#   score_Wald_scatter(wald_out = wald_out,
#                      score_results = score_results,
#                      phenotype_path = phenotype_path,
#                      outfile = scorevswald_outfile)
  
#   GMMAT_manhattan(wald_out, score_results, phenotype_path = phenotype_path)
  
}

GMMAT_workflow_no_cov <- function(phenotype_path, kinship_path, gds_path, ncores,
                           alpha = 0.05, n_SNPs = 100){
    
  phenotype <- fread(phenotype_path)
  colnames(phenotype)[3] <- "V3" # So formula works easily
  kinship <- load_kinship(kinship_path, labels = phenotype$V1)
  
  output_prefix <- paste0(outdir,
                          file_path_sans_ext(file_path_sans_ext(basename(phenotype_path))),
                          "_logitlink")
  
  score_outfile <- paste0(output_prefix,
                          ".score",
                          ".glmm")
  
  Wald_outfile <- paste0(output_prefix,
                         ".Wald",
                         ".glmm")
  
  scorevswald_outfile <- paste0(output_prefix,
                                ".scorevsWald",
                                ".png")

# Perform score test on all SNPs ------------------------------------------

  
  score_test(phenotype, kinship, Covar = NULL, gds_path, score_outfile, ncores)

# Filter to top SNPs from score results -----------------------------------

  
  score_results <- fread(score_outfile)
  score_results_filtered <- filter_SNPs(x = score_results,
                                        missing_threshold = 0.1,
                                        contig = TRUE,
                                        n_SNPs = n_SNPs)
  snps <- score_results_filtered$SNP
  

# Perform Wald test on top SNPs from score test ---------------------------


  Sys.setenv(MKL_NUM_THREADS = 1)
  
  phenotype <- na.omit(phenotype)
  
  wald_out <- wald_test(phenotype = phenotype,
                        Covar = NULL,
                        kinship = kinship,
                        gds_path = gds_path,
                        snps = snps,
                        outfile = Wald_outfile)

# Make plots --------------------------------------------------------------

  
  colnames(score_results)[4:ncol(score_results)] <- 
    c("REF", "ALT", "N", "MISSRATE", "AF_score_test", "SCORE", "VAR_score_test", "p_score_test")
  
  colnames(wald_out)[7:11] <- c("AF_Wald", "Beta_Wald", "SE_Wald", "p_Wald", "converged_Wald")
  
  
  
  score_Wald_scatter(wald_out = wald_out,
                     score_results = score_results,
                     phenotype_path = phenotype_path,
                     outfile = scorevswald_outfile)
  
  GMMAT_manhattan(wald_out, score_results, phenotype_path = phenotype_path)
  
}
