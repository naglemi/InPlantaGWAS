find_genes <- function(annotations, GWAS_hits, distance_threshold){
  all_annotations <- foreach(chr = levels(factor(GWAS_hits$CHR)),
                             .combine = rbind) %do% {

                               annotations_this_chr <- annotations[which(annotations$seqid == chr), ]
      
                               GWAS_hits_this_chr <- GWAS_hits[which(
                                 GWAS_hits$CHR == chr), ]
      
                               GWAS_hits_this_chr <- 
                                 GWAS_hits_this_chr[order(
                                   GWAS_hits_this_chr$PVAL, decreasing = TRUE), ]
      
                               for(i in 1:nrow(annotations_this_chr)){
                                 #for(i in 1:5){
                                 
                                 #print(paste("Start position of this annotation is",
                                 #            annotations_this_chr$start[i]))
                                 
                                 distance_from_hits <- 
                                   annotations_this_chr$start[i] -
                                   GWAS_hits_this_chr$POS
                                 
                                 closest_hit <- 
                                   GWAS_hits_this_chr[which(
                                     abs(distance_from_hits) == min((abs(distance_from_hits)))), ]
                                   
                                 # This message floods the console
#                                  if(length(closest_hit$POS) > 1){
#                                    print(paste("On chr", chr))
#                                    print(annotations_this_chr$start[i])
#                                    print(closest_hit)
#                                    message("Note in this situation the start site of a gene is exactly in the middle of two trait-associated loci. We will select the one with the lower p-value.")
#                                  }
                                 
                                 annotations_this_chr$closest_hit_position[i] <- closest_hit$POS[1]
                                 annotations_this_chr$closest_hit_pval[i] <- closest_hit$PVAL[1]
                                 
                                 #print(paste("Hit closest to this annotation is", closest_hit_position))
                               }
                               
                               annotations_this_chr
                             }
  
  all_genes <- all_annotations[which(all_annotations$sequence_type == "gene"), ]
  all_genes$distance <- abs(all_genes$start - all_genes$closest_hit_position)
  all_genes_near_hits <- all_genes[which(all_genes$distance < distance_threshold), ]
  
  attributes <- stringr::str_split_fixed(all_genes_near_hits$attributes,
                                         ";",
                                         3)
  
  attributes <- as.data.frame(attributes)
  colnames(attributes) <- c("ID", "Name", "ancestorIdentifier")
  attributes$Name <- gsub("Name=", "", attributes$Name)
  attributes$link <- paste0("http://browser.planteome.org/api/search/annotation?bioentity=",
                            attributes$Name)
  
  all_genes_near_hits <- cbind(all_genes_near_hits,
                               attributes$Name,
                               attributes$link)
  
  all_genes_near_hits
}

pull_gene_info <- function(GWAS_genes, im, verbose=TRUE){
    
    GWAS_genes$gene_info <- GWAS_genes$latest_transcript <- rep(NA, nrow(GWAS_genes))

    queryNames = getTemplateQuery(
      im = im, 
      name = "Gene_info"
    )

    for(i in 1:nrow(GWAS_genes)){

        # Build and run query to Phytomine
        this_name <- as.character(GWAS_genes$`attributes$Name`[i])
        
        if(verbose == TRUE){
            message(paste("Pulling gene info data for", this_name))
        }
        
        queryNames$where[[1]]$value <- this_name
        name_query_results <- runQuery(im, queryNames)

        # Parse data and add it to GWAS gene dataframe
        name_query_results$version <- stringr::str_split_fixed(name_query_results$`Gene.transcripts.primaryIdentifier`, "\\.", 3)[,3]
        name_query_results <- name_query_results[which(name_query_results$Gene.organism.shortName == "P. trichocarpa"), ]

        if(nrow(name_query_results) == 0){
            name_latest_transcript_annotation <- NA
            if(verbose == TRUE){
                message(paste("No gene info data for", this_name))
            }

            name_latest_transcript <- NA
        }
        if(nrow(name_query_results) > 0){
            if(verbose == TRUE){
                message("Success.")
            }
            name_latest_transcript_annotation <-
              name_query_results$Gene.briefDescription[which(name_query_results$version == max(name_query_results$version))]

            name_latest_transcript <-
              name_query_results$Gene.transcripts.primaryIdentifier[which(name_query_results$version == max(name_query_results$version))]

        }

        GWAS_genes$gene_info[i] <- name_latest_transcript_annotation
        GWAS_genes$latest_transcript[i] <- name_latest_transcript
    }
    
    GWAS_genes
}

pull_homologs <- function(GWAS_genes, im, database = "Pull"){
    # Build and run query

    if(database == "Pull"){
        queryHomologs = getTemplateQuery(
            im = im, 
            name = "Two_Species_Homologs")
        
        queryHomologs$where[[1]]$value <- "P. trichocarpa"
        queryHomologs$where[[2]]$value <- "A. thaliana"
        
        all_homologs <- runQuery(im, queryHomologs)
    }


    if(database != "Pull"){
        all_homologs <- read.csv(database)
    }

    # Add homolog data for each gene to our database of GWAS hits

    GWAS_genes$homologs <- rep(NA, nrow(GWAS_genes))
    
    all_homologs$Homolog.gene.primaryIdentifier <- as.character(all_homologs$Homolog.gene.primaryIdentifier)
    GWAS_genes$`attributes$Name` <- as.character(GWAS_genes$`attributes$Name`)

    for(i in 1:nrow(GWAS_genes)){
      #print(i)
      all_homologs_this_gene <- all_homologs$Homolog.ortholog_gene.primaryIdentifier[which(
        all_homologs$Homolog.gene.primaryIdentifier == GWAS_genes$`attributes$Name`[i])]

      GWAS_genes$homologs[i] <- paste(all_homologs_this_gene,
                                      collapse = ";")
    }

    GWAS_genes
}

pull_descriptions <- function(GWAS_genes, im){
    queryGO = getTemplateQuery(
      im = im, 
      name = "GO_terms_for_transcript"
    )

    GWAS_genes$GO_descriptions <- rep(NA, nrow(GWAS_genes))

    for(i in 1:nrow(GWAS_genes)){
        this_transcript <- as.character(GWAS_genes$`latest_transcript`[i])
        queryGO$where[[1]]$value <- this_transcript

        GO_query_results <- runQuery(im, queryGO)

        # Parse data
        descriptions <- GO_query_results$Gene.ontologyAnnotations.ontologyTerm.description
        descriptions <- descriptions[descriptions != ""]
        descriptions <- paste(descriptions,
                              collapse = "; ")

        # Make more readable
        descriptions <- gsub(";", "; ", descriptions)

        # Add to df
        GWAS_genes$GO_descriptions[i] <- descriptions
    }
    
    GWAS_genes
    
}

annotate_genes <- function(GWAS_genes, im){
    message("Pulling gene info")
    GWAS_genes <- pull_gene_info(GWAS_genes, im, verbose = FALSE)
    message("Pulling gene homologs")
    GWAS_genes <- pull_homologs(GWAS_genes, im, database = "Potri_Attha_homologs.csv")
    message("Pulling gene descriptions")
    GWAS_genes <- pull_descriptions(GWAS_genes, im)
    
    GWAS_genes
}