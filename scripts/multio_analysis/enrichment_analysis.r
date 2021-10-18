# enrichment_analysis.r
# Ashley Mae Conard
# Last Modified: 03/21
# Performs gene ontology and gene set enrichment analysis with clusterProfiler.         
# Ref: altered from https://github.com/ashleymaeconard/TIMEOR
# TODO: update R ts 4.0, use argparse for R

# Checking input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Type: Rscript enrichment_analysis.r 
            1) /PATH/TO/INPUT_FILE
            2) EXPERIMENT_NAME (e.g. test_results, write 'test') 
            3) OUTDIR
            4) ORGANISM (dme (Drosophila melanogaster), hsa (Homo sapiens) or mmu (Mus musculus)
            5) ADJ_PVAL (recommend 0.05)", call.=FALSE)
} else if (length(args) == 5) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 1) /PATH/TO/INPUT_FILE
                  2) EXPERIMENT_NAME (e.g. test_results, write 'test') 
                  3) OUTDIR
                  4) ORGANISM (dme (Drosophila melanogaster), hsa (Homo sapiens) or mmu (Mus musculus)
                  5) ADJ_PVAL (recommend 0.05))")
}

# Assigning input arguments 
IN_OUTPUT <- args[1]
EXPERIMENT_NAME <- args[2]
OUTDIR <- args[3] # 0 or 1 
ORGANISM <- args[4] # dme (Drosophila melanogaster), hsa (Homo sapiens), mmu (Mus musculus)
options(digits=5)
ADJ_PVAL <- as.double(args[5]) # recommend 0.05
geneNumClust <- "1" # TODO: update to remove this

# Assigning organism library
if(ORGANISM=="dme"){
	    ORG_DB="org.Dm.eg.db"
} else if(ORGANISM=="hse"){
	    ORGANISM = "hsa"
    ORG_DB="org.Hs.eg.db"
}else if(ORGANISM=="mus"){
	    ORGANISM = "mmu"
    ORG_DB="org.Mm.eg.db"
} else{
	    stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
}

# Loading packages
library(clusterProfiler)
library(ORG_DB, character.only = TRUE) # organism database library

# Create output directory as needed
dir.create(file.path(dirname(OUTDIR), basename(OUTDIR)), showWarnings = FALSE)

assess_enrichment <- function(geneENS, type_enr){
    ### Calculates GO over-representation (enrichement) test for the input gene list. ###
    ego <- enrichGO(gene          = geneENS,
                    OrgDb         = ORG_DB,
                    keyType       = "ENSEMBL",
                    ont           = type_enr,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = ADJ_PVAL,
                    qvalueCutoff  = ADJ_PVAL)                
    ego <- setReadable(ego, OrgDb = ORG_DB)
    return(ego)
}

assess_gsea <- function(geneENS, type_enr){
    ### Performs gene set enrichment analysis for the input gene list. ###
    gse <- gseGO(gene         = geneENS,
                ont           = type_enr, 
                keyType       = "ENSEMBL", 
                nPerm         = 10000, 
                minGSSize     = 3, 
                maxGSSize     = 800, 
                pvalueCutoff  = ADJ_PVAL, 
                verbose       = TRUE, 
                OrgDb         = ORG_DB, 
                pAdjustMethod = "BH") 
    gse <- setReadable(gse, OrgDb = ORG_DB)
    return(gse)
}

generate_plots <- function(egoo, gl, ex, type_enrich, outdir, num_gene_cluster){
    ### Returns multiple plots of enrichment analysis for the input gene list. ###
    # Checking if the enrichment object is empty, and if so, exit function and script
    if(nrow(egoo)>0){

        # Generating subfolder within specific cluster 
        subdirec <- file.path(outdir, type_enrich) 
        if (!dir.exists(subdirec)){
            dir.create(subdirec)
        } else {
            cat(subdirec, "subdirectory exists.")
        }

        # Plotting GSEA 
        svg(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "gsea_dotplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")))
        plt_dotplot <- dotplot(egsea)
        print(plt_dotplot)
        dev.off()

        # Plotting enrichment dot plot
        svg(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_dotplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")))
        plt_dotplot <- dotplot(egoo)
        print(plt_dotplot)
        dev.off()
        
        # Plotting relationship graph plot
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_emaplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
        plt_emaplot <- emapplot(egoo) 
        print(plt_emaplot)
        dev.off()

        # Plotting barplot 
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_barplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
        plt_bar<- barplot(egoo, showCategory=10)
        print(plt_bar)
        dev.off()
        print("here here")

        # Plotting gene concept network
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_conceptplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 30, height = 15)
        p2 <- cnetplot(egoo, categorySize="pvalue", foldChange=gl)
        plt_cnet <- cnetplot(egoo, foldChange=gl, categorySize="pvalue", colorEdge = TRUE) 
        plt_cnet_cir <- cnetplot(egoo, foldChange=gl, circular = TRUE, categorySize="pvalue", colorEdge = TRUE)
        plt_concept <- plot_grid(plt_cnet, plt_cnet_cir, ncol=2)
        print(plt_concept)
        dev.off()
        print("here here, here")

        # Plotting phylogeny plot (if there is more than 1 enriched GO term)
        if(dim(egoo)[1]>1){
            pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_goplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
            plt_goplot <- goplot(egoo)
            print(plt_goplot)
            dev.off()
        }
    } else{
        print(paste(type_enrich," is empty", sep=""))
    }
}

main <- function(){
    ### Formats input gene list and runs GO analysis and GSEA.###
    dir_file <- IN_OUTPUT
    
    if (! length(dir_file)>1){ # check to make sure only one geneList file per cluster
        #geneNumClust=(sub(".*_([^.]+)\\.csv.*", "\\1", dir_file))
        OUTPUT_DIR <- OUTDIR
        
	# Generating ranked gene list for gene set enrichment analysis (GSEA)
        f <- read.csv(dir_file, sep="\n", header = FALSE)
        geneList<-f[,1]
        names(geneList)<- as.character(f[,0])
        geneList<-sort(geneList, decreasing=TRUE)
        
        # Reading the gene file for each cluster    
        a <- dput(as.character(f[1][,1]))
        
        # Removing ID and/or transcript ID if present at end of gene name
        if (grepl("_", a[1])){
            cat("Removing values after and including '_' potentially at end of gene/transcript name\n")
            xs <- gsub("\\-R[A-Z][_|\\>].*","",a)
        } else if (grepl("-R[A-Z]\\>", a[1])) {
            cat("Removing ID at end of gene name\n")
            xs <- gsub("\\-R[A-Z]\\>.*","",a)
        } else{
            xs <- unique(a)
        }

        # Converting gene names to ensemble and entrezids
        gene = bitr(geneID=xs, fromType="SYMBOL", toType=c("ENSEMBL","ENTREZID"), OrgDb=ORG_DB)
        cat("gene----")
        print(gene)

        # Writing an intermediate file with gene names
        type_enrichment <- c("BP","MF","CC") # Biological Process (BP), Molecular Function (MF), Cellular Component (CC)

        # Plotting GO and gene set enrichment per type (BP, MF, CC)
        for(types in type_enrichment){
            print(paste("Processing", geneNumClust, sep=" "))
            print(paste("Determining",types,"enrichment", sep=" "))
            ego <- assess_enrichment(gene$ENSEMBL, types)
            egsea <- assess_gsea(gene$ENSEMBL, types)
            gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
            
            print(paste("Plotting",types,"enrichment for", EXPERIMENT_NAME, sep=" "))
            generate_plots(ego, geneList, EXPERIMENT_NAME, types, OUTPUT_DIR, geneNumClust)
        }
        
    }
}

if(!interactive()) {
    main()
}
