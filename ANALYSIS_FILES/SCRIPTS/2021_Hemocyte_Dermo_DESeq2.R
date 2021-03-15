
## 2021 Dermo and Hemocyte Transcriptome Analysis ##

# Erin Roberts
# University of Rhode Island

#### load libraries ####
library(DESeq2)  
library(ggplot2)
library(magrittr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(questionr)
library(apeglm)
library(genefilter)
library(fission)
library(tidyr)
library(stringr)
library(rtracklayer)
library(UpSetR)
library(ComplexHeatmap)
library(reshape2)
library(plyr)
library(Repitools)
library(purrr)
library(tibble)
library(ggfortify)
library(ggpubr)
library(viridis)
library(extrafont)
library(limma)
library(data.table)

#### LOADING SAVED GENOME, APOPTOSIS NAMES, IAP XP LISTS ####
Apoptosis_frames <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_apoptosis_products.RData")
annotations <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_annotations.RData")

# Full IAP list with domain type
load(file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/IAP_domain_structure_no_dup_rm.RData")
# load DEG apop list joined with type from IAP script
load(file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/C_vir_C_gig_apop_LFC_IAP_OG_domain_structure")
# Load IAP pathway list 
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/combined_gene_name_org_yes_no_table_unique_pathway_joined.RData")

# load data frames with IAP and GIMAP XM and XP information with haplotigs already collapsed (no domain information)
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_CG_uniq_XP_XM.Rdata")
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_CV_uniq_XP_XM.Rdata")

#### HEMOCYTE TRANSCRIPTOME ANALYSIS ####

## LOAD DATA
Zhang_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Zhang_transcript_count_matrix.csv", header=TRUE,
                         row.names = "transcript_id")
head(Zhang_counts)
colnames(Zhang_counts)

# remove MSTRG novel transcript lines (can assess these later if necessary)
Zhang_counts <- Zhang_counts[!grepl("MSTRG", row.names(Zhang_counts)),]

# Cute the "rna-" from the beginning of rownames
remove_rna = function(x){
  return(gsub("rna-","",x))
}
row.names(Zhang_counts) <- remove_rna(row.names(Zhang_counts))
head(Zhang_counts)

#Load in sample metadata
Zhang_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Zhang_coldata2.csv", row.names = 1 )
View(Zhang_coldata)  
nrow(Zhang_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(Zhang_coldata) %in% colnames(Zhang_counts))  #Should return TRUE
# returns TRUE
all(colnames(Zhang_counts) %in% rownames(Zhang_coldata))  
# returns TRUE
all(rownames(Zhang_coldata) == colnames(Zhang_counts))    # should return TRUE
# returns TRUE


### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
Zhang_counts_matrix <- as.matrix(Zhang_counts)
Zhangrlogcounts <- rlog(Zhang_counts_matrix, blind =TRUE)

# run PCA
pcZhang <- prcomp(t(Zhangrlogcounts))
# plot PCA
autoplot(pcZhang)

# Lets add colour to look at the clustering for Status
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour="condition", 
         size=5) # Vibrio tubiashi and LPS cluster together (somewhat)
# with MSTRG removed, V aes and V alg2 cluster, Mlut, LPS and PBS cluster
# control and PBS clustered together, LPS and M lut clustered together (see this in downstream PCAs too)
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour="group_by_sim", 
         size=5) 

# Is there a time effect? 
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour="time", 
         size=5)  # time may have some effect but because there are no replicates at this time point it is difficult to tell. 
#Going to remove this from the DESeq2 formula
# Where do the different conditions cluster?
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour="condition", 
         size=5)
# Plot PCA 2 and 3 for comparison
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour = "condition", 
         size = 5,
         x = 2,
         y = 3) # control and PBS clustering here, LPS and V. tubiashi still clustering
# with MSTRG removed, V tub and V aes are closest cluster, they also cluster with LPS
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour = "group_by_sim", 
         size = 5,
         x = 2,
         y = 3)
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour = "condition", 
         size = 5,
         x = 3,
         y = 4) # LPS and V. tubiashi still clustering
# with MSTRG removed, V ang and V alg1 cluster the closest 

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Design here is to combine control and PBS as control and the others as challenge
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Creating three here so I can compare the results
Zhang_dds <- DESeqDataSetFromMatrix(countData = Zhang_counts,
                                    colData = Zhang_coldata,
                                    design = ~time + group_by_sim) # add time to control for injection and time effect

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
Zhang_dds <- Zhang_dds[ rowSums(counts(Zhang_dds)) > 10, ]

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(Zhang_coldata$condition) #"Control" "LPS"     "M_lut"   "PBS"     "V_aes"   "V_alg_1" "V_alg_2" "V_ang"   "V_tub"  
levels(Zhang_coldata$group_by_sim) # "control"             "LPS_M_lut"           "V_aes_V_alg1_V_alg2" "V_tub_V_ang"    
levels(Zhang_coldata$time) # 12hr, no injection
Zhang_dds$time <- factor(Zhang_dds$time , levels = c("No_injection","12h"))

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
Zhang_dds_rlog <- rlog(Zhang_dds, blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(Zhang_dds_rlog, intgroup=c("group_by_sim", "condition"))
# control and PBS still clustering, LPS M. lut and V. aes clustering
# with MSTRG removed, the LPS and M. lut cluster most closely 

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene, 
# and fitting a generalized linear model.

Zhang_dds_deseq <- DESeq(Zhang_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Zhang_dds_deseq) # [1] "Intercept"                                   "time_12h_vs_No_injection"                   
#[3] "group_by_sim_LPS_M_lut_vs_control"           "group_by_sim_V_aes_V_alg1_V_alg2_vs_control"
#[5] "group_by_sim_V_tub_V_ang_vs_control" 

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
Zhang_dds_deseq_res_V_alg1 <- results(Zhang_dds_deseq, alpha=0.05, name = "group_by_sim_V_aes_V_alg1_V_alg2_vs_control"  )
Zhang_dds_deseq_res_V_tub <- results(Zhang_dds_deseq, alpha=0.05, name= "group_by_sim_V_tub_V_ang_vs_control" )
Zhang_dds_deseq_res_LPS <- results(Zhang_dds_deseq, alpha=0.05, name= "group_by_sim_LPS_M_lut_vs_control")

head(Zhang_dds_deseq_res_V_alg1) #  group by sim Vibrio vs control 
head(Zhang_dds_deseq_res_V_tub) # group by sim LPS M lut Vtub vs control 
head(Zhang_dds_deseq_res_LPS) # group by sim V aes v alg2 vs control 

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.

# NOTES from Michael love on Lfcshrinkage (https://support.bioconductor.org/p/77461/): 
# https://support.bioconductor.org/p/110307/ # very helpful distinction between lfcestimate and lfc shrinkage
# The difference between results() and lfcShrink() is that the former does not provide fold change shrinkage. 
# The latter function calls results() internally to create the p-value and adjusted p-value columns, 
# which provide inference on the maximum likelihood LFC. The shrunken fold changes are useful for ranking genes by 
# effect size and for visualization.
# The shrinkage is generally useful, which is why it is enabled by default. Full methods are described in the DESeq2 paper (see DESeq2 citation),
# but in short, it looks at the largest fold changes that are not due to low counts and uses these to inform a prior distribution. 
# So the large fold changes from genes with lots of statistical information are not shrunk, while the imprecise fold changes are shrunk. 
# This allows you to compare all estimated LFC across experiments, for example, which is not really feasible without the use of a prior.
# THE lfcshrinkage is not Affecting the p values at all, but its just shrinking the log2 fold change and calculating a new standard error for it 
# https://support.bioconductor.org/p/95695/

# Notes on setting up coefficients for apeglm, https://support.bioconductor.org/p/115435/ , https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
# Although apeglm cannot be used with contrast, we note that many designs can be easily rearranged such that what was a contrast becomes its own coefficient.
# In this case, the dispersion does not have to be estimated again, as the designs are equivalent, up to the meaning of the coefficients. 
# Instead, one need only run nbinomWaldTest to re-estimate MLE coefficients – these are necessary for apeglm – and then run lfcShrink specifying 
# the coefficient of interest in resultsNames(dds)
# The user would for example, either change the levels of dds$condition or replace the design using design(dds)<-, then run nbinomWaldTest followed by lfcShrink

# For each LFCshrink I can pass to it my res object for each so that I can keep my alpha setting at 0.05. Doing this procedure will 
# keep the p-values and padj from the results() call, and simply update the LFCs so they are posterior estimates.

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm

Zhang_dds_deseq_res_V_alg1_LFC <- lfcShrink(Zhang_dds_deseq, coef="group_by_sim_V_aes_V_alg1_V_alg2_vs_control" , type= "apeglm", res=Zhang_dds_deseq_res_V_alg1)
# Review results object summary
summary(Zhang_dds_deseq_res_V_alg1_LFC)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 1149, 3.4%
#LFC < 0 (down)     : 173, 0.52%
#outliers [1]       : 0, 0%
#low counts [2]     : 5831, 17%
#(mean count < 5)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

Zhang_dds_deseq_res_V_tub_LFC <- lfcShrink(Zhang_dds_deseq, coef="group_by_sim_V_tub_V_ang_vs_control" , type= "apeglm", res=Zhang_dds_deseq_res_V_tub)
summary(Zhang_dds_deseq_res_V_tub_LFC)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 665, 2%
#LFC < 0 (down)     : 217, 0.65%
#outliers [1]       : 0, 0%
#low counts [2]     : 9719, 29%
#(mean count < 10)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

Zhang_dds_deseq_res_LFC_LPS <- lfcShrink(Zhang_dds_deseq, coef="group_by_sim_LPS_M_lut_vs_control", type= "apeglm", res=Zhang_dds_deseq_res_LPS)
summary(Zhang_dds_deseq_res_LFC_LPS)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 620, 1.9%
#LFC < 0 (down)     : 213, 0.64%
#outliers [1]       : 0, 0%
#low counts [2]     : 11662, 35%
#(mean count < 13)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
#plotMA(Zhang_dds_deseq_res_LFC, ylim = c(-5, 5))
## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
#hist(Zhang_dds_deseq_res_LFC$padj[Zhang_dds_deseq_res_LFC$baseMean > 1], breaks = 0:20/20,
#   col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Zhang_dds_deseq_res_V_alg1_LFC_sig <- subset(Zhang_dds_deseq_res_V_alg1_LFC, padj < 0.05)
Zhang_dds_deseq_res_V_alg1_LFC_sig $transcript_id <- row.names(Zhang_dds_deseq_res_V_alg1_LFC_sig )
Zhang_dds_deseq_res_V_alg1_LFC_sig  <- as.data.frame(Zhang_dds_deseq_res_V_alg1_LFC_sig )
nrow(Zhang_dds_deseq_res_V_alg1_LFC_sig )  #1322

Zhang_dds_deseq_res_V_tub_LFC_sig <- subset(Zhang_dds_deseq_res_V_tub_LFC, padj < 0.05)
Zhang_dds_deseq_res_V_tub_LFC_sig  $transcript_id <- row.names(Zhang_dds_deseq_res_V_tub_LFC_sig  )
Zhang_dds_deseq_res_V_tub_LFC_sig   <- as.data.frame(Zhang_dds_deseq_res_V_tub_LFC_sig )
nrow(Zhang_dds_deseq_res_V_tub_LFC_sig  )  #882

Zhang_dds_deseq_res_LFC_LPS_sig <- subset(Zhang_dds_deseq_res_LFC_LPS, padj < 0.05)
Zhang_dds_deseq_res_LFC_LPS_sig $transcript_id <- row.names(Zhang_dds_deseq_res_LFC_LPS_sig )
Zhang_dds_deseq_res_LFC_LPS_sig  <- as.data.frame(Zhang_dds_deseq_res_LFC_LPS_sig)
nrow(Zhang_dds_deseq_res_LFC_LPS_sig)  #833

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons
# topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 40)
# mat  <- assay(vsd)[ topVarGenes, ]
# mat  <- mat - rowMeans(mat)
# anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
# pheatmap(mat, annotation_col = anno)

topVarGenes_Zhang_dds_rlog <-  head(order(rowVars(assay(Zhang_dds_rlog )), decreasing = TRUE), 100)
family_Zhang_mat <- assay(Zhang_dds_rlog)[topVarGenes_Zhang_dds_rlog,]
family_Zhang_mat <- family_Zhang_mat - rowMeans(family_Zhang_mat)
family_Zhang_anno <- as.data.frame(colData(Zhang_dds_rlog)[, c("condition","group_by_sim")])
family_Zhang_heatmap <- pheatmap(family_Zhang_mat , annotation_col = family_Zhang_anno)
head(family_Zhang_mat)
# With MSTRG removed, the control PBS and control control cluster together, then the LPS, M lut and V alg1 cluster
# and the other Vibrio challenges cluster 

# reorder annotation table to match ordering in heatmap 
family_Zhang_heatmap_reorder <-rownames(family_Zhang_mat[family_Zhang_heatmap$tree_row[["order"]],])
# annotate the row.names
family_Zhang_mat_prot <- as.data.frame(family_Zhang_heatmap_reorder)
colnames(family_Zhang_mat_prot)[1] <- "transcript_id"
family_Zhang_mat_prot_annot <- left_join(family_Zhang_mat_prot, select(C_gig_rtracklayer_transcripts, transcript_id, product, gene), by = "transcript_id")

# Gene clustering heatmap with only apoptosis genes #
# vector C_gig_apop transcript IDs
C_gig_rtracklayer_apop_product_final_transcript_id <- as.vector(C_gig_rtracklayer_apop_product_final$transcript_id)
# Search original Zhang_counts for apoptosis genes and do rlog on just these
Zhang_counts_apop <- Zhang_counts[row.names(Zhang_counts) %in% C_gig_rtracklayer_apop_product_final_transcript_id,]
nrow(Zhang_counts_apop) #833
head(Zhang_counts_apop)

Zhang_counts_apop_dds <- DESeqDataSetFromMatrix(countData = Zhang_counts_apop,
                                                colData = Zhang_coldata,
                                                design = ~time + group_by_sim) # add time to control for injection and time effect
# Prefiltering the data and running rlog
Zhang_counts_apop_dds <- Zhang_counts_apop_dds[ rowSums(counts(Zhang_counts_apop_dds)) > 10, ]
Zhang_counts_apop_rlog <- rlog(Zhang_counts_apop_dds, blind=TRUE)

## PCA plot of rlog transformed counts for apoptosis
plotPCA(Zhang_counts_apop_rlog , intgroup=c("group_by_sim", "condition"))

# heatmap of all apoptosis genes 
Zhang_counts_apop_assay <-  assay(Zhang_counts_apop_rlog)[,]
Zhang_counts_apop_assay_mat <- Zhang_counts_apop_assay - rowMeans(Zhang_counts_apop_assay)
Zhang_counts_apop_assay_anno <- as.data.frame(colData(Zhang_counts_apop_rlog )[, c("condition","group_by_sim")])
Zhang_counts_apop_assay_heatmap <- pheatmap(Zhang_counts_apop_assay_mat  , annotation_col = Zhang_counts_apop_assay_anno)
head(Zhang_counts_apop_assay_mat )
# control and PBS grouping, V_aes and V_ alg2 grouping, V_ang V alg 1 and V tub clustering

# heatmap of most variable apoptosis genes (this selects genes with the greatest variance in the sample)
topVarGenes_Zhang_counts_apop_assay <-  head(order(rowVars(assay(Zhang_counts_apop_rlog)), decreasing = TRUE), 100) 
top_Var_Zhang_counts_apop_assay_mat<- assay(Zhang_counts_apop_rlog )[topVarGenes_Zhang_counts_apop_assay,]
top_Var_Zhang_counts_apop_assay_mat <- top_Var_Zhang_counts_apop_assay_mat - rowMeans(top_Var_Zhang_counts_apop_assay_mat)
top_Var_Zhang_counts_apop_assay_anno <- as.data.frame(colData(Zhang_counts_apop_rlog )[, c("condition","group_by_sim")])
top_Var_Zhang_counts_apop_assay_heatmap <- pheatmap(top_Var_Zhang_counts_apop_assay_mat  , annotation_col = top_Var_Zhang_counts_apop_assay_anno)
head(top_Var_Zhang_counts_apop_assay_mat )
# same grouping as above with top 200, top 100 changes grouping M lut LPS V alg1 and V alg2 group together in a cluster, while V tub, V ang and V aes group with control

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_Zhang_counts_apop_assay_heatmap_reorder <-rownames(top_Var_Zhang_counts_apop_assay_mat[top_Var_Zhang_counts_apop_assay_heatmap $tree_row[["order"]],])
# annotate the row.names
top_Var_Zhang_counts_apop_assay_prot <- as.data.frame(top_Var_Zhang_counts_apop_assay_heatmap_reorder)
colnames(top_Var_Zhang_counts_apop_assay_prot)[1] <- "transcript_id"
top_Var_Zhang_counts_apop_assay_prot_annot <- left_join(top_Var_Zhang_counts_apop_assay_prot, select(C_gig_rtracklayer_transcripts, transcript_id, product, gene), by = "transcript_id")
#isolate interesting clusters

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP <- merge(Zhang_dds_deseq_res_V_alg1_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP_arranged <- arrange(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP , -log2FoldChange) 
nrow(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP) #26

Zhang_dds_deseq_res_V_tub_LFC_sig_APOP <- merge(Zhang_dds_deseq_res_V_tub_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_deseq_res_V_tub_LFC_sig_APOP_arranged <- arrange(Zhang_dds_deseq_res_V_tub_LFC_sig_APOP  , -log2FoldChange) 
nrow(Zhang_dds_deseq_res_V_tub_LFC_sig_APOP) # 18

Zhang_dds_deseq_res_LFC_LPS_sig_APOP <- merge(Zhang_dds_deseq_res_LFC_LPS_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_deseq_res_LFC_LPS_sig_APOP_arranged <- arrange(Zhang_dds_deseq_res_LFC_LPS_sig_APOP , -log2FoldChange) 
nrow(Zhang_dds_deseq_res_LFC_LPS_sig_APOP )  #20

# Compare apoptosis genes between group_by_sim groups
Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP$group_by_sim <- "V_aes_V_alg1_V_alg2"
Zhang_dds_deseq_res_V_tub_LFC_sig_APOP$group_by_sim <- "V_tub_V_ang"
Zhang_dds_deseq_res_LFC_LPS_sig_APOP$group_by_sim <- "LPS_M_lut"

# combine data frames 
Zhang_upset_all_sig_APOP <- rbind(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                                  Zhang_dds_deseq_res_V_tub_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                                  Zhang_dds_deseq_res_LFC_LPS_sig_APOP[,c("product","group_by_sim","log2FoldChange")])

Zhang_upset_all_sig_APOP_joined <- full_join( Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP[,c("product","group_by_sim")], Zhang_dds_deseq_res_V_tub_LFC_sig_APOP[,c("product","group_by_sim")], by ="product")
Zhang_upset_all_sig_APOP_joined <- full_join(Zhang_upset_all_sig_APOP_joined , Zhang_dds_deseq_res_LFC_LPS_sig_APOP[,c("product","group_by_sim")], by ="product")
Zhang_upset_all_sig_APOP_joined <- na.omit(Zhang_upset_all_sig_APOP_joined)
nrow(Zhang_upset_all_sig_APOP_joined) # 13 shared between all

# Convert into wide format using reshape
Zhang_upset_all_sig_APOP_tally <- Zhang_upset_all_sig_APOP %>% group_by(product) %>% tally() 
Zhang_upset_all_sig_APOP_upset <- Zhang_upset_all_sig_APOP %>% group_by(product) %>% mutate(value=1) %>% spread(group_by_sim, value, fill =0 )
Zhang_upset_all_sig_APOP_upset <- as.matrix(Zhang_upset_all_sig_APOP_upset)

# Make plot
Zhang_full_LFC_plot <- ggplot(Zhang_upset_all_sig_APOP, aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

Zhang_dds_deseq_res_V_alg1_APOP_plot <- ggplot(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Zhang V. alg, V. aes vs Control") +
  ylab("Log2 Fold Change")
Zhang_dds_deseq_res_V_tub_APOP_plot <- ggplot(Zhang_dds_deseq_res_V_tub_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("V. tub V. ang vs Control") +
  ylab("Log2 Fold Change")
Zhang_dds_deseq_res_LPS_APOP_plot <- ggplot(Zhang_dds_deseq_res_LFC_LPS_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Zhang LPS and M Lut LFC vs Control") +
  ylab("Log2 Fold Change")

ggarrange(Zhang_dds_deseq_res_V_alg1_APOP_plot , Zhang_dds_deseq_res_V_tub_APOP_plot,Zhang_dds_deseq_res_LPS_APOP_plot)

## Extract IAP and GIMAP specific manually curated protein lists
Zhang_dds_deseq_res_V_alg1_LFC_sig_IAP <- merge(Zhang_dds_deseq_res_V_alg1_LFC_sig, BIR_XP_gff_CG_uniq_XP_XM, by = "transcript_id")
Zhang_dds_deseq_res_V_alg1_LFC_sig_GIMAP <- merge(Zhang_dds_deseq_res_V_alg1_LFC_sig, AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM, by = "transcript_id")

Zhang_dds_deseq_res_V_tub_LFC_sig_IAP <- merge(Zhang_dds_deseq_res_V_tub_LFC_sig , BIR_XP_gff_CG_uniq_XP_XM, by = "transcript_id")
Zhang_dds_deseq_res_V_tub_LFC_sig_GIMAP <- merge(Zhang_dds_deseq_res_V_tub_LFC_sig , AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM, by = "transcript_id")

Zhang_dds_deseq_res_LFC_LPS_sig_IAP <- merge(Zhang_dds_deseq_res_LFC_LPS_sig, BIR_XP_gff_CG_uniq_XP_XM, by = "transcript_id")
Zhang_dds_deseq_res_LFC_LPS_sig_GIMAP <- merge(Zhang_dds_deseq_res_LFC_LPS_sig, AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM, by = "transcript_id")

