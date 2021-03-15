
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

# C_vir_rtracklayer_transcripts
C_vir_rtracklayer_transcripts <- C_vir_rtracklayer %>% filter(grepl("rna",ID))

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
hemo_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/Hemocyte_2021_transcript_count_matrix.csv", header=TRUE,
                         row.names = "transcript_id")
colnames(hemo_counts)[1:12] <-  c("1_Dermo_GDC_R1_001",  "1_Dermo_ZVAD_R1_001", "1_Dermo_R1_001"     , "1_control_R1_001" ,   "2_Dermo_GDC_R1_001" , "2_Dermo_ZVAD_R1_001" ,"2_Dermo_R1_001",     
                                                                    "2_control_R1_001"  ,  "3_Dermo_GDC_R1_001" , "3_Dermo_ZVAD_R1_001", "3_Dermo_R1_001"   ,   "3_control_R1_001"   )
head(hemo_counts)
colnames(hemo_counts)

# remove MSTRG novel transcript lines (can assess these later)
hemo_counts <- hemo_counts[!grepl("MSTRG", row.names(hemo_counts)),]

# Cute the "rna-" from the beginning of rownames
remove_rna = function(x){
  return(gsub("rna-","",x))
}
row.names(hemo_counts) <- remove_rna(row.names(hemo_counts))
head(hemo_counts)

#Load in sample metadata
hemo_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/Hemo_Pmar_coldata.csv", row.names = 1 )
View(hemo_coldata)  
nrow(hemo_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
hemo_coldata <- hemo_coldata[colnames(hemo_counts),]  
# make the Pmar group the reference level

all(rownames(hemo_coldata) %in% colnames(hemo_counts))  #Should return TRUE
# returns TRUE
all(colnames(hemo_counts) %in% rownames(hemo_coldata))  
# returns TRUE
all(rownames(hemo_coldata) == colnames(hemo_counts))    # should return TRUE
# returns TRUE

### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
hemo_counts_matrix <- as.matrix(hemo_counts)
hemo_rlog_counts <- rlog(hemo_counts_matrix, blind =TRUE)

# run PCA
pchemo <- prcomp(t(hemo_rlog_counts))
# plot PCA
autoplot(pchemo)

# Lets add colour to look at the clustering for Status
autoplot(pchemo,
         data = hemo_coldata, 
         colour="condition", 
         size=5) 
# little clustering by condition

autoplot(pchemo,
         data = as.data.frame(hemo_coldata), 
         colour="pool", # clustering by pool 
         size=5) 
# PC1 and PC2 don't explain much of the variation, which means that the Pool effect is not too large 

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Creating three here so I can compare the results
hemo_dds <- DESeqDataSetFromMatrix(countData = hemo_counts,
                                    colData = hemo_coldata,
                                    design = ~condition) # only compare by condition

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
hemo_dds <- hemo_dds[ rowSums(counts(hemo_dds)) > 10, ]

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(hemo_coldata$condition)  # control is currently listed first, so this looks good 

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
hemo_dds_rlog <- rlog(hemo_dds, blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(hemo_dds_rlog, intgroup=c("condition")) # a bit more clustering by condition now 

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene, 
# and fitting a generalized linear model.
hemo_dds_deseq <- DESeq(hemo_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(hemo_dds_deseq)

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
hemo_dds_deseq_res_Pmar <- results(hemo_dds_deseq, alpha=0.05, name = "condition_Pmar_vs_control"   )
hemo_dds_deseq_res_Pmar_GDC <- results(hemo_dds_deseq, alpha=0.05, name= "condition_Pmar_GDC_vs_control" )
hemo_dds_deseq_res_Pmar_ZVAD <- results(hemo_dds_deseq, alpha=0.05, name= "condition_Pmar_ZVAD_vs_control")

head(hemo_dds_deseq_res_Pmar) #  condition Pmar vs control
head(hemo_dds_deseq_res_Pmar_GDC) # condition Pmar GDC vs control
head(hemo_dds_deseq_res_Pmar_ZVAD) # condition Pmar ZVAD vs control 

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

hemo_dds_deseq_res_Pmar_LFC <- lfcShrink(hemo_dds_deseq, coef="condition_Pmar_vs_control" , type= "apeglm", res=hemo_dds_deseq_res_Pmar)
# Review results object summary
summary(hemo_dds_deseq_res_Pmar_LFC)
#out of 45323 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 219, 0.48%
#LFC < 0 (down)     : 299, 0.66%
#outliers [1]       : 3499, 7.7%
#low counts [2]     : 3515, 7.8%
#(mean count < 2)

hemo_dds_deseq_res_Pmar_GDC_LFC <- lfcShrink(hemo_dds_deseq, coef="condition_Pmar_GDC_vs_control" , type= "apeglm", res=hemo_dds_deseq_res_Pmar_GDC)
summary(hemo_dds_deseq_res_Pmar_GDC_LFC)
#out of 45323 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 808, 1.8%
#LFC < 0 (down)     : 769, 1.7%
#outliers [1]       : 3499, 7.7%
#low counts [2]     : 10077, 22%
#(mean count < 6)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

hemo_dds_deseq_res_Pmar_ZVAD_LFC <- lfcShrink(hemo_dds_deseq, coef="condition_Pmar_ZVAD_vs_control", type= "apeglm", res=hemo_dds_deseq_res_Pmar_ZVAD)
summary(hemo_dds_deseq_res_Pmar_ZVAD_LFC)
# out of 45323 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 440, 0.97%
# LFC < 0 (down)     : 382, 0.84%
# outliers [1]       : 3499, 7.7%
# low counts [2]     : 11569, 26%
# (mean count < 7)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

### REPEAT WITH P.MAR AS THE COMPARISON GROUP 

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

hemo_counts_Pmar <- hemo_counts[,c(1:3,5:7,9:11)]
hemo_coldata_Pmar <- hemo_coldata %>% filter(condition != "control")
hemo_coldata_Pmar$condition <- droplevels(hemo_coldata_Pmar$condition)
levels(hemo_coldata_Pmar$condition)

## Creating three here so I can compare the results
hemo_Pmar_dds <- DESeqDataSetFromMatrix(countData = hemo_counts_Pmar,
                                   colData = hemo_coldata_Pmar,
                                   design = ~condition) # only compare by condition

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
hemo_Pmar_dds <- hemo_Pmar_dds[ rowSums(counts(hemo_Pmar_dds)) > 10, ]

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
# set Pmarinus to be the reference level
levels(hemo_coldata_Pmar$condition) 

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
hemo_Pmar_dds_rlog <- rlog(hemo_Pmar_dds, blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(hemo_Pmar_dds_rlog, intgroup=c("condition")) # a bit more clustering by condition now 

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene, 
# and fitting a generalized linear model.
hemo_Pmar_dds_deseq <- DESeq(hemo_Pmar_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(hemo_Pmar_dds_deseq)

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar <- results(hemo_Pmar_dds_deseq, alpha=0.05, name= "condition_Pmar_GDC_vs_Pmar" )
hemo_Pmar_dds_deseq_res_Pmar_ZVAD_Pmar <- results(hemo_Pmar_dds_deseq, alpha=0.05, name= "condition_Pmar_ZVAD_vs_Pmar")

head(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar) # condition Pmar GDC vs Pmar
head(hemo_Pmar_dds_deseq_res_Pmar_ZVAD_Pmar) # condition Pmar ZVAD vs Pmar

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

hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC <- lfcShrink(hemo_Pmar_dds_deseq, coef="condition_Pmar_GDC_vs_Pmar" , type= "apeglm", res=hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar)
summary(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC)
  #out of 42050 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 695, 1.7%
  #LFC < 0 (down)     : 501, 1.2%
  #outliers [1]       : 3222, 7.7%
  #low counts [2]     : 6515, 15%
  #(mean count < 4)
  #[1] see 'cooksCutoff' argument of ?results
  #[2] see 'independentFiltering' argument of ?results

hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC <- lfcShrink(hemo_Pmar_dds_deseq, coef="condition_Pmar_ZVAD_vs_Pmar", type= "apeglm", res=hemo_Pmar_dds_deseq_res_Pmar_ZVAD_Pmar)
summary(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC)
  #out of 42050 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 344, 0.82%
  #LFC < 0 (down)     : 220, 0.52%
  #outliers [1]       : 3222, 7.7%
  #low counts [2]     : 7319, 17%
  #(mean count < 4)
  #[1] see 'cooksCutoff' argument of ?results
  #[2] see 'independentFiltering' argument of ?results

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(hemo_dds_deseq_res_Pmar_LFC, ylim = c(-5, 5))
plotMA(hemo_dds_deseq_res_Pmar_ZVAD_LFC, ylim = c(-5, 5))
plotMA(hemo_dds_deseq_res_Pmar_GDC_LFC, ylim = c(-5, 5))
plotMA(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC, ylim = c(-5, 5))
plotMA(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC, ylim = c(-5, 5))

# these plots all have genes with the same level of expression

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
hemo_dds_deseq_res_Pmar_LFC_sig <- subset(hemo_dds_deseq_res_Pmar_LFC, padj < 0.05)
hemo_dds_deseq_res_Pmar_LFC_sig $ID<- row.names(hemo_dds_deseq_res_Pmar_LFC_sig )
hemo_dds_deseq_res_Pmar_LFC_sig  <- as.data.frame(hemo_dds_deseq_res_Pmar_LFC_sig )
nrow(hemo_dds_deseq_res_Pmar_LFC_sig )  #518

hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig <- subset(hemo_dds_deseq_res_Pmar_ZVAD_LFC, padj < 0.05)
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig$ID <- row.names(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig  )
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig   <- as.data.frame(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig )
nrow(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig  )  #822

hemo_dds_deseq_res_Pmar_GDC_LFC_sig <- subset(hemo_dds_deseq_res_Pmar_GDC_LFC, padj < 0.05)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig $ID <- row.names(hemo_dds_deseq_res_Pmar_GDC_LFC_sig )
hemo_dds_deseq_res_Pmar_GDC_LFC_sig  <- as.data.frame(hemo_dds_deseq_res_Pmar_GDC_LFC_sig)
nrow(hemo_dds_deseq_res_Pmar_GDC_LFC_sig)  #1577

hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig <- subset(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC, padj < 0.05)
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig $ID <- row.names(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig )
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig  <- as.data.frame(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig)
nrow(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig)  #564

hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig <- subset(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC, padj < 0.05)
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig $ID <- row.names(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig )
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig  <- as.data.frame(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig)
nrow(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig)  #1196

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons
# topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 40)
# mat  <- assay(vsd)[ topVarGenes, ]
# mat  <- mat - rowMeans(mat)
# anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
# pheatmap(mat, annotation_col = anno)

topVarGenes_hemo_dds_rlog <-  head(order(rowVars(assay(hemo_dds_rlog )), decreasing = TRUE), 100)
family_hemo_mat <- assay(hemo_dds_rlog)[topVarGenes_hemo_dds_rlog,]
family_hemo_mat <- family_hemo_mat - rowMeans(family_hemo_mat)
family_hemo_anno <- as.data.frame(colData(hemo_dds_rlog)[, c("condition")])
rownames(family_hemo_anno) <- colnames(family_hemo_mat)
family_hemo_heatmap <- pheatmap(family_hemo_mat , annotation_col = family_hemo_anno)
head(family_hemo_mat)

# reorder annotation table to match ordering in heatmap 
family_hemo_heatmap_reorder <-rownames(family_hemo_mat[family_hemo_heatmap$tree_row[["order"]],])
# annotate the row.names
family_hemo_mat_prot <- as.data.frame(family_hemo_heatmap_reorder)
colnames(family_hemo_mat_prot)[1] <- "ID"
family_hemo_mat_prot_annot <- left_join(family_hemo_mat_prot, select(C_vir_rtracklayer_transcripts, ID, product, gene), by = "ID")

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge

# Pmar vs control
hemo_dds_deseq_res_Pmar_LFC_sig_APOP <- merge(hemo_dds_deseq_res_Pmar_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
hemo_dds_deseq_res_Pmar_LFC_sig_APOP_arranged <- arrange(hemo_dds_deseq_res_Pmar_LFC_sig_APOP , -log2FoldChange) 
nrow(hemo_dds_deseq_res_Pmar_LFC_sig_APOP) #15

# pmar ZVAD vs control
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP <- merge(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig , C_vir_rtracklayer_apop_product_final, by = "ID")
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_arranged <- arrange(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP  , -log2FoldChange) 
nrow(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP) # 21

# pmar GDC vs control
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP <- merge(hemo_dds_deseq_res_Pmar_GDC_LFC_sig , C_vir_rtracklayer_apop_product_final, by = "ID")
hemo_dds_deseq_res_Pmar_GDC_LFC_LFC_sig_APOP_arranged <- arrange(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP  , -log2FoldChange) 
nrow(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP) # 62

# pmar GDC vs Pmar
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP <- merge(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig , C_vir_rtracklayer_apop_product_final, by = "ID")
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_arranged <- arrange(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP  , -log2FoldChange) 
nrow(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP) # 53

# pmar ZVAD vs Pmar
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP <- merge(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig , C_vir_rtracklayer_apop_product_final, by = "ID")
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_arranged <- arrange(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP  , -log2FoldChange) 
nrow(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP) # 15

# Compare apoptosis genes between group_by_sim groups
hemo_dds_deseq_res_Pmar_LFC_sig_APOP$condition <- "Pmar_vs_control"
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP$condition <- "Pmar_ZVAD_vs_control"
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP$condition <- "Pmar_GDC_vs_control"
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP$condition <- "Pmar_ZVAD_vs_Pmar"
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP$condition <- "Pmar_GDC_vs_Pmar"

# combine data frames 
hemo_upset_all_sig_APOP <- rbind(hemo_dds_deseq_res_Pmar_LFC_sig_APOP[,c("product","condition","log2FoldChange")],
                                 hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP[,c("product","condition","log2FoldChange")],
                                 hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP[,c("product","condition","log2FoldChange")],
                                 hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP[,c("product","condition","log2FoldChange")],
                                 hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP[,c("product","condition","log2FoldChange")])

# Convert into wide format using reshape
hemo_upset_all_sig_APOP_tally <- hemo_upset_all_sig_APOP %>% group_by(product) %>% tally() 
hemo_upset_all_sig_APOP_upset <- hemo_upset_all_sig_APOP %>% group_by(product) %>% mutate(value=1) %>% spread(condition, value, fill =0 )
hemo_upset_all_sig_APOP_upset <- as.matrix(hemo_upset_all_sig_APOP_upset)

# Make plot
hemo_full_LFC_plot <- ggplot(hemo_upset_all_sig_APOP, aes(x=product,y=log2FoldChange, fill=condition )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

# plot of control vs P. mar plots
hemo_full_LFC_plot_control_Pmar <- hemo_upset_all_sig_APOP %>% filter(condition == "Pmar_vs_control" | condition == "Pmar_ZVAD_vs_Pmar" | condition =="Pmar_GDC_vs_Pmar") %>%
  ggplot(., aes(x=product,y=log2FoldChange, fill=condition )) + geom_col(position="dodge") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

# plot of control vs P. mar
hemo_dds_deseq_res_Pmar_LFC_sig_APOP_plot <- ggplot(hemo_dds_deseq_res_Pmar_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") +  ggtitle("P.mar. vs control")  +
  ylab("Log2 Fold Change")

# plot of Pmar ZVAD vs control
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_plot <- ggplot(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("P.mar + ZVAD vs control") +
  ylab("Log2 Fold Change")

# plot of Pmar GDC vs control 
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_plot <- ggplot(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("P.mar + GDC vs control") +
  ylab("Log2 Fold Change")

# plot of Pmar ZVAD vs P.mar
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_plot <- ggplot(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("P.mar + ZVAD vs P.mar") +
  ylab("Log2 Fold Change")

# plot of Pmar GDC vs P.mar 
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_plot <- ggplot(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("P.mar + GDC vs P.mar") +
  ylab("Log2 Fold Change")

ggarrange(hemo_dds_deseq_res_Pmar_LFC_sig_APOP_plot, hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_plot , hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_plot )
ggarrange(hemo_dds_deseq_res_Pmar_LFC_sig_APOP_plot, hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_plot, hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_plot)

## Extract IAP specific manually curated protein lists
hemo_dds_deseq_res_Pmar_LFC_sig_APOP_IAP <- merge(hemo_dds_deseq_res_Pmar_LFC_sig, BIR_XP_gff_CV_uniq_XP_XM, by = "ID")
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_IAP <- merge(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig, BIR_XP_gff_CV_uniq_XP_XM, by = "ID")
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_IAP <- merge(hemo_dds_deseq_res_Pmar_GDC_LFC_sig, BIR_XP_gff_CV_uniq_XP_XM, by = "ID")

hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_IAP <- merge(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig, BIR_XP_gff_CV_uniq_XP_XM, by = "ID")
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_IAP <- merge(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig, BIR_XP_gff_CV_uniq_XP_XM, by = "ID")

# Join with domain info by protein id
hemo_dds_deseq_res_Pmar_LFC_sig_APOP_IAP_dm <- left_join(hemo_dds_deseq_res_Pmar_LFC_sig_APOP_IAP, IAP_domain_structure_no_dup_rm[,-6])
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_IAP_dm <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_IAP, IAP_domain_structure_no_dup_rm[,-6])
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_IAP_dm <- left_join(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_IAP, IAP_domain_structure_no_dup_rm[,-6])

hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_IAP_dm <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_IAP, IAP_domain_structure_no_dup_rm[,-6])
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_IAP_dm <- left_join(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_IAP, IAP_domain_structure_no_dup_rm[,-6])

