
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
library(topGO)

#### LOADING SAVED GENOME, APOPTOSIS NAMES, IAP XP LISTS ####
Apoptosis_frames <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_apoptosis_products.RData")
annotations <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_annotations.RData")

Perkinsus_rtracklayer <- readGFF("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/GCF_000006405.1_JCVI_PMG_1.0_genomic.gff")
Perkinsus_rtracklayer <- as.data.frame(Perkinsus_rtracklayer)

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

### Volcano plots of significant genes
# compute significance 
hemo_dds_deseq_res_Pmar_LFC_sig_volcano <- hemo_dds_deseq_res_Pmar_LFC_sig
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano <- hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig

hemo_dds_deseq_res_Pmar_LFC_sig_volcano$log10 <- -log10(hemo_dds_deseq_res_Pmar_LFC_sig$padj)
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano$log10 <- -log10(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig$padj)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano$log10 <- -log10(hemo_dds_deseq_res_Pmar_GDC_LFC_sig$padj)

hemo_dds_deseq_res_Pmar_LFC_sig_volcano_plot <- ggplot(data = as.data.frame(hemo_dds_deseq_res_Pmar_LFC_sig_volcano),
                                                       aes(x=log2FoldChange, y=log10)) + geom_point() + theme_bw() + 
  labs(y = "-log10(adjusted p-value)")

hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano_plot <- ggplot(data = as.data.frame(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano), 
                                                            aes(x=log2FoldChange, y=log10)) + geom_point() + theme_bw() +
  labs(y = "-log10(adjusted p-value)")
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano_plot <- ggplot(data = as.data.frame(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano), 
                                                           aes(x=log2FoldChange, y=log10)) + geom_point() + theme_bw() +
  labs(y = "-log10(adjusted p-value)")

# annot all 
hemo_dds_deseq_res_Pmar_LFC_sig_volcano_annot <- hemo_dds_deseq_res_Pmar_LFC_sig_volcano %>% mutate(ID = rownames(.)) %>% 
  left_join(., dplyr::select(C_vir_rtracklayer_transcripts, ID, product, gene, transcript_id), by = "ID")

hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano_annot <- hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano %>% mutate(ID = rownames(.)) %>% 
  left_join(., dplyr::select(C_vir_rtracklayer_transcripts, ID, product, gene,transcript_id), by = "ID")

hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano_annot <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano %>% mutate(ID = rownames(.)) %>% 
  left_join(., dplyr::select(C_vir_rtracklayer_transcripts, ID, product, gene,transcript_id), by = "ID")

# annote those greater than 5
hemo_dds_deseq_res_Pmar_LFC_sig_volcano_5_annot <- hemo_dds_deseq_res_Pmar_LFC_sig_volcano_annot %>% filter(log2FoldChange >= 5.0 | log2FoldChange <= -5.0)

hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano_5_annot <- hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano_annot %>% filter(log2FoldChange >= 5.0 | log2FoldChange <= -5.0)

hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano_5_annot <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano_annot %>% filter(log2FoldChange >= 5.0 | log2FoldChange <= -5.0)

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
hemo_dds_deseq_res_Pmar_LFC_sig_APOP_IAP <- merge(hemo_dds_deseq_res_Pmar_LFC_sig, BIR_XP_gff_CV_uniq_XP_XM, by = "ID") %>% distinct(ID, .keep_all = TRUE)
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_IAP <- merge(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig, BIR_XP_gff_CV_uniq_XP_XM, by = "ID") %>% distinct(ID, .keep_all = TRUE)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_IAP <- merge(hemo_dds_deseq_res_Pmar_GDC_LFC_sig, BIR_XP_gff_CV_uniq_XP_XM, by = "ID") %>% distinct(ID, .keep_all = TRUE)
 
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_IAP <- merge(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig, BIR_XP_gff_CV_uniq_XP_XM, by = "ID") %>% distinct(ID, .keep_all = TRUE)
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_IAP <- merge(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig, BIR_XP_gff_CV_uniq_XP_XM, by = "ID") %>% distinct(ID, .keep_all = TRUE)

# Join with domain info by protein id
hemo_dds_deseq_res_Pmar_LFC_sig_APOP_IAP_dm <- left_join(hemo_dds_deseq_res_Pmar_LFC_sig_APOP_IAP, IAP_domain_structure_no_dup_rm[,-6]) %>% distinct(ID, .keep_all = TRUE)
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_IAP_dm <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_IAP, IAP_domain_structure_no_dup_rm[,-6]) %>% distinct(ID, .keep_all = TRUE)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_IAP_dm <- left_join(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_IAP, IAP_domain_structure_no_dup_rm[,-6]) %>% distinct(ID, .keep_all = TRUE)

hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_IAP_dm <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_IAP, IAP_domain_structure_no_dup_rm[,-6]) %>% distinct(ID, .keep_all = TRUE)
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_IAP_dm <- left_join(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_IAP, IAP_domain_structure_no_dup_rm[,-6]) %>% distinct(ID, .keep_all = TRUE)

## Plot with IAP domain information
hemo_dds_deseq_res_Pmar_LFC_sig_APOP_IAP_dm
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_IAP_dm
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_IAP_dm
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_IAP_dm
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_IAP_dm

# plot of control vs P. mar
hemo_dds_deseq_res_Pmar_LFC_sig_APOP_plot_IAP <- hemo_dds_deseq_res_Pmar_LFC_sig_APOP %>% left_join(., hemo_dds_deseq_res_Pmar_LFC_sig_APOP_IAP_dm[,c("ID","Domain_Name")]) %>%
ggplot(., aes(x=product, y = log2FoldChange, fill=Domain_Name)) + 
  geom_col(position="dodge") + 
  #facet_grid(.~Domain_Name) +
  coord_flip() + 
 # scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") +  
  ggtitle("P.mar. vs control")  +
  ylab("Log2 Fold Change")

ggsave(hemo_dds_deseq_res_Pmar_LFC_sig_APOP_plot_IAP,  file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES/hemo_dds_deseq_res_Pmar_LFC_sig_APOP_plot_IAP_PLOT",
       device = "tiff",
       height = 8, width = 8)

# plot of Pmar ZVAD vs control
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_plot_IAP <- hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP %>% left_join(., hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_IAP_dm[,c("ID","Domain_Name")]) %>%
  ggplot(., aes(x=product, y = log2FoldChange, fill=Domain_Name)) + geom_col(position="dodge") +
  #facet_grid(.~Domain_Name) +
  coord_flip() +
  #scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + 
  ggtitle("P.mar + ZVAD vs control") +
  ylab("Log2 Fold Change")

ggsave(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_plot_IAP,  file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES/hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_plot_IAP_PLOT",
       device = "tiff",
       height = 8, width = 8)

# plot of Pmar GDC vs control 
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_plot_IAP <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP %>% left_join(., hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_IAP_dm[,c("ID","Domain_Name")]) %>%
  ggplot(., aes(x=product, y = log2FoldChange, fill=Domain_Name)) + 
  geom_col(position="dodge") +
  #facet_grid(.~Domain_Name) +
  coord_flip() + 
  #scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + 
  ggtitle("P.mar + GDC vs control") +
  ylab("Log2 Fold Change")

ggsave(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_plot_IAP,  file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES/hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_plot_IAP_PLOT",
       device = "tiff",
       height = 8, width = 10)

# plot of Pmar ZVAD vs P.mar
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_plot_IAP <- hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP %>% left_join(., hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_IAP_dm[,c("ID","Domain_Name")]) %>%
  ggplot(., aes(x=product, y = log2FoldChange, fill=Domain_Name)) + 
  geom_col(position="dodge") +
  coord_flip() + 
  #scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + 
  ggtitle("P.mar + ZVAD vs P.mar") +
  ylab("Log2 Fold Change")

ggsave(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_plot_IAP,  file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES/hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_plot_IAP_PLOT",
       device = "tiff",
       height = 8, width = 8)

# plot of Pmar GDC vs P.mar 
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_plot_IAP <- hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP %>% left_join(., hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_IAP_dm[,c("ID","Domain_Name")]) %>%
  ggplot(., aes(x=product, y = log2FoldChange, fill=Domain_Name)) + geom_col(position="dodge") +
  coord_flip() + 
  #scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + 
  ggtitle("P.mar + GDC vs P.mar") +
  ylab("Log2 Fold Change")

ggsave(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_plot_IAP,  file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES/hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_plot_IAP_PLOT",
       device = "tiff",
       height = 8, width = 10)

## Upset plot heatmap of significant apoptosis expression across all treatments 

# combine all dataframes
C_vir_hemo_comb <- rbind(hemo_dds_deseq_res_Pmar_LFC_sig_APOP,
                         hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP,
                         hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP) %>% mutate(transcript_product = paste(product, transcript_id)) %>%
  dplyr::select(transcript_product, condition, log2FoldChange)
C_vir_hemo_comb_spread <- spread(C_vir_hemo_comb, condition, log2FoldChange, fill = 0)
C_vir_hemo_comb_spread <- column_to_rownames(C_vir_hemo_comb_spread , var = "transcript_product") 
C_vir_hemo_comb_spread_mat <- as.matrix(C_vir_hemo_comb_spread)

C_vir_labels =c( "P. mar + GDC vs\nControl", "P. mar vs\nControl", "P. mar + ZVAD vs\nControl")
# create named vector to hold column names
C_vir_column_labels = structure(paste0(C_vir_labels), names = paste0(colnames(C_vir_hemo_comb_spread_mat)))

pdf("./FIGURES/C_vir_hemo_comb_spread_mat.pdf", width = 12, height = 10)
C_vir_heatmap <- ComplexHeatmap::Heatmap(C_vir_hemo_comb_spread_mat, border = TRUE, 
                        #column_title = ComplexHeatmap::gt_render("*C. virginica* Experimental Group"), 
                        column_title_side = "bottom", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                        row_title = "Apoptosis Transcript and Product Name", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                        row_dend_width = unit(2, "cm"),
                        column_labels = C_vir_column_labels[colnames(C_vir_hemo_comb_spread_mat)],
                        # apply split by k-meams clustering to highlight groups
                        row_km = 3, column_km = 2, row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 10),
                        heatmap_legend_param = list(title = "Log2 Fold Change"))
ComplexHeatmap::draw(C_vir_heatmap, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 100), "mm")) #bottom, left, top, right paddings
dev.off()

### Upset plot of significant LFCs > 1

hemo_dds_deseq_res_Pmar_LFC_sig_volcano_annot_1 <- hemo_dds_deseq_res_Pmar_LFC_sig_volcano_annot %>% filter(log2FoldChange >= 1.0 | log2FoldChange <= -1.0) %>%
  filter(!is.na(product)) %>% mutate(condition = "control")

hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano_annot_1 <- hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano_annot %>% 
  filter(log2FoldChange >= 1.0 | log2FoldChange <= -1.0) %>% 
  filter(!is.na(product)) %>% mutate(condition = "Pmar_ZVAD")

hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano_annot_1 <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano_annot %>% 
  filter(log2FoldChange >= 1.0 | log2FoldChange <= -1.0) %>%
  filter(!is.na(product)) %>% mutate(condition = "Pmar_GDC")


C_vir_hemo_1_comb <- rbind(hemo_dds_deseq_res_Pmar_LFC_sig_volcano_annot_1,
                           hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano_annot_1,
                           hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano_annot_1) %>% mutate(transcript_product = paste(product, transcript_id)) %>%
  dplyr::select(transcript_product, condition, log2FoldChange)
C_vir_hemo_1_comb_spread <- spread(C_vir_hemo_1_comb, condition, log2FoldChange, fill = 0)
C_vir_hemo_1_comb_spread <- column_to_rownames(C_vir_hemo_1_comb_spread , var = "transcript_product") 
C_vir_hemo_1_comb_spread_mat <- as.matrix(C_vir_hemo_1_comb_spread)

C_vir_1_labels =c( "P. mar + GDC vs\nControl", "P. mar vs\nControl", "P. mar + ZVAD vs\nControl")
# create named vector to hold column names
C_vir_1_column_labels = structure(paste0(C_vir_1_labels), names = paste0(colnames(C_vir_hemo_1_comb_spread_mat)))

pdf("./FIGURES/C_vir_hemo_1_comb_spread_mat.pdf", width = 12, height = 10)
C_vir_heatmap <- ComplexHeatmap::Heatmap(C_vir_hemo_1_comb_spread_mat, border = TRUE, 
                                         #column_title = ComplexHeatmap::gt_render("*C. virginica* Experimental Group"), 
                                         column_title_side = "bottom", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                         row_title = "Transcript and Product Name LFC > abs(1)", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                         row_dend_width = unit(2, "cm"),
                                         column_labels = C_vir_1_column_labels[colnames(C_vir_hemo_1_comb_spread_mat)],
                                         # apply split by k-meams clustering to highlight groups
                                         row_km = 4, column_km = 2, 
                                         row_names_gp = gpar(fontsize = 2),
                                         column_names_gp = gpar(fontsize = 10),
                                         heatmap_legend_param = list(title = "Log2 Fold Change"))
ComplexHeatmap::draw(C_vir_heatmap, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 100), "mm")) #bottom, left, top, right paddings
dev.off()

pdf("./FIGURES/C_vir_hemo_1_comb_spread_mat_tall.pdf", width = 12, height = 20)
C_vir_heatmap <- ComplexHeatmap::Heatmap(C_vir_hemo_1_comb_spread_mat, border = TRUE, 
                                         #column_title = ComplexHeatmap::gt_render("*C. virginica* Experimental Group"), 
                                         column_title_side = "bottom", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                         row_title = "Transcript and Product Name LFC > abs(1)", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                         row_dend_width = unit(2, "cm"),
                                         column_labels = C_vir_1_column_labels[colnames(C_vir_hemo_1_comb_spread_mat)],
                                         # apply split by k-meams clustering to highlight groups
                                         row_km = 4, column_km = 2, 
                                         row_names_gp = gpar(fontsize = 2),
                                         column_names_gp = gpar(fontsize = 10),
                                         heatmap_legend_param = list(title = "Log2 Fold Change"))
ComplexHeatmap::draw(C_vir_heatmap, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 100), "mm")) #bottom, left, top, right paddings
dev.off()


## Compare lists of apoptosis transcripts between the inhibitor vs control and inhibitor vs parasite lists 
hemo_dds_deseq_res_Pmar_LFC_sig_APOP$DEG_comp <- "vs_control"
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP$DEG_comp <- "vs_control"
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP$DEG_comp <- "vs_control"
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP$DEG_comp <- "vs_pmar"
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP$DEG_comp <- "vs_pmar"

# compare GDC 
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared <- full_join(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                    "transcript_id","gene","product","DEG_comp")],
                                                                     hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                     "transcript_id","gene","product","DEG_comp")], 
                                                                    by = c("ID", "transcript_id","gene","product")) %>% filter(!is.na(DEG_comp.x) & !is.na(DEG_comp.y))

hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_unique_vs_control <- left_join(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                    "transcript_id","gene","product","DEG_comp")],
                                                                    hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                    "transcript_id","gene","product","DEG_comp")], 
                                                                    by = c("ID", "transcript_id","gene","product")) %>%
                                                                    filter(is.na(log2FoldChange.y)) %>% distinct(ID, .keep_all = TRUE)
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_unique_vs_Pmar <- left_join( hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                              "transcript_id","gene","product","DEG_comp")],
                                                                          hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                              "transcript_id","gene","product","DEG_comp")], 
                                                                               by = c("ID", "transcript_id","gene","product")) %>%
                                                                              filter(is.na(log2FoldChange.y)) %>% distinct(ID, .keep_all = TRUE)
# repeat for ZVAD
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_shared <- full_join(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                    "transcript_id","gene","product","DEG_comp")],
                                                                    hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                    "transcript_id","gene","product","DEG_comp")], 
                                                                     by = c("ID", "transcript_id","gene","product")) %>% filter(!is.na(DEG_comp.x) & !is.na(DEG_comp.y))

hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_unique_vs_control <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                               "transcript_id","gene","product","DEG_comp")],
                                                                               hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                                "transcript_id","gene","product","DEG_comp")], 
                                                                               by = c("ID", "transcript_id","gene","product")) %>%
                                                                              filter(is.na(log2FoldChange.y)) %>% distinct(ID, .keep_all = TRUE)

hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP_unique_vs_Pmar <- left_join( hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                                       "transcript_id","gene","product","DEG_comp")],
                                                                                       hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP[,c("ID","log2FoldChange",
                                                                                      "transcript_id","gene","product","DEG_comp")], 
                                                                                       by = c("ID", "transcript_id","gene","product")) %>%
                                                                                      filter(is.na(log2FoldChange.y)) %>% distinct(ID, .keep_all = TRUE)



# Plot these shared or unique DEGs
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared_vs_control <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared[,1:6] %>% 
  dplyr::rename(log2FoldChange = log2FoldChange.x, DEG_comp = DEG_comp.x)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared_vs_pmar <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared[,c(1,3:5,7,8)] %>% 
  dplyr::rename(log2FoldChange = log2FoldChange.y, DEG_comp = DEG_comp.y)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared_join <- rbind(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared_vs_control,hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared_vs_pmar)

hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared_plot <-  hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared_join %>%
  left_join(., hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_APOP_IAP_dm[,c("ID","Domain_Name")]) %>%
  ggplot(., aes(x=product, y = log2FoldChange, fill=Domain_Name)) + 
  geom_col(position="dodge") +
  facet_grid(.~DEG_comp) +
  coord_flip() + 
  #scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + 
  ggtitle("P.mar + GDC vs control or P.mar + GDC vs P.mar") +
  ylab("Log2 Fold Change")

ggsave(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared_plot,  file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES/hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_shared_plot_PLOT",
       device = "tiff",
       height = 8, width = 12)

#GDC unique vs control not including the parasite 
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_unique_vs_control_no_pm <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_unique_vs_control[!(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_unique_vs_control$ID %in% hemo_dds_deseq_res_Pmar_LFC_sig_APOP$ID),]

hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_unique_vs_control_no_pm_plot <-  hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_unique_vs_control_no_pm %>%
  left_join(., hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_IAP_dm[,c("ID","Domain_Name")]) %>%
  ggplot(., aes(x=product, y = log2FoldChange.x, fill=Domain_Name)) + 
  geom_col(position="dodge") +
  coord_flip() + 
  #scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + 
  ggtitle("P.mar + GDC vs control without parasite alone shared") +
  ylab("Log2 Fold Change")

ggsave(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_unique_vs_control_no_pm_plot ,  file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES/hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_unique_vs_control_no_pm_PLOT",
       device = "tiff",
       height = 8, width = 12)

### Assess non-apoptotic differentially expressed genes in each and over 1 
hemo_dds_deseq_res_Pmar_LFC_sig_ID  <- merge(hemo_dds_deseq_res_Pmar_LFC_sig , C_vir_rtracklayer_transcripts, by = "ID") 
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID <- merge(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig, C_vir_rtracklayer_transcripts, by = "ID")
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID <- merge(hemo_dds_deseq_res_Pmar_GDC_LFC_sig, C_vir_rtracklayer_transcripts, by = "ID")
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID <- merge(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig, C_vir_rtracklayer_transcripts, by = "ID")
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID <- merge(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig, C_vir_rtracklayer_transcripts, by = "ID")

# filter to be those greater than 1 or less than -1
hemo_dds_deseq_res_Pmar_LFC_sig_ID_1  <- merge(hemo_dds_deseq_res_Pmar_LFC_sig , C_vir_rtracklayer_transcripts, by = "ID") %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1)
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID_1 <- merge(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig, C_vir_rtracklayer_transcripts, by = "ID") %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID_1 <- merge(hemo_dds_deseq_res_Pmar_GDC_LFC_sig, C_vir_rtracklayer_transcripts, by = "ID") %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1)
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID_1 <- merge(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig, C_vir_rtracklayer_transcripts, by = "ID") %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1)
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID_1 <- merge(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig, C_vir_rtracklayer_transcripts, by = "ID") %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1)

# explore genes changing more than 1 fold
ggplot(hemo_dds_deseq_res_Pmar_LFC_sig_ID_1, aes(x = product, y = log2FoldChange)) + geom_col(position = "dodge") + coord_flip()

ggplot(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID_1, aes(x = product, y = log2FoldChange)) + geom_col(position = "dodge") + coord_flip()

ggplot(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID_1, aes(x = product, y = log2FoldChange)) + geom_col(position = "dodge") + coord_flip()

ggplot(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID_1, aes(x = product, y = log2FoldChange)) + geom_col(position = "dodge") + coord_flip()

ggplot(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID_1, aes(x = product, y = log2FoldChange)) + geom_col(position = "dodge") + coord_flip()

### Run Interproscan to get GO terms for all those significant transcript IDs

# concatenate all transcript ID lists so I can create a lookup list in terminal
hemo_dds_deseq_sig_XP <- rbind(hemo_dds_deseq_res_Pmar_LFC_sig_ID,
                                hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID,
                                hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID,
                                hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID,
                                hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID) %>% distinct(ID) %>%
                                dplyr::rename(Parent = ID)
hemo_dds_deseq_sig_XP <- as.data.frame(hemo_dds_deseq_sig_XP)
# join protein_id by searching for the transcript ID in the parent column
C_vir_rtracklayer_XP <- C_vir_rtracklayer %>% filter(!is.na(protein_id))
C_vir_rtracklayer_XP$Parent <- as.character(C_vir_rtracklayer_XP$Parent)
hemo_dds_deseq_sig_XP_df <- left_join(hemo_dds_deseq_sig_XP, unique(C_vir_rtracklayer_XP[,c("protein_id","Parent")])) %>% filter(!is.na(protein_id))

# NA's due to XR proteins from non-coding RNA

# export protein IDs to lookup in bluewaves
write.table(hemo_dds_deseq_sig_XP_df$protein_id, file = "hemo_dds_deseq_sig_XP_df_lookup.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

## Get gene universe of all the expressed genes in each experiment! Need to have the GO in this full list in order to do enrichment 

hemo_counts_universe <- as.data.frame(row.names(hemo_dds)) # use this data frame with the extras already removed
colnames(hemo_counts_universe)[1] <- "Parent" 
hemo_counts_universe_XP <- left_join(hemo_counts_universe, unique(C_vir_rtracklayer_XP[,c("protein_id","Parent")]))  %>% filter(!is.na(protein_id))
nrow(hemo_counts_universe_XP) # 41475

# export into files with 5000 in each 
write.table(hemo_counts_universe_XP[1:5000, "protein_id"], file = "hemo_counts_universe_XP_5000.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_counts_universe_XP[5001:10000, "protein_id"], file = "hemo_counts_universe_XP_10000.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_counts_universe_XP[10001:15000, "protein_id"], file = "hemo_counts_universe_XP_15000.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_counts_universe_XP[15001:20000, "protein_id"], file = "hemo_counts_universe_XP_20000.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_counts_universe_XP[20001:25000, "protein_id"], file = "hemo_counts_universe_XP_25000.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_counts_universe_XP[25001:30000, "protein_id"], file = "hemo_counts_universe_XP_30000.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_counts_universe_XP[30001:35000, "protein_id"], file = "hemo_counts_universe_XP_35000.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_counts_universe_XP[35001:40000, "protein_id"], file = "hemo_counts_universe_XP_40000.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_counts_universe_XP[40001:41475, "protein_id"], file = "hemo_counts_universe_XP_41475.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)


#### HEMOCYTE GO ANALYSIS ####

GO_sig_terms <- rtracklayer::readGFF("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/hemo_dds_deseq_sig_XP_df_lookup_seq.fa_1.gff3")
GO_sig_terms <- as.data.frame(GO_sig_terms)

GO_sig_terms_found <- GO_sig_terms %>% filter(Ontology_term != "character(0)") %>% dplyr::rename(protein_id = seqid) 

# join with GO terms for each list
hemo_dds_deseq_res_Pmar_LFC_sig_ID <- hemo_dds_deseq_res_Pmar_LFC_sig_ID %>% dplyr::rename(Parent_gene = Parent, Parent = ID) %>% dplyr::select(-protein_id)
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID <- hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID %>% dplyr::rename(Parent_gene = Parent, Parent = ID) %>% dplyr::select(-protein_id)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID %>% dplyr::rename(Parent_gene = Parent, Parent = ID) %>% dplyr::select(-protein_id)
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID <- hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID %>% dplyr::rename(Parent_gene = Parent, Parent = ID) %>% dplyr::select(-protein_id)
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID <- hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID %>% dplyr::rename(Parent_gene = Parent, Parent = ID) %>% dplyr::select(-protein_id)

hemo_dds_deseq_res_Pmar_LFC_sig_ID$Parent <-  as.character(hemo_dds_deseq_res_Pmar_LFC_sig_ID$Parent)
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID$Parent <-  as.character(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID$Parent)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID$Parent <-  as.character(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID$Parent)
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID$Parent <-  as.character(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID$Parent)
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID$Parent <-  as.character(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID$Parent)

hemo_dds_deseq_res_Pmar_LFC_sig_ID_GO <- left_join(hemo_dds_deseq_res_Pmar_LFC_sig_ID, unique(C_vir_rtracklayer_XP[,c("protein_id","Parent")])) 
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID_GO <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID, unique(C_vir_rtracklayer_XP[,c("protein_id","Parent")])) 
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID_GO <- left_join(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID, unique(C_vir_rtracklayer_XP[,c("protein_id","Parent")])) 
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID_GO <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID, unique(C_vir_rtracklayer_XP[,c("protein_id","Parent")])) 
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID_GO <- left_join(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID, unique(C_vir_rtracklayer_XP[,c("protein_id","Parent")])) 

hemo_dds_deseq_res_Pmar_LFC_sig_ID_GO  <- left_join(hemo_dds_deseq_res_Pmar_LFC_sig_ID_GO, GO_sig_terms_found[,c("protein_id","Ontology_term")])
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID_GO <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID_GO, GO_sig_terms_found[,c("protein_id","Ontology_term")])
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID_GO <- left_join(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID_GO, GO_sig_terms_found[,c("protein_id","Ontology_term")])
hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID_GO <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_Pmar_LFC_sig_ID_GO, GO_sig_terms_found[,c("protein_id","Ontology_term")])
hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID_GO <- left_join(hemo_Pmar_dds_deseq_res_Pmar_GDC_Pmar_LFC_sig_ID_GO, GO_sig_terms_found[,c("protein_id","Ontology_term")])

## Perform GO enrichment 

#### PERKINSUS TRANSCRIPTOME ANALYSIS ####

## LOAD DATA
perk_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/Dermo_2021_transcript_count_matrix.csv", header=TRUE,
                        row.names = "transcript_id")
colnames(perk_counts)[1:12] <-  c("1_Dermo_GDC_R1_001",  "1_Dermo_ZVAD_R1_001", "1_Dermo_R1_001"     , "1_control_R1_001" ,   "2_Dermo_GDC_R1_001" , "2_Dermo_ZVAD_R1_001" ,"2_Dermo_R1_001",     
                                  "2_control_R1_001"  ,  "3_Dermo_GDC_R1_001" , "3_Dermo_ZVAD_R1_001", "3_Dermo_R1_001"   ,   "3_control_R1_001"   )
head(perk_counts)
colnames(perk_counts)

# remove MSTRG novel transcript lines (can assess these later)
perk_counts <- perk_counts[!grepl("MSTRG", row.names(perk_counts)),]

# Cute the "rna-" from the beginning of rownames
remove_rna = function(x){
  return(gsub("rna-","",x))
}
row.names(perk_counts) <- remove_rna(row.names(perk_counts))
head(perk_counts)

#Load in sample metadata
perk_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/Hemo_Pmar_coldata.csv", row.names = 1 )
View(perk_coldata)  
nrow(perk_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
perk_coldata <- perk_coldata[colnames(perk_counts),]  

all(rownames(perk_coldata) %in% colnames(perk_counts))  #Should return TRUE
# returns TRUE
all(colnames(perk_counts) %in% rownames(perk_coldata))  
# returns TRUE
all(rownames(perk_coldata) == colnames(perk_counts))    # should return TRUE
# returns TRUE

# remove the control samples (since no Perkinsus was added here!)
perk_coldata <- perk_coldata %>% filter(condition != "control")
perk_coldata$condition <- droplevels(perk_coldata$condition)
perk_counts <- perk_counts[,c(1:3,5:7,9:11)]

### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
perk_counts_matrix <- as.matrix(perk_counts)
perk_rlog_counts <- rlog(perk_counts_matrix, blind =TRUE)

# run PCA
pcperk <- prcomp(t(perk_rlog_counts))
# plot PCA
autoplot(pcperk)

# Lets add colour to look at the clustering for Status
autoplot(pcperk,
         data = perk_coldata, 
         colour="condition", 
         size=5) 
# little clustering by condition

autoplot(pcperk,
         data = as.data.frame(perk_coldata), 
         colour="pool", # clustering by pool 
         size=5) 

# clustering mostly by pool, but PC1 and PC2 only explain about 30% of the total sample variation

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Creating three here so I can compare the results
perk_dds <- DESeqDataSetFromMatrix(countData = perk_counts,
                                   colData = perk_coldata,
                                   design = ~condition) # only compare by condition

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
perk_dds <- perk_dds[ rowSums(counts(perk_dds)) > 10, ]

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(perk_coldata$condition)  # Pmar is currently listed first, so this looks good 

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
perk_dds_rlog <- rlog(perk_dds, blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(perk_dds_rlog, intgroup=c("condition")) # a bit more clustering by condition now 

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene, 
# and fitting a generalized linear model.
perk_dds_deseq <- DESeq(perk_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(perk_dds_deseq)

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
perk_dds_deseq_res_Pmar_GDC <- results(perk_dds_deseq, alpha=0.05, name= "condition_Pmar_GDC_vs_Pmar" )
perk_dds_deseq_res_Pmar_ZVAD <- results(perk_dds_deseq, alpha=0.05, name= "condition_Pmar_ZVAD_vs_Pmar")

head(perk_dds_deseq_res_Pmar_GDC) # condition Pmar GDC vs control
head(perk_dds_deseq_res_Pmar_ZVAD) # condition Pmar ZVAD vs control 

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

perk_dds_deseq_res_Pmar_GDC_LFC <- lfcShrink(perk_dds_deseq, coef="condition_Pmar_GDC_vs_Pmar" , type= "apeglm", res=perk_dds_deseq_res_Pmar_GDC)
summary(perk_dds_deseq_res_Pmar_GDC_LFC)
#out of 16842 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 20, 0.12%
#LFC < 0 (down)     : 19, 0.11%
#outliers [1]       : 252, 1.5%
#low counts [2]     : 5537, 33%
#(mean count < 8)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

perk_dds_deseq_res_Pmar_ZVAD_LFC <- lfcShrink(perk_dds_deseq, coef="condition_Pmar_ZVAD_vs_Pmar", type= "apeglm", res=perk_dds_deseq_res_Pmar_ZVAD)
summary(perk_dds_deseq_res_Pmar_ZVAD_LFC)
#out of 16842 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 19, 0.11%
#LFC < 0 (down)     : 27, 0.16%
#outliers [1]       : 252, 1.5%
#low counts [2]     : 2939, 17%
#(mean count < 3)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering

perk_dds_deseq_res_Pmar_ZVAD_LFC_sig <- subset(perk_dds_deseq_res_Pmar_ZVAD_LFC, padj < 0.05)
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig$transcript_id <- row.names(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig  )
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig   <- as.data.frame(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig )
nrow(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig  )  #46

perk_dds_deseq_res_Pmar_GDC_LFC_sig <- subset(perk_dds_deseq_res_Pmar_GDC_LFC, padj < 0.05)
perk_dds_deseq_res_Pmar_GDC_LFC_sig $transcript_id <- row.names(perk_dds_deseq_res_Pmar_GDC_LFC_sig )
perk_dds_deseq_res_Pmar_GDC_LFC_sig  <- as.data.frame(perk_dds_deseq_res_Pmar_GDC_LFC_sig)
nrow(perk_dds_deseq_res_Pmar_GDC_LFC_sig)  #39

# Annotate these genes
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot <- merge(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig, unique(dplyr::select(Perkinsus_rtracklayer, transcript_id, product))) %>% 
  mutate(product = case_when(product == "conserved hypothetical protein" | product == "hypothetical protein" ~ paste(product, transcript_id, sep = "-"),
                             product != "conserved hypothetical protein" | product != "hypothetical protein" ~ product))

perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot <- merge(perk_dds_deseq_res_Pmar_GDC_LFC_sig, unique(dplyr::select(Perkinsus_rtracklayer,transcript_id, product))) %>% 
  mutate(product = case_when(product == "conserved hypothetical protein" | product == "hypothetical protein" ~ paste(product, transcript_id, sep = "-"),
                             product != "conserved hypothetical protein" | product != "hypothetical protein" ~ product))

perk_dds_deseq_res_Pmar_GDC_ZVAD_LFC_sig_annot_comb <- intersect(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot[,c("transcript_id","product")], perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot[,c("transcript_id","product")])

# overall, little change from control in terms of Perkinsus expression. Which means that these inhibitors didn't really affect the parasite
# both treatments are significantly differentially expressing the same genes. 
    # This could indicate perhaps PCR bias in the Perkinsus genes that were sequenced, or that treatments caused little change
    # (which is good considering we only wanted the inhibitors to affect the hemocytes)

# plot LFC 
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_plot <- ggplot(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot, aes(x= product, y = log2FoldChange)) + 
  geom_col() + coord_flip() + ggtitle("P. mar. ZVAD\n vs P. mar control")

perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_plot <- ggplot(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot, aes(x= product, y = log2FoldChange)) +
  geom_col() + coord_flip() + ggtitle("P. mar. GDC\n vs P. mar control")

combined_Pmar <- ggarrange(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_plot, perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_plot)

ggsave(combined_Pmar, device = "tiff",
       width = 12, height = 7,
       file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/combined_Pmar_plot")

## Perkinsus expression upset plot 
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot$condition <- "P_mar_ZVAD_vs_Pmar"
perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot$condition <- "P_mar_GDC_vs_Pmar"

# combine all dataframes
Perk_comb <- rbind(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot,
                   perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot) %>%
  dplyr::select(product, condition, log2FoldChange)
Perk_comb_spread <- spread(Perk_comb, condition, log2FoldChange, fill = 0)
Perk_comb_spread <- column_to_rownames(Perk_comb_spread , var = "product") 
Perk_comb_spread_mat <- as.matrix(Perk_comb_spread)

Perk_labels =c( "P. mar + GDC vs\nP. mar Control", "P. mar + ZVAD vs\nP. mar Control")
# create named vector to hold column names
Perk_column_labels = structure(paste0(Perk_labels), names = paste0(colnames(Perk_comb_spread_mat)))

pdf("./FIGURES/Perk_spread_mat.pdf", width = 8, height = 8)
Perk_heatmap <- ComplexHeatmap::Heatmap(Perk_comb_spread_mat, border = TRUE, 
                        column_title_side = "bottom", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                        row_title = "Transcript and Product Name", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                        row_dend_width = unit(2, "cm"),
                        column_labels = Perk_column_labels[colnames(Perk_comb_spread_mat)],
                        # apply split by k-meams clustering to highlight groups
                        #row_km = 3, column_km = 2, 
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        heatmap_legend_param = list(title = "Log2 Fold Change"))
ComplexHeatmap::draw(Perk_heatmap, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 40), "mm")) #bottom, left, top, right paddings

dev.off()
