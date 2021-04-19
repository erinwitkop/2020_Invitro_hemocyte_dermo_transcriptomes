
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
library(GOSim)
library(GO.db)
library(PCAtools)


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
pdf("./FIGURES/pchemo_condition.pdf", height = 5, width = 5)
autoplot(pchemo,
         data = hemo_coldata, 
         colour="condition", 
         size=5) + ggtitle("Hemocyte expression")
dev.off()
# little clustering by condition

pdf("./FIGURES/pchemo_pool.pdf", height = 5, width = 5)
autoplot(pchemo,
         data = hemo_coldata, 
         colour="pool", 
         size=5) + ggtitle("Hemocyte expression")
dev.off()
# PC1 and PC2 don't explain much of the variation, which means that the Pool effect is not too large 

## Plot total reads in each sample
hemo_counts_total <- colSums(hemo_counts)
hemo_counts_total <- as.data.frame(hemo_counts_total)
hemo_counts_total <- hemo_counts_total %>% mutate(sample_name = rownames(.)) %>% 
  mutate(condition= case_when(
    grepl("GDC", sample_name) ~"GDC",
    grepl("ZVAD", sample_name) ~"ZVAD",
    grepl("control",sample_name)~"control",
    TRUE~"Dermo"))
hemo_counts_total$sample_name <- factor(hemo_counts_total$sample_name, 
                                        levels = c("1_control_R1_001","2_control_R1_001","3_control_R1_001",
                                             "1_Dermo_R1_001", "2_Dermo_R1_001","3_Dermo_R1_001",
                                                   "1_Dermo_GDC_R1_001",  "2_Dermo_GDC_R1_001","3_Dermo_GDC_R1_001" ,
                                                   "1_Dermo_ZVAD_R1_001","2_Dermo_ZVAD_R1_001",  "3_Dermo_ZVAD_R1_001" ))

hemo_counts_total_plot <- ggplot(hemo_counts_total, aes(x=sample_name, y = hemo_counts_total, fill = condition)) + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(hemo_counts_total_plot, file = "./FIGURES/hemo_counts_total_plot.tiff", device = "tiff",
       height=5, width = 5)


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

head(assay(hemo_dds_rlog),3)

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

# Annotate all the significant proteins in each 

hemo_dds_deseq_res_Pmar_LFC_sig_all_annot <- left_join(hemo_dds_deseq_res_Pmar_LFC_sig,dplyr::select(C_vir_rtracklayer_transcripts, ID, product, gene), by = "ID")

hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_all_annot <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig,dplyr::select(C_vir_rtracklayer_transcripts, ID, product, gene), by = "ID")

hemo_dds_deseq_res_Pmar_GDC_LFC_sig_all_annot <- left_join(hemo_dds_deseq_res_Pmar_GDC_LFC_sig,dplyr::select(C_vir_rtracklayer_transcripts, ID, product, gene), by = "ID")

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
family_hemo_mat_prot_annot <- left_join(family_hemo_mat_prot, dplyr::select(C_vir_rtracklayer_transcripts, ID, product, gene), by = "ID")

### Volcano plots of significant genes
# compute significance 
hemo_dds_deseq_res_Pmar_LFC_sig_volcano <- hemo_dds_deseq_res_Pmar_LFC_sig
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano <- hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig

hemo_dds_deseq_res_Pmar_LFC_sig_volcano$log10 <- -log10(hemo_dds_deseq_res_Pmar_LFC_sig$padj)
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano$log10 <- -log10(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig$padj)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano$log10 <- -log10(hemo_dds_deseq_res_Pmar_GDC_LFC_sig$padj)

# plot the volcano plots
hemo_dds_deseq_res_Pmar_LFC_sig_volcano_plot <- ggplot(data = as.data.frame(hemo_dds_deseq_res_Pmar_LFC_sig_volcano),
                                                       aes(x=log2FoldChange, y=log10)) + geom_point() + theme_bw() + 
  labs(y = "-log10(adjusted p-value)", title = "P.mar vs. control hemocytes")

hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano_plot <- ggplot(data = as.data.frame(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano), 
                                                            aes(x=log2FoldChange, y=log10)) + geom_point() + theme_bw() +
  labs(y = "-log10(adjusted p-value)", title = "P.mar + ZVAD vs. control hemocytes")

hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano_plot <- ggplot(data = as.data.frame(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano), 
                                                           aes(x=log2FoldChange, y=log10)) + geom_point() + theme_bw() +
  labs(y = "-log10(adjusted p-value)", title = "P.mar + GDC vs. control hemocytes")

hemo_volcano <- ggarrange(hemo_dds_deseq_res_Pmar_LFC_sig_volcano_plot, 
          hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_volcano_plot,
          hemo_dds_deseq_res_Pmar_GDC_LFC_sig_volcano_plot)

ggsave(hemo_volcano, file = "./FIGURES/hemo_volcano_plot", device = "tiff", height = 10, width = 10)


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

### Upset plot heatmap of significant apoptosis expression across all treatments  ####

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

## Plot apoptosis transcformed counts across each sample ####
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons

C_vir_hemo_comb_ID <-  rbind(hemo_dds_deseq_res_Pmar_LFC_sig_APOP,
                             hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP,
                             hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP) %>% dplyr::pull(., ID)
class(C_vir_hemo_comb_ID)

# find rownames matching these IDs
apop_hemo_mat <- assay(hemo_dds_rlog)[C_vir_hemo_comb_ID,]
apop_hemo_mat <- as.data.frame(apop_hemo_mat) %>% mutate(ID = rownames(.)) %>% 
  left_join(., C_vir_rtracklayer_apop_product_final[,c("ID","product","transcript_id")], by = "ID") %>% filter(!is.na(product)) %>%
  mutate(transcript_product = paste(product, transcript_id, sep = "-")) %>% dplyr::select(-ID,-product, -transcript_id)
rownames(apop_hemo_mat) <- apop_hemo_mat$transcript_product
apop_hemo_mat <- apop_hemo_mat[,-13]
apop_hemo_mat <- as.matrix(apop_hemo_mat)  

apop_hemo_anno <- as.data.frame(colData(hemo_dds_rlog)[, c("condition")])
rownames(apop_hemo_anno) <- colnames(apop_hemo_mat)
colnames(apop_hemo_anno)[1] <- "Condition"

pdf("./FIGURES/C_vir_apop_hemo_mat1.pdf", width = 12, height = 12)
pheatmap(apop_hemo_mat , annotation_col = apop_hemo_anno)
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

### Run Interproscan to get GO terms for all those significant transcript IDs ####

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

# topGO tutorial: https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

GO_sig_terms <- rtracklayer::readGFF("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/hemo_dds_deseq_sig_XP_df_lookup_seq.fa_1.gff3")
GO_sig_terms <- as.data.frame(GO_sig_terms)

GO_sig_terms_found <- GO_sig_terms %>% filter(Ontology_term != "character(0)") %>% dplyr::rename(protein_id = seqid) 

# upload GO universe, every Interproscan line with GO term
GO_universe <-  rtracklayer::readGFF("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/hemo_universe_all.GO.gff3")
GO_universe <- as.data.frame(GO_universe)
nrow(GO_universe) # 164675 lines 

# get list of gene identifiers and gene significance values, focusing just on the GDC vs control and ZVAD vs control 
hemo_dds_deseq_res_Pmar_LFC_sig_gene_list <- hemo_dds_deseq_res_Pmar_LFC_sig_ID %>% dplyr::select(ID, padj) %>% dplyr::rename(transcript_id = ID)
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_gene_list <- hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_ID %>% dplyr::select(ID, padj) %>% dplyr::rename(transcript_id = ID)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_gene_list <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_ID %>% dplyr::select(ID, padj) %>% dplyr::rename(transcript_id = ID)

## get gene universe list of GO terms and GO identifiers 
# first join with the rna ID
GO_universe_rna <- GO_universe %>% dplyr::rename(protein_id = seqid) %>% left_join(., unique(C_vir_rtracklayer_XP[,c("protein_id","Parent")])) %>% dplyr::rename(transcript_id = Parent)

# find transcripts with GO terms and correctly format
GO_universe_rna_found <- GO_universe_rna %>% filter(Ontology_term != "character(0)") %>% distinct(Ontology_term, protein_id, .keep_all = TRUE)
nrow(GO_universe_rna_found)
class(GO_universe_rna_found$Ontology_term) # AsIs

GO_universe_rna_found$Ontology_term <- unlist(as.character(GO_universe_rna_found$Ontology_term))
class(GO_universe_rna_found$Ontology_term)
GO_universe_rna_found$Ontology_term <- gsub('\\\"', "", GO_universe_rna_found$Ontology_term, fixed = TRUE)
GO_universe_rna_found$Ontology_term <- gsub('\"', "", GO_universe_rna_found$Ontology_term, fixed = TRUE)
GO_universe_rna_found$Ontology_term <- gsub('c(', "", GO_universe_rna_found$Ontology_term, fixed = TRUE)
GO_universe_rna_found$Ontology_term <- gsub(')', "", GO_universe_rna_found$Ontology_term, fixed = TRUE)

View(GO_universe_rna_found)
# save this
save(GO_universe_rna_found, file = "GO_universe_rna_found.RData")

# format dataframe for use in topGO for custom annotations
GO_universe_rna_found_geneID2GO <- GO_universe_rna_found %>% dplyr::select(transcript_id, Ontology_term)
# export as text file to get tab separate file and then re-load using read mappings function
write.table(GO_universe_rna_found_geneID2GO, file = "GO_universe_rna_found_geneID2GO.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
GO_universe_rna_found_geneID2GO_mapping <-  readMappings(file= "GO_universe_rna_found_geneID2GO.txt")
# save as Rdata for use in WGCNA workspace
save(GO_universe_rna_found_geneID2GO_mapping, file = "GO_universe_rna_found_geneID2GO_mapping.RData")

# join GO terms for reference later and only keep the list of sig genes in each group that have a GO term
hemo_dds_deseq_res_Pmar_LFC_sig_gene_list_GO <- left_join(hemo_dds_deseq_res_Pmar_LFC_sig_gene_list, GO_universe_rna_found) %>% filter(!is.na(Ontology_term)) %>% distinct(transcript_id, .keep_all = TRUE)
nrow(hemo_dds_deseq_res_Pmar_LFC_sig_gene_list_GO) # 328
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_gene_list_GO <- left_join(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_gene_list, GO_universe_rna_found) %>% filter(!is.na(Ontology_term)) %>% distinct(transcript_id, .keep_all = TRUE)
nrow(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_gene_list_GO) # 537
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_gene_list_GO <- left_join(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_gene_list, GO_universe_rna_found) %>% filter(!is.na(Ontology_term)) %>% distinct(transcript_id, .keep_all = TRUE)
nrow(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_gene_list_GO) # 925

# get those sig LFCs with a GO term and put in correct format for topGO
geneNames <- names(GO_universe_rna_found_geneID2GO_mapping)
head(geneNames)

hemo_dds_deseq_res_Pmar_LFC_sig_gene_list_GO_only <- as.factor(hemo_dds_deseq_res_Pmar_LFC_sig_gene_list_GO$transcript_id)
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_gene_list_GO_only <-  as.factor(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_gene_list_GO$transcript_id)
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_gene_list_GO_only <-  as.factor(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_gene_list_GO$transcript_id)

hemo_dds_deseq_res_Pmar_LFC_sig_gene_list_factor <- factor(as.integer(geneNames %in% hemo_dds_deseq_res_Pmar_LFC_sig_gene_list_GO_only))
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_gene_list_factor <- factor(as.integer(geneNames %in% hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_gene_list_GO_only))
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_gene_list_factor <- factor(as.integer(geneNames %in% hemo_dds_deseq_res_Pmar_GDC_LFC_sig_gene_list_GO_only))

names(hemo_dds_deseq_res_Pmar_LFC_sig_gene_list_factor)   <- geneNames
names(hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_gene_list_factor) <- geneNames
names(hemo_dds_deseq_res_Pmar_GDC_LFC_sig_gene_list_factor)  <- geneNames

### Make topGO data object 
Pmar_control_GOdata <- new("topGOdata", description = "Pmar vs control Gene Enrichment", 
                           # I want to test MF
                           ontology = "MF",
                           # define here the genes of interest
                    allGenes = hemo_dds_deseq_res_Pmar_LFC_sig_gene_list_factor,
                    nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)

Pmar_ZVAD_control_GOdata <- new("topGOdata", description = "Pmar vs control Gene Enrichment", 
                           # I want to test MF
                           ontology = "MF",
                           # define here the genes of interest
                           allGenes = hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_gene_list_factor,
                           nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)

Pmar_GDC_control_GOdata <- new("topGOdata", description = "Pmar vs control Gene Enrichment", 
                                # I want to test MF
                                ontology = "MF",
                                # define here the genes of interest
                                allGenes = hemo_dds_deseq_res_Pmar_GDC_LFC_sig_gene_list_factor,
                                nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)

#nodeSize=used to prune the GO hierarchy from the terms which have less than 1 annotated genes
#annFUN.gene2GO = this function is used when the annotations are provided as a gene-to-GOs mapping.

### Perform Encrichment tests 
Pmar_control_Fisher_Weight <- runTest(Pmar_control_GOdata, algorithm = "weight01", statistic = "fisher")
Pmar_ZVAD_control_Fisher_Weight <- runTest(Pmar_ZVAD_control_GOdata, algorithm = "weight01", statistic = "fisher")
Pmar_GDC_control_Fisher_Weight <- runTest(Pmar_GDC_control_GOdata, algorithm = "weight01", statistic = "fisher")


## Analyze enrichment test results 

# see how many results we get where weight01 gives a P-value <= 0.05:
Pmar_control_summary <- summary(attributes(Pmar_control_Fisher_Weight)$score <= 0.05)
# 18 significant

Pmar_ZVAD_control_summary <- summary(attributes(Pmar_ZVAD_control_Fisher_Weight)$score <= 0.05)
# 29 significant

Pmar_GDC_control_summary <- summary(attributes(Pmar_GDC_control_Fisher_Weight)$score <= 0.05)
# 44 significant
  
#print out the top significant results
Pmar_control_Res <- GenTable(Pmar_control_GOdata, topgoFisher = Pmar_control_Fisher_Weight, orderBy = "topgoFisher", topNodes = 18)
Pmar_ZVAD_control_Res <- GenTable(Pmar_ZVAD_control_GOdata, topgoFisher = Pmar_ZVAD_control_Fisher_Weight, orderBy = "topgoFisher", topNodes = 29)
Pmar_GDC_control_Res <- GenTable(Pmar_GDC_control_GOdata, topgoFisher = Pmar_GDC_control_Fisher_Weight, orderBy = "topgoFisher", topNodes = 44)

write.csv(Pmar_control_Res , file = "Pmar_control_Res_GO.csv")
write.csv(Pmar_ZVAD_control_Res , file = "Pmar_ZVAD_control_Res_GO.csv")
write.csv(Pmar_GDC_control_Res , file = "Pmar_GDC_control_Res_GO.csv")

# Plot results in REVIGO

# are there common significant terms across all? 

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
pdf("./FIGURES/pcperk_condition.pdf",height = 5, width = 5)
autoplot(pcperk,
         data = perk_coldata, 
         colour="condition", 
         size=5) + ggtitle("Perkinsus expression")
dev.off()
# little clustering by condition

pdf("./FIGURES/pcperk_pool.pdf",height = 5, width = 5)
autoplot(pcperk,
         data = as.data.frame(perk_coldata), 
         colour="pool", # clustering by pool 
         size=5) + ggtitle("Perkinsus expression")
dev.off()

# clustering mostly by pool, but PC1 and PC2 only explain about 30% of the total sample variation

## Plot total reads in each sample
perk_counts_total <- colSums(perk_counts)
perk_counts_total <- as.data.frame(perk_counts_total)
perk_counts_total <- perk_counts_total %>% mutate(sample_name = rownames(.)) %>% 
  mutate(condition= case_when(
    grepl("GDC", sample_name) ~"GDC",
    grepl("ZVAD", sample_name) ~"ZVAD",
    TRUE~"Dermo"))
perk_counts_total$sample_name <- factor(perk_counts_total$sample_name, 
                                        levels = c("1_Dermo_R1_001", "2_Dermo_R1_001","3_Dermo_R1_001",
                                                   "1_Dermo_GDC_R1_001",  "2_Dermo_GDC_R1_001","3_Dermo_GDC_R1_001" ,
                                                    "1_Dermo_ZVAD_R1_001","2_Dermo_ZVAD_R1_001",  "3_Dermo_ZVAD_R1_001" ))

perk_counts_total_plot <- ggplot(perk_counts_total, aes(x=sample_name, y = perk_counts_total, fill = condition)) + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(perk_counts_total_plot, file = "./FIGURES/perk_counts_total_plot.tiff", device = "tiff",
       height=5, width = 5)

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

# Join with these their Interproscan results
perk_dds_deseq_res_Pmar_GDC_ZVAD_LFC_sig_annot_comb_Interpro <-  left_join(perk_dds_deseq_res_Pmar_GDC_ZVAD_LFC_sig_annot_comb, Perk_Interpro_GO_terms_XP, by = "transcript_id")

# plot LFC 
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_plot <- ggplot(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot, aes(x= product, y = log2FoldChange)) + 
  geom_col() + coord_flip() + ggtitle("P. mar. ZVAD\n vs P. mar control")

perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_plot <- ggplot(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot, aes(x= product, y = log2FoldChange)) +
  geom_col() + coord_flip() + ggtitle("P. mar. GDC\n vs P. mar control")

combined_Pmar <- ggarrange(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_plot, perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_plot)

ggsave(combined_Pmar, device = "tiff",
       width = 12, height = 7,
       file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/combined_Pmar_plot")

## Perkinsus expression upset plot ####
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

## Plot transcformed counts across each sample ####
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons

perk_comb_ID <-  rbind(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot,
                       perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot) %>% dplyr::pull(., transcript_id)
class(perk_comb_ID)

# find rownames matching these IDs
perk_mat <- assay(perk_dds_rlog)[perk_comb_ID,]
perk_mat <- as.data.frame(perk_mat) %>% mutate(transcript_id = rownames(.)) %>% 
  left_join(., unique(Perkinsus_rtracklayer[,c("product","transcript_id")]), by = "transcript_id") %>% filter(!is.na(product)) %>%
  mutate(transcript_product = paste(product, transcript_id, sep = "-")) %>% dplyr::select(-product, -transcript_id)
rownames(perk_mat) <- perk_mat$transcript_product
perk_mat <- perk_mat[,-10]
perk_mat <- as.matrix(perk_mat)  

perk_anno <- as.data.frame(colData(perk_dds_rlog)[, c("condition")])
rownames(perk_anno) <- colnames(perk_mat)
colnames(perk_anno)[1] <- "Condition"

pdf("./FIGURES/perk_mat1.pdf", width = 12, height = 12)
pheatmap(perk_mat , annotation_col = perk_anno)
dev.off()

### Volcano plots of significant genes ####
# compute significance 
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot$log10 <- -log10(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot$padj)
perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot$log10 <- -log10(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot$padj)

# plot the volcano plots
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_volcano_plot <- ggplot(data = as.data.frame(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot),
                                                       aes(x=log2FoldChange, y=log10)) + geom_point() + theme_bw() + 
  labs(y = "-log10(adjusted p-value)", title = "P.mar + ZVAD vs.\ncontrol P. mar") + xlim(-10,15)

perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_volcano_plot <- ggplot(data = as.data.frame(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot), 
                                                            aes(x=log2FoldChange, y=log10)) + geom_point() + theme_bw() +
  labs(y = "-log10(adjusted p-value)", title = "P.mar + GDC vs.\ncontrol P. mar") + xlim(-10,15)


perk_volcano <- ggarrange(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_volcano_plot, 
                          perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_volcano_plot)

ggsave(perk_volcano, file = "./FIGURES/perk_volcano_plot", device = "tiff", height = 3, width = 5)


#### Heatmaps of the genes with the most variable expression patterns ####

topVarGenes_perk_dds_rlog <-  head(order(rowVars(assay(perk_dds_rlog )), decreasing = TRUE), 100)
family_perk_mat <- assay(perk_dds_rlog)[topVarGenes_perk_dds_rlog,]
family_perk_mat <- family_perk_mat - rowMeans(family_perk_mat)
family_perk_anno <- as.data.frame(colData(perk_dds_rlog)[, c("condition")])
rownames(family_perk_anno) <- colnames(family_perk_mat)
family_perk_heatmap <- pheatmap(family_perk_mat , annotation_col = family_perk_anno)
head(family_perk_mat)

# reorder annotation table to match ordering in heatmap 
family_perk_heatmap_reorder <-rownames(family_perk_mat[family_perk_heatmap$tree_row[["order"]],])
# annotate the row.names
family_perk_mat_prot <- as.data.frame(family_perk_heatmap_reorder)
colnames(family_perk_mat_prot)[1] <- "Parent"
family_perk_mat_prot_annot <- left_join(family_perk_mat_prot, dplyr::select(Perkinsus_rtracklayer_XP, Parent, product), by = "Parent") %>% distinct(Parent,product)

# the results of this are not very informative!

#### Export protein lists to run Interproscan ####

# concatenate all transcript ID lists so I can create a lookup list in terminal
perk_comb_ID_df <- as.data.frame(perk_comb_ID) %>% dplyr::rename(Parent = perk_comb_ID)

Perkinsus_rtracklayer_XP <- Perkinsus_rtracklayer %>% filter(!is.na(protein_id)) %>% tidyr::separate(Parent, into = c("rna","Parent"),sep = "-")
perk_comb_ID_sig_XP_df <- left_join(perk_comb_ID_df, unique(Perkinsus_rtracklayer_XP[,c("Parent","protein_id")])) %>% filter(!is.na(protein_id))

# NA's due to XR proteins from non-coding RNA

# export protein IDs to lookup in bluewaves
write.table(perk_comb_ID_sig_XP_df$protein_id, file = "perk_comb_ID_sig_XP_df_lookup.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Also export the full list of transcripts to run in Interproscan 

perk_dds_rlog_mat_XP <- as.data.frame(assay(perk_dds_rlog)) %>% mutate(Parent = rownames(.)) %>%
  left_join(., unique(Perkinsus_rtracklayer_XP[,c("Parent","Name")])) %>% dplyr::select(Name) %>% filter(!is.na(Name))

# export protein IDs to lookup in bluewaves
write.table(perk_dds_rlog_mat_XP$Name, file = "perk_dds_rlog_mat_XP_lookup.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


#### Perkinsus GO term analysis ####

# Download all Interproscan and GO terms for the full Perkinsus transcriptome 
# removed the lines with sequence and no Interproscan information before loading
Perk_Interpro_GO_terms <- rtracklayer::readGFF("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/perk_dds_rlog_mat_XP_lookup_split_all.XP.gff3")
Perk_Interpro_GO_terms <- as.data.frame(Perk_Interpro_GO_terms)
nrow(Perk_Interpro_GO_terms)

# Join with XM terms
Perk_Interpro_GO_terms_XP <- Perk_Interpro_GO_terms %>% dplyr::rename(protein_id = seqid) %>%
  left_join(., unique(Perkinsus_rtracklayer_XP[,c("Parent","protein_id","product")])) %>% dplyr::rename(transcript_id = Parent) %>% mutate(transcript_id_product = paste(product, transcript_id, sep = "-"))
# save this data
save(Perk_Interpro_GO_terms_XP, file="Perk_Interpro_GO_terms_XP.RData")

# Download all GO and Interproscan terms for the significant DEGs
Perk_GO_sig_terms <- rtracklayer::readGFF("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/perk_comb_ID_sig_XP_df_lookup.fa_1.gff3")
Perk_Interpro_sig_terms <- rtracklayer::readGFF("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/perk_comb_ID_sig_XP_df_lookup.fa.gff3")
Perk_GO_sig_terms <- as.data.frame(Perk_GO_sig_terms)
Perk_Interpro_sig_terms <- as.data.frame(Perk_Interpro_sig_terms)
  
# save this
save(Perk_GO_sig_terms_found, file = "Perk_GO_sig_terms_found.RData")

### Format Full GO term list for export and topGO later
# format dataframe for use in topGO for custom annotations

# get the GO terms for each protein 
Perk_GO_terms_found <- Perk_Interpro_GO_terms_XP %>% filter(Ontology_term != "character(0)") %>% 
  distinct(Ontology_term, protein_id, .keep_all = TRUE) 
class(Perk_GO_terms_found$Ontology_term) # AsIs

Perk_GO_terms_found$Ontology_term <- unlist(as.character(Perk_GO_terms_found$Ontology_term))
class(Perk_GO_terms_found$Ontology_term)
Perk_GO_terms_found$Ontology_term <- gsub('\\\"', "", Perk_GO_terms_found$Ontology_term, fixed = TRUE)
Perk_GO_terms_found$Ontology_term <- gsub('\"', "", Perk_GO_terms_found$Ontology_term, fixed = TRUE)
Perk_GO_terms_found$Ontology_term <- gsub('c(', "", Perk_GO_terms_found$Ontology_term, fixed = TRUE)
Perk_GO_terms_found$Ontology_term <- gsub(')', "", Perk_GO_terms_found$Ontology_term, fixed = TRUE)
Perk_GO_terms_found$Ontology_term <- gsub(' ', "", Perk_GO_terms_found$Ontology_term, fixed = TRUE)

Perk_GO_terms_found_geneID2GO <- Perk_GO_terms_found %>% dplyr::select(transcript_id, Ontology_term)
# export as text file to get tab separate file and then re-load using read mappings function
write.table(Perk_GO_terms_found_geneID2GO, file = "Perk_GO_terms_found_geneID2GO.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
Perk_GO_terms_found_geneID2GO_mapping <-  readMappings(file= "Perk_GO_terms_found_geneID2GO.txt")
# save as Rdata for use in WGCNA workspace
save(Perk_GO_terms_found_geneID2GO_mapping, file = "Perk_GO_terms_found_geneID2GO_mapping.RData")


### Make REVIGO plots of only the significant DEG perkinsus transcripts
# get the GO terms for each protein 
Perk_GO_sig_terms_found <- Perk_GO_sig_terms %>% filter(Ontology_term != "character(0)") %>% dplyr::rename(protein_id = seqid) %>% distinct(Ontology_term, protein_id, .keep_all = TRUE)
class(Perk_GO_sig_terms_found$Ontology_term) # AsIs

Perk_GO_sig_terms_found$Ontology_term <- unlist(as.character(Perk_GO_sig_terms_found$Ontology_term))
class(Perk_GO_sig_terms_found$Ontology_term)
Perk_GO_sig_terms_found$Ontology_term <- gsub('\\\"', "", Perk_GO_sig_terms_found$Ontology_term, fixed = TRUE)
Perk_GO_sig_terms_found$Ontology_term <- gsub('\"', "", Perk_GO_sig_terms_found$Ontology_term, fixed = TRUE)
Perk_GO_sig_terms_found$Ontology_term <- gsub('c(', "", Perk_GO_sig_terms_found$Ontology_term, fixed = TRUE)
Perk_GO_sig_terms_found$Ontology_term <- gsub(')', "", Perk_GO_sig_terms_found$Ontology_term, fixed = TRUE)
Perk_GO_sig_terms_found$Ontology_term <- gsub(' ', "", Perk_GO_sig_terms_found$Ontology_term, fixed = TRUE)

View(Perk_GO_sig_terms_found)

# separate each GO term
Perk_GO_sig_terms_found_sep <- Perk_GO_sig_terms_found %>%  separate_rows(Ontology_term, 1, sep = ",") 

# Look up GO terms and use REVIGO to plot by frequency (note this is not an enrichment analysis..just looking at overall frequency)
Perk_GO_sig_terms_found_sep_GO <- Perk_GO_sig_terms_found_sep %>% dplyr::select(Ontology_term)
Perk_GO_sig_terms_found_sep_GO_count <- Perk_GO_sig_terms_found_sep_GO %>% group_by(Ontology_term) %>% count()
# ran REVIGO saying "higher value is better"
# visualized plot by value

# Match Molecular process GO terms to their proteins 
P_marinus_MF <- read.csv("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/MF_REVIGO_P_marinus.csv")
P_marinus_MF <- P_marinus_MF %>% dplyr::rename(Ontology_term = TermID) %>% dplyr::rename(Ontology_desc = Name)
P_marinus_MF$Ontology_desc <- as.character(P_marinus_MF$Ontology_desc)

Perk_GO_sig_terms_found_sep_MF <- left_join(Perk_GO_sig_terms_found_sep,P_marinus_MF[, c("Ontology_desc", "Ontology_term")])

# get list of distinct protein and GO MF IDs 
Perk_GO_sig_terms_found_sep_MF_uniq <- Perk_GO_sig_terms_found_sep_MF %>% distinct(protein_id, Ontology_term, Ontology_desc) %>% filter(!is.na(Ontology_desc))

# join with parent to get XMs
Perk_GO_sig_terms_found_sep_MF_uniq_XM <- Perk_GO_sig_terms_found_sep_MF_uniq %>% dplyr::rename(Name = protein_id) %>% 
  left_join(., unique(Perkinsus_rtracklayer_XP[,c("Name","Parent")])) %>% dplyr::rename(protein_id = Name, transcript_id = Parent)

# Join back with LFC and plot 
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_GO <- left_join(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot,Perk_GO_sig_terms_found_sep_MF_uniq_XM)
perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_GO <- left_join(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot , Perk_GO_sig_terms_found_sep_MF_uniq_XM)

# Plot each list in REVIGO
View(unique(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_GO$Ontology_term))
View(unique(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_GO$Ontology_term))

# plot LFC 
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_GO_plot <- ggplot(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_GO, aes(x= Ontology_desc, y = log2FoldChange)) + 
  geom_col() + coord_flip() + ggtitle("P. mar. ZVAD\n vs P. mar control")

perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_GO_plot <- ggplot(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_GO, aes(x= Ontology_desc, y = log2FoldChange)) + 
  geom_col() + coord_flip() + ggtitle("P. mar. GDC\n vs P. mar control")

# Not very informative! 

# Join instead with Interproscan annotations 
Perk_GO_sig_terms_found_sep_MF_XM <-  Perk_GO_sig_terms_found_sep_MF %>% distinct(protein_id, as.character(Dbxref), as.character(signature_desc)) %>% dplyr::rename(Name = protein_id) %>% 
  left_join(., unique(Perkinsus_rtracklayer_XP[,c("Name","Parent")])) %>% dplyr::rename(protein_id = Name, transcript_id = Parent)

perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_Interpro <- left_join(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot,Perk_GO_sig_terms_found_sep_MF_XM) %>% left_join(., Perk_GO_sig_terms_found_sep_MF_uniq_XM)
perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_Interpro <- left_join(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot, Perk_GO_sig_terms_found_sep_MF_XM)  %>% left_join(., Perk_GO_sig_terms_found_sep_MF_uniq_XM)

# simplify to one line per protein
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_Interpro_uniq <- perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_Interpro %>% distinct(protein_id, Ontology_term, .keep_all = TRUE) %>% distinct(protein_id, .keep_all = TRUE )
perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_Interpro_uniq <- perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_Interpro %>% distinct(protein_id, Ontology_term, .keep_all = TRUE)%>% distinct(protein_id, .keep_all = TRUE )

# rename columns
colnames(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_Interpro_uniq)[11:12] <- c("Dbxref","signature_desc")
colnames(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_Interpro_uniq )[11:12] <- c("Dbxref","signature_desc")

# create combined annotation with product name, Interpro, and GO term for proteins
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_Interpro_uniq_comb <- perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_Interpro_uniq  %>% mutate(comb_product = paste(product,signature_desc, Ontology_desc, sep = "-"))
perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_Interpro_uniq_comb <- perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_Interpro_uniq %>% mutate(comb_product = paste(product,signature_desc,Ontology_desc, sep = "-"))

## Remake expression upset plot but replacing hypothetical proteins with any interesting Interproscan terms

# combine dataframes
Perk_comb_GO <- rbind(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_Interpro_uniq_comb ,
                      perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_Interpro_uniq_comb) %>%
  dplyr::select(comb_product, condition, log2FoldChange)
Perk_comb_spread_GO <- spread(Perk_comb_GO, condition, log2FoldChange, fill = 0)
Perk_comb_spread_GO <- column_to_rownames(Perk_comb_spread_GO , var = "comb_product") 
Perk_comb_spread_GO_mat <- as.matrix(Perk_comb_spread_GO)

Perk_labels =c( "P. mar + GDC vs\nP. mar Control", "P. mar + ZVAD vs\nP. mar Control")
# create named vector to hold column names
Perk_column_labels = structure(paste0(Perk_labels), names = paste0(colnames(Perk_comb_spread_mat)))

pdf("./FIGURES/Perk_spread_GO_mat.pdf", width = 10, height = 8)
Perk_heatmap_GO <- ComplexHeatmap::Heatmap(Perk_comb_spread_GO_mat, border = TRUE, 
                                        column_title_side = "bottom", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                        row_title = "Transcript and Product Name", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                        row_dend_width = unit(2, "cm"),
                                        column_labels = Perk_column_labels[colnames(Perk_comb_spread_mat)],
                                        # apply split by k-meams clustering to highlight groups
                                        #row_km = 3, column_km = 2, 
                                        row_names_gp = gpar(fontsize = 7),
                                        column_names_gp = gpar(fontsize = 8),
                                        heatmap_legend_param = list(title = "Log2 Fold Change"))
ComplexHeatmap::draw(Perk_heatmap_GO, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 100), "mm")) #bottom, left, top, right paddings

dev.off()

## Analyze Interproscan resuls of the significant GDC and ZVAD results

perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_Interpro_all <- left_join(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot,Perk_Interpro_GO_terms_XP,by = "transcript_id") 
perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_Interpro_all  <- left_join(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot,Perk_Interpro_GO_terms_XP, by = "transcript_id") 

# export to csv so I can copy and paste in Interpro descriptions
write.table(perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_Interpro_all, file = "perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot_Interpro_all.txt", sep= "\t", row.names = FALSE)
write.table(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_Interpro_all, file = "perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot_Interpro_all.txt", sep= "\t", row.names = FALSE)


#### HEMOCYTE PCA ANALYSIS ####

# helpful PCA plot tutorial
#https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

# Load phenotype data

load("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/PCA_pheno_2020.RData")
PCA_pheno_2020

# set the sample names to be the same as the count data names for each pool and treatment 
colnames(hemo_counts)
#[1] "1_Dermo_GDC_R1_001"  "1_Dermo_ZVAD_R1_001" "1_Dermo_R1_001"      "1_control_R1_001"    "2_Dermo_GDC_R1_001"  "2_Dermo_ZVAD_R1_001"
#[7] "2_Dermo_R1_001"      "2_control_R1_001"    "3_Dermo_GDC_R1_001"  "3_Dermo_ZVAD_R1_001" "3_Dermo_R1_001"      "3_control_R1_001"   

# mutate treatment and Pool ID
PCA_pheno_2020_all <- PCA_pheno_2020 %>% 
  # remove the beads data
  filter(Treat != "BEADS_LPS") %>%
  mutate(ID = case_when(
  ID == "Pool1" ~ "1",
  ID == "Pool2" ~ "2",
  ID == "Pool3" ~ "3",
  TRUE ~ NA_character_)) %>%
  mutate(Treat = case_when(
    Treat == "Control_hemo" ~ "control_R1_001",
    Treat == "Dermo_GDC" ~ "Dermo_GDC_R1_001",
    Treat == "Dermo_ZVAD" ~ "Dermo_ZVAD_R1_001",
    Treat == "Dermo" ~ "Dermo_R1_001",
    TRUE ~ NA_character_)) %>% mutate(Sample_Name = paste(ID,Treat, sep = "_"))

#### HEMOCYTE APOPTOSIS EXPRESSION PCA ###

# Use the matrix generated in the hemocyte section to plot apoptosis transcript expression
colnames(apop_hemo_mat)

# separate phenotype data into two groups, one with control data and one without
PCA_pheno_2020_control <- PCA_pheno_2020_all %>% dplyr::select(-contains("hemo_perk")) 
PCA_pheno_2020_perk <- PCA_pheno_2020_all %>% filter(!grepl("control",Treat)) 

# make metadata table with just ID and treat
PCA_pheno_2020_control_metadata <- PCA_pheno_2020_control %>% # set rownames to sample name and then remove
  column_to_rownames(., var = "Sample_Name") %>% dplyr::select(ID,Treat) 
PCA_pheno_2020_perk_metadata <- PCA_pheno_2020_control %>% # set rownames to sample name and then remove
  column_to_rownames(., var = "Sample_Name") %>% dplyr::select(ID,Treat) 
  
## Hemocyte Apoptosis PCA including Control samples

# put metatdata in same order as expression matrix 
PCA_pheno_2020_control_metadata <- PCA_pheno_2020_control_metadata[colnames(apop_hemo_mat),]

# check sample name match between metadata and expression data
all(colnames(apop_hemo_mat) == rownames(PCA_pheno_2020_control_metadata)) # TRUE

## Join together expression data and phenotype data all into one dataframe
# transpose the metadata table so that the ID column is the 
PCA_pheno_2020_control_trans <- PCA_pheno_2020_control %>% column_to_rownames(., var = "Sample_Name") %>% dplyr::select(-ID)
class(PCA_pheno_2020_control_trans$Percent_of_this_plot_arcsine_APOP_hemo_alone)
PCA_pheno_2020_control_transpose <- transpose(PCA_pheno_2020_control_trans)
rownames(PCA_pheno_2020_control_transpose) <- colnames(PCA_pheno_2020_control_trans)
colnames(PCA_pheno_2020_control_transpose) <- rownames(PCA_pheno_2020_control_trans)
#remove top row
PCA_pheno_2020_control_transpose <- PCA_pheno_2020_control_transpose[-1,]
# put in correct order
PCA_pheno_2020_control_transpose <- PCA_pheno_2020_control_transpose[,colnames(apop_hemo_mat)]
PCA_pheno_2020_control_transpose_mat <- data.matrix(PCA_pheno_2020_control_transpose)

# bind together the apop_hemo_mat and this matrix for all the samples
all(colnames(apop_hemo_mat) == colnames(PCA_pheno_2020_control_transpose)) # TRUE first make sure samples are in the same order
apop_hemo_mat_pheno <- rbind(apop_hemo_mat,PCA_pheno_2020_control_transpose)
class(apop_hemo_mat_pheno)
apop_hemo_mat_pheno <- data.matrix(apop_hemo_mat_pheno)
class(apop_hemo_mat_pheno)

## compute PCAs, remove lower 10% of variables based on variance

# PCA hemo apop expression plus phenotype
hemo_apop_control_pca <- pca(apop_hemo_mat_pheno , metadata = PCA_pheno_2020_control_metadata) 
# PCA phenotype only
pheno_pca <- pca(PCA_pheno_2020_control_transpose_mat, metadata = PCA_pheno_2020_control_metadata) 
# PCA with just the expression data for hemocytes
hemo_apop_control_pca_no_pheno <- pca(apop_hemo_mat , metadata = PCA_pheno_2020_control_metadata) 

## Plot hemocyte apoptosis plus the phenotype data
# scree plot
hemo_apop_control_pca_scree <- screeplot(hemo_apop_control_pca,
          components = getComponents(hemo_apop_control_pca, 1:12))

# biplot
hemo_apop_control_pca_biplot <- biplot(hemo_apop_control_pca, 
       showLoadings = TRUE, 
       ntopLoadings = 20, title = "Biplot of Top 20 Loadings") 

hemo_apop_control_pca_biplot <- hemo_apop_control_pca_biplot + theme(plot.margin = unit(c(0, 0, 0, 0), "null"))
       
# about 70% of the variation explained when apoptosis expression and cell death phenotypes 
  # PC1 explains the difference between the GDC plot and the control and Dermo/ZVAD
  # PC2 explains the difference between control samples and Dermo2ZVAD which is an outlier, and the dermo and ZVAD samples
  # ZVAD and Dermo samples always cluster and the GDC samples always cluster, and the controls cluster

# plot loadings
# For each PC of interest, ‘plotloadings’ determines the variables falling within the top/bottom 5% of the loadings range, 
  # and then creates a final consensus list of these. These variables are then plotted.
  # loadings describe how much each variable contributes to a particular principal component. 
  # Large loadings (positive or negative) indicate that a particular variable has a strong relationship to a particular principal component. 
  # The sign of a loading indicates whether a variable and a principal component are positively or negatively correlated.
hemo_apop_control_pca_loadings <- plotloadings(hemo_apop_control_pca,
             components = getComponents(hemo_apop_control_pca, c(1,2)), # makes point sizes proportional to the loadings
             rangeRetain = 0.1,
             labSize = 3.0,
             absolute = FALSE,
             title = 'Loadings plot of Top 10% variables',
             shape = 23, shapeSizeRange = c(1, 5),
             col = c('limegreen', 'black', 'red'),
             drawConnectors = TRUE)

# Plot all hemo plus pheno data together
hemo_apop_control_pca_multiplot <- cowplot::plot_grid(hemo_apop_control_pca_biplot,  hemo_apop_control_pca_loadings,
                                                      nrow=2, labels = "AUTO", axis = "h",align = "h")


ggsave(plot = hemo_apop_control_pca_multiplot, device = "tiff", filename = "hemo_apop_control_pca_multiplot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES",
       height = 20, width = 20)


# plot just phenotypes no expression data
biplot(pheno_pca, showLoadings = TRUE) 
  # GDC phenotypes cluster closely, while the Dermo and ZVAD are mostly mixed, and then the control samples are more distinguished


# plot just the apoptosis expression and no phenotypes
biplot(hemo_apop_control_pca_no_pheno, 
       #showLoadings = TRUE, 
       ntopLoadings = 20)
  # about 50% of variation explained..having the phenotype data added explains more of the variance in the model
  # again we have the GDC clustering, control clustering, and the Dermo and ZVAD clustering

plotloadings(hemo_apop_control_pca_no_pheno,
             components = getComponents(hemo_apop_control_pca_no_pheno, c(1,2)), # makes point sizes proportional to the loadings
             rangeRetain = 0.1,
             labSize = 2.0,
             absolute = FALSE,
             title = 'Loadings plot',
             subtitle = 'PCs 1 and 2',
             caption = 'Top 10% variables',
             shape = 23, shapeSizeRange = c(1, 5),
             col = c('limegreen', 'black', 'red'),
             drawConnectors = TRUE)


## Hemocyte apoptosis with just Dermo samples
# put metatdata in same order as expression matrix 
# remove control samples from apop_hemo_mat
apop_hemo_mat_perk <- as.data.frame(apop_hemo_mat) %>% dplyr::select(-contains("control"))
PCA_pheno_2020_perk_metadata <- PCA_pheno_2020_perk_metadata[colnames(apop_hemo_mat_perk),]

# check sample name match between metadata and expression data
all(colnames(apop_hemo_mat_perk) == rownames(PCA_pheno_2020_perk_metadata)) # TRUE

## Join together expression data and phenotype data all into one dataframe
# transpose the metadata table so that the ID column is the 
PCA_pheno_2020_perk_trans <- PCA_pheno_2020_perk %>% column_to_rownames(., var = "Sample_Name") %>% dplyr::select(-ID)
class(PCA_pheno_2020_perk_trans$Percent_of_this_plot_arcsine_APOP_hemo_alone)

PCA_pheno_2020_perk_transpose <- transpose(PCA_pheno_2020_perk_trans)
rownames(PCA_pheno_2020_perk_transpose) <- colnames(PCA_pheno_2020_perk_trans)
colnames(PCA_pheno_2020_perk_transpose) <- rownames(PCA_pheno_2020_perk_trans)
#remove top row
PCA_pheno_2020_perk_transpose <- PCA_pheno_2020_perk_transpose[-1,]
# put in correct order
PCA_pheno_2020_perk_transpose <- PCA_pheno_2020_perk_transpose[,colnames(apop_hemo_mat_perk)]
PCA_pheno_2020_perk_transpose_mat <- data.matrix(PCA_pheno_2020_perk_transpose)

# bind together the apop_hemo_mat and this matrix for all the samples
all(colnames(apop_hemo_mat_perk) == colnames(PCA_pheno_2020_perk_transpose)) # TRUE first make sure samples are in the same order
apop_hemo_mat_pheno_perk <- rbind(apop_hemo_mat_perk,PCA_pheno_2020_perk_transpose)
class(apop_hemo_mat_pheno_perk)
apop_hemo_mat_pheno_perk <- data.matrix(apop_hemo_mat_pheno)
class(apop_hemo_mat_pheno_perk)

## compute PCAs
# PCA with apoptosis phenotype and gene expression
hemo_apop_perk_pca <- pca(apop_hemo_mat_pheno_perk , metadata = PCA_pheno_2020_perk_metadata) 
# only phenotype
pheno_pca_perk <- pca(PCA_pheno_2020_perk_transpose_mat, metadata = PCA_pheno_2020_perk_metadata) 
# only apoptosis expression
hemo_apop_perk_pca_no_pheno <- pca(apop_hemo_mat_perk, metadata = PCA_pheno_2020_perk_metadata)

## plot Hemocyte expression plus phenotype

screeplot(hemo_apop_perk_pca, components = getComponents(hemo_apop_perk_pca, 1:20))

# biplot
hemo_apop_perk_pca_bipplot <- biplot(hemo_apop_perk_pca, ntopLoadings = 20, showLoadings = TRUE,
                             title = "Top 20 Loadings") # about 70% of the variation explained when apoptosis expression and cell death phenotypes 
# with controls removed, PC1 explains more of the variation while PC2 explains slightly less
# PC1 explains the difference between ZVAD and control , and the GDC samples.. so it really segregates the response to GDC
# PC2 really only separates the GDC and Dermo/ZVAD from one ZVAD sample

# plot loadings
# For each PC of interest, ‘plotloadings’ determines the variables falling within the top/bottom 5% of the loadings range, 
# and then creates a final consensus list of these. These variables are then plotted.
# loadings describe how much each variable contributes to a particular principal component. 
# Large loadings (positive or negative) indicate that a particular variable has a strong relationship to a particular principal component. 
# The sign of a loading indicates whether a variable and a principal component are positively or negatively correlated.
hemo_apop_perk_pca_loadings <- plotloadings(hemo_apop_perk_pca,
             components = getComponents(hemo_apop_perk_pca, c(1,2)), # makes point sizes proportional to the loadings
             rangeRetain = 0.1,
             labSize = 3.0,
             absolute = FALSE,
             title = 'Loadings plot of Top 10% variables',
             shape = 23, shapeSizeRange = c(1, 5),
             col = c('limegreen', 'black', 'red'),
             drawConnectors = TRUE)

# only need to pay attention to the PC1 transcripts that segregate
#for PC1 the mitochondrial response is strongly negatively correlated with PC1
# anything in the negatives here are more correlated with the GDC response
  # only caspase 8 and the mitochondrial phenotype are correlated with the GDC at a level of top 10% of variables

# Plot all hemo plus pheno data together
hemo_apop_perk_pca_multiplot <- cowplot::plot_grid(hemo_apop_perk_pca_bipplot,  hemo_apop_perk_pca_loadings,
                                                      nrow=2, labels = "AUTO", axis = "h",align = "h")

ggsave(plot = hemo_apop_perk_pca_multiplot, device = "tiff", filename = "hemo_apop_perk_pca_multiplot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES",
       height = 20, width = 20)


## plot just the phenotypes

biplot(pheno_pca_perk, showLoadings = TRUE)
# the firs two principal components explain almost all the variance..with the mitochondrial assay and the apoptosis assay results being
# strongly negatively correlated. GDC sample still cluster strongly and the dermo only and ZVAD cluster together

#### HEMOCYTE PCA WITH MEAN EXPRESSION ####

# helpful PCA plot tutorial
#https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

# Load phenotype data

load("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/COMBINED_ANALYSIS/R_ANALYSIS/PCA_pheno_2020.RData")
PCA_pheno_2020

# set the sample names to be the same as the count data names for each pool and treatment 
colnames(hemo_counts)
#[1] "1_Dermo_GDC_R1_001"  "1_Dermo_ZVAD_R1_001" "1_Dermo_R1_001"      "1_control_R1_001"    "2_Dermo_GDC_R1_001"  "2_Dermo_ZVAD_R1_001"
#[7] "2_Dermo_R1_001"      "2_control_R1_001"    "3_Dermo_GDC_R1_001"  "3_Dermo_ZVAD_R1_001" "3_Dermo_R1_001"      "3_control_R1_001"   

# mutate treatment and Pool ID
PCA_pheno_2020_all <- PCA_pheno_2020 %>% 
  # remove the beads data
  filter(Treat != "BEADS_LPS") %>%
  mutate(ID = case_when(
    ID == "Pool1" ~ "1",
    ID == "Pool2" ~ "2",
    ID == "Pool3" ~ "3",
    TRUE ~ NA_character_)) %>%
  mutate(Treat = case_when(
    Treat == "Control_hemo" ~ "control_R1_001",
    Treat == "Dermo_GDC" ~ "Dermo_GDC_R1_001",
    Treat == "Dermo_ZVAD" ~ "Dermo_ZVAD_R1_001",
    Treat == "Dermo" ~ "Dermo_R1_001",
    TRUE ~ NA_character_)) %>% mutate(Sample_Name = paste(ID,Treat, sep = "_"))

#### HEMOCYTE APOPTOSIS EXPRESSION PCA ###

# Use the matrix generated in the hemocyte section to plot apoptosis transcript expression
colnames(apop_hemo_mat)

# separate phenotype data into two groups, one with control data and one without
PCA_pheno_2020_control <- PCA_pheno_2020_all %>% dplyr::select(-contains("hemo_perk")) 
PCA_pheno_2020_perk <- PCA_pheno_2020_all %>% filter(!grepl("control",Treat)) 

# make metadata table with just ID and treat
PCA_pheno_2020_control_metadata <- PCA_pheno_2020_control %>% # set rownames to sample name and then remove
  column_to_rownames(., var = "Sample_Name") %>% dplyr::select(ID,Treat) 
PCA_pheno_2020_perk_metadata <- PCA_pheno_2020_control %>% # set rownames to sample name and then remove
  column_to_rownames(., var = "Sample_Name") %>% dplyr::select(ID,Treat) 

## Hemocyte Apoptosis PCA including Control samples

# put metatdata in same order as expression matrix 
PCA_pheno_2020_control_metadata <- PCA_pheno_2020_control_metadata[colnames(apop_hemo_mat),]

# check sample name match between metadata and expression data
all(colnames(apop_hemo_mat) == rownames(PCA_pheno_2020_control_metadata)) # TRUE

## Join together expression data and phenotype data all into one dataframe
# transpose the metadata table so that the ID column is the 
PCA_pheno_2020_control_trans <- PCA_pheno_2020_control %>% column_to_rownames(., var = "Sample_Name") %>% dplyr::select(-ID)
class(PCA_pheno_2020_control_trans$Percent_of_this_plot_arcsine_APOP_hemo_alone)
PCA_pheno_2020_control_transpose <- transpose(PCA_pheno_2020_control_trans)
rownames(PCA_pheno_2020_control_transpose) <- colnames(PCA_pheno_2020_control_trans)
colnames(PCA_pheno_2020_control_transpose) <- rownames(PCA_pheno_2020_control_trans)
#remove top row
PCA_pheno_2020_control_transpose <- PCA_pheno_2020_control_transpose[-1,]
# put in correct order
PCA_pheno_2020_control_transpose <- PCA_pheno_2020_control_transpose[,colnames(apop_hemo_mat)]
PCA_pheno_2020_control_transpose_mat <- data.matrix(PCA_pheno_2020_control_transpose)

# bind together the apop_hemo_mat and this matrix for all the samples
all(colnames(apop_hemo_mat) == colnames(PCA_pheno_2020_control_transpose)) # TRUE first make sure samples are in the same order
apop_hemo_mat_pheno <- rbind(apop_hemo_mat,PCA_pheno_2020_control_transpose)
class(apop_hemo_mat_pheno)
apop_hemo_mat_pheno <- data.matrix(apop_hemo_mat_pheno)
class(apop_hemo_mat_pheno)

## compute PCAs, remove lower 10% of variables based on variance

# PCA hemo apop expression plus phenotype
hemo_apop_control_pca <- pca(apop_hemo_mat_pheno , metadata = PCA_pheno_2020_control_metadata) 
# PCA phenotype only
pheno_pca <- pca(PCA_pheno_2020_control_transpose_mat, metadata = PCA_pheno_2020_control_metadata) 
# PCA with just the expression data for hemocytes
hemo_apop_control_pca_no_pheno <- pca(apop_hemo_mat , metadata = PCA_pheno_2020_control_metadata) 

## Plot hemocyte apoptosis plus the phenotype data
# scree plot
hemo_apop_control_pca_scree <- screeplot(hemo_apop_control_pca,
                                         components = getComponents(hemo_apop_control_pca, 1:12))

# biplot
hemo_apop_control_pca_biplot <- biplot(hemo_apop_control_pca, 
                                       showLoadings = TRUE, 
                                       ntopLoadings = 20, title = "Biplot of Top 20 Loadings") 

hemo_apop_control_pca_biplot <- hemo_apop_control_pca_biplot + theme(plot.margin = unit(c(0, 0, 0, 0), "null"))

# about 70% of the variation explained when apoptosis expression and cell death phenotypes 
# PC1 explains the difference between the GDC plot and the control and Dermo/ZVAD
# PC2 explains the difference between control samples and Dermo2ZVAD which is an outlier, and the dermo and ZVAD samples
# ZVAD and Dermo samples always cluster and the GDC samples always cluster, and the controls cluster

# plot loadings
# For each PC of interest, ‘plotloadings’ determines the variables falling within the top/bottom 5% of the loadings range, 
# and then creates a final consensus list of these. These variables are then plotted.
# loadings describe how much each variable contributes to a particular principal component. 
# Large loadings (positive or negative) indicate that a particular variable has a strong relationship to a particular principal component. 
# The sign of a loading indicates whether a variable and a principal component are positively or negatively correlated.
hemo_apop_control_pca_loadings <- plotloadings(hemo_apop_control_pca,
                                               components = getComponents(hemo_apop_control_pca, c(1,2)), # makes point sizes proportional to the loadings
                                               rangeRetain = 0.1,
                                               labSize = 3.0,
                                               absolute = FALSE,
                                               title = 'Loadings plot of Top 10% variables',
                                               shape = 23, shapeSizeRange = c(1, 5),
                                               col = c('limegreen', 'black', 'red'),
                                               drawConnectors = TRUE)

# Plot all hemo plus pheno data together
hemo_apop_control_pca_multiplot <- cowplot::plot_grid(hemo_apop_control_pca_biplot,  hemo_apop_control_pca_loadings,
                                                      nrow=2, labels = "AUTO", axis = "h",align = "h")


ggsave(plot = hemo_apop_control_pca_multiplot, device = "tiff", filename = "hemo_apop_control_pca_multiplot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES",
       height = 20, width = 20)


# plot just phenotypes no expression data
biplot(pheno_pca, showLoadings = TRUE) 
# GDC phenotypes cluster closely, while the Dermo and ZVAD are mostly mixed, and then the control samples are more distinguished


# plot just the apoptosis expression and no phenotypes
biplot(hemo_apop_control_pca_no_pheno, 
       #showLoadings = TRUE, 
       ntopLoadings = 20)
# about 50% of variation explained..having the phenotype data added explains more of the variance in the model
# again we have the GDC clustering, control clustering, and the Dermo and ZVAD clustering

plotloadings(hemo_apop_control_pca_no_pheno,
             components = getComponents(hemo_apop_control_pca_no_pheno, c(1,2)), # makes point sizes proportional to the loadings
             rangeRetain = 0.1,
             labSize = 2.0,
             absolute = FALSE,
             title = 'Loadings plot',
             subtitle = 'PCs 1 and 2',
             caption = 'Top 10% variables',
             shape = 23, shapeSizeRange = c(1, 5),
             col = c('limegreen', 'black', 'red'),
             drawConnectors = TRUE)


## Hemocyte apoptosis with just Dermo samples
# put metatdata in same order as expression matrix 
# remove control samples from apop_hemo_mat
apop_hemo_mat_perk <- as.data.frame(apop_hemo_mat) %>% dplyr::select(-contains("control"))
PCA_pheno_2020_perk_metadata <- PCA_pheno_2020_perk_metadata[colnames(apop_hemo_mat_perk),]

# check sample name match between metadata and expression data
all(colnames(apop_hemo_mat_perk) == rownames(PCA_pheno_2020_perk_metadata)) # TRUE

## Join together expression data and phenotype data all into one dataframe
# transpose the metadata table so that the ID column is the 
PCA_pheno_2020_perk_trans <- PCA_pheno_2020_perk %>% column_to_rownames(., var = "Sample_Name") %>% dplyr::select(-ID)
class(PCA_pheno_2020_perk_trans$Percent_of_this_plot_arcsine_APOP_hemo_alone)

PCA_pheno_2020_perk_transpose <- transpose(PCA_pheno_2020_perk_trans)
rownames(PCA_pheno_2020_perk_transpose) <- colnames(PCA_pheno_2020_perk_trans)
colnames(PCA_pheno_2020_perk_transpose) <- rownames(PCA_pheno_2020_perk_trans)
#remove top row
PCA_pheno_2020_perk_transpose <- PCA_pheno_2020_perk_transpose[-1,]
# put in correct order
PCA_pheno_2020_perk_transpose <- PCA_pheno_2020_perk_transpose[,colnames(apop_hemo_mat_perk)]
PCA_pheno_2020_perk_transpose_mat <- data.matrix(PCA_pheno_2020_perk_transpose)

# bind together the apop_hemo_mat and this matrix for all the samples
all(colnames(apop_hemo_mat_perk) == colnames(PCA_pheno_2020_perk_transpose)) # TRUE first make sure samples are in the same order
apop_hemo_mat_pheno_perk <- rbind(apop_hemo_mat_perk,PCA_pheno_2020_perk_transpose)
class(apop_hemo_mat_pheno_perk)
apop_hemo_mat_pheno_perk <- data.matrix(apop_hemo_mat_pheno)
class(apop_hemo_mat_pheno_perk)

## compute PCAs
# PCA with apoptosis phenotype and gene expression
hemo_apop_perk_pca <- pca(apop_hemo_mat_pheno_perk , metadata = PCA_pheno_2020_perk_metadata) 
# only phenotype
pheno_pca_perk <- pca(PCA_pheno_2020_perk_transpose_mat, metadata = PCA_pheno_2020_perk_metadata) 
# only apoptosis expression
hemo_apop_perk_pca_no_pheno <- pca(apop_hemo_mat_perk, metadata = PCA_pheno_2020_perk_metadata)

## plot Hemocyte expression plus phenotype

screeplot(hemo_apop_perk_pca, components = getComponents(hemo_apop_perk_pca, 1:20))

# biplot
hemo_apop_perk_pca_bipplot <- biplot(hemo_apop_perk_pca, ntopLoadings = 20, showLoadings = TRUE,
                                     title = "Top 20 Loadings") # about 70% of the variation explained when apoptosis expression and cell death phenotypes 
# with controls removed, PC1 explains more of the variation while PC2 explains slightly less
# PC1 explains the difference between ZVAD and control , and the GDC samples.. so it really segregates the response to GDC
# PC2 really only separates the GDC and Dermo/ZVAD from one ZVAD sample

# plot loadings
# For each PC of interest, ‘plotloadings’ determines the variables falling within the top/bottom 5% of the loadings range, 
# and then creates a final consensus list of these. These variables are then plotted.
# loadings describe how much each variable contributes to a particular principal component. 
# Large loadings (positive or negative) indicate that a particular variable has a strong relationship to a particular principal component. 
# The sign of a loading indicates whether a variable and a principal component are positively or negatively correlated.
hemo_apop_perk_pca_loadings <- plotloadings(hemo_apop_perk_pca,
                                            components = getComponents(hemo_apop_perk_pca, c(1,2)), # makes point sizes proportional to the loadings
                                            rangeRetain = 0.1,
                                            labSize = 3.0,
                                            absolute = FALSE,
                                            title = 'Loadings plot of Top 10% variables',
                                            shape = 23, shapeSizeRange = c(1, 5),
                                            col = c('limegreen', 'black', 'red'),
                                            drawConnectors = TRUE)

# only need to pay attention to the PC1 transcripts that segregate
#for PC1 the mitochondrial response is strongly negatively correlated with PC1
# anything in the negatives here are more correlated with the GDC response
# only caspase 8 and the mitochondrial phenotype are correlated with the GDC at a level of top 10% of variables

# Plot all hemo plus pheno data together
hemo_apop_perk_pca_multiplot <- cowplot::plot_grid(hemo_apop_perk_pca_bipplot,  hemo_apop_perk_pca_loadings,
                                                   nrow=2, labels = "AUTO", axis = "h",align = "h")

ggsave(plot = hemo_apop_perk_pca_multiplot, device = "tiff", filename = "hemo_apop_perk_pca_multiplot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES",
       height = 20, width = 20)


## plot just the phenotypes

biplot(pheno_pca_perk, showLoadings = TRUE)
# the firs two principal components explain almost all the variance..with the mitochondrial assay and the apoptosis assay results being
# strongly negatively correlated. GDC sample still cluster strongly and the dermo only and ZVAD cluster together



#### PERKKINSUS EXPRESSION PCA ####

# Use the phenotype data and metadata tables that were supset for just the Perkinsus data above 
PCA_pheno_2020_perk
PCA_pheno_2020_perk_metadata

# get full counts matrix from all perkinsus data
perk_dds_rlog_mat <- assay(perk_dds_rlog)

# make the product name the rownames so I can interpret the biplots better
perk_dds_rlog_mat <- as.data.frame(perk_dds_rlog_mat) %>% mutate(Parent = rownames(.)) %>%
  left_join(., unique(Perkinsus_rtracklayer_XP[,c("product","Parent")])) %>% mutate(transcript_id_product = paste(product, Parent, sep = "-")) %>%
    column_to_rownames(., var = "transcript_id_product") %>% dplyr::select(-product, - Parent)
perk_dds_rlog_mat <- as.matrix(perk_dds_rlog_mat)
colnames(perk_dds_rlog_mat)

# check sample name match between metadata and expression data
all(colnames(perk_dds_rlog_mat) == rownames(PCA_pheno_2020_perk_metadata)) # TRUE

# bind together the perk_dds_rlog_mat and the previously transposed data for just the perkinsus samples
all(colnames(perk_dds_rlog_mat) == colnames(PCA_pheno_2020_perk_transpose)) # TRUE first make sure samples are in the same order

perk_dds_rlog_mat_pheno_perk <- rbind(perk_dds_rlog_mat,PCA_pheno_2020_perk_transpose)
class(perk_dds_rlog_mat_pheno_perk)
perk_dds_rlog_mat_pheno_perk <- data.matrix(perk_dds_rlog_mat_pheno_perk)
class(perk_dds_rlog_mat_pheno_perk)

## compute PCAs
# PCA with apoptosis phenotype and gene expression
perk_dds_rlog_mat_pheno_perk_pca <- pca(perk_dds_rlog_mat_pheno_perk, metadata = PCA_pheno_2020_perk_metadata) 

# only perkinsus expression
perk_dds_rlog_mat_pheno_perk_no_pheno <- pca(perk_dds_rlog_mat, metadata = PCA_pheno_2020_perk_metadata)

## plot Perkinsus expression plus phenotype

#plot PCs as screeplot
perk_dds_rlog_mat_pheno_perk_pca_scree <- screeplot(perk_dds_rlog_mat_pheno_perk_pca, components = getComponents(perk_dds_rlog_mat_pheno_perk_pca, 1:20))
  # variation is spread on a lot of PCs

# biplot
perk_dds_rlog_mat_pheno_perk_pca_bipplot_1_2 <- biplot(perk_dds_rlog_mat_pheno_perk_pca, ntopLoadings = 25, showLoadings = TRUE,
                                     title = "Top 25 Loadings", x="PC1",y="PC2") 
# PCs 1 and 2 explain only about 30% of the total variation. Not great clustering 
perk_dds_rlog_mat_pheno_perk_pca_bipplot_1_3 <- biplot(perk_dds_rlog_mat_pheno_perk_pca, ntopLoadings = 20, showLoadings = TRUE,
                                                   title = "Top 20 Loadings", x="PC1",y="PC3") 

# plot as pairings plot to view all PCs
perk_dds_rlog_mat_pheno_perk_pca_pairs_Treat <- pairsplot(perk_dds_rlog_mat_pheno_perk_pca, colby = "Treat")
  # PCs 1,2, and 3 explain most of the variation and have some grouping by treatment
perk_dds_rlog_mat_pheno_perk_pca_pairs_ID <- pairsplot(perk_dds_rlog_mat_pheno_perk_pca, colby = "ID")

# plot loadings
# For each PC of interest, ‘plotloadings’ determines the variables falling within the top/bottom 5% of the loadings range, 
# and then creates a final consensus list of these. These variables are then plotted.
# loadings describe how much each variable contributes to a particular principal component. 
# Large loadings (positive or negative) indicate that a particular variable has a strong relationship to a particular principal component. 
# The sign of a loading indicates whether a variable and a principal component are positively or negatively correlated.
perk_dds_rlog_mat_pheno_perk_pca_loadings <- plotloadings(perk_dds_rlog_mat_pheno_perk_pca,
                                            components = getComponents(perk_dds_rlog_mat_pheno_perk_pca, c(1,2)), # makes point sizes proportional to the loadings
                                            rangeRetain = 0.05,
                                            labSize = 3.0,
                                            absolute = FALSE,
                                            title = 'Loadings plot of Top 10% variables',
                                            shape = 23, shapeSizeRange = c(1, 5),
                                            col = c('limegreen', 'black', 'red'),
                                            drawConnectors = TRUE)

# assess the variables from the 1% loading cutoff
perk_PC1_2_1 <- data.frame("transcript_id_product" = c(
                              "UDP-N-acteylglucosamine pyrophosphorylase, putative-XM_002786674.1" ,
                              "hypothetical protein-XM_002773609.1",
                              "calmodulin-domain protein kinase, putative-XM_002767007.1")) %>%
                                left_join(., Perk_Interpro_GO_terms_XP)

# Plot all hemo plus pheno data together
perk_dds_rlog_mat_pheno_perk_pca_multiplot <- cowplot::plot_grid(perk_dds_rlog_mat_pheno_perk_pca_scree,  perk_dds_rlog_mat_pheno_perk_pca_bipplot_1_2,
                                                                 nrow=2, labels = "AUTO", axis = "h",align = "h", rel_heights = c(0.3,1), rel_widths = c(0.5,1))

ggsave(plot = perk_dds_rlog_mat_pheno_perk_pca_multiplot, device = "tiff", filename = "perk_dds_rlog_mat_pheno_perk_pca_multiplot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES",
       height = 20, width = 20)

### Examine loadings from PC1 and PC2 
perk_dds_rlog_mat_pheno_perk_pca_loadings <- perk_dds_rlog_mat_pheno_perk_pca$loadings %>% rownames_to_column(., var = "transcript_id_product")

## Examine PC1 loadings
# sort by the top absolute value loadings in PC1
perk_dds_rlog_mat_pheno_perk_pca_loadings_PC1 <- perk_dds_rlog_mat_pheno_perk_pca_loadings[order(abs(perk_dds_rlog_mat_pheno_perk_pca_loadings$PC1), decreasing = TRUE),]

# analyze top 25 
perk_dds_rlog_mat_pheno_perk_pca_loadings_PC1_50 <- perk_dds_rlog_mat_pheno_perk_pca_loadings_PC1[1:50,]

# join with the Interproscan information
perk_dds_rlog_mat_pheno_perk_pca_loadings_PC1_50_Interpro <- perk_dds_rlog_mat_pheno_perk_pca_loadings_PC1_50 %>% dplyr::select(transcript_id_product, PC1,PC2) %>% 
  left_join(., Perk_Interpro_GO_terms_XP)

## Examine PC2 loadings
# sort by the top absolute value loadings in PC2
perk_dds_rlog_mat_pheno_perk_pca_loadings_PC2 <- perk_dds_rlog_mat_pheno_perk_pca_loadings[order(abs(perk_dds_rlog_mat_pheno_perk_pca_loadings$PC2), decreasing = TRUE),]

# analyze top 25 
perk_dds_rlog_mat_pheno_perk_pca_loadings_PC2_25 <- perk_dds_rlog_mat_pheno_perk_pca_loadings_PC2[1:25,]

# join with the Interproscan information
perk_dds_rlog_mat_pheno_perk_pca_loadings_PC2_25_Interpro <- perk_dds_rlog_mat_pheno_perk_pca_loadings_PC2_25 %>% dplyr::select(transcript_id_product, PC1,PC2) %>% 
                                                                         left_join(., Perk_Interpro_GO_terms_XP)

## review particular transcripts of interest in the raw data - why weren't some of them differentially expressed?

# WAG22 antigen precursor, putative-XM_002780500.1 
assay(perk_dds_rlog)["XM_002780500.1",] # much higher expression in GDC-1 likely just due to random chance and not due to a real treatment effect
# ABC transporter, putative-XM_002773540.1
assay(perk_dds_rlog)["XM_002773540.1",] # really is much higher expression across all GDC samples
#ATP-dependent RNA helicase, putative-XM_002788234.1 
assay(perk_dds_rlog)["XM_002788234.1",] # this is similarly scewed due to GDC 1

#### EXPRESSION PCAS WITH THE MEAN EXPRESSION ####

# Performing the PCAs with mean expression should help reduce skewing based on expression in a single treatment

# get full counts matrix from all perkinsus data
perk_dds_rlog_mat_df <- as.data.frame(perk_dds_rlog_mat)

# get mean across only the rows for each group
perk_dds_rlog_mat_df$Dermo_mean <- rowMeans(subset(perk_dds_rlog_mat_df, select = c("1_Dermo_R1_001", 
                                                                                    "2_Dermo_R1_001",
                                                                                    "3_Dermo_R1_001")), na.rm = TRUE)
perk_dds_rlog_mat_df$Dermo_GDC_mean <- rowMeans(subset(perk_dds_rlog_mat_df, select = c("1_Dermo_GDC_R1_001", 
                                                                                    "2_Dermo_GDC_R1_001",
                                                                                    "3_Dermo_GDC_R1_001")), na.rm = TRUE)

perk_dds_rlog_mat_df$Dermo_ZVAD_mean <- rowMeans(subset(perk_dds_rlog_mat_df, select =c("1_Dermo_ZVAD_R1_001", 
                                                                                        "2_Dermo_ZVAD_R1_001",
                                                                                        "3_Dermo_ZVAD_R1_001")), na.rm = TRUE)

# Subset just these mean counts columns
colnames(perk_dds_rlog_mat_df)
perk_dds_rlog_mat_df_rowmeans <- perk_dds_rlog_mat_df[,c("Dermo_mean","Dermo_GDC_mean","Dermo_ZVAD_mean")]
perk_dds_rlog_mat_df_rowmeans <- as.matrix(perk_dds_rlog_mat_df_rowmeans)

# Format metadata to now match the rowmeans dataframe
PCA_pheno_2020_perk_metadata_rowmeans <- PCA_pheno_2020_perk_metadata %>% mutate(Treat = case_when(
      Treat == "Dermo_R1_001" ~ "Dermo_mean",
      Treat == "Dermo_GDC_R1_001" ~ "Dermo_GDC_mean",
      Treat == "Dermo_ZVAD_R1_001" ~ "Dermo_ZVAD_mean")) %>% distinct(Treat) %>% column_to_rownames(., var = "Treat")
# fix order
PCA_pheno_2020_perk_metadata_rowmeans <- PCA_pheno_2020_perk_metadata_rowmeans[colnames(perk_dds_rlog_mat_df_rowmeans),]

# check sample name match between metadata and expression data
all(colnames(perk_dds_rlog_mat_df_rowmeans) == rownames(PCA_pheno_2020_perk_metadata_rowmeans)) # TRUE

# get the mean of the phenotype columns
class(PCA_pheno_2020_perk_transpose)
PCA_pheno_2020_perk_transpose_rowmeans <- PCA_pheno_2020_perk_transpose %>% rownames_to_column(., var = "ID")
PCA_pheno_2020_perk_transpose_rowmeans$ID <- as.factor(PCA_pheno_2020_perk_transpose_rowmeans$ID )
PCA_pheno_2020_perk_transpose_rowmeans <- PCA_pheno_2020_perk_transpose_rowmeans %>%  
  mutate_if(is.character,as.numeric) %>% column_to_rownames(., var="ID")
PCA_pheno_2020_perk_transpose_rowmeans$Dermo_mean <-   rowMeans(subset(PCA_pheno_2020_perk_transpose_rowmeans, 
                                                                       select = c("1_Dermo_R1_001", 
                                                                                  "2_Dermo_R1_001",
                                                                                  "3_Dermo_R1_001")), na.rm = TRUE)
PCA_pheno_2020_perk_transpose_rowmeans$Dermo_GDC_mean <-   rowMeans(subset(PCA_pheno_2020_perk_transpose_rowmeans, 
                                                                       select = c("1_Dermo_GDC_R1_001", 
                                                                                  "2_Dermo_GDC_R1_001",
                                                                                  "3_Dermo_GDC_R1_001")), na.rm = TRUE)
PCA_pheno_2020_perk_transpose_rowmeans$Dermo_ZVAD_mean <-   rowMeans(subset(PCA_pheno_2020_perk_transpose_rowmeans, 
                                                                           select = c("1_Dermo_ZVAD_R1_001", 
                                                                                      "2_Dermo_ZVAD_R1_001",
                                                                                      "3_Dermo_ZVAD_R1_001")), na.rm = TRUE)
# subset for just mean rows 
PCA_pheno_2020_perk_transpose_rowmeans <- PCA_pheno_2020_perk_transpose_rowmeans %>% dplyr::select(Dermo_mean, Dermo_GDC_mean, Dermo_ZVAD_mean)
  
# bind together the rowmeans counts and the phenotype data
perk_dds_rlog_mat_pheno_perk_rowmeans <- rbind(perk_dds_rlog_mat_df_rowmeans,PCA_pheno_2020_perk_transpose_rowmeans)
class(perk_dds_rlog_mat_pheno_perk_rowmeans)
perk_dds_rlog_mat_pheno_perk_rowmeans <- data.matrix(perk_dds_rlog_mat_pheno_perk_rowmeans)
class(perk_dds_rlog_mat_pheno_perk_rowmeans)

## compute PCAs
# PCA with apoptosis phenotype and gene expression
perk_dds_rlog_mat_pheno_perk_rowmeans_pca <- pca(perk_dds_rlog_mat_pheno_perk_rowmeans, metadata = PCA_pheno_2020_perk_metadata_rowmeans) 

# only perkinsus expression
perk_dds_rlog_mat_df_rowmeans_no_pheno <- pca(perk_dds_rlog_mat_df_rowmeans, metadata = PCA_pheno_2020_perk_metadata_rowmeans)

## plot Perkinsus expression plus phenotype

#plot PCs as screeplot
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_scree <- screeplot(perk_dds_rlog_mat_pheno_perk_rowmeans_pca, components = getComponents(perk_dds_rlog_mat_pheno_perk_rowmeans_pca, 1:3))
# variation is entirely in the first two PCs

# biplot
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_bipplot_1_2 <- biplot(perk_dds_rlog_mat_pheno_perk_rowmeans_pca, ntopLoadings = 20, showLoadings = TRUE,
                                                       title = "Top 20 Loadings", x="PC1",y="PC2") 

# plot loadings
# For each PC of interest, ‘plotloadings’ determines the variables falling within the top/bottom 5% of the loadings range, 
# and then creates a final consensus list of these. These variables are then plotted.
# loadings describe how much each variable contributes to a particular principal component. 
# Large loadings (positive or negative) indicate that a particular variable has a strong relationship to a particular principal component. 
# The sign of a loading indicates whether a variable and a principal component are positively or negatively correlated.
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_plot <- plotloadings(perk_dds_rlog_mat_pheno_perk_rowmeans_pca,
                                                          components = getComponents(perk_dds_rlog_mat_pheno_perk_rowmeans_pca, c(1,2)), # makes point sizes proportional to the loadings
                                                          rangeRetain = 0.01,
                                                          labSize = 3.0,
                                                          absolute = FALSE,
                                                          title = 'Loadings plot of Top 10% variables',
                                                          shape = 23, shapeSizeRange = c(1, 5),
                                                          col = c('limegreen', 'black', 'red'),
                                                          drawConnectors = TRUE)

# Plot all hemo plus pheno data together
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_multiplot <- cowplot::plot_grid(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_scree,  perk_dds_rlog_mat_pheno_perk_rowmeans_pca_bipplot_1_2,
                                                                 nrow=2, labels = "AUTO", axis = "h",align = "h", rel_heights = c(0.3,1), rel_widths = c(0.5,1))

ggsave(plot = perk_dds_rlog_mat_pheno_perk_rowmeans_pca_multiplot, device = "tiff", filename = "perk_dds_rlog_mat_pheno_perk_rowmeans_pca_multiplot.tiff",
       path = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/FIGURES",
       height = 20, width = 20)

### Examine loadings from PC1 and PC2 
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings <- perk_dds_rlog_mat_pheno_perk_rowmeans_pca$loadings %>% rownames_to_column(., var = "transcript_id_product")

## Examine PC1 loadings
# sort by the top absolute value loadings in PC1
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1 <- perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings[order(abs(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings$PC1), decreasing = TRUE),]

# analyze top 25 
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50 <- perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1[1:50,] %>% filter(!grepl("Percent",transcript_id_product))

# join with the Interproscan information
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro <- perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50 %>% dplyr::select(transcript_id_product, PC1,PC2) %>% 
  left_join(., Perk_Interpro_GO_terms_XP)

# Join with LFC information from ZVAD and GDC experiments
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_ZVAD <- left_join(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro, 
                                                                                    perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot[,c("transcript_id","log2FoldChange","condition")]) %>%
                                                                                  filter(!is.na(log2FoldChange))
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_GDC <- left_join(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro,
                                                                                    perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot[,c("transcript_id","log2FoldChange","condition")]) %>%
                                                                                    filter(!is.na(log2FoldChange))
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_comb <- rbind(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_ZVAD,
                                                                                     perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_GDC)                                                                                    
# how many overlaps in the top 50 loadings?
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_comb %>% filter(grepl("GDC", condition)) %>% distinct(transcript_id) # only 7 overlaps in the top 50 loadings...., 9 in top 100
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_comb %>% filter(grepl("ZVAD", condition)) %>% distinct(transcript_id) # only two ZVAD overlaps

#whats in those top 50 loadings, export to review in excel
write.table(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_GDC, "perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_GDC.txt",row.names = FALSE, sep = "\t")
write.table(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_ZVAD, "perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_Interpro_LFC_ZVAD.txt",row.names = FALSE, sep = "\t")

## Plot the rlog transformed counts of the top 50 loadings as a heatmap
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_rownames <- perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50 %>% column_to_rownames(., "transcript_id_product") 
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_heatmap <- perk_dds_rlog_mat_df[rownames(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_rownames),]
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_heatmap_mat <- as.matrix(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_heatmap)
perk_anno_mean <- data.frame(Condition = c("Dermo_mean","Dermo_GDC_mean","Dermo_ZVAD_mean"))
rownames(perk_anno_mean) <- c("Dermo_mean","Dermo_GDC_mean","Dermo_ZVAD_mean")
perk_anno_mean_comb <- rbind(perk_anno, perk_anno_mean) %>% rownames_to_column(., var = "sample") %>% mutate(ID = case_when(
  grepl("1_D", sample) ~ "1",
  grepl("2_D", sample) ~ "2",
  grepl("3_D", sample) ~ "3",
  grepl("mean", sample) ~ "mean")) %>% column_to_rownames(.,var = "sample")

pdf("./FIGURES/perk_mat1_top_50_loadings1.pdf", width = 12, height = 12)
pheatmap(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC1_50_heatmap_mat , annotation_col = perk_anno_mean_comb)
dev.off()


## Examine PC2 loadings
# sort by the top absolute value loadings
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2 <- perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings[order(abs(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings$PC2), decreasing = TRUE),]

# analyze top 25 
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50 <- perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2[1:50,] %>% filter(!grepl("Percent",transcript_id_product))

# join with the Interproscan information
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro <- perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50 %>% dplyr::select(transcript_id_product, PC2,PC2) %>% 
  left_join(., Perk_Interpro_GO_terms_XP)

# Join with LFC information from ZVAD and GDC experiments
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro_LFC_ZVAD <- left_join(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro, 
                                                                                         perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot[,c("transcript_id","log2FoldChange","condition")]) %>%
  filter(!is.na(log2FoldChange))
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro_LFC_GDC <- left_join(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro,
                                                                                        perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot[,c("transcript_id","log2FoldChange","condition")]) %>%
  filter(!is.na(log2FoldChange))
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro_LFC_comb <- rbind(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro_LFC_ZVAD,
                                                                                     perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro_LFC_GDC)                                                                                    
# how many overlaps in the top 50 loadings?
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro_LFC_comb %>% filter(grepl("GDC", condition)) %>% distinct(transcript_id) # 0 overlaps in the top 50 loadings...., 9 in top 100
perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro_LFC_comb %>% filter(grepl("ZVAD", condition)) %>% distinct(transcript_id) # only 10 ZVAD overlaps

#whats in those top 50 loadings, export to review in excel
write.table(perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro_LFC_ZVAD, "perk_dds_rlog_mat_pheno_perk_rowmeans_pca_loadings_PC2_50_Interpro_LFC_ZVAD.txt",row.names = FALSE, sep = "\t")



## review particular transcripts of interest in the raw data - why weren't some of them differentially expressed?
assay(perk_dds_rlog)["XM_002775492.1",] # very similar expression across all challenges
assay(perk_dds_rlog)["XM_002775509.1",] # protein phosphatase - is higher expression in GDC

#### EXPORT TRANSCRIPTOME DATA AND METADATA FOR WGCNA ANALYSIS ####
save(perk_dds_rlog, hemo_dds_rlog, file = "2021_Hemocyte_Dermo_expression_rlog_matrices.RData") 
save(perk_coldata, hemo_coldata, file = "2021_Hemocyte_Dermo_expression_coldata.RData")


#### BIPLOT SOURCE CODE ####
#' Draw a bi-plot, comparing 2 selected principal components / eigenvectors.
#'
#' @param pcaobj Object of class 'pca' created by pca().
#' @param x A principal component to plot on x-axis. All principal component
#'   names are stored in pcaobj$label.
#' @param y A principal component to plot on y-axis. All principal component
#'   names are stored in pcaobj$label.
#' @param showLoadings Logical, indicating whether or not to overlay
#'   variable loadings.
#' @param ntopLoadings If showLoadings == TRUE, select this many variables
#'   based on absolute ordered variable loading for each PC in the biplot.
#'   As a result of looking across 2 PCs, it can occur whereby greater than
#'   this number are actually displayed.
#' @param showLoadingsNames Logical, indicating to show variable loadings names
#'   or not.
#' @param colLoadingsNames If 'showLoadings == TRUE', colour of text labels.
#' @param sizeLoadingsNames If 'showLoadings == TRUE', size of text labels.
#' @param boxedLoadingsNames Logical, if 'showLoadings == TRUE', draw text
#'   labels in boxes.
#' @param fillBoxedLoadings When 'boxedLoadingsNames == TRUE', this controls
#'   the background fill of the boxes. To control both the fill and
#'   transparency, user can specify a value of the form
#'   'alpha(<colour>, <alpha>)'.
#' @param drawConnectorsLoadings If 'showLoadings == TRUE', draw line connectors
#'   to the variable loadings arrows in order to fit more labels in the plot
#'   space.
#' @param widthConnectorsLoadings If 'showLoadings == TRUE', width of the line
#'   connectors drawn to the variable loadings arrows.
#' @param colConnectorsLoadings If 'showLoadings == TRUE', colour of the line
#'   connectors drawn to the variable loadings arrows.
#' @param lengthLoadingsArrowsFactor If 'showLoadings == TRUE', multiply the
#'   internally-determined length of the variable loadings arrows by this
#'   factor.
#' @param colLoadingsArrows If showLoadings == TRUE, colour of the variable
#'   loadings arrows.
#' @param widthLoadingsArrows If showLoadings == TRUE, width of the variable
#'   loadings arrows.
#' @param alphaLoadingsArrow If showLoadings == TRUE, colour transparency of
#'   the variable loadings arrows.
#' @param colby If NULL, all points will be coloured differently. If not NULL,
#'   value is assumed to be a column name in pcaobj$metadata relating to some
#'   grouping/categorical variable.
#' @param colkey Vector of name-value pairs relating to value passed to 'col',
#'   e.g., c(A='forestgreen', B='gold').
#' @param colLegendTitle Title of the legend for the variable specified
#'   by 'colby'.
#' @param singlecol If specified, all points will be shaded by this colour.
#'   Overrides 'col'.
#' @param shape If NULL, all points will be have the same shape. If not NULL,
#'   value is assumed to be a column name in pcaobj$metadata relating to some
#'   grouping/categorical variable.
#' @param shapekey Vector of name-value pairs relating to value passed to
#'   'shape', e.g., c(A=10, B=21).
#' @param shapeLegendTitle Title of the legend for the variable specified
#'   by 'shape'.
#' @param pointSize Size of plotted points.
#' @param legendPosition Position of legend ('top', 'bottom', 'left', 'right',
#'   'none').
#' @param legendLabSize Size of plot legend text.
#' @param legendTitleSize Size of plot legend title text.
#' @param legendIconSize Size of plot legend icons / symbols.
#' @param encircle Logical, indicating whether to draw a polygon around
#'   the groups specified by 'colby'.
#' @param encircleFill Logical, if 'encircle == TRUE', this determines
#'   whether to fill the encircled region or not.
#' @param encircleFillKey Vector of name-value pairs relating to value passed to
#'   'encircleFill', e.g., c(A='forestgreen', B='gold'). If NULL, the fill
#'   is controlled by whatever has already been used for 'colby' / 'colkey'.
#' @param encircleAlpha Alpha for purposes of controlling colour transparency of
#'   the encircled region. Used when 'encircle == TRUE'.
#' @param encircleLineSize Line width of the encircled line when
#'   'encircle == TRUE'.
#' @param encircleLineCol Colour of the encircled line when
#'   'encircle == TRUE'.
#' @param ellipse Logical, indicating whether to draw a stat ellipse around
#'   the groups specified by 'colby'.
#' @param ellipseConf Confidence intervals of the stat ellipses when
#'   ellipse == TRUE.
#' @param ellipseFill Logical, if 'ellipse == TRUE', this determines
#'   whether to fill the region or not.
#' @param ellipseFillKey Vector of name-value pairs relating to value passed to
#'   'ellipseFill', e.g., c(A='forestgreen', B='gold'). If NULL, the fill
#'   is controlled by whatever has already been used for 'colby' / 'colkey'.
#' @param ellipseAlpha Alpha for purposes of controlling colour transparency of
#'   the ellipse region. Used when 'ellipse == TRUE'.
#' @param ellipseLineSize Line width of the ellipse line when 'ellipse == TRUE'.
#' @param ellipseLineCol Colour of the ellipse line when 'ellipse == TRUE'.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#' @param lab A vector containing labels to add to the plot. 
#' @param labSize Size of labels.
#' @param labhjust Horizontal adjustment of label.
#' @param labvjust Vertical adjustment of label.
#' @param boxedLabels Logical, draw text labels in boxes.
#' @param selectLab A vector containing a subset of lab to plot.
#' @param drawConnectors Logical, indicating whether or not to connect plot
#'   labels to their corresponding points by line connectors.
#' @param widthConnectors Line width of connectors.
#' @param colConnectors Line colour of connectors.
#' @param xlab Label for x-axis.
#' @param xlabAngle Rotation angle of x-axis labels.
#' @param xlabhjust Horizontal adjustment of x-axis labels.
#' @param xlabvjust Vertical adjustment of x-axis labels.
#' @param ylab Label for y-axis.
#' @param ylabAngle Rotation angle of y-axis labels.
#' @param ylabhjust Horizontal adjustment of y-axis labels.
#' @param ylabvjust Vertical adjustment of y-axis labels.
#' @param axisLabSize Size of x- and y-axis labels.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param caption Plot caption.
#' @param titleLabSize Size of plot title.
#' @param subtitleLabSize Size of plot subtitle.
#' @param captionLabSize Size of plot caption.
#' @param hline Draw one or more horizontal lines passing through this/these
#'   values on y-axis. For single values, only a single numerical value is
#'   necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
#' @param hlineType Line type for hline ('blank', 'solid', 'dashed', 'dotted',
#'   'dotdash', 'longdash', 'twodash').
#' @param hlineCol Colour of hline.
#' @param hlineWidth Width of hline.
#' @param vline Draw one or more vertical lines passing through this/these
#'   values on x-axis. For single values, only a single numerical value is
#'   necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
#' @param vlineType Line type for vline ('blank', 'solid', 'dashed', 'dotted',
#'   'dotdash', 'longdash', 'twodash').
#' @param vlineCol Colour of vline.
#' @param vlineWidth Width of vline.
#' @param gridlines.major Logical, indicating whether or not to draw major
#'   gridlines.
#' @param gridlines.minor Logical, indicating whether or not to draw minor
#'   gridlines.
#' @param borderWidth Width of the border on the x and y axes.
#' @param borderColour Colour of the border on the x and y axes.
#' @param returnPlot Logical, indicating whether or not to return the plot
#'   object.
#'
#' @details Draw a bi-plot, comparing 2 selected principal components / eigenvectors.
#'
#' @return A \code{\link{ggplot2}} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#'   options(scipen=10)
#'   options(digits=6)
#'
#'   col <- 20
#'   row <- 20000
#'   mat1 <- matrix(
#'     rexp(col*row, rate = 0.1),
#'     ncol = col)
#'   rownames(mat1) <- paste0('gene', 1:nrow(mat1))
#'   colnames(mat1) <- paste0('sample', 1:ncol(mat1))
#'
#'   mat2 <- matrix(
#'   rexp(col*row, rate = 0.1),
#'     ncol = col)
#'   rownames(mat2) <- paste0('gene', 1:nrow(mat2))
#'   colnames(mat2) <- paste0('sample', (ncol(mat1)+1):(ncol(mat1)+ncol(mat2)))
#'
#'   mat <- cbind(mat1, mat2)
#'
#'   metadata <- data.frame(row.names = colnames(mat))
#'   metadata$Group <- rep(NA, ncol(mat))
#'   metadata$Group[seq(1,40,2)] <- 'A'
#'   metadata$Group[seq(2,40,2)] <- 'B'
#'   metadata$CRP <- sample.int(100, size=ncol(mat), replace=TRUE)
#'   metadata$ESR <- sample.int(100, size=ncol(mat), replace=TRUE)
#'
#'   p <- pca(mat, metadata = metadata, removeVar = 0.1)
#'
#'   biplot(p)
#'
#'   biplot(p, colby = 'Group', shape = 'Group')
#'
#'   biplot(p, colby = 'Group', colkey = c(A = 'forestgreen', B = 'gold'),
#'     legendPosition = 'right')
#'
#'   biplot(p, colby = 'Group', colkey = c(A='forestgreen', B='gold'),
#'     shape = 'Group', shapekey = c(A=10, B=21), legendPosition = 'bottom')
#'
#' @import ggplot2
#' @import ggrepel
#' 
#' @export
biplot <- function(
  pcaobj,
  x = 'PC1',
  y = 'PC2',
  showLoadings = FALSE,
  ntopLoadings = 5,
  showLoadingsNames = if (showLoadings) TRUE else FALSE,
  colLoadingsNames = 'black',
  sizeLoadingsNames = 3,
  boxedLoadingsNames = TRUE,
  fillBoxedLoadings = alpha('white', 1/4),
  drawConnectorsLoadings = TRUE,
  widthConnectorsLoadings = 0.5,
  colConnectorsLoadings = 'grey50',
  lengthLoadingsArrowsFactor = 1.5,
  colLoadingsArrows = 'black',
  widthLoadingsArrows = 0.5,
  alphaLoadingsArrow = 1.0,
  colby = NULL,
  colkey = NULL,
  colLegendTitle = if (!is.null(colby)) colby else NULL,
  singlecol = NULL,
  shape = NULL,
  shapekey = NULL,
  shapeLegendTitle = if (!is.null(shape)) shape else NULL,
  pointSize = 3.0,
  legendPosition = 'none',
  legendLabSize = 12,
  legendTitleSize = 14,
  legendIconSize = 5.0,
  encircle = FALSE,
  encircleFill = TRUE,
  encircleFillKey = NULL,
  encircleAlpha = 1/4,
  encircleLineSize = 0.25,
  encircleLineCol = NULL,
  ellipse = FALSE,
  ellipseConf = 0.95,
  ellipseFill = TRUE,
  ellipseFillKey = NULL,
  ellipseAlpha = 1/4,
  ellipseLineSize = 0.25,
  ellipseLineCol = NULL,
  xlim = if(showLoadings) c(min(pcaobj$rotated[,x]) - 5, max(pcaobj$rotated[,x]) + 5) else
    c(min(pcaobj$rotated[,x]) - 1, max(pcaobj$rotated[,x]) + 1),
  ylim = if(showLoadings) c(min(pcaobj$rotated[,y]) - 5, max(pcaobj$rotated[,y]) + 5) else
    c(min(pcaobj$rotated[,y]) - 1, max(pcaobj$rotated[,y]) + 1),
  lab = rownames(pcaobj$metadata),
  labSize = 3.0,
  labhjust = 1.5,
  labvjust = 0,
  boxedLabels = FALSE,
  selectLab = NULL,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  xlab = paste0(x, ', ', round(pcaobj$variance[x], digits = 2), '% variation'),
  xlabAngle = 0,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = paste0(y, ', ', round(pcaobj$variance[y], digits = 2), '% variation'),
  ylabAngle = 0,
  ylabhjust = 0.5,
  ylabvjust = 0.5,
  axisLabSize = 16,
  title = '',
  subtitle = '',
  caption = '',
  titleLabSize = 16,
  subtitleLabSize = 12,
  captionLabSize = 12,
  hline = NULL,
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  borderWidth = 0.8,
  borderColour = 'black',
  returnPlot = TRUE)
{
  
  labFun <- xidx <- yidx <- NULL
  
  # create a base theme that will later be modified
  th <- theme_bw(base_size = 24) +
    
    theme(
      legend.background = element_rect(),
      
      plot.title = element_text(angle = 0, size = titleLabSize,
                                face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = subtitleLabSize,
                                   face = 'plain', vjust = 1),
      plot.caption = element_text(angle = 0, size = captionLabSize,
                                  face = 'plain', vjust = 1),
      
      axis.text.x = element_text(angle = xlabAngle, size = axisLabSize,
                                 hjust = xlabhjust, vjust = xlabvjust),
      axis.text.y = element_text(angle = ylabAngle, size = axisLabSize,
                                 hjust = ylabhjust, vjust = ylabvjust),
      axis.title = element_text(size=axisLabSize),
      
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(size = legendLabSize),
      
      title = element_text(size = legendLabSize),
      legend.title = element_text(size = legendTitleSize))
  
  # set plot data labels (e.g. sample names)
  plotobj <- NULL
  plotobj$x <- pcaobj$rotated[,x]
  plotobj$y <- pcaobj$rotated[,y]
  if (!is.null(lab)) {
    plotobj$lab <- lab
  }
  plotobj <- as.data.frame(plotobj, stringsAsFactors = FALSE)
  
  # If user has supplied values in selectLab, convert labels to
  # NA and then re-set with those in selectLab
  if (!is.null(selectLab)) {
    if (is.null(lab)) {
      stop(paste0('You have specified lab as NULL ',
                  '- no labels can be selected!'))
    } else {
      names.new <- rep(NA, length(plotobj$lab))
      indices <- which(plotobj$lab %in% selectLab)
      names.new[indices] <- plotobj$lab[indices]
      plotobj$lab <- names.new
    }
  }
  
  # decide on how to colour the points, and specify the shape of these
  if (is.null(colby)) {
    if (!is.null(lab)) {
      plotobj$col <- lab
    } else {
      plotobj$col <- seq_len(length(pcaobj$yvars))
    }
  } else {
    plotobj$col <- pcaobj$metadata[,colby]
  }
  if (!is.null(shape)) {
    plotobj$shape <- pcaobj$metadata[,shape]
  }
  
  # create the plot object
  plot <- ggplot(plotobj, aes(x = x, y = y)) + th +
    
    guides(fill = guide_legend(),
           shape = guide_legend(),
           colour = guide_legend(override.aes = list(size = legendIconSize)))
  
  # if user specified a colour with 'singlecol', colour all points by this
  # otherwise, colour all points differently using ggplot engine.
  # shape of points remains independent of colouring
  if (is.null(singlecol)) {
    if (!is.null(shape)) {
      plot <- plot + geom_point(aes(color = col, shape = shape),
                                size = pointSize)
    } else {
      plot <- plot + geom_point(aes(color = col),
                                size = pointSize)
    }
  } else if (!is.null(singlecol)) {
    if (!is.null(shape)) {
      plot <- plot + geom_point(aes(color = singlecol, shape = shape),
                                size = pointSize)
    } else {
      plot <- plot + geom_point(aes(color = singlecol),
                                size = pointSize)
    }
  }
  
  # sort out custom colour pairing, and custom shapes
  if (!is.null(colkey)) {
    plot <- plot + scale_colour_discrete('') +
      scale_color_manual(values = colkey)
  }
  if (!is.null(shapekey)) {
    plot <- plot + scale_shape_manual(values = shapekey)
  }
  
  # plot loadings arrows?
  if (showLoadings) {
    # get top ntopLoadings to display
    xidx <- order(abs(pcaobj$loadings[,x]), decreasing = TRUE)
    yidx <- order(abs(pcaobj$loadings[,y]), decreasing = TRUE)
    vars <- unique(c(
      rownames(pcaobj$loadings)[xidx][seq_len(ntopLoadings)],
      rownames(pcaobj$loadings)[yidx][seq_len(ntopLoadings)]))
    
    # get scaling parameter to match between variable loadings and rotated loadings
    r <- min(
      (max(pcaobj$rotated[,x]) - min(pcaobj$rotated[,x]) /
         (max(pcaobj$loadings[,x]) - min(pcaobj$loadings[,x]))),
      (max(pcaobj$rotated[,y]) - min(pcaobj$rotated[,y]) /
         (max(pcaobj$loadings[,y]) - min(pcaobj$loadings[,y]))))
    
    plot <- plot +
      geom_segment(data = pcaobj$loadings[vars,],
                   aes(x = 0, y = 0,
                       xend = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                       yend = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor),
                   arrow = arrow(length = unit(1/2, 'picas'), ends = 'last'), 
                   color = colLoadingsArrows,
                   size = widthLoadingsArrows,
                   alpha = alphaLoadingsArrow,
                   show.legend = NA)
    
    if (showLoadingsNames) {
      if (drawConnectorsLoadings) {
        if (boxedLoadingsNames) {
          plot <- plot + coord_equal() +
            geom_label_repel(data = pcaobj$loadings[vars,], 
                             aes(label = vars,
                                 x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                                 y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor,
                                 hjust = 0),
                             color = colLoadingsNames,
                             size = sizeLoadingsNames,
                             fill = fillBoxedLoadings,
                             segment.color = colConnectorsLoadings,
                             segment.size = widthConnectorsLoadings)
        } else {
          plot <- plot + coord_equal() +
            geom_text_repel(data = pcaobj$loadings[vars,], 
                            aes(label = vars,
                                x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                                y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor,
                                hjust = 0),
                            color = colLoadingsNames,
                            size = sizeLoadingsNames,
                            segment.color = colConnectorsLoadings,
                            segment.size = widthConnectorsLoadings)
        }
      } else {
        if (boxedLoadingsNames) {
          plot <- plot + coord_equal() +
            geom_label(data = pcaobj$loadings[vars,], 
                       aes(label = vars,
                           x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                           y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor,
                           hjust = 0),
                       color = colLoadingsNames,
                       size = sizeLoadingsNames,
                       fill = NA)
        } else {
          plot <- plot + coord_equal() +
            geom_text(data = pcaobj$loadings[vars,], 
                      aes(label = vars,
                          x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                          y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor,
                          hjust = 0),
                      color = colLoadingsNames,
                      size = sizeLoadingsNames,
                      check_overlap = TRUE)
        }
      }
    }
  }
  
  # add elements to the plot for xy labeling and axis limits
  plot <- plot + xlab(xlab) + ylab(ylab)
  if (!is.null(xlim)) {
    plot <- plot + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    plot <- plot + ylim(ylim[1], ylim[2])
  }
  
  # add elements to the plot for title, subtitle, caption, and legend titles
  plot <- plot + labs(title = title, 
                      subtitle = subtitle, caption = caption,
                      fill = '', colour = colLegendTitle, shape = shapeLegendTitle)
  
  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline,
                              linetype = vlineType,
                              colour = vlineCol,
                              size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = hline,
                              linetype = hlineType,
                              colour = hlineCol,
                              size = hlineWidth)
  }
  
  # border around plot
  plot <- plot +
    theme(panel.border = element_rect(
      colour = borderColour,
      fill = NA,
      size = borderWidth))
  
  # gridlines
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }
  
  # labeling
  if (boxedLabels) {
    if (drawConnectors) {
      labFun <- function(...) geom_label_repel(...)
    } else {
      labFun <- function(...) geom_label(...)
    }
  } else {
    if (drawConnectors) {
      labFun <- function(...) geom_text_repel(...)
    } else {
      labFun <- function(...) geom_text(...)
    }
  }
  
  if (!is.null(lab)) {
    if (drawConnectors && is.null(selectLab)) {
      plot <- plot + labFun(
        data = plotobj,
        aes(label = lab),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        hjust = labhjust,
        vjust = labvjust)
    } else if (drawConnectors && !is.null(selectLab)) {
      plot <- plot + labFun(
        data=subset(plotobj,
                    !is.na(plotobj[,'lab'])),
        aes(label = lab),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        hjust = labhjust,
        vjust = labvjust)
    } else if (!drawConnectors && !is.null(selectLab)) {
      if (boxedLabels) {
        plot <- plot + labFun(
          data=subset(plotobj,
                      !is.na(plotobj[,'lab'])),
          aes(label = lab),
          size = labSize,
          hjust = labhjust,
          vjust = labvjust)
      } else {
        plot <- plot + labFun(
          data=subset(plotobj,
                      !is.na(plotobj[,'lab'])),
          aes(label = lab),
          size = labSize,
          check_overlap = TRUE,
          hjust = labhjust,
          vjust = labvjust)
      }
    } else if (!drawConnectors && is.null(selectLab)) {
      if (boxedLabels) {
        plot <- plot + labFun(
          data = plotobj,
          aes(label = lab),
          size = labSize,
          check_overlap = TRUE,
          hjust = labhjust,
          vjust = labvjust)
      } else {
        plot <- plot + labFun(
          data = plotobj,
          aes(label = lab),
          size = labSize,
          check_overlap = TRUE,
          hjust = labhjust,
          vjust = labvjust)
      }
    }
  }
  
  # encircle
  if (encircle) {
    if (encircleFill) {
      if (is.null(encircleLineCol)) {
        plot <- plot +
          ggalt::geom_encircle(
            aes(group = col,
                fill = col,
                colour = col),
            alpha = encircleAlpha,
            size = encircleLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      } else {
        plot <- plot +
          ggalt::geom_encircle(
            aes(group = col,
                fill = col),
            colour = encircleLineCol,
            alpha = encircleAlpha,
            size = encircleLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      }
    } else {
      if (is.null(encircleLineCol)) {
        plot <- plot +
          ggalt::geom_encircle(
            aes(group = col,
                colour = col),
            fill = NA,
            alpha = encircleAlpha,
            size = encircleLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      } else {
        plot <- plot +
          ggalt::geom_encircle(
            aes(group = col),
            colour = encircleLineCol,
            fill = NA,
            alpha = encircleAlpha,
            size = encircleLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      }
    }
    
    if (encircleFill) {
      if (is.null(encircleFillKey)) {
        if (!is.null(colkey)) {
          plot <- plot + scale_fill_manual(values = colkey)
        }
      } else {
        plot <- plot + scale_fill_manual(values = encircleFillKey)
      }
    }
  }
  
  # ellipse
  if (ellipse) {
    if (ellipseFill) {
      if (is.null(ellipseLineCol)) {
        plot <- plot +
          stat_ellipse(
            aes(group = col,
                fill = col,
                colour = col),
            geom = 'polygon',
            level = ellipseConf,
            alpha = ellipseAlpha,
            size = ellipseLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      } else {
        plot <- plot +
          stat_ellipse(
            aes(group = col,
                fill = col),
            colour = ellipseLineCol,
            geom = 'polygon',
            level = ellipseConf,
            alpha = ellipseAlpha,
            size = ellipseLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      }
    } else {
      if (is.null(ellipseLineCol)) {
        plot <- plot +
          stat_ellipse(
            aes(group = col,
                colour = col),
            fill = NA,
            geom = 'polygon',
            level = ellipseConf,
            alpha = ellipseAlpha,
            size = ellipseLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      } else {
        plot <- plot +
          stat_ellipse(
            aes(group = col),
            colour = ellipseLineCol,
            fill = NA,
            geom = 'polygon',
            level = ellipseConf,
            alpha = ellipseAlpha,
            size = ellipseLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      }
    }
    
    if (ellipseFill) {
      if (is.null(ellipseFillKey)) {
        if (!is.null(colkey)) {
          plot <- plot + scale_fill_manual(values = colkey)
        }
      } else {
        plot <- plot + scale_fill_manual(values = ellipseFillKey)
      }
    }
  }
  
  # return plot?
  if (returnPlot) {
    return(plot)
  } else if (!returnPlot) {
    plot
  }
}
