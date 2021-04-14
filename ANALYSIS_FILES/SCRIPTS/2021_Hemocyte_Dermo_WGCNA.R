# WGCNA Analysis of Hemocyte and Perkinsus transcriptome from 2020 In vitro experiment
# Erin Roberts

### LOAD PACKAGES ####
library(tidyverse)
library(limma)
library(WGCNA) # v WGCNA_1.68
options(stringsAsFactors = FALSE) # run every time
allowWGCNAThreads()
library(cluster)
library(anRichment)
library(anRichmentMethods)
library(plyr)
library(dplyr)
library(magicfor)
cor <- WGCNA::cor # run every time
library(UpSetR)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(VennDiagram)
library(ComplexHeatmap)
library(pheatmap)
library(gt)
library(paletteer)
library(rtracklayer)
# Using R version 3.6.1

# helpful WGCNA tutorials and FAQs
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

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
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_CV_uniq_XP_XM.Rdata")

# Load dataframes with Interproscan terms
load(file="/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/Perk_Interpro_GO_terms_XP.RData")
View(Perk_Interpro_GO_terms_XP)
load(file= "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/GO_universe_rna_found.RData")

## LOAD TRANSFORMED EXPRESSION DATA AND COLDATA ###

load(file='/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/2021_Hemocyte_Dermo_expression_rlog_matrices.RData')
load(file="/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/2021_Hemocyte_Dermo_expression_coldata.RData")

#####  DATA FORMATTING ####

hemo_dds_rlog_matrix <- assay(hemo_dds_rlog)
hemo_dds_rlog_matrix <- t(hemo_dds_rlog_matrix )
head(hemo_dds_rlog_matrix)
ncol(hemo_dds_rlog_matrix) 
colnames(hemo_dds_rlog_matrix)

perk_dds_rlog_matrix <- assay(perk_dds_rlog)
perk_dds_rlog_matrix <- t(perk_dds_rlog_matrix )
head(perk_dds_rlog_matrix)
ncol(perk_dds_rlog_matrix) 
colnames(perk_dds_rlog_matrix)


### PICK SOFT THRESHOLD ###

# Pick soft threshold
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# save for running TOM from cluster
#write.table(Dermo_Tolerant_dds_vst_matrix, file="/Volumes/My Passport for Mac/2021_Hemocyte_Dermo_WGCNA/Dermo_Tolerant_dds_vst_matrix.table")

# Hemocyte data 
# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
hemo_dds_rlog_matrix_sft <- pickSoftThreshold(hemo_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
save(hemo_dds_rlog_matrix_sft, file="/Volumes/My Passport for Mac/2021_Hemocyte_Dermo_WGCNA/hemo_dds_rlog_matrix_sft")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(hemo_dds_rlog_matrix_sft$fitIndices[,1], -sign(hemo_dds_rlog_matrix_sft$fitIndices[,3])*hemo_dds_rlog_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(hemo_dds_rlog_matrix_sft$fitIndices[,1], -sign(hemo_dds_rlog_matrix_sft$fitIndices[,3])*hemo_dds_rlog_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(hemo_dds_rlog_matrix_sft$fitIndices[,1], hemo_dds_rlog_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(hemo_dds_rlog_matrix_sft$fitIndices[,1], hemo_dds_rlog_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Selecting Softthreshold of 7 since this is lowest value past 0.9 we start to see flattening 

## Perkinsus data 
# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
perk_dds_rlog_matrix_sft <- pickSoftThreshold(perk_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
save(perk_dds_rlog_matrix_sft, file="/Volumes/My Passport for Mac/2021_Hemocyte_Dermo_WGCNA/perk_dds_rlog_matrix_sft")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(perk_dds_rlog_matrix_sft$fitIndices[,1], -sign(perk_dds_rlog_matrix_sft$fitIndices[,3])*perk_dds_rlog_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(perk_dds_rlog_matrix_sft$fitIndices[,1], -sign(perk_dds_rlog_matrix_sft$fitIndices[,3])*perk_dds_rlog_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(perk_dds_rlog_matrix_sft$fitIndices[,1], perk_dds_rlog_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(perk_dds_rlog_matrix_sft$fitIndices[,1], perk_dds_rlog_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Selecting Softthreshold of 7 since this is lowest value past 0.9 we start to see flattening 

#### ONE STEP NETWORK CONSTRUCTION, MODULE DETECTION, MODULE DENDROGRAM INSPECTION ####
hemo_full_net = blockwiseModules(hemo_dds_rlog_matrix, power = 7, # picked suitable power in the code above 
                               TOMType = "signed", # use signed TOM type
                               networkType= "signed hybrid", # use signed hybrid network type
                               corType = "bicor", # use suggested bicor
                               TminModuleSize = 30, # recommended default
                               reassignThreshold = 0, # recommended default
                               mergeCutHeight = 0.25, # recommended default
                               numericLabels = TRUE, # recommended default
                               pamRespectsDendro = FALSE,# recommended default
                               verbose = 3, 
                               maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer

# How many modules identified
table(hemo_full_net$colors) # 121 total modules
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
hemo_full_mergedColors = labels2colors(hemo_full_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(hemo_full_net$dendrograms[[1]], hemo_full_mergedColors[hemo_full_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
hemo_full_moduleLabels = hemo_full_net$colors
hemo_full_moduleColors = labels2colors(hemo_full_net$colors)
hemo_full_MEs = hemo_full_net$MEs
hemo_full_geneTree = hemo_full_net$dendrograms[[1]]
 #save network
save(hemo_full_net, hemo_full_mergedColors, hemo_full_moduleLabels, hemo_full_moduleColors, hemo_full_MEs, hemo_full_geneTree, 
     file = "/Volumes/My Passport for Mac/2021_Hemocyte_Dermo_WGCNA/hemo_full_network.RData")


## Repeat for Perkinsus expression

perk_full_net = blockwiseModules(perk_dds_rlog_matrix, power = 7, # picked suitable power in the code above 
                                 TOMType = "signed", # use signed TOM type
                                 networkType= "signed hybrid", # use signed hybrid network type
                                 corType = "bicor", # use suggested bicor
                                 TminModuleSize = 30, # recommended default
                                 reassignThreshold = 0, # recommended default
                                 mergeCutHeight = 0.25, # recommended default
                                 numericLabels = TRUE, # recommended default
                                 pamRespectsDendro = FALSE,# recommended default
                                 verbose = 3, 
                                 maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer

# How many modules identified
table(perk_full_net$colors) # 134 total modules
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
perk_full_mergedColors = labels2colors(perk_full_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(perk_full_net$dendrograms[[1]], perk_full_mergedColors[perk_full_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
perk_full_moduleLabels = perk_full_net$colors
perk_full_moduleColors = labels2colors(perk_full_net$colors)
perk_full_MEs = perk_full_net$MEs
perk_full_geneTree = perk_full_net$dendrograms[[1]]
#save network
save(perk_full_net, perk_full_mergedColors, perk_full_moduleLabels, perk_full_moduleColors, perk_full_MEs, perk_full_geneTree, 
     file = "/Volumes/My Passport for Mac/2021_Hemocyte_Dermo_WGCNA/perk_full_network.RData")

#### BINARIZE CATEGORICAL VARIABLES TO TEST MODULE ASSOCIATIONS WITH CHALLENGE ####
# out = binarizeCategoricalVariable(x,
# includePairwise = TRUE,
# includeLevelVsAll = FALSE);

hemo_coldata_collapse <- hemo_coldata %>% dplyr::select(condition)
perk_coldata_collapse <- perk_coldata %>% dplyr::select(condition)

all(row.names(hemo_coldata_collapse) == row.names(hemo_dds_rlog_matrix)) # TRUE
all(row.names(perk_coldata_collapse) == row.names(perk_dds_rlog_matrix) ) # TRUE

# binarize
hemo_coldata_collapse_binarize <- binarizeCategoricalColumns.pairwise(hemo_coldata_collapse)
row.names(hemo_coldata_collapse_binarize) <- row.names(hemo_coldata_collapse)
colnames(hemo_coldata_collapse_binarize)

perk_coldata_collapse_binarize <- binarizeCategoricalColumns.pairwise(perk_coldata_collapse)
row.names(perk_coldata_collapse_binarize) <- row.names(perk_coldata_collapse)
colnames(perk_coldata_collapse_binarize)

#### QUANTIFY MODULE ASSOCIATIONS WITH CHALLENGE ####

## Going to quantify for each set any significant associations with treatment 

# tutorial for this section: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf
# Define numbers of genes and samples
hemo_full_nGenes = ncol(hemo_dds_rlog_matrix)
hemo_full_nSamples = nrow(hemo_dds_rlog_matrix)

# Recalculate MEs with color labels
hemo_full_MEs0 = moduleEigengenes(hemo_dds_rlog_matrix, hemo_full_moduleColors)$eigengenes
hemo_full_MEs = orderMEs(hemo_full_MEs0)
hemo_full_moduleTraitCor = cor(hemo_full_MEs, hemo_coldata_collapse_binarize, use = "p");
hemo_full_moduleTraitPvalue = corPvalueStudent(hemo_full_moduleTraitCor, hemo_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
hemo_full_textMatrix = paste(signif(hemo_full_moduleTraitCor, 2), "\n(",
                                  signif(hemo_full_moduleTraitPvalue, 1), ")", sep = "");
dim(hemo_full_textMatrix) = dim(hemo_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = hemo_full_moduleTraitCor,
               xLabels = names(hemo_coldata_collapse_binarize),
               yLabels = names(hemo_full_MEs),
               ySymbols = names(hemo_full_MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = hemo_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
hemo_full_moduleTraitCor_df <- as.data.frame(hemo_full_moduleTraitCor) %>% dplyr::select(condition.Pmar.vs.control,condition.Pmar_GDC.vs.control,condition.Pmar_ZVAD.vs.control)
colnames(hemo_full_moduleTraitCor_df) <- c("condition.Pmar.vs.control.moduleTraitCor", "condition.Pmar_GDC.vs.control.moduleTraitCor", "condition.Pmar_ZVAD.vs.control.moduleTraitCor")
hemo_full_moduleTraitPvalue_df <- as.data.frame(hemo_full_moduleTraitPvalue) %>% dplyr::select(condition.Pmar.vs.control,condition.Pmar_GDC.vs.control,condition.Pmar_ZVAD.vs.control)
colnames(hemo_full_moduleTraitPvalue_df) <- c("condition.Pmar.vs.control.moduleTraitPvalue", "condition.Pmar_GDC.vs.control.moduleTraitPvalue", "condition.Pmar_ZVAD.vs.control.moduleTraitPvalue")

hemo_full_moduleTraitCor_Pval_df <- cbind(hemo_full_moduleTraitCor_df, hemo_full_moduleTraitPvalue_df) 

hemo_full_moduleTraitCor_Pval_df_Pmar_vs_control <- hemo_full_moduleTraitCor_Pval_df %>% dplyr::select(contains("condition.Pmar.vs.control"))
hemo_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_control <- hemo_full_moduleTraitCor_Pval_df %>% dplyr::select(contains("condition.Pmar_GDC.vs.control"))
hemo_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_control <- hemo_full_moduleTraitCor_Pval_df %>% dplyr::select(contains("condition.Pmar_ZVAD.vs.control"))

# Significantly correlated modules
hemo_full_moduleTraitCor_Pval_df_Pmar_vs_control_sig <- hemo_full_moduleTraitCor_Pval_df_Pmar_vs_control %>% filter(condition.Pmar.vs.control.moduleTraitPvalue <= 0.05)  %>% rownames_to_column(., "mod_names") 
nrow(hemo_full_moduleTraitCor_Pval_df_Pmar_vs_control_sig) #26 significant modules

hemo_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_control_sig <- hemo_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_control %>% filter(condition.Pmar_GDC.vs.control.moduleTraitPvalue<= 0.05)  %>% rownames_to_column(., "mod_names") 
nrow(hemo_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_control_sig) #41 significant modules, only 10 have positive correlations with GDC

hemo_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_control_sig <- hemo_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_control %>% filter(condition.Pmar_ZVAD.vs.control.moduleTraitPvalue<= 0.05)  %>% rownames_to_column(., "mod_names") 
nrow(hemo_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_control_sig) #32 significant modules

# compare modules
hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare <- hemo_full_moduleTraitCor_Pval_df_Pmar_vs_control_sig %>% 
  full_join(.,hemo_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_control_sig) %>% 
  full_join(.,hemo_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_control_sig) %>% dplyr::select(!contains("Pvalue"))
  
# find those shared between all 
hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_shared <- drop_na(hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare)
# 8 shared between all 
# interesting modules that have an opposite sign between treatments
   # MEyellow
   # MEdarkseagreen4
   # MEplum4

### Repeat for P_marinus counts 

# tutorial for this section: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf
# Define numbers of genes and samples
perk_full_nGenes = ncol(perk_dds_rlog_matrix)
perk_full_nSamples = nrow(perk_dds_rlog_matrix)

# Recalculate MEs with color labels
perk_full_MEs0 = moduleEigengenes(perk_dds_rlog_matrix, perk_full_moduleColors)$eigengenes
perk_full_MEs = orderMEs(perk_full_MEs0)
perk_full_moduleTraitCor = cor(perk_full_MEs, perk_coldata_collapse_binarize, use = "p");
perk_full_moduleTraitPvalue = corPvalueStudent(perk_full_moduleTraitCor, perk_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
perk_full_textMatrix = paste(signif(perk_full_moduleTraitCor, 2), "\n(",
                             signif(perk_full_moduleTraitPvalue, 1), ")", sep = "");
dim(perk_full_textMatrix) = dim(perk_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = perk_full_moduleTraitCor,
               xLabels = names(perk_coldata_collapse_binarize),
               yLabels = names(perk_full_MEs),
               ySymbols = names(perk_full_MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = perk_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
perk_full_moduleTraitCor_df <- as.data.frame(perk_full_moduleTraitCor) %>% dplyr::select(condition.Pmar_GDC.vs.Pmar,condition.Pmar_ZVAD.vs.Pmar)
colnames(perk_full_moduleTraitCor_df) <- c( "condition.Pmar_GDC.vs.Pmar.moduleTraitCor", "condition.Pmar_ZVAD.vs.Pmar.moduleTraitCor")
perk_full_moduleTraitPvalue_df <- as.data.frame(perk_full_moduleTraitPvalue) %>% dplyr::select(condition.Pmar_GDC.vs.Pmar,condition.Pmar_ZVAD.vs.Pmar)
colnames(perk_full_moduleTraitPvalue_df) <- c( "condition.Pmar_GDC.vs.Pmar.moduleTraitPvalue", "condition.Pmar_ZVAD.vs.Pmar.moduleTraitPvalue")

perk_full_moduleTraitCor_Pval_df <- cbind(perk_full_moduleTraitCor_df, perk_full_moduleTraitPvalue_df) 

perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar <- perk_full_moduleTraitCor_Pval_df %>% dplyr::select(contains("condition.Pmar_GDC.vs.Pmar"))
perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar <- perk_full_moduleTraitCor_Pval_df %>% dplyr::select(contains("condition.Pmar_ZVAD.vs.Pmar"))

# Significantly correlated modules

perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar_sig <- perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar %>% filter(condition.Pmar_GDC.vs.Pmar.moduleTraitPvalue<= 0.05)  %>% rownames_to_column(., "mod_names") 
nrow(perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar_sig) #19 significant modules

perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar_sig <- perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar %>% filter(condition.Pmar_ZVAD.vs.Pmar.moduleTraitPvalue<= 0.05)  %>% rownames_to_column(., "mod_names") 
nrow(perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar_sig) #24 significant modules

# compare modules
perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare <- 
  full_join(perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar_sig,perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar_sig) %>% dplyr::select(!contains("Pvalue"))

# find those shared between all 
perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare_shared <- drop_na(perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare)
# 5 shared between all, all move in the same direction 

#### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES ####

# Perform first for hemocyte samples
hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare

hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_list <- hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare$mod_names
hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_list_rm <- str_remove(hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_list, "ME")

# Use function to lookup all apop names for each significant module
matrix= hemo_dds_rlog_matrix
moduleColors=hemo_full_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_list_rm) <- c("red"            ,"yellow"        ,  "skyblue3"       , "blueviolet"    ,  "darkorange2"    , "orangered4"      ,"orangered3"     , "purple"        ,  "lightcyan"     , 
                                                                      "deeppink"       ,"salmon2"       ,  "darkolivegreen4", "darkviolet"    ,  "bisque4"        , "darkseagreen4"   ,"navajowhite"    , "navajowhite1"  ,  "antiquewhite2" , 
                                                                      "green4"         ,"mediumpurple1" ,  "honeydew"       , "blue"          ,  "plum4"          , "lightsteelblue"  ,"lightpink2"     , "mistyrose"     ,  "magenta3"      , 
                                                                      "chocolate4"     ,"thistle"       ,  "antiquewhite1"  , "navajowhite2"  ,  "royalblue"      , "skyblue1"        ,"skyblue2"       , "lavenderblush1",  "darkgreen"     , 
                                                                      "lightslateblue" ,"green"         ,  "salmon4"        , "plum2"         ,  "lightcoral"     , "yellow3"         ,"tan4"           , "lightpink3"    ,  "mediumpurple3" , 
                                                                      "plum"           ,"cyan"          ,  "darkred"        , "paleturquoise" ,  "black"          , "darkorange"      ,"yellowgreen"    , "skyblue4"      ,  "plum1"         , 
                                                                      "blue4"          ,"coral1"        ,  "darkslateblue"  , "darkolivegreen",  "palevioletred2" , "lightsteelblue1")

hemo_full_module_apop <- lapply(hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_list_rm,  lookup_mod_apop)
hemo_full_module_apop_df <- do.call(rbind,hemo_full_module_apop)
hemo_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(hemo_full_module_apop_df))
hemo_full_module_apop_df$mod_names <- gsub("^","ME",hemo_full_module_apop_df$mod_names)
# add module significance
hemo_full_module_apop_df <- left_join(hemo_full_module_apop_df,hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare)
hemo_full_module_apop_df$exp <- "hemo"

# which modules have apoptosis gene names
hemo_full_module_apop_df %>% distinct(mod_names) %>% View() # 57

# which modules have more than five apoptosis gene names
hemo_full_module_apop_df_5_greater <- hemo_full_module_apop_df %>% group_by(mod_names) %>% filter(n() >= 5) %>% distinct(mod_names)
hemo_full_module_apop_df_5_greater_annot <- hemo_full_module_apop_df[(hemo_full_module_apop_df$mod_names %in% hemo_full_module_apop_df_5_greater$mod_names),]

hemo_full_module_apop_df_5_greater_control_Pmar <- hemo_full_module_apop_df_5_greater_annot %>% filter(!is.na(condition.Pmar.vs.control.moduleTraitCor)) %>% 
  dplyr::select(mod_names) %>% distinct()
hemo_full_module_apop_df_5_greater_control_Pmar_GDC <- hemo_full_module_apop_df_5_greater_annot %>% filter(!is.na(condition.Pmar_GDC.vs.control.moduleTraitCor)) %>% 
  dplyr::select(mod_names) %>% distinct()
hemo_full_module_apop_df_5_greater_control_Pmar_ZVAD <- hemo_full_module_apop_df_5_greater_annot %>% filter(!is.na(condition.Pmar_ZVAD.vs.control.moduleTraitCor)) %>% 
  dplyr::select(mod_names) %>% distinct()

### can't repeat for P. marinus because I haven't systematically annotated their apoptosis genes 

### Gene relationship to trait and important modules: Gene Significance and Module Membership ####
# Interpretation notes from Langfelder and Horvath 2008
  # The sign of module membership encodes whether the gene has a positive or a negative relationship with the blue module eigengene. 
  # The module membership measure can be defined for all input genes (irrespective of their original module membership). 
  # It turns out that the module membership measure is highly related to the intramodular connectivity kIM. Highly connected 
  # intramodular hub genes tend to have high module membership values to the respective module

  #The higher the mean gene significance in a module, the more significantly related the module is to the clinical trait of interest. B.

## Hemocyte MM for all treatments
# names (colors) of the modules
hemo_full_modNames = substring(names(hemo_full_MEs), 3)
hemo_full_geneModuleMembership = as.data.frame(cor(hemo_dds_rlog_matrix, hemo_full_MEs, use = "p"))

# calculate the module membership. this will be the same across every treatment, the treatment really only matters for gene trait significance
hemo_full_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(hemo_full_geneModuleMembership), hemo_full_nSamples))

names(hemo_full_geneModuleMembership) = paste("MM", hemo_full_modNames, sep="")
names(hemo_full_MMPvalue) = paste("p.MM", hemo_full_modNames, sep="")

## calculate gene trait significance for each treatment
# Control vs. Pmar
# Define variable injected 
hemo_full_control_Pmar = as.data.frame(hemo_coldata_collapse_binarize$condition.Pmar.vs.control);
names(hemo_full_control_Pmar) = "control_Pmar"
hemo_full_geneTraitSignificance_control_Pmar = as.data.frame(cor(hemo_dds_rlog_matrix,hemo_full_control_Pmar, use = "p"))
hemo_full_GSPvalue_control_Pmar = as.data.frame(corPvalueStudent(as.matrix(hemo_full_geneTraitSignificance_control_Pmar), hemo_full_nSamples))

names(hemo_full_geneTraitSignificance_control_Pmar) = paste("GS.", names(hemo_full_control_Pmar), sep="")
names(hemo_full_GSPvalue_control_Pmar) = paste("p.GS.", names(hemo_full_control_Pmar), sep="")

# Control vs Pmar_GDC
hemo_full_control_Pmar_GDC = as.data.frame(hemo_coldata_collapse_binarize$condition.Pmar_GDC.vs.control);
names(hemo_full_control_Pmar_GDC) = "control_Pmar_GDC"
hemo_full_geneTraitSignificance_control_Pmar_GDC = as.data.frame(cor(hemo_dds_rlog_matrix,hemo_full_control_Pmar_GDC, use = "p"))
hemo_full_GSPvalue_control_Pmar_GDC = as.data.frame(corPvalueStudent(as.matrix(hemo_full_geneTraitSignificance_control_Pmar_GDC), hemo_full_nSamples))

names(hemo_full_geneTraitSignificance_control_Pmar_GDC) = paste("GS.", names(hemo_full_control_Pmar_GDC), sep="")
names(hemo_full_GSPvalue_control_Pmar_GDC) = paste("p.GS.", names(hemo_full_control_Pmar_GDC), sep="")

# Control vs Pmar_ZVAD
hemo_full_control_Pmar_ZVAD = as.data.frame(hemo_coldata_collapse_binarize$condition.Pmar_ZVAD.vs.control);
names(hemo_full_control_Pmar_ZVAD) = "control_Pmar_ZVAD"
hemo_full_geneTraitSignificance_control_Pmar_ZVAD = as.data.frame(cor(hemo_dds_rlog_matrix,hemo_full_control_Pmar_ZVAD, use = "p"))
hemo_full_GSPvalue_control_Pmar_ZVAD = as.data.frame(corPvalueStudent(as.matrix(hemo_full_geneTraitSignificance_control_Pmar_ZVAD), hemo_full_nSamples))

names(hemo_full_geneTraitSignificance_control_Pmar_ZVAD) = paste("GS.", names(hemo_full_control_Pmar_ZVAD), sep="")
names(hemo_full_GSPvalue_control_Pmar_ZVAD) = paste("p.GS.", names(hemo_full_control_Pmar_ZVAD), sep="")

## Perk MM for all treatments
# names (colors) of the modules
perk_full_modNames = substring(names(perk_full_MEs), 3)
perk_full_geneModuleMembership = as.data.frame(cor(perk_dds_rlog_matrix, perk_full_MEs, use = "p"))

# calculate the module membership. this will be the same across every treatment, the treatment really only matters for gene trait significance
perk_full_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(perk_full_geneModuleMembership), perk_full_nSamples))

names(perk_full_geneModuleMembership) = paste("MM", perk_full_modNames, sep="")
names(perk_full_MMPvalue) = paste("p.MM", perk_full_modNames, sep="")

## calculate gene trait significance for each treatment
# Pmar vs Pmar_GDC
perk_full_Pmar_Pmar_GDC = as.data.frame(perk_coldata_collapse_binarize$condition.Pmar_GDC.vs.Pmar);
names(perk_full_Pmar_Pmar_GDC) = "Pmar_Pmar_GDC"
perk_full_geneTraitSignificance_Pmar_Pmar_GDC = as.data.frame(cor(perk_dds_rlog_matrix,perk_full_Pmar_Pmar_GDC, use = "p"))
perk_full_GSPvalue_Pmar_Pmar_GDC = as.data.frame(corPvalueStudent(as.matrix(perk_full_geneTraitSignificance_Pmar_Pmar_GDC), perk_full_nSamples))

names(perk_full_geneTraitSignificance_Pmar_Pmar_GDC) = paste("GS.", names(perk_full_Pmar_Pmar_GDC), sep="")
names(perk_full_GSPvalue_Pmar_Pmar_GDC) = paste("p.GS.", names(perk_full_Pmar_Pmar_GDC), sep="")

# Pmar vs Pmar_ZVAD
perk_full_Pmar_Pmar_ZVAD = as.data.frame(perk_coldata_collapse_binarize$condition.Pmar_ZVAD.vs.Pmar);
names(perk_full_Pmar_Pmar_ZVAD) = "Pmar_Pmar_ZVAD"
perk_full_geneTraitSignificance_Pmar_Pmar_ZVAD = as.data.frame(cor(perk_dds_rlog_matrix,perk_full_Pmar_Pmar_ZVAD, use = "p"))
perk_full_GSPvalue_Pmar_Pmar_ZVAD = as.data.frame(corPvalueStudent(as.matrix(perk_full_geneTraitSignificance_Pmar_Pmar_ZVAD), perk_full_nSamples))

names(perk_full_geneTraitSignificance_Pmar_Pmar_ZVAD) = paste("GS.", names(perk_full_Pmar_Pmar_ZVAD), sep="")
names(perk_full_GSPvalue_Pmar_Pmar_ZVAD) = paste("p.GS.", names(perk_full_Pmar_Pmar_ZVAD), sep="")

### Intramodular analysis: identifying genes with high GS and MM

## Hemocyte intramodular analysis - perform for each treatment 
# Hemocyte Control vs P. mar
# only view the modules that are interesting for apoptosis (>=5 apoptosis transcripts found)
hemo_full_module_apop_df_5_greater_control_Pmar_list <- as.character(unlist(hemo_full_module_apop_df_5_greater_control_Pmar))
hemo_full_module_apop_df_5_greater_control_Pmar_list <- str_remove(hemo_full_module_apop_df_5_greater_control_Pmar_list, "ME")

GS_MM_plot <- function(list) {
hemo_full_module = list 
hemo_full_column = match(hemo_full_module, hemo_full_modNames)
hemo_full_moduleGenes = hemo_full_moduleColors==hemo_full_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(hemo_full_geneModuleMembership [hemo_full_moduleGenes, hemo_full_column]),
                   abs(hemo_full_geneTraitSignificance_control_Pmar[hemo_full_moduleGenes, 1]),
                   xlab = paste("Module Membership in", hemo_full_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = hemo_full_module)
quartz.save(paste("./FIGURES/control_Pmar",list, sep ="_"), type = "png", device = dev.cur(), dpi = 100)
}
lapply(hemo_full_module_apop_df_5_greater_control_Pmar_list,  GS_MM_plot)

## modules with a good correlation between GS and module membership (looking at those with higher than 0.4 correlation)
# high correlation of GS and MM illustrates that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait
  # red = cor = 0.44
  # lightcyan = 0.6
  # darkseagreen4 = 0.4 
  # antiquewhite2 = 0.44
  # blue = 0.4

# Hemocyte Control vs P. mar_GDC
# only view the modules that are interesting for apoptosis (>=5 apoptosis transcripts found)
hemo_full_module_apop_df_5_greater_control_Pmar_GDC_list <- as.character(unlist(hemo_full_module_apop_df_5_greater_control_Pmar_GDC))
hemo_full_module_apop_df_5_greater_control_Pmar_GDC_list <- str_remove(hemo_full_module_apop_df_5_greater_control_Pmar_GDC_list, "ME")

GS_MM_plot <- function(list) {
  hemo_full_module = list 
  hemo_full_column = match(hemo_full_module, hemo_full_modNames)
  hemo_full_moduleGenes = hemo_full_moduleColors==hemo_full_module
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(hemo_full_geneModuleMembership [hemo_full_moduleGenes, hemo_full_column]),
                     abs(hemo_full_geneTraitSignificance_control_Pmar_GDC[hemo_full_moduleGenes, 1]),
                     xlab = paste("Module Membership in", hemo_full_module, "module"),
                     ylab = "Gene significance for challenge",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = hemo_full_module)
  quartz.save(paste("./FIGURES/control_Pmar_GDC",list, sep ="_"), type = "png", device = dev.cur(), dpi = 100)
}
lapply(hemo_full_module_apop_df_5_greater_control_Pmar_GDC_list,  GS_MM_plot)

## modules with a good correlation between GS and module membership (looking at those with higher than 0.4 correlation)
# high correlation of GS and MM illustrates that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait
  # yellow = 0.71
  # lightcyan = 0.55
  # antique white2 = 0.75
  # navajowhite2 = 0.59
  # darkgreen = 0.4
  # lighpink3 = 0.62
  # mediumpurple3 = 0.53
  # plum = 0.73
  # cyan = 0.46
  # darkred = 0.64 
  # paleturquoise = 0.6
  # black = 0.61
  # darkorange = 0.62

# Hemocyte Control vs P. mar_ZVAD
# only view the modules that are interesting for apoptosis (>=5 apoptosis transcripts found)
hemo_full_module_apop_df_5_greater_control_Pmar_ZVAD_list <- as.character(unlist(hemo_full_module_apop_df_5_greater_control_Pmar_ZVAD))
hemo_full_module_apop_df_5_greater_control_Pmar_ZVAD_list <- str_remove(hemo_full_module_apop_df_5_greater_control_Pmar_ZVAD_list, "ME")

GS_MM_plot <- function(list) {
  hemo_full_module = list 
  hemo_full_column = match(hemo_full_module, hemo_full_modNames)
  hemo_full_moduleGenes = hemo_full_moduleColors==hemo_full_module
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(hemo_full_geneModuleMembership [hemo_full_moduleGenes, hemo_full_column]),
                     abs(hemo_full_geneTraitSignificance_control_Pmar_ZVAD[hemo_full_moduleGenes, 1]),
                     xlab = paste("Module Membership in", hemo_full_module, "module"),
                     ylab = "Gene significance for challenge",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = hemo_full_module)
  quartz.save(paste("./FIGURES/control_Pmar_ZVAD",list, sep ="_"), type = "png", device = dev.cur(), dpi = 100)
}
lapply(hemo_full_module_apop_df_5_greater_control_Pmar_ZVAD_list,  GS_MM_plot)

## modules with a good correlation between GS and module membership (looking at those with higher than 0.4 correlation)
# high correlation of GS and MM illustrates that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait
  # yellow = 0.48
  # orangered4 = 0.52
  # lightcyan = 0.56
  # plum = 0.61
  # darkslateblue = 0.57

## Parasite intramodular analysis - perform for each treatment 
# P. mar Control vs P. mar_GDC
# view the modules significant with GDC treatment 
perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar_sig_list <-  as.character(unlist(perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar_sig$mod_names))
perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar_sig_list <- str_remove(perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar_sig_list, "ME")

GS_MM_plot_perk <- function(list) {
  perk_full_module = list 
  perk_full_column = match(perk_full_module, perk_full_modNames)
  perk_full_moduleGenes = perk_full_moduleColors==perk_full_module
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(perk_full_geneModuleMembership [perk_full_moduleGenes, perk_full_column]),
                     abs(perk_full_geneTraitSignificance_Pmar_Pmar_GDC[perk_full_moduleGenes, 1]),
                     xlab = paste("Module Membership in", perk_full_module, "module"),
                     ylab = "Gene significance for challenge",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = perk_full_module)
  quartz.save(paste("./FIGURES/perk_Pmar_Pmar_GDC",list, sep ="_"), type = "png", device = dev.cur(), dpi = 100)
}
lapply(perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar_sig_list,  GS_MM_plot_perk)

## modules with a good correlation between GS and module membership (looking at those with higher than 0.4 correlation)
# high correlation of GS and MM illustrates that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait
  #darkseagreen2 = 0.67
  #orange = 0.52
  # indianred3 = 0.48
  # yellow4 = 0.66
  #steelblue = 0.42
  #darkorange2 = 0.57
  #darkturqoise = 0.49
  # indianred4 = 0.47
  #blue4 = 0.78
  #lightblue4 = 0.78
# P. mar Control vs P. mar_GDC
# view the modules significant with GDC treatment 
perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar_sig_list <-  as.character(unlist(perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar_sig$mod_names))
perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar_sig_list <- str_remove(perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar_sig_list, "ME")

GS_MM_plot_perk <- function(list) {
  perk_full_module = list 
  perk_full_column = match(perk_full_module, perk_full_modNames)
  perk_full_moduleGenes = perk_full_moduleColors==perk_full_module
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(perk_full_geneModuleMembership [perk_full_moduleGenes, perk_full_column]),
                     abs(perk_full_geneTraitSignificance_Pmar_Pmar_ZVAD[perk_full_moduleGenes, 1]),
                     xlab = paste("Module Membership in", perk_full_module, "module"),
                     ylab = "Gene significance for challenge",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = perk_full_module)
  quartz.save(paste("./FIGURES/perk_Pmar_Pmar_ZVAD",list, sep ="_"), type = "png", device = dev.cur(), dpi = 100)
}
lapply(perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar_sig_list,  GS_MM_plot_perk)

## modules with a good correlation between GS and module membership (looking at those with higher than 0.4 correlation)
# high correlation of GS and MM illustrates that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait
  #darkseagreen2 = 0.73
  # darkseagreen4 = 0.42
  #lightpink3 = 0.57
  #pink3 = 0.78
  #navajowhite2 = 0.7
  # darkturqiouse = 0.61
  #blue4 = 0.53

#### IDENTIFY HUB GENES IN EACH SIG MODULE ####

# Hemocytes control vs Pmar
hemo_full_colorh_control_Pmar = hemo_full_module_apop_df_5_greater_control_Pmar

perk_full_Module_hub_genes_control_Pmar <- chooseTopHubInEachModule(
  hemo_dds_rlog_matrix, 
  hemo_full_colorh_control_Pmar, 
  power = 7,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
class(hemo_full_Module_hub_genes_control_Pmar)
hemo_full_Module_hub_genes_control_Pmar_df <- as.data.frame(hemo_full_Module_hub_genes_control_Pmar) %>% rownames_to_column(., var = "mod_names")
colnames(hemo_full_Module_hub_genes_control_Pmar_df)[2] <- "ID"
hemo_full_Module_hub_genes_control_Pmar_all_control_Pmar <- merge(hemo_full_Module_hub_genes_control_Pmar_df, C_vir_rtracklayer, by = "ID")
hemo_full_Module_hub_genes_control_Pmar_all_control_Pmar$product
  # [1] "titin homolog, transcript variant X8"                        "uncharacterized LOC111125453, transcript variant X2"        
  # [3] "protocadherin-9-like, transcript variant X8"                 "uncharacterized LOC111122522, transcript variant X2"        
  # [5] "ribosome-binding protein 1-like, transcript variant X18"     "retina and anterior neural fold homeobox protein 2-like"    
  # [7] "uncharacterized LOC111128504, transcript variant X2"         "neural cell adhesion molecule 1-like, transcript variant X3"
  # [9] "uncharacterized LOC111119964, transcript variant X3" 

# Hemocytes control vs Pmar_GDC
hemo_full_colorh_control_Pmar_GDC = hemo_full_module_apop_df_5_greater_control_Pmar_GDC_list

hemo_full_Module_hub_genes_control_Pmar_GDC <- chooseTopHubInEachModule(
  hemo_dds_rlog_matrix, 
  hemo_full_colorh_control_Pmar_GDC, 
  power = 7,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)

hemo_full_Module_hub_genes_control_Pmar_GDC_df <- as.data.frame(hemo_full_Module_hub_genes_control_Pmar_GDC) %>% rownames_to_column(., var = "mod_names")
colnames(hemo_full_Module_hub_genes_control_Pmar_GDC_df)[2] <- "ID"
hemo_full_Module_hub_genes_control_Pmar_GDC_all_control_Pmar_GDC <- merge(hemo_full_Module_hub_genes_control_Pmar_GDC_df, C_vir_rtracklayer, by = "ID")
hemo_full_Module_hub_genes_control_Pmar_GDC_all_control_Pmar_GDC$product

hemo_full_Module_hub_genes_control_Pmar_GDC_all_control_Pmar_GDC$ID %in% BIR_XP_gff_CV_uniq_XP_XM$ID # none are IAP


# Hemocytes control vs Pmar_ZVAD
hemo_full_colorh_control_Pmar_ZVAD = hemo_full_module_apop_df_5_greater_control_Pmar_ZVAD_list

hemo_full_Module_hub_genes_control_Pmar_ZVAD <- chooseTopHubInEachModule(
  hemo_dds_rlog_matrix, 
  hemo_full_colorh_control_Pmar_ZVAD, 
  power = 7,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)

hemo_full_Module_hub_genes_control_Pmar_ZVAD_df <- as.data.frame(hemo_full_Module_hub_genes_control_Pmar_ZVAD) %>% rownames_to_column(., var = "mod_names")
colnames(hemo_full_Module_hub_genes_control_Pmar_ZVAD_df)[2] <- "ID"
hemo_full_Module_hub_genes_control_Pmar_ZVAD_all_control_Pmar_ZVAD <- merge(hemo_full_Module_hub_genes_control_Pmar_ZVAD_df, C_vir_rtracklayer, by = "ID")
hemo_full_Module_hub_genes_control_Pmar_ZVAD_all_control_Pmar_ZVAD$product

hemo_full_Module_hub_genes_control_Pmar_ZVAD_all_control_Pmar_ZVAD$ID %in% BIR_XP_gff_CV_uniq_XP_XM$ID # none are IAP

# Perkinsus hub genes 
# GDC 
perk_full_colorh_control_Pmar_GDC = perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar_sig_list

perk_full_Module_hub_genes_control_Pmar_GDC <- chooseTopHubInEachModule(
  perk_dds_rlog_matrix, 
  perk_full_colorh_control_Pmar_GDC, 
  power = 7,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)

perk_full_Module_hub_genes_control_Pmar_GDC_df <- as.data.frame(perk_full_Module_hub_genes_control_Pmar_GDC) %>% rownames_to_column(., var = "mod_names")
colnames(perk_full_Module_hub_genes_control_Pmar_GDC_df)[2] <- "Name"
perk_full_Module_hub_genes_control_Pmar_GDC_all_control_Pmar_GDC <- merge(perk_full_Module_hub_genes_control_Pmar_GDC_df, Perkinsus_rtracklayer, by = "Name")
perk_full_Module_hub_genes_control_Pmar_GDC_all_control_Pmar_GDC$product

# ZVAD
perk_full_colorh_control_Pmar_ZVAD = perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar_sig_list

perk_full_Module_hub_genes_control_Pmar_ZVAD <- chooseTopHubInEachModule(
  perk_dds_rlog_matrix, 
  perk_full_colorh_control_Pmar_ZVAD, 
  power = 7,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)

perk_full_Module_hub_genes_control_Pmar_ZVAD_df <- as.data.frame(perk_full_Module_hub_genes_control_Pmar_ZVAD) %>% rownames_to_column(., var = "mod_names")
colnames(perk_full_Module_hub_genes_control_Pmar_ZVAD_df)[2] <- "Name"
perk_full_Module_hub_genes_control_Pmar_ZVAD_all_control_Pmar_ZVAD <- merge(perk_full_Module_hub_genes_control_Pmar_ZVAD_df, Perkinsus_rtracklayer, by = "Name")
perk_full_Module_hub_genes_control_Pmar_ZVAD_all_control_Pmar_ZVAD$product

### EXPORT DATAFRAMES FOR ANALYSIS ####

# DATA TO INTERGRATE FOR VISUALIZATION IN EXCEL
  # GOALS: examine the pathways that are in modules and significantly positively correlated with challenge
  # in the paper, the WGCNA data just provides further evidence of pathways that are working together, and is complementary to the
    # DESeq2 paper. Additionally, identification of modules with a high correlation between gene significance and module membership shows
    # where genes are highly significantly associated with challenge and are also highly important members of the module.. these are modules
    # with potentially interesting candidates that are driving the response to challenge
    # finally identified hub genes may also be interesting candidates for study as they are highly connected and may be most important
  
  # HEMOCYTES 
      # 1. Apoptosis Product information for significant modules with apoptosis members 
      # 2. module significance with challenge (GS) and module membership (MM)
      # 3. Correlation between GS and MM for these modules
      # 4. Hub genes in these modules 
      # 5. Module membership for each gene in each module - correlation between the expression values and the module membership
            # - Highly connected intramodular hub genes tend to have high module membership values to the respective module. 

#### Calculate module membership value for each gene in each module ####
hemo_datKME=signedKME(hemo_dds_rlog_matrix, hemo_full_MEs, outputColumnName="MM.")
head(hemo_datKME)

perk_datKME=signedKME(perk_dds_rlog_matrix, perk_full_MEs, outputColumnName="MM.")
head(perk_datKME)


#### Find genes with high gene significance and high intramodular connectivity in interesting modules ####
# Remember that module membership is highly related to intramodular connectivity

# Which modules should be looked at?
  # For hemocytes: modules that are significant for treatment, have >5 apoptosis transcripts, and have high correlation between GS and MM
  # For hemocytes: same as above, except no filtering of modules for apoptosis related transcripts
# filter genes for each module of interest that have correlation between 

# annotate intramodular hub genes in each module
# Use function to lookup all apop names for each significant module
matrix = hemo_dds_rlog_matrix
lookup =   C_vir_rtracklayer_transcripts
datKME =  hemo_datKME

lookup_annot_intramodular <- function(list) {
FilterGenes = abs(GS1)> .2 & abs(datKME[list])>.8
FilterGenes_annot <- data.frame("ID" = dimnames(data.frame(matrix))[[2]][FilterGenes]) %>% left_join(., lookup) %>% mutate(mod_names = list)
}

## Hemocytes first
# control_Pmar
  # modules of interest fitting criteria above
  # red = cor = 0.44
  # lightcyan = 0.6
  # darkseagreen4 = 0.4 
  # antiquewhite2 = 0.44
  # blue = 0.4

GS1 = hemo_full_geneTraitSignificance_control_Pmar
hemo_control_Pmar_intramodular <- c("MM.red", "MM.lightcyan", "MM.darkseagreen4", "MM.antiquewhite2", "MM.blue")
names(hemo_control_Pmar_intramodular) <- c("MM.red", "MM.lightcyan", "MM.darkseagreen4", "MM.antiquewhite2", "MM.blue")

FilterGenes_control_Pmar <- lapply(hemo_control_Pmar_intramodular, lookup_annot_intramodular)
FilterGenes_control_Pmar_df <- do.call(rbind,FilterGenes_control_Pmar) %>% mutate(group = "control_Pmar")

# control_Pmar_GDC 
  # modules of interest fitting criteria above  
    # yellow = 0.71
    # lightcyan = 0.55
    # antique white2 = 0.75
    # navajowhite2 = 0.59
    # darkgreen = 0.4
    # lightpink3 = 0.62
    # mediumpurple3 = 0.53
    # plum = 0.73
    # cyan = 0.46
    # darkred = 0.64 
    # paleturquoise = 0.6
    # black = 0.61
    # darkorange = 0.62

GS1 = hemo_full_geneTraitSignificance_control_Pmar_GDC
hemo_control_Pmar_GDC_intramodular <- c("MM.yellow","MM.lightcyan","MM.antiquewhite2","MM.navajowhite2","MM.darkgreen","MM.lightpink3",
                                        "MM.mediumpurple3" ,"MM.plum","MM.cyan","MM.darkred","MM.paleturquoise" ,"MM.black","MM.darkorange")
names(hemo_control_Pmar_GDC_intramodular) <- c("MM.yellow","MM.lightcyan","MM.antiquewhite2","MM.navajowhite2","MM.darkgreen","MM.lighypink3",
                                               "MM.mediumpurple3" ,"MM.plum","MM.cyan","MM.darkred","MM.paleturquoise" ,"MM.black","MM.darkorange")

FilterGenes_control_Pmar_GDC <- lapply(hemo_control_Pmar_GDC_intramodular, lookup_annot_intramodular)
FilterGenes_control_Pmar_GDC_df <- do.call(rbind,FilterGenes_control_Pmar_GDC) %>% mutate(group = "control_Pmar_GDC")

# control_Pmar_ZVAD 
  # modules of interest fitting criteria above
      # yellow = 0.48
      # orangered4 = 0.52
      # lightcyan = 0.56
      # plum = 0.61
      # darkslateblue = 0.57

GS1 = hemo_full_geneTraitSignificance_control_Pmar_ZVAD
hemo_control_Pmar_ZVAD_intramodular <- c("MM.yellow", "MM.orangered4", "MM.lightcyan", "MM.plum","MM.darkslateblue")
names(hemo_control_Pmar_ZVAD_intramodular) <- c("MM.yellow", "MM.orangered4", "MM.lightcyan", "MM.plum","MM.darkslateblue")

FilterGenes_control_Pmar_ZVAD <- lapply(hemo_control_Pmar_ZVAD_intramodular, lookup_annot_intramodular)
FilterGenes_control_Pmar_ZVAD_df <- do.call(rbind,FilterGenes_control_Pmar_ZVAD) %>% mutate(group = "control_Pmar_ZVAD")

## Combine hemocyte data frames
FilterGenes_comb <- rbind(FilterGenes_control_Pmar_df,
                          FilterGenes_control_Pmar_GDC_df,
                          FilterGenes_control_Pmar_ZVAD_df)


### Repeat finding of intramodular hub genes, but for P. marinus samples now 
# annotate intramodular hub genes in each module
# Use function to lookup all apop names for each significant module
matrix = perk_dds_rlog_matrix
lookup =   Perkinsus_rtracklayer
GS1 = hemo_full_geneTraitSignificance_control_Pmar
datKME =  hemo_datKME

lookup_annot_intramodular <- function(list) {
  FilterGenes = abs(GS1)> .2 & abs(datKME[list])>.8
  FilterGenes_annot <- data.frame("ID" = dimnames(data.frame(matrix))[[2]][FilterGenes]) %>% left_join(., lookup) %>% mutate(mod_names = list)
}


# control_Pmar
# modules of interest fitting criteria above
# red = cor = 0.44
# lightcyan = 0.6
# darkseagreen4 = 0.4 
# antiquewhite2 = 0.44
# blue = 0.4
hemo_control_Pmar_intramodular <- c("MM.red", "MM.lightcyan", "MM.darkseagreen4", "MM.antiquewhite2", "MM.blue")
names(hemo_control_Pmar_intramodular) <- c("MM.red", "MM.lightcyan", "MM.darkseagreen4", "MM.antiquewhite2", "MM.blue")

FilterGenes_control_Pmar <- lapply(hemo_control_Pmar_intramodular, lookup_annot_intramodular)
FilterGenes_control_Pmar_df <- do.call(rbind,FilterGenes_control_Pmar) %>% mutate(group = "control_Pmar")





#### EXPORT WGNCA MATRIX TO CALCULATE INTRAMODULAR CONNECTIVITY IN BLUEWAVES ####

#Export full Matrices as tables for each experiment so that I can import the matrices in bluewaves 
write.table(hemo_dds_rlog_matrix, file="/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/WGCNA/hemo_dds_rlog_matrix.table")
write.table(perk_dds_rlog_matrix, file="/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/WGCNA/perk_dds_rlog_matrix.table")

# Export moduleColors file for each
save(hemo_full_moduleColors, file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/WGCNA/hemo_full_moduleColors.RData")
save(perk_full_moduleColors, file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/WGCNA/perk_full_moduleColors.RData")

# Export MEs for each
save(hemo_full_MEs, file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/WGCNA/hemo_full_MEs.RData")
save(perk_full_MEs, file = "/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/WGCNA/perk_full_MEs.RData")

