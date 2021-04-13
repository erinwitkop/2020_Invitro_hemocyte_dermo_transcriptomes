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

  #The higher the mean gene significance in a module, the more significantly related the mod- ule is to the clinical trait of interest. B.

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
ggsave(plot = last_plot(), device = "png", 
       file=paste("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/WGCNA/FIGURES/control_Pmar",list,".png",sep=""))
}

lapply(hemo_full_module_apop_df_5_greater_control_Pmar_list,  GS_MM_plot)

#### CALCULATE INTRAMODULAR CONNECTIVITY ####

hemo_ADJ1=abs(cor(hemo_dds,use="p"))^6 # this command caused my RStudio to crash..
hemo_Alldegrees1=intramodularConnectivity(hemo_ADJ1, hemo_colorh1)
head(hemo_Alldegrees1)

## Plot gene significance and intramodular connectivity for each module

colorlevels=unique(colorh1)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (colorh1==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=colorh1[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}

#### IDENTIFY HUB GENES IN EACH SIG MODULE ####
hemo_full_colorh = c("darkslateblue", "turquoise",     "greenyellow",   "skyblue3" ,     "cyan"  ,        "red"  ,         "tan" )

hemo_full_Module_hub_genes <- chooseTopHubInEachModule(
  hemo_dds_rlog_matrix_common, 
  hemo_full_colorh, 
  power = 3,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
class(hemo_full_Module_hub_genes)
hemo_full_Module_hub_genes_df <- as.data.frame(hemo_full_Module_hub_genes)
colnames(hemo_full_Module_hub_genes_df)[1] <- "ID"
hemo_full_Module_hub_genes_apop <- merge(hemo_full_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(hemo_full_Module_hub_genes_apop) # 7, none involed in apoptosis

## Compare turquoise module from consensus network to the full network...doesn't make sense to do this since I haven't confirmed that they are preserved yet 

hemo_full_module_apop_df_turq <- hemo_full_module_apop_df %>% filter(mod_names == "MEturquoise")
hemo_full_module_apop_df_turq$type <- "full"
Dermo_Tol_module_apop_df_turq <- Dermo_Tol_module_apop_df %>% filter(mod_names == "MEturquoise")
Dermo_Tol_module_apop_df_turq$type <- "consensus"

Dermo_Tol_turq_comparison <- full_join(hemo_full_module_apop_df_turq[,c("product","transcript_id","type")], Dermo_Tol_module_apop_df_turq[,c("product","transcript_id","type")], by ="transcript_id")
# few shared genes 

## Export modules to cytoscape for visualization ###
# Recalculate topological overlap if needed
#hemo_full_TOM = TOMsimilarityFromExpr(hemo_dds_rlog_matrix,
#                                           power = 3, # picked suitable power in the code above 
#                                           TOMType = "signed", # use signed TOM type
#                                           networkType= "signed hybrid", # use signed hybrid network type
#                                           corType = "bicor") # use suggested bicor
#
#save(hemo_full_TOM, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/hemo_full_TOM.RData" )

# Select modules
hemo_full_modules = c("darkslateblue", "turquoise", "greenyellow",   "skyblue3" , "cyan"  ,"red"  , "tan" )
# Select module probes
hemo_full_probes = colnames(hemo_dds_rlog_matrix)
# export moduleColors file for use in cluster
write.table(hemo_full_moduleColors, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/hemo_full_moduleColors.table")

hemo_full_inModule = is.finite(match(hemo_full_moduleColors, hemo_full_modules))
hemo_full_modProbes = hemo_full_probes[hemo_full_inModule]
hemo_full_modGenes = C_vir_rtracklayer$ID[match(hemo_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
hemo_full_modTOM = hemo_full_TOM[hemo_full_inModule, hemo_full_inModule]

dimnames(hemo_full_modTOM) = list(hemo_full_modProbes, hemo_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
#hemo_full_cyt = exportNetworkToCytoscape(hemo_full_modTOM,
#                               edgeFile = paste("CytoscapeInput-edges-", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
#                               weighted = TRUE,
#                               threshold = 0.02,
#                               nodeNames = hemo_full_modProbes,
#                               altNodeNames = hemo_full_modGenes,
#                               nodeAttr = hemo_full_moduleColors[hemo_full_inModule])
#

## Upload finished cytoscape network node file so that annotation information can be added
# The below code may not be necessary because I can add the apoptosis data table as a separate data table cytoscape will automatically attach 
# it to the network 
#hemo_full_cyt_node <- read.table(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/CytoscapeInput-nodes-Dermo_Tull_fulldarkslateblue-turquoise-greenyellow-skyblue3-cyan-red-tan.txt",
#                                      sep="\t", skip=1) # skip the first line because it has the original column names
#
## original col names were : nodeName  altName nodeAttr[nodesPresent, ]
#
#colnames(hemo_full_cyt_node)[c(1:3)] <- c("ID","altName","nodesPresent")
#hemo_full_cyt_node <- left_join(hemo_full_cyt_node, C_vir_rtracklayer_apop_product_final[,c("product","gene","ID","transcript_id")], by ="ID")
#
#write.table(hemo_full_cyt_node, sep = " ", quote= FALSE, row.names=FALSE, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/CytoscapeInput-nodes-Dermo_Tull_fulldarkslateblue-turquoise-greenyellow-skyblue3-cyan-red-tan_ANNOT.txt")
## Compare IAPs and GIMAPs between Consensus set and Full set

# Consensus set in significant modules
Dermo_Tol_module_apop_df_IAP <- Dermo_Tol_module_apop_df[grepl("IAP", Dermo_Tol_module_apop_df$product, ignore.case = TRUE),]
Dermo_Tol_module_apop_df_IAP <- Dermo_Tol_module_apop_df_IAP[order(Dermo_Tol_module_apop_df_IAP$product),]
Dermo_Tol_module_apop_df_GIMAP <- Dermo_Tol_module_apop_df[grepl("IMAP", Dermo_Tol_module_apop_df$product, ignore.case = TRUE),]
Dermo_Tol_module_apop_df_GIMAP <- Dermo_Tol_module_apop_df_GIMAP[order(Dermo_Tol_module_apop_df_GIMAP$product),]

# full set in significant modules
hemo_full_module_apop_df_IAP <-   hemo_full_module_apop_df[grepl("IAP", hemo_full_module_apop_df$product, ignore.case = TRUE),]
hemo_full_module_apop_df_IAP <-   hemo_full_module_apop_df_IAP[order(hemo_full_module_apop_df_IAP$product),]
hemo_full_module_apop_df_GIMAP <- hemo_full_module_apop_df[grepl("IMAP", hemo_full_module_apop_df$product, ignore.case = TRUE),]
hemo_full_module_apop_df_GIMAP <- hemo_full_module_apop_df_GIMAP[order(hemo_full_module_apop_df_GIMAP$product),]

setdiff(Dermo_Tol_module_apop_df_IAP$transcript_id,hemo_full_module_apop_df_IAP$transcript_id)
setdiff(hemo_full_module_apop_df_IAP$transcript_id,Dermo_Tol_module_apop_df_IAP$transcript_id)



