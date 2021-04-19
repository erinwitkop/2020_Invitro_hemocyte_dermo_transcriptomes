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
library(topGO)
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

# Import GO annotation data 
load(file = "../GO_universe_rna_found_geneID2GO_mapping.RData")
load(file = "../Perk_GO_terms_found_geneID2GO_mapping.RData")

## LOAD TRANSFORMED EXPRESSION DATA AND COLDATA ###

load(file='/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/2021_Hemocyte_Dermo_expression_rlog_matrices.RData')
load(file="/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/2021_Hemocyte_Dermo_expression_coldata.RData")

#####  DATA FORMATTING ####

hemo_dds_rlog_matrix <- assay(hemo_dds_rlog)
hemo_dds_rlog_matrix <- t(hemo_dds_rlog_matrix )
head(hemo_dds_rlog_matrix)
ncol(hemo_dds_rlog_matrix) 
colnames(hemo_dds_rlog_matrix)
names(hemo_dds_rlog_matrix) <- names(hemo_dds_rlog)

perk_dds_rlog_matrix <- assay(perk_dds_rlog)
perk_dds_rlog_matrix <- t(perk_dds_rlog_matrix )
head(perk_dds_rlog_matrix)
ncol(perk_dds_rlog_matrix) 
colnames(perk_dds_rlog_matrix)
names(perk_dds_rlog_matrix) <- names(perk_dds_rlog)

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

# all mod_names significant in either challenge
perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare_mod_names <- perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare$mod_names

# find those shared between all 
perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare_shared <- drop_na(perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare)
# 5 shared between all, all move in the same direction 

## Heatmap of only the significantly correlated modules 
# Graph and color code each the strength of association (correlation) of module eigengenes and trait

# subset perk_full_moduleTraitCor, perk_full_moduleTraitPvalue, perk_full_MEs for only those modules significant in either challenge
perk_full_moduleTraitCor_sig <- perk_full_moduleTraitCor[rownames(perk_full_moduleTraitCor) %in% perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare_mod_names,]
perk_full_moduleTraitCor_sig <- perk_full_moduleTraitCor_sig[,-3]
perk_full_moduleTraitPvalue_sig <- perk_full_moduleTraitPvalue[rownames(perk_full_moduleTraitPvalue) %in% perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare_mod_names,]
perk_full_moduleTraitPvalue_sig <- perk_full_moduleTraitPvalue_sig[,-3]
perk_full_MEs_sig <- perk_full_MEs[,colnames(perk_full_MEs) %in% perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare_mod_names]
perk_coldata_collapse_binarize_sig <- perk_coldata_collapse_binarize[,-3]

# Will display correlations and their p-values
perk_full_textMatrix_sig = paste(signif(perk_full_moduleTraitCor_sig, 2), "\n(",
                             signif(perk_full_moduleTraitPvalue_sig, 1), ")", sep = "");
dim(perk_full_textMatrix_sig) = dim(perk_full_moduleTraitCor_sig)

# make plot
sizeGrWindow(10,6)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = perk_full_moduleTraitCor_sig,
               xLabels = c("P.mar. and GDC-0152", "P. mar. and ZVAD-fmk"),
               yLabels = names(perk_full_MEs_sig),
               ySymbols = names(perk_full_MEs_sig),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = perk_full_textMatrix_sig,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.7,
               #zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("P. marinus Module-Challenge Relationships"))
# have to save manually...weird!


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

# P. mar Control vs P. mar_ZVAD
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

#### IDENTIFY HUB GENE IN EACH SIG MODULE ####

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

### DATA ANALYSIS NOTES ####

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


#### Identify all intramodular hub genes: high gene significance and high intramodular connectivity in interesting modules ####
# Remember that module membership is highly related to intramodular connectivity

# Which modules should be looked at?
  # For hemocytes: modules that are significant for treatment, have >5 apoptosis transcripts, and have high (>0.4) correlation between GS and MM
  # For hemocytes: same as above, except no filtering of modules for apoptosis related transcripts
# filter genes for each module of interest that have correlation between 

# annotate intramodular hub genes in each module
matrix = hemo_dds_rlog_matrix
lookup =   C_vir_rtracklayer_transcripts
datKME =  hemo_datKME

lookup_annot_intramodular <- function(list) {
FilterGenes = abs(GS1)> .6 & abs(datKME[list])>.8
FilterGenes_annot <- data.frame("ID" = dimnames(data.frame(matrix))[[2]][FilterGenes]) %>% left_join(., lookup) %>% mutate(mod_names = list)
}

## Hemocytes first
# control_Pmar
  # >0.6
    # lightcyan = 0.6 (highest association)
  # 0.4-0.6
    # red = cor = 0.44
    # darkseagreen4 = 0.4 
    # antiquewhite2 = 0.44
    # blue = 0.4

GS1 = hemo_full_geneTraitSignificance_control_Pmar
hemo_control_Pmar_intramodular <- c("MM.red", "MM.lightcyan", "MM.darkseagreen4", "MM.antiquewhite2", "MM.blue")
names(hemo_control_Pmar_intramodular) <- c("MM.red", "MM.lightcyan", "MM.darkseagreen4", "MM.antiquewhite2", "MM.blue")

FilterGenes_control_Pmar <- lapply(hemo_control_Pmar_intramodular, lookup_annot_intramodular)
FilterGenes_control_Pmar_df <- do.call(rbind,FilterGenes_control_Pmar) %>% mutate(group = "control_Pmar")

# Join with trait correlation for each module
FilterGenes_control_Pmar_df$mod_names <- str_replace(FilterGenes_control_Pmar_df$mod_names, "MM.","ME")
FilterGenes_control_Pmar_df <- left_join(FilterGenes_control_Pmar_df, hemo_full_moduleTraitCor_Pval_df_Pmar_vs_control_sig) %>% 
        dplyr::rename(moduleTraitCor = condition.Pmar.vs.control.moduleTraitCor, moduleTraitPvalue = condition.Pmar.vs.control.moduleTraitPvalue)
# join with gene trait significance 
hemo_full_geneTraitSignificance_control_Pmar_hub <- hemo_full_geneTraitSignificance_control_Pmar %>% rownames_to_column(., var = "ID") %>% filter(ID %in% unique(FilterGenes_control_Pmar_df$ID))
FilterGenes_control_Pmar_df <- left_join(FilterGenes_control_Pmar_df,hemo_full_geneTraitSignificance_control_Pmar_hub) %>% dplyr::rename(GS = GS.control_Pmar)

# control_Pmar_GDC 
  # modules of interest fitting criteria above  
  # >0.6
    # yellow = 0.71
    # antique white2 = 0.75
    # lightpink3 = 0.62
    # plum = 0.73 (highest association)
    # darkred = 0.64 
    # paleturquoise = 0.6
    # black = 0.61
    # darkorange = 0.62
  # 0.4 - 0.6
    # lightcyan = 0.55
    # navajowhite2 = 0.59
    # darkgreen = 0.4
    # mediumpurple3 = 0.53
    # cyan = 0.46

GS1 = hemo_full_geneTraitSignificance_control_Pmar_GDC
hemo_control_Pmar_GDC_intramodular <- c("MM.yellow","MM.lightcyan","MM.antiquewhite2","MM.navajowhite2","MM.darkgreen","MM.lightpink3",
                                        "MM.mediumpurple3" ,"MM.plum","MM.cyan","MM.darkred","MM.paleturquoise" ,"MM.black","MM.darkorange")
names(hemo_control_Pmar_GDC_intramodular) <- c("MM.yellow","MM.lightcyan","MM.antiquewhite2","MM.navajowhite2","MM.darkgreen","MM.lightpink3",
                                               "MM.mediumpurple3" ,"MM.plum","MM.cyan","MM.darkred","MM.paleturquoise" ,"MM.black","MM.darkorange")

FilterGenes_control_Pmar_GDC <- lapply(hemo_control_Pmar_GDC_intramodular, lookup_annot_intramodular)
FilterGenes_control_Pmar_GDC_df <- do.call(rbind,FilterGenes_control_Pmar_GDC) %>% mutate(group = "control_Pmar_GDC")

# Join with trait correlation for each module
FilterGenes_control_Pmar_GDC_df$mod_names <- str_replace(FilterGenes_control_Pmar_GDC_df$mod_names, "MM.","ME")
FilterGenes_control_Pmar_GDC_df <- left_join(FilterGenes_control_Pmar_GDC_df, hemo_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_control_sig) %>% 
  dplyr::rename(moduleTraitCor = condition.Pmar_GDC.vs.control.moduleTraitCor, moduleTraitPvalue = condition.Pmar_GDC.vs.control.moduleTraitPvalue)
# join with gene trait significance 
hemo_full_geneTraitSignificance_control_Pmar_GDC_hub <- hemo_full_geneTraitSignificance_control_Pmar_GDC %>% rownames_to_column(., var = "ID") %>% filter(ID %in% unique(FilterGenes_control_Pmar_GDC_df$ID))
FilterGenes_control_Pmar_GDC_df <- left_join(FilterGenes_control_Pmar_GDC_df,hemo_full_geneTraitSignificance_control_Pmar_GDC_hub) %>% dplyr::rename(GS = GS.control_Pmar_GDC)
 
# control_Pmar_ZVAD 
  # modules of interest fitting criteria above
  # > 0.6 
    # plum = 0.61
  # 0.4-0.6
      # yellow = 0.48
      # orangered4 = 0.52
      # lightcyan = 0.56
      # darkslateblue = 0.57

GS1 = hemo_full_geneTraitSignificance_control_Pmar_ZVAD
hemo_control_Pmar_ZVAD_intramodular <- c("MM.yellow", "MM.orangered4", "MM.lightcyan", "MM.plum","MM.darkslateblue")
names(hemo_control_Pmar_ZVAD_intramodular) <- c("MM.yellow", "MM.orangered4", "MM.lightcyan", "MM.plum","MM.darkslateblue")

FilterGenes_control_Pmar_ZVAD <- lapply(hemo_control_Pmar_ZVAD_intramodular, lookup_annot_intramodular)
FilterGenes_control_Pmar_ZVAD_df <- do.call(rbind,FilterGenes_control_Pmar_ZVAD) %>% mutate(group = "control_Pmar_ZVAD")

# Join with trait correlation for each module
FilterGenes_control_Pmar_ZVAD_df$mod_names <- str_replace(FilterGenes_control_Pmar_ZVAD_df$mod_names, "MM.","ME")
FilterGenes_control_Pmar_ZVAD_df <- left_join(FilterGenes_control_Pmar_ZVAD_df, hemo_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_control_sig) %>% 
  dplyr::rename(moduleTraitCor = condition.Pmar_ZVAD.vs.control.moduleTraitCor, moduleTraitPvalue = condition.Pmar_ZVAD.vs.control.moduleTraitPvalue)
# join with gene trait significance 
hemo_full_geneTraitSignificance_control_Pmar_ZVAD_hub <- hemo_full_geneTraitSignificance_control_Pmar_ZVAD %>% rownames_to_column(., var = "ID") %>% filter(ID %in% unique(FilterGenes_control_Pmar_ZVAD_df$ID))
FilterGenes_control_Pmar_ZVAD_df <- left_join(FilterGenes_control_Pmar_ZVAD_df,hemo_full_geneTraitSignificance_control_Pmar_ZVAD_hub) %>% dplyr::rename(GS = GS.control_Pmar_ZVAD)

## Combine hemocyte data frames
FilterGenes_comb <- rbind(FilterGenes_control_Pmar_df,
                          FilterGenes_control_Pmar_GDC_df,
                          FilterGenes_control_Pmar_ZVAD_df)

# how many in each module per treatment?
FilterGenes_comb %>% group_by(mod_names, group) %>% dplyr::count() %>% View()

# are any of these apoptotic?
FilterGenes_comb_apop <- FilterGenes_comb[FilterGenes_comb$ID %in% C_vir_rtracklayer_apop_product_final$ID,] %>% 
  dplyr::select(ID, gene,product, transcript_id, mod_names, group, moduleTraitCor,moduleTraitPvalue, GS)

# how many apop intramodular hub genes in each module per treatment?
FilterGenes_comb_apop %>% group_by(mod_names, group) %>% dplyr::count() %>% View()

### Repeat finding of intramodular hub genes, but for P. marinus samples now 
# annotate intramodular hub genes in each module
# Use function to lookup all apop names for each significant module
matrix = perk_dds_rlog_matrix
lookup =   Perkinsus_rtracklayer
datKME =  perk_datKME

lookup_annot_intramodular_perk <- function(list) {
  FilterGenes = abs(GS1)> .6 & abs(datKME[list])>.8
  FilterGenes_annot <- data.frame("Name" = dimnames(data.frame(matrix))[[2]][FilterGenes]) %>% left_join(., lookup) %>% mutate(mod_names = list)
}

# Pmar_vs_Pmar_GDC
# modules of interest fitting criteria above
# >0.6
    #darkseagreen2 = 0.67
    #yellow4 = 0.66   
    #blue4 = 0.78
    #lightblue4 = 0.78
# 0.4-0.6
    #orange = 0.52
    #indianred3 = 0.48
    #steelblue = 0.42
    #darkorange2 = 0.57
    #darkturquoise = 0.49
    #indianred4 = 0.47

GS1 = perk_full_geneTraitSignificance_Pmar_Pmar_GDC
Pmar_intramodular <- c("MM.darkseagreen2","MM.orange","MM.indianred3","MM.yellow4","MM.steelblue","MM.darkorange2","MM.darkturquoise","MM.indianred4","MM.blue4","MM.lightblue4")
names(Pmar_intramodular) <- c("MM.darkseagreen2","MM.orange","MM.indianred3","MM.yellow4","MM.steelblue","MM.darkorange2","MM.darkturquoise","MM.indianred4","MM.blue4","MM.lightblue4")

FilterGenes_Pmar_Pmar_GDC <- lapply(Pmar_intramodular, lookup_annot_intramodular_perk)
FilterGenes_Pmar_Pmar_GDC_df <- do.call(rbind,FilterGenes_Pmar_Pmar_GDC) %>% mutate(group = "Pmar_Pmar_GDC")

# Join with trait correlation for each module
FilterGenes_Pmar_Pmar_GDC_df$mod_names <- str_replace(FilterGenes_Pmar_Pmar_GDC_df$mod_names, "MM.","ME")
FilterGenes_Pmar_Pmar_GDC_df <- perk_full_moduleTraitCor_Pval_df_Pmar_GDC_vs_Pmar %>% rownames_to_column(., "mod_names") %>%
  left_join(FilterGenes_Pmar_Pmar_GDC_df, .) %>%
  dplyr::rename(moduleTraitCor = condition.Pmar_GDC.vs.Pmar.moduleTraitCor, moduleTraitPvalue = condition.Pmar_GDC.vs.Pmar.moduleTraitPvalue)
# join with gene trait significance 
perk_full_geneTraitSignificance_Pmar_Pmar_GDC_hub <- perk_full_geneTraitSignificance_Pmar_Pmar_GDC %>% rownames_to_column(., var = "transcript_id") %>% filter(transcript_id %in% unique(FilterGenes_Pmar_Pmar_GDC_df$transcript_id))
FilterGenes_Pmar_Pmar_GDC_df <- left_join(FilterGenes_Pmar_Pmar_GDC_df,perk_full_geneTraitSignificance_Pmar_Pmar_GDC_hub) %>% dplyr::rename(GS = GS.Pmar_Pmar_GDC)

## Pmar_Pmar_ZVAD
# >0.6
  #darkseagreen2 = 0.73
  #pink3 = 0.78
  #navajowhite2 = 0.7
  # darkturqiouse = 0.61
#0.4-0.6
  #darkseagreen4 = 0.42
  #lightpink3 = 0.57
  #blue4 = 0.53

GS1 = perk_full_geneTraitSignificance_Pmar_Pmar_ZVAD
Pmar_intramodular <- c("MM.darkseagreen2","MM.darkseagreen4","MM.lightpink3","MM.pink3","MM.navajowhite2","MM.darkturquoise","MM.blue4")
names(Pmar_intramodular) <- c("MM.darkseagreen2","MM.darkseagreen4","MM.lightpink3","MM.pink3","MM.navajowhite2","MM.darkturquoise","MM.blue4")

FilterGenes_Pmar_Pmar_ZVAD <- lapply(Pmar_intramodular, lookup_annot_intramodular_perk)
FilterGenes_Pmar_Pmar_ZVAD_df <- do.call(rbind,FilterGenes_Pmar_Pmar_ZVAD) %>% mutate(group = "Pmar_Pmar_ZVAD")

# Join with trait correlation for each module
FilterGenes_Pmar_Pmar_ZVAD_df$mod_names <- str_replace(FilterGenes_Pmar_Pmar_ZVAD_df$mod_names, "MM.","ME")
FilterGenes_Pmar_Pmar_ZVAD_df <- perk_full_moduleTraitCor_Pval_df_Pmar_ZVAD_vs_Pmar %>% rownames_to_column(., "mod_names") %>%
  left_join(FilterGenes_Pmar_Pmar_ZVAD_df, .) %>%
  dplyr::rename(moduleTraitCor = condition.Pmar_ZVAD.vs.Pmar.moduleTraitCor, moduleTraitPvalue = condition.Pmar_ZVAD.vs.Pmar.moduleTraitPvalue)
# join with gene trait significance 
perk_full_geneTraitSignificance_Pmar_Pmar_ZVAD_hub <- perk_full_geneTraitSignificance_Pmar_Pmar_ZVAD %>% rownames_to_column(., var = "transcript_id") %>% filter(transcript_id %in% unique(FilterGenes_Pmar_Pmar_ZVAD_df$transcript_id))
FilterGenes_Pmar_Pmar_ZVAD_df <- left_join(FilterGenes_Pmar_Pmar_ZVAD_df,perk_full_geneTraitSignificance_Pmar_Pmar_ZVAD_hub) %>% dplyr::rename(GS = GS.Pmar_Pmar_ZVAD)

# Combine P.mar dataframes
FilterGenes_Pmar_comb <- rbind(FilterGenes_Pmar_Pmar_ZVAD_df, FilterGenes_Pmar_Pmar_GDC_df)

# how many in each module per treatment?
FilterGenes_Pmar_comb %>% group_by(mod_names, group) %>% dplyr::mutate(count = n()) %>% dplyr::distinct(mod_names, group, count, moduleTraitCor, moduleTraitPvalue) %>% View()
  # many intramodular hub

# Combine with Interproscan 
FilterGenes_Pmar_comb_Interpro <- left_join(FilterGenes_Pmar_comb[,c("Name","product","transcript_id","moduleTraitCor","moduleTraitPvalue","mod_names","GS")], 
                                            Perk_Interpro_GO_terms_XP[,c("protein_id","source","transcript_id", "Ontology_term",
                                                                         "Dbxref","signature_desc")], by = "transcript_id")

#### Pmar intramodular hub gene analysis for interesting modules ####

# Interesting modules
  # What is criteria for an interesting module?
      # 1. Module is significantly correlated with challenge group
      # 2. Module has high (>0.6) correlation between gene significance and module membership
      # 3. Module is either uniquely and highly (>0.8) significant to the particular challenge, or lowly (<0.4) correlated with the other treatment    
  # Based on this the following modules are interest 
  # 1. Pmar vs GDC: the steelblue, darkorange2, lightblue4
  # 2. Pmar vs ZVAD: the linkpink3, pink3, and navajowhite2 are of interest

# filter out GDC modules from the list of intramodular hub genes (high GS and high modulemembership)
FilterGenes_Pmar_comb_Interpro_GDC_steelblue <- FilterGenes_Pmar_comb_Interpro %>% filter(mod_names == "MEsteelblue")
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2 <- FilterGenes_Pmar_comb_Interpro %>% filter(mod_names == "MEdarkorange2")
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4 <- FilterGenes_Pmar_comb_Interpro %>% filter(mod_names == "MElightblue4")

# filter out ZVAD modules from the list of intramodular hub genes (high GS and high modulemembership)
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3 <- FilterGenes_Pmar_comb_Interpro %>% filter(mod_names == "MElightpink3")
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3 <- FilterGenes_Pmar_comb_Interpro %>% filter(mod_names == "MEpink3")
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2 <- FilterGenes_Pmar_comb_Interpro %>% filter(mod_names == "MEnavajowhite2")

# How many intramodular hub genes in each with very high positive GS? Perform for interesting modules lightblue4, pink3, navajowhite2
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4 %>% distinct(transcript_id, GS) %>% filter( GS >= 0.6) %>% dplyr::count() # 137
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3  %>% distinct(transcript_id, GS) %>% filter( GS >= 0.6) %>% dplyr::count() # 51
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3  %>% distinct(transcript_id, GS) %>% filter( GS >= 0.6) %>% dplyr::count() # 52

FilterGenes_Pmar_comb_Interpro_GDC_lightblue4.8 <- FilterGenes_Pmar_comb_Interpro_GDC_lightblue4 %>% filter( GS >= 0.8)
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3.8 <- FilterGenes_Pmar_comb_Interpro_ZVAD_pink3 %>% filter( GS >= 0.8)
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3.8 <- FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3 %>% filter( GS >= 0.8)

## assess GO terms from important module intramodular hub genes 
# first format the GO terms so I can get unique terms for each protein

## GDC steelblue
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO <- FilterGenes_Pmar_comb_Interpro_GDC_steelblue %>% 
  filter(Ontology_term != "character(0)" & !is.na(transcript_id)) %>% distinct(transcript_id,Ontology_term)
# format GO columns
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term <- unlist(as.character(FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term))
class(FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term)
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term <- gsub('\\\"', "", FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term <- gsub('\"', "", FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term <- gsub('c(', "", FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term <- gsub(')', "", FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term <- gsub(' ', "", FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO$Ontology_term, fixed = TRUE)
# separate rows now into separate rows and get unique GO terms for each protein
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO <- separate_rows(FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO,2,sep = ",") %>% distinct(transcript_id, Ontology_term)
# get term frequency
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO_freq <- FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO %>% ungroup() %>% 
  group_by(Ontology_term) %>% dplyr::count()

## GDC darkorange2
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO <- FilterGenes_Pmar_comb_Interpro_GDC_darkorange2 %>% 
  filter(Ontology_term != "character(0)" & !is.na(transcript_id)) %>% distinct(transcript_id,Ontology_term)
# format GO columns
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term <- unlist(as.character(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term))
class(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term)
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term <- gsub('\\\"', "", FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term <- gsub('\"', "", FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term <- gsub('c(', "", FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term <- gsub(')', "", FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term <- gsub(' ', "", FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO$Ontology_term, fixed = TRUE)
# separate rows now into separate rows and get unique GO terms for each protein
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO <- separate_rows(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO,2,sep = ",") %>% distinct(transcript_id, Ontology_term)
# get term frequency
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO_freq <- FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO %>% ungroup() %>% 
  group_by(Ontology_term) %>% dplyr::count()

## GDC lightblue4
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO <- FilterGenes_Pmar_comb_Interpro_GDC_lightblue4 %>% 
  filter(Ontology_term != "character(0)" & !is.na(transcript_id)) %>% distinct(transcript_id,Ontology_term)
# format GO columns
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term <- unlist(as.character(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term))
class(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term)
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term <- gsub('\\\"', "", FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term <- gsub('\"', "", FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term <- gsub('c(', "", FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term <- gsub(')', "", FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term <- gsub(' ', "", FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO$Ontology_term, fixed = TRUE)
# separate rows now into separate rows and get unique GO terms for each protein
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO <- separate_rows(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO,2,sep = ",") %>% distinct(transcript_id, Ontology_term)
# get term frequency
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO_freq <- FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO %>% ungroup() %>% 
  group_by(Ontology_term) %>% dplyr::count()

## ZVAD lightpink3
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO <- FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3 %>% 
  filter(Ontology_term != "character(0)" & !is.na(transcript_id)) %>% distinct(transcript_id,Ontology_term)
# format GO columns
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term <- unlist(as.character(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term))
class(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term)
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term <- gsub('\\\"', "", FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term <- gsub('\"', "", FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term <- gsub('c(', "", FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term <- gsub(')', "", FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term <- gsub(' ', "", FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO$Ontology_term, fixed = TRUE)
# separate rows now into separate rows and get unique GO terms for each protein
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO <- separate_rows(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO,2,sep = ",") %>% distinct(transcript_id, Ontology_term)
# get term frequency
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO_freq <- FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO %>% ungroup() %>% 
  group_by(Ontology_term) %>% dplyr::count()

## ZVAD pink3
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO <- FilterGenes_Pmar_comb_Interpro_ZVAD_pink3 %>% 
  filter(Ontology_term != "character(0)" & !is.na(transcript_id)) %>% distinct(transcript_id,Ontology_term)
# format GO columns
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term <- unlist(as.character(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term))
class(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term)
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term <- gsub('\\\"', "", FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term <- gsub('\"', "", FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term <- gsub('c(', "", FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term <- gsub(')', "", FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term <- gsub(' ', "", FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO$Ontology_term, fixed = TRUE)
# separate rows now into separate rows and get unique GO terms for each protein
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO <- separate_rows(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO,2,sep = ",") %>% distinct(transcript_id, Ontology_term)
# get term frequency
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO_freq <- FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO %>% ungroup() %>% 
  group_by(Ontology_term) %>% dplyr::count()

## ZVAD navajowhite2
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO <- FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2 %>% 
  filter(Ontology_term != "character(0)" & !is.na(transcript_id)) %>% distinct(transcript_id,Ontology_term)
# format GO columns
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term <- unlist(as.character(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term))
class(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term)
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term <- gsub('\\\"', "", FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term <- gsub('\"', "", FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term <- gsub('c(', "", FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term <- gsub(')', "", FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term, fixed = TRUE)
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term <- gsub(' ', "", FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO$Ontology_term, fixed = TRUE)
# separate rows now into separate rows and get unique GO terms for each protein
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO <- separate_rows(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO,2,sep = ",") %>%
  distinct(transcript_id, Ontology_term) %>% filter(!is.na(Ontology_term))
# get term frequency
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO_freq <- FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO %>% ungroup() %>% 
  group_by(Ontology_term) %>% dplyr::count()

# export all to REVIGO to get pictures 
write.table(FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO_freq, "FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO.txt",row.names = FALSE,quote = FALSE)
write.table(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO_freq, "FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO.txt",row.names = FALSE,quote = FALSE)
write.table(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO_freq, "FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO.txt",row.names = FALSE,quote = FALSE)
write.table(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO_freq, "FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO.txt",row.names = FALSE,quote = FALSE)
write.table(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO_freq, "FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO.txt",row.names = FALSE,quote = FALSE)
write.table(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO_freq, "FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO.txt",row.names = FALSE,quote = FALSE)

### PMAR GO ENRICHMENT FOR INTERESTING MODULE INTRAMODULAR HUB GENES ####
Perk_geneNames <- names(Perk_GO_terms_found_geneID2GO_mapping)
# topGO tutorial: https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
head(Perk_geneNames)

## GO enrichment for just the Intramodular hub genes 
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO_only     <-    as.factor(na.omit(unique(FilterGenes_Pmar_comb_Interpro_GDC_steelblue$transcript_id)))
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO_only   <-  as.factor(na.omit(unique(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2$transcript_id)))
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO_only    <-   as.factor(na.omit(unique(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4$transcript_id)))
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO_only   <-  as.factor(na.omit(unique(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3$transcript_id)))
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO_only        <-       as.factor(na.omit(unique(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3$transcript_id)))
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO_only <- as.factor(na.omit(unique(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2$transcript_id)))

FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO_only_factor <- factor(as.integer(Perk_geneNames %in% FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO_only))
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO_only_factor <- factor(as.integer(Perk_geneNames %in% FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO_only))
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO_only_factor <- factor(as.integer(Perk_geneNames %in% FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO_only))
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO_only_factor <- factor(as.integer(Perk_geneNames %in% FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO_only))
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO_only_factor <- factor(as.integer(Perk_geneNames %in% FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO_only))
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO_only_factor <- factor(as.integer(Perk_geneNames %in% FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO_only))

names(FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO_only_factor) <- Perk_geneNames
names(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO_only_factor) <- Perk_geneNames
names(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO_only_factor) <- Perk_geneNames
names(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO_only_factor) <- Perk_geneNames
names(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO_only_factor) <- Perk_geneNames
names(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO_only_factor) <- Perk_geneNames

### Make topGO data object 
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata <- new("topGOdata", description = "GDC_steelblue Gene Enrichment", 
                           # I want to test MF
                           ontology = "MF",
                           # define here the genes of interest
                           allGenes = FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GO_only_factor,
                           nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata <- new("topGOdata", description = "GDC_darkorange2 Gene Enrichment", 
                                                           # I want to test MF
                                                           ontology = "MF",
                                                           # define here the genes of interest
                                                           allGenes = FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GO_only_factor,
                                                           nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata <- new("topGOdata", description = "GDC_lightblue4 Gene Enrichment", 
                                                           # I want to test MF
                                                           ontology = "MF",
                                                           # define here the genes of interest
                                                           allGenes = FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO_only_factor,
                                                           nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata <- new("topGOdata", description = "ZVAD_lightpink3 Gene Enrichment", 
                                                            # I want to test MF
                                                            ontology = "MF",
                                                            # define here the genes of interest
                                                            allGenes = FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GO_only_factor,
                                                            nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata <- new("topGOdata", description = "ZVAD_pink3 Gene Enrichment", 
                                                             # I want to test MF
                                                             ontology = "MF",
                                                             # define here the genes of interest
                                                             allGenes = FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GO_only_factor,
                                                             nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata <- new("topGOdata", description = "ZVAD_navajowhite2 Gene Enrichment", 
                                                        # I want to test MF
                                                        ontology = "MF",
                                                        # define here the genes of interest
                                                        allGenes = FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GO_only_factor,
                                                        nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

#nodeSize=used to prune the GO hierarchy from the terms which have less than 1 annotated genes
#annFUN.gene2GO = this function is used when the annotations are provided as a gene-to-GOs mapping.

### Perform Encrichment tests 
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata_Fisher_Weight <- runTest(FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata, algorithm = "weight01", statistic = "fisher")
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata_Fisher_Weight <- runTest(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata, algorithm = "weight01", statistic = "fisher")
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Fisher_Weight <- runTest(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata, algorithm = "weight01", statistic = "fisher")
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata_Fisher_Weight <- runTest(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata, algorithm = "weight01", statistic = "fisher")
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_Fisher_Weight <- runTest(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata, algorithm = "weight01", statistic = "fisher")
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Fisher_Weight <- runTest(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata, algorithm = "weight01", statistic = "fisher")

## Analyze enrichment test results 
# see how many results we get where weight01 gives a P-value <= 0.05
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata_summary <- summary(attributes(FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata_Fisher_Weight)$score <= 0.05)
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata_summary <- summary(attributes(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata_Fisher_Weight)$score <= 0.05)
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_summary <- summary(attributes(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Fisher_Weight)$score <= 0.05)
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata_summary <- summary(attributes(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata_Fisher_Weight)$score <= 0.05)
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_summary <- summary(attributes(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_Fisher_Weight)$score <= 0.05)
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_summary <- summary(attributes(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Fisher_Weight)$score <= 0.05)

#print out the top results
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)
    #GO.ID                                        Term Annotated Significant Expected topgoFisher
    #1  GO:0005200      structural constituent of cytoskeleton        15           2     0.13      0.0045
    #2  GO:0005506                            iron ion binding        50           3     0.45      0.0046
    #3  GO:0003824                          catalytic activity      3307          40    29.58      0.0059
    #4  GO:0008410                    CoA-transferase activity         5           1     0.04      0.0333
    #5  GO:0003934               GTP cyclohydrolase I activity         5           1     0.04      0.0333
    #6  GO:0004857                   enzyme inhibitor activity         5           1     0.04      0.0333
    #7  GO:0004807         triose-phosphate isomerase activity         5           1     0.04      0.0333
    #8  GO:0016972                      thiol oxidase activity         7           1     0.06      0.0464
    #9  GO:0016684 oxidoreductase activity, acting on perox...         7           1     0.06      0.0464

FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)

FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)
  #GO.ID                                        Term Annotated Significant Expected topgoFisher
  #1  GO:0004190        aspartic-type endopeptidase activity        68           3     0.25     0.00068
  #2  GO:0016646 oxidoreductase activity, acting on the C...         6           1     0.02     0.01532
  #3  GO:0003924                             GTPase activity        84           2     0.31     0.01947
  #4  GO:0015299           solute:proton antiporter activity         9           1     0.03     0.02289
  #5  GO:0046933 proton-transporting ATP synthase activit...        11           1     0.04     0.02791
  #6  GO:0008536                          Ran GTPase binding        14           1     0.05     0.03539
  #7  GO:0016765 transferase activity, transferring alkyl...        20           1     0.07     0.05019 

FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)
  #GO.ID                                        Term Annotated Significant Expected topgoFisher
  #1  GO:0008168                  methyltransferase activity       109           5     1.48      0.0026
  #2  GO:0009678 hydrogen-translocating pyrophosphatase a...        11           2     0.15      0.0040
  #3  GO:0010181                                 FMN binding        12           2     0.16      0.0047
  #4  GO:0004427            inorganic diphosphatase activity        16           2     0.22      0.0084
  #5  GO:0019825                              oxygen binding        24           2     0.33      0.0185
  #6  GO:0020037                                heme binding        66           3     0.90      0.0200
  #7  GO:0003824                          catalytic activity      3307          47    45.00      0.0219
  #8  GO:0043169                              cation binding       641          10     8.72      0.0332
  #9  GO:0017076                   purine nucleotide binding       921          11    12.53      0.0338
  #10 GO:0004198 calcium-dependent cysteine-type endopept...         5           1     0.07      0.0430
  #11 GO:0008234            cysteine-type peptidase activity       168           5     2.29      0.0464
  #12 GO:0008270                            zinc ion binding       228           5     3.10      0.0495


#### PMAR GO ENRICHMENT FOR INTERESTING MODULE ALL TRANSCRIPTS ####
## GO enrichment for all genes in the module
# get full list of genes for each module 
GDC_steelblue_all_genes <- names(perk_dds_rlog_matrix)[perk_full_moduleColors =="steelblue"] # only 127 in the module 
GDC_steelblue_all_genes <- as.factor(GDC_steelblue_all_genes[!is.na(GDC_steelblue_all_genes)])

GDC_darkorange2_all_genes <- names(perk_dds_rlog_matrix)[perk_full_moduleColors =="darkorange2"] 
GDC_darkorange2_all_genes <- as.factor(GDC_darkorange2_all_genes[!is.na(GDC_darkorange2_all_genes)])

GDC_lightblue4_all_genes <- names(perk_dds_rlog_matrix)[perk_full_moduleColors =="lightblue4"]  
GDC_lightblue4_all_genes <- as.factor(GDC_lightblue4_all_genes[!is.na(GDC_lightblue4_all_genes)])

ZVAD_lightpink3_all_genes <- names(perk_dds_rlog_matrix)[perk_full_moduleColors =="lightpink3"]  
ZVAD_lightpink3_all_genes <- as.factor(ZVAD_lightpink3_all_genes[!is.na(ZVAD_lightpink3_all_genes)])

ZVAD_pink3_all_genes <- names(perk_dds_rlog_matrix)[perk_full_moduleColors =="pink3"]  
ZVAD_pink3_all_genes <- as.factor(ZVAD_pink3_all_genes[!is.na(ZVAD_pink3_all_genes)])

ZVAD_navajowhite2_all_genes <- names(perk_dds_rlog_matrix)[perk_full_moduleColors =="navajowhite2"]  
ZVAD_navajowhite2_all_genes <- ZVAD_navajowhite2_all_genes[!is.na(ZVAD_navajowhite2_all_genes)]

GDC_steelblue_all_genes_factor <- factor(as.integer(Perk_geneNames %in% GDC_steelblue_all_genes))
GDC_darkorange2_all_genes_factor <- factor(as.integer(Perk_geneNames %in% GDC_darkorange2_all_genes))
GDC_lightblue4_all_genes_factor <- factor(as.integer(Perk_geneNames %in% GDC_lightblue4_all_genes))
ZVAD_lightpink3_all_genes_factor <- factor(as.integer(Perk_geneNames %in% ZVAD_lightpink3_all_genes))
ZVAD_pink3_all_genes_factor <- factor(as.integer(Perk_geneNames %in% ZVAD_pink3_all_genes))
ZVAD_navajowhite2_all_genes_factor <- factor(as.integer(Perk_geneNames %in% ZVAD_navajowhite2_all_genes))

names(GDC_steelblue_all_genes_factor) <- Perk_geneNames
names(GDC_darkorange2_all_genes_factor) <- Perk_geneNames
names(GDC_lightblue4_all_genes_factor) <- Perk_geneNames
names(ZVAD_lightpink3_all_genes_factor) <- Perk_geneNames
names(ZVAD_pink3_all_genes_factor) <- Perk_geneNames
names(ZVAD_navajowhite2_all_genes_factor) <- Perk_geneNames

### Make topGO data object 
GDC_steelblue_GOdata <- new("topGOdata", description = "GDC_steelblue Gene Enrichment", 
                                                           # I want to test MF
                                                           ontology = "MF",
                                                           # define here the genes of interest
                                                           allGenes = GDC_steelblue_all_genes_factor,
                                                           nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

GDC_darkorange2_GOdata <- new("topGOdata", description = "GDC_darkorange2 Gene Enrichment", 
                                                             # I want to test MF
                                                             ontology = "MF",
                                                             # define here the genes of interest
                                                             allGenes = GDC_darkorange2_all_genes_factor,
                                                             nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

GDC_lightblue4_GOdata <- new("topGOdata", description = "GDC_lightblue4 Gene Enrichment", 
                                                            # I want to test MF
                                                            ontology = "MF",
                                                            # define here the genes of interest
                                                            allGenes = GDC_lightblue4_all_genes_factor,
                                                            nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

ZVAD_lightpink3_GOdata <- new("topGOdata", description = "ZVAD_lightpink3 Gene Enrichment", 
                                                             # I want to test MF
                                                             ontology = "MF",
                                                             # define here the genes of interest
                                                             allGenes = ZVAD_lightpink3_all_genes_factor,
                                                             nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

ZVAD_pink3_GOdata <- new("topGOdata", description = "ZVAD_pink3 Gene Enrichment", 
                                                        # I want to test MF
                                                        ontology = "MF",
                                                        # define here the genes of interest
                                                        allGenes = ZVAD_pink3_all_genes_factor,
                                                        nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

ZVAD_navajowhite2_GOdata <- new("topGOdata", description = "ZVAD_navajowhite2 Gene Enrichment", 
                                                               # I want to test MF
                                                               ontology = "MF",
                                                               # define here the genes of interest
                                                               allGenes = ZVAD_navajowhite2_all_genes_factor,
                                                               nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

#nodeSize=used to prune the GO hierarchy from the terms which have less than 1 annotated genes
#annFUN.gene2GO = this function is used when the annotations are provided as a gene-to-GOs mapping.

### Perform Encrichment tests 
GDC_steelblue_GOdata_Fisher_Weight <- runTest(GDC_steelblue_GOdata, algorithm = "weight01", statistic = "fisher")
GDC_darkorange2_GOdata_Fisher_Weight <- runTest(GDC_darkorange2_GOdata, algorithm = "weight01", statistic = "fisher")
GDC_lightblue4_GOdata_Fisher_Weight <- runTest(GDC_lightblue4_GOdata, algorithm = "weight01", statistic = "fisher")
ZVAD_lightpink3_GOdata_Fisher_Weight <- runTest(ZVAD_lightpink3_GOdata, algorithm = "weight01", statistic = "fisher")
ZVAD_pink3_GOdata_Fisher_Weight <- runTest(ZVAD_pink3_GOdata, algorithm = "weight01", statistic = "fisher")
ZVAD_navajowhite2_GOdata_Fisher_Weight <- runTest(ZVAD_navajowhite2_GOdata, algorithm = "weight01", statistic = "fisher")

## Analyze enrichment test results 
# see how many results we get where weight01 gives a P-value <= 0.05
GDC_steelblue_GOdata_summary <- summary(attributes(GDC_steelblue_GOdata_Fisher_Weight)$score <= 0.05) # 9 sig
GDC_darkorange2_GOdata_summary <- summary(attributes(GDC_darkorange2_GOdata_Fisher_Weight)$score <= 0.05) # 8 sig
GDC_lightblue4_GOdata_summary <- summary(attributes(GDC_lightblue4_GOdata_Fisher_Weight)$score <= 0.05) # 4 sig
ZVAD_lightpink3_GOdata_summary <- summary(attributes(ZVAD_lightpink3_GOdata_Fisher_Weight)$score <= 0.05) # 9 sig
ZVAD_pink3_GOdata_summary <- summary(attributes(ZVAD_pink3_GOdata_Fisher_Weight)$score <= 0.05) # 3 sig
ZVAD_navajowhite2_GOdata_summary <- summary(attributes(ZVAD_navajowhite2_GOdata_Fisher_Weight)$score <= 0.05) # 11

#print out the top results, though only GDC_lightblue4 is sig
GDC_steelblue_GOdata_Res <- GenTable(GDC_steelblue_GOdata, topgoFisher = GDC_steelblue_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
GDC_darkorange2_GOdata_Res <- GenTable(GDC_darkorange2_GOdata, topgoFisher = GDC_darkorange2_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
GDC_lightblue4_GOdata_Res <- GenTable(GDC_lightblue4_GOdata, topgoFisher = GDC_lightblue4_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
      #GO.ID                                        Term Annotated Significant Expected topgoFisher
      #1  GO:0005506                            iron ion binding        50           2     0.12      0.0044
      #2  GO:0008375      acetylglucosaminyltransferase activity        11           1     0.03      0.0218

ZVAD_lightpink3_GOdata_Res <- GenTable(ZVAD_lightpink3_GOdata, topgoFisher = ZVAD_lightpink3_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
ZVAD_pink3_GOdata_Res <- GenTable(ZVAD_pink3_GOdata, topgoFisher = ZVAD_pink3_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
      #GO.ID                                        Term Annotated Significant Expected topgoFisher
      #1  GO:0004190        aspartic-type endopeptidase activity        68           2     0.07     0.00085
      #2  GO:0046933 proton-transporting ATP synthase activit...        11           1     0.01     0.00731
      #3  GO:0016887                             ATPase activity       116           2     0.12     0.03838

ZVAD_navajowhite2_GOdata_Res <- GenTable(ZVAD_navajowhite2_GOdata, topgoFisher = ZVAD_navajowhite2_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
      #GO.ID                                        Term Annotated Significant Expected topgoFisher
      #1  GO:0043169                              cation binding       641           7     2.87      0.0097
      #2  GO:0003824                          catalytic activity      3307          17    14.79      0.0099
      #3  GO:0008483                       transaminase activity         6           1     0.03      0.0170
      #4  GO:0004376     glycolipid mannosyltransferase activity         6           1     0.03      0.0170
      #5  GO:0005507                          copper ion binding        10           1     0.04      0.0282
      #6  GO:0004129               cytochrome-c oxidase activity        10           1     0.04      0.0282
      #7  GO:0030976              thiamine pyrophosphate binding        11           1     0.05      0.0310
      #8  GO:0008168                  methyltransferase activity       109           2     0.49      0.0384
      #9  GO:0008138 protein tyrosine/serine/threonine phosph...        14           1     0.06      0.0393

#### PMAR TOP INTERESTING MODULES GO VISUALIZATION ####

## Dotplot of top modules GO terms
  # GDC_lightblue4, ZVAD navajowhite2 and ZVAD pink 3

FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_Res
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Res
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Res

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

