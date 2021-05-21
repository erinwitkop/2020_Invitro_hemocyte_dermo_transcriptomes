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
library(clusterProfiler)
library(DOSE)
library(viridis)
library(metafor)
# Using R version 3.6.1

# helpful WGCNA tutorials and FAQs
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

#### LOADING SAVED GENOME, APOPTOSIS NAMES, IAP XP LISTS ####
Apoptosis_frames <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_apoptosis_products.RData")
annotations <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_annotations.RData")

Perkinsus_rtracklayer <- readGFF("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/GCF_000006405.1_JCVI_PMG_1.0_genomic.gff")
Perkinsus_rtracklayer <- as.data.frame(Perkinsus_rtracklayer)
Perkinsus_rtracklayer_transcripts <- Perkinsus_rtracklayer %>% filter(!is.na(transcript_id))
# save
save(Perkinsus_rtracklayer_transcripts, file = "Perkinsus_rtracklayer_transcripts.RData")

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

# Import DEG lists
# perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot, perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot 
load(file = "../2021_Dermo_DEG_annot.RData")
# hemo_dds_deseq_res_Pmar_LFC_sig_APOP, hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP, hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP
load(file = "../2021_Hemocyte_DEG_apop_annot.RData")

## Import hemocyte phenotype data from gene expression data
# PCA_pheno_2020_all
load(file = "../PCA_pheno_2020_all.RData")

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

pdf("./FIGURES/hemo_scale_free_plot.pdf", width = 4, height = 5 )
# Scale-free topology fit index as a function of the soft-thresholding power
plot(hemo_dds_rlog_matrix_sft$fitIndices[,1], -sign(hemo_dds_rlog_matrix_sft$fitIndices[,3])*hemo_dds_rlog_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Hemocyte Scale Independence"));
text(hemo_dds_rlog_matrix_sft$fitIndices[,1], -sign(hemo_dds_rlog_matrix_sft$fitIndices[,3])*hemo_dds_rlog_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()


sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
pdf("./FIGURES/hemo_mean_connectivity_plot.pdf", width = 4, height = 5 )
# Mean connectivity as a function of the soft-thresholding power
plot(hemo_dds_rlog_matrix_sft$fitIndices[,1], hemo_dds_rlog_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Hemocyte Mean Connectivity"))
text(hemo_dds_rlog_matrix_sft$fitIndices[,1], hemo_dds_rlog_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
dev.off()

# Selecting Softthreshold of 7 since this is lowest value past 0.9 we start to see flattening 

## Perkinsus data 
# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
perk_dds_rlog_matrix_sft <- pickSoftThreshold(perk_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
save(perk_dds_rlog_matrix_sft, file="/Volumes/My Passport for Mac/2021_Hemocyte_Dermo_WGCNA/perk_dds_rlog_matrix_sft")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

pdf("./FIGURES/perk_scale_free_plot.pdf", width = 4, height = 5 )
# Scale-free topology fit index as a function of the soft-thresholding power
plot(perk_dds_rlog_matrix_sft$fitIndices[,1], -sign(perk_dds_rlog_matrix_sft$fitIndices[,3])*perk_dds_rlog_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("P. marinus Scale Independence"));
text(perk_dds_rlog_matrix_sft$fitIndices[,1], -sign(perk_dds_rlog_matrix_sft$fitIndices[,3])*perk_dds_rlog_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()
dev.off()

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
pdf("./FIGURES/perk_mean_connectivity_plot.pdf", width = 4, height = 5 )
# Mean connectivity as a function of the soft-thresholding power
plot(perk_dds_rlog_matrix_sft$fitIndices[,1], perk_dds_rlog_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("P. marinus Mean Connectivity"))
text(perk_dds_rlog_matrix_sft$fitIndices[,1], perk_dds_rlog_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
dev.off()
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
pdf("./FIGURES/hemo_dendrogram.pdf", width = 10, height = 5 )
plotDendroAndColors(hemo_full_net$dendrograms[[1]], hemo_full_mergedColors[hemo_full_net$blockGenes[[1]]],
                    "Hemocyte\nModule Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

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
pdf("./FIGURES/perk_dendrogram.pdf", width = 10, height = 5 )
plotDendroAndColors(perk_full_net$dendrograms[[1]], perk_full_mergedColors[perk_full_net$blockGenes[[1]]],
                    "P. marinus\nModule Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

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

#### QUANTIFY MODULE ASSOCIATIONS WITH CHALLENGE ONLY ####

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
  
hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_sig_modnames <- hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare$mod_names

# find those shared between all 
hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_shared <- drop_na(hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare)
# 8 shared between all 
# interesting modules that have an opposite sign between treatments
   # MEyellow
   # MEdarkseagreen4
   # MEplum4

## Heatmap of only the significantly correlated modules 
# Graph and color code each the strength of association (correlation) of module eigengenes and trait

# subset hemo_full_moduleTraitCor, hemo_full_moduleTraitPvalue, hemo_full_MEs for only those modules significant in either challenge
hemo_full_moduleTraitCor_sig <- hemo_full_moduleTraitCor[rownames(hemo_full_moduleTraitCor) %in% hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_sig_modnames,]
hemo_full_moduleTraitCor_sig <- hemo_full_moduleTraitCor_sig[,c(1:3)]
hemo_full_moduleTraitPvalue_sig <- hemo_full_moduleTraitPvalue[rownames(hemo_full_moduleTraitPvalue) %in% hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_sig_modnames,]
hemo_full_moduleTraitPvalue_sig <- hemo_full_moduleTraitPvalue_sig[,c(1:3)]
hemo_full_MEs_sig <- hemo_full_MEs[,colnames(hemo_full_MEs) %in% hemo_full_moduleTraitCor_Pval_df_Pmar_sig_compare_sig_modnames]
hemo_coldata_collapse_binarize_sig <- hemo_coldata_collapse_binarize[,c(1:3)]

# Will display correlations and their p-values
hemo_full_textMatrix_sig = paste(signif(hemo_full_moduleTraitCor_sig, 2), "\n(",
                                 signif(hemo_full_moduleTraitPvalue_sig, 1), ")", sep = "");
dim(hemo_full_textMatrix_sig) = dim(hemo_full_moduleTraitCor_sig)

# make plot
sizeGrWindow(10,6)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = hemo_full_moduleTraitCor_sig,
               xLabels = c("Con vs. Pmar", "P.mar. GDC vs. Con", "P. mar. ZVAD vs. Con"),
               yLabels = names(hemo_full_MEs_sig),
               ySymbols = names(hemo_full_MEs_sig),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = hemo_full_textMatrix_sig,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.7,
               #zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Hemocyte Module-Challenge Relationships"))
# have to save manually...weird!


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

#### ANALYSIS OF DATA FROM MODULE ASSOCIATION WITH CHALLENGE ONLY ####
#### annotate apoptosis genes in significant modules ####

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
verboseScatterplot(abs(hemo_full_geneModuleMembership[hemo_full_moduleGenes, hemo_full_column]),
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

#### Identify hub genes in each module  ####

# Hemocytes control vs Pmar
hemo_full_colorh_control_Pmar = hemo_full_module_apop_df_5_greater_control_Pmar

hemo_full_Module_hub_genes_control_Pmar <- chooseTopHubInEachModule(
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

### Data analysis notes ####

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
moduleColors = hemo_full_moduleColors
lookup =   C_vir_rtracklayer_transcripts
datKME =  hemo_datKME

lookup_annot_intramodular <- function(list) {
module_name = str_remove(list, "MM.")
mod_list <- names(matrix)[moduleColors == module_name]
FilterGenes = abs(GS1)> .6 & abs(datKME[list])>.8
FilterGenes_annot <- data.frame("ID" = dimnames(data.frame(matrix))[[2]][FilterGenes]) %>% left_join(., lookup) %>% mutate(mod_names = list) %>% filter(ID %in% mod_list)
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
FilterGenes_comb_hub_count <- FilterGenes_comb %>% group_by(mod_names, group) %>% dplyr::mutate(total_hub = n()) %>% dplyr::distinct(mod_names, group, total_hub, moduleTraitCor, moduleTraitPvalue)

# are any of these apoptotic?
FilterGenes_comb_apop <- FilterGenes_comb[FilterGenes_comb$ID %in% C_vir_rtracklayer_apop_product_final$ID,] %>% 
  dplyr::select(ID, gene,product, transcript_id, mod_names, group, moduleTraitCor,moduleTraitPvalue, GS)

# how many apop intramodular hub genes in each module per treatment?
FilterGenes_comb_apop_count <- FilterGenes_comb_apop %>% group_by(mod_names, group) %>% dplyr::count() %>% dplyr::rename(apop_hub = n)
FilterGenes_comb_apop_count_mod_names <- FilterGenes_comb_apop_count %>% ungroup() %>% distinct(mod_names) 
FilterGenes_comb_apop_count_mod_names <- FilterGenes_comb_apop_count_mod_names$mod_names

FilterGenes_comb_hub_count_apop <- left_join(FilterGenes_comb_hub_count, FilterGenes_comb_apop_count)

View(FilterGenes_comb_hub_count_apop )

# combine with Interproscan results
FilterGenes_comb_Interpro <- left_join(FilterGenes_comb[,c("ID","product","transcript_id","moduleTraitCor","moduleTraitPvalue","mod_names","GS")], 
          GO_universe_rna_found[,c("ID","source","transcript_id", "Ontology_term",
                                       "Dbxref","signature_desc")], by =  c("ID", "transcript_id"))

## Heatmap of only the significantly correlated modules that have apoptosis hub genes
# Graph and color code each the strength of association (correlation) of module eigengenes and trait

# subset hemo_full_moduleTraitCor, hemo_full_moduleTraitPvalue, hemo_full_MEs for only those modules significant in either challenge
hemo_full_moduleTraitCor_sig_apop <- hemo_full_moduleTraitCor[rownames(hemo_full_moduleTraitCor) %in% FilterGenes_comb_apop_count_mod_names,]
hemo_full_moduleTraitCor_sig_apop <- hemo_full_moduleTraitCor_sig_apop[,c(1:3)]
hemo_full_moduleTraitPvalue_sig_apop <- hemo_full_moduleTraitPvalue[rownames(hemo_full_moduleTraitPvalue) %in% FilterGenes_comb_apop_count_mod_names,]
hemo_full_moduleTraitPvalue_sig_apop <- hemo_full_moduleTraitPvalue_sig_apop[,c(1:3)]
hemo_full_MEs_sig_apop <- hemo_full_MEs[,colnames(hemo_full_MEs) %in% FilterGenes_comb_apop_count_mod_names]
hemo_coldata_collapse_binarize_sig_apop <- hemo_coldata_collapse_binarize[,c(1:3)]

# Will display correlations and their p-values
hemo_full_textMatrix_sig_apop = paste(signif(hemo_full_moduleTraitCor_sig_apop, 2), "\n(",
                                 signif(hemo_full_moduleTraitPvalue_sig_apop, 1), ")", sep = "");
dim(hemo_full_textMatrix_sig_apop) = dim(hemo_full_moduleTraitCor_sig_apop)

# make plot
sizeGrWindow(10,6)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = hemo_full_moduleTraitCor_sig_apop,
               xLabels = c("Con vs. Pmar", "P.mar. GDC vs. Con", "P. mar. ZVAD vs. Con"),
               yLabels = names(hemo_full_MEs_sig_apop),
               ySymbols = names(hemo_full_MEs_sig_apop),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = hemo_full_textMatrix_sig_apop,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.7,
               #zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Hemocyte Module-Challenge Relationships"))
# have to save manually...weird!

### Repeat finding of intramodular hub genes, but for P. marinus samples now 
# annotate intramodular hub genes in each module
# Use function to lookup all apop names for each significant module
matrix = perk_dds_rlog_matrix
lookup =   Perkinsus_rtracklayer
datKME =  perk_datKME
moduleColors = perk_full_moduleColors

lookup_annot_intramodular_perk <- function(list) {
  module_name = str_remove(list, "MM.")
  mod_list <- names(matrix)[moduleColors == module_name]
  FilterGenes = abs(GS1)> .6 & abs(datKME[list])>.8
  FilterGenes_annot <- data.frame("Name" = dimnames(data.frame(matrix))[[2]][FilterGenes]) %>% left_join(., lookup) %>% mutate(mod_names = list) %>% filter(Name %in% mod_list)
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

#### Hemocyte intramodular hub gene analysis for interesting modules ####


### ASSESS OVERLAPS WITH DEG RESULTS

# find the hits for each treatment group and then count the number of hits in each modules
hemo_dds_deseq_res_Pmar_LFC_sig_APOP_sig_modules <- FilterGenes_comb[FilterGenes_comb$transcript_id %in% hemo_dds_deseq_res_Pmar_LFC_sig_APOP$transcript_id,] %>% 
  filter(group == "control_Pmar") %>% 
  group_by(mod_names) %>% dplyr::mutate(module_count = n()) %>% dplyr::select(product, transcript_id, mod_names, group, moduleTraitCor, moduleTraitPvalue, GS, module_count)

hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_sig_modules <- FilterGenes_comb[FilterGenes_comb$transcript_id %in% hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP$transcript_id,] %>% 
  filter(group == "control_Pmar_GDC") %>% 
  group_by(mod_names) %>% dplyr::mutate(module_count = n()) %>% dplyr::select(product, transcript_id, mod_names, group, moduleTraitCor, moduleTraitPvalue, GS, module_count)

hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_sig_modules <- FilterGenes_comb[FilterGenes_comb$transcript_id %in% hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP$transcript_id,] %>% 
  filter(group == "control_Pmar_ZVAD") %>% 
  group_by(mod_names) %>% dplyr::mutate(module_count = n()) %>% dplyr::select(product, transcript_id, mod_names, group, moduleTraitCor, moduleTraitPvalue, GS, module_count)


# any apoptosis transcripts found in every GDC responsive module ? (black, navajowhite2, darkred, lightpink3)
mod_list <- c("MEblack", "MEnavajowhite2", "MEdarkred", "MElightpink3")
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_sig_modules_sig_shared_important <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_sig_modules %>% filter(mod_names %in% mod_list) %>% 
distinct(product, mod_names, transcript_id) %>% ungroup() %>% 
group_by(transcript_id) %>% dplyr::mutate(count = n()) %>% filter(count == 4) %>% distinct(product, transcript_id)
# 0

hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_sig_modules_sig_shared <- hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_sig_modules %>% filter(mod_names %in% mod_list) %>% 
  distinct(product, mod_names, transcript_id) %>% ungroup() %>% 
  group_by(transcript_id) %>% dplyr::mutate(count = n()) 

# how many overlaps?
hemo_dds_deseq_res_Pmar_LFC_sig_APOP_sig_modules %>% distinct(module_count, mod_names, group) %>% arrange(desc(module_count))
# 0
hemo_dds_deseq_res_Pmar_GDC_LFC_sig_APOP_sig_modules %>% distinct(module_count, mod_names, group) %>% arrange(desc(module_count)) # the GDC 
 #mod_names       group            module_count
 #<chr>           <chr>                   <int>
 #  1 MEyellow        control_Pmar_GDC           14
 #2 MEnavajowhite2  control_Pmar_GDC            5
 #3 MElightpink3    control_Pmar_GDC            2
 #4 MEcyan          control_Pmar_GDC            1
 #5 MEdarkred       control_Pmar_GDC            1
 #6 MEpaleturquoise control_Pmar_GDC            1
 #7 MEblack         control_Pmar_GDC            1
 #8 MEdarkorange    control_Pmar_GDC            1
hemo_dds_deseq_res_Pmar_ZVAD_LFC_sig_APOP_sig_modules %>% distinct(module_count, mod_names, group) %>% arrange(desc(module_count))
    #mod_names    group             module_count
    #<chr>        <chr>                    <int>
    #1 MEyellow     control_Pmar_ZVAD            2
    #2 MEorangered4 control_Pmar_ZVAD            1


#### Hemocyte GO enrichment for interesting module all transcripts  ####
## GO enrichment for all genes in the module
# get full list of genes for each module 
hemo_geneNames <- names(GO_universe_rna_found_geneID2GO_mapping)

FilterGenes_comb_apop_count_mod_names

hemo_MEantiquewhite2 <- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors =="steelblue"]
hemo_MEblack   <- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "black"]
hemo_MEblue   <- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "blue"]
hemo_MEcyan   <- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "cyan"]
hemo_MEdarkgreen  <- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "darkgreen"]
hemo_MEdarkorange  <- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "darkorange"]
hemo_MEdarkred  <- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "darkred"]
hemo_MEdarkseagreen4<- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "darkseagreen4"]
hemo_MEdarkslateblue<- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "darkslateblue"]
hemo_MElightcyan<- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "lightcyan"]
hemo_MElightpink3<- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "lightpink3"]
hemo_MEmediumpurple3 <- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "mediumpurple3"]
hemo_MEnavajowhite2 <- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "navajowhite2"]
hemo_MEorangered4<- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "orangered4"]
hemo_MEpaleturquoise<- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "paleturquoise"]
hemo_MEplum<- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "plum"]
hemo_MEred<- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "red"]
hemo_MEyellow<- names(hemo_dds_rlog_matrix)[hemo_full_moduleColors == "yellow"]

hemo_MEantiquewhite2 <- as.factor(hemo_MEantiquewhite2[!is.na(hemo_MEantiquewhite2)])
hemo_MEblack <- as.factor(hemo_MEblack[!is.na(hemo_MEblack)])
hemo_MEblue  <- as.factor(hemo_MEblue [!is.na(hemo_MEblue )])
hemo_MEcyan  <- as.factor(hemo_MEcyan [!is.na(hemo_MEcyan )])
hemo_MEdarkgreen <- as.factor(hemo_MEdarkgreen[!is.na(hemo_MEdarkgreen)])
hemo_MEdarkorange <- as.factor(hemo_MEdarkorange[!is.na(hemo_MEdarkorange)])
hemo_MEdarkred <- as.factor(hemo_MEdarkred[!is.na(hemo_MEdarkred)])
hemo_MEdarkseagreen4 <- as.factor(hemo_MEdarkseagreen4[!is.na(hemo_MEdarkseagreen4)])
hemo_MEdarkslateblue <- as.factor(hemo_MEdarkslateblue[!is.na(hemo_MEdarkslateblue)])
hemo_MElightcyan <- as.factor(hemo_MElightcyan[!is.na(hemo_MElightcyan)])
hemo_MElightpink3 <- as.factor(hemo_MElightpink3[!is.na(hemo_MElightpink3)])
hemo_MEmediumpurple3 <- as.factor(hemo_MEmediumpurple3[!is.na(hemo_MEmediumpurple3)])
hemo_MEnavajowhite2 <- as.factor(hemo_MEnavajowhite2[!is.na(hemo_MEnavajowhite2)])
hemo_MEorangered4 <- as.factor(hemo_MEorangered4[!is.na(hemo_MEorangered4)])
hemo_MEpaleturquoise <- as.factor(hemo_MEpaleturquoise[!is.na(hemo_MEpaleturquoise)])
hemo_MEplum <- as.factor(hemo_MEplum[!is.na(hemo_MEplum)])
hemo_MEred <- as.factor(hemo_MEred[!is.na(hemo_MEred)])
hemo_MEyellow <- as.factor(hemo_MEyellow[!is.na(hemo_MEyellow)])

hemo_MEantiquewhite2_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEantiquewhite2))
hemo_MEblack_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEblack))
hemo_MEblue_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEblue))
hemo_MEcyan_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEcyan))
hemo_MEdarkgreen_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEdarkgreen))
hemo_MEdarkorange_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEdarkorange))
hemo_MEdarkred_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEdarkred))
hemo_MEdarkseagreen4_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEdarkseagreen4))
hemo_MEdarkslateblue_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEdarkslateblue))
hemo_MElightcyan_factor <- factor(as.integer(hemo_geneNames %in% hemo_MElightcyan))
hemo_MElightpink3_factor <- factor(as.integer(hemo_geneNames %in% hemo_MElightpink3))
hemo_MEmediumpurple3_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEmediumpurple3))
hemo_MEnavajowhite2_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEnavajowhite2))
hemo_MEorangered4_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEorangered4))
hemo_MEpaleturquoise_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEpaleturquoise))
hemo_MEplum_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEplum))
hemo_MEred_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEred))
hemo_MEyellow_factor <- factor(as.integer(hemo_geneNames %in% hemo_MEyellow))

names(hemo_MEantiquewhite2_factor) <- hemo_geneNames
names(hemo_MEblack_factor) <- hemo_geneNames
names(hemo_MEblue_factor) <- hemo_geneNames
names(hemo_MEcyan_factor) <- hemo_geneNames
names(hemo_MEdarkgreen_factor) <- hemo_geneNames
names(hemo_MEdarkorange_factor) <- hemo_geneNames
names(hemo_MEdarkred_factor) <- hemo_geneNames
names(hemo_MEdarkseagreen4_factor) <- hemo_geneNames
names(hemo_MEdarkslateblue_factor) <- hemo_geneNames
names(hemo_MElightcyan_factor) <- hemo_geneNames
names(hemo_MElightpink3_factor) <- hemo_geneNames
names(hemo_MEmediumpurple3_factor) <- hemo_geneNames
names(hemo_MEnavajowhite2_factor) <- hemo_geneNames
names(hemo_MEorangered4_factor) <- hemo_geneNames
names(hemo_MEpaleturquoise_factor) <- hemo_geneNames
names(hemo_MEplum_factor) <- hemo_geneNames
names(hemo_MEred_factor) <- hemo_geneNames
names(hemo_MEyellow_factor) <- hemo_geneNames

### Make topGO data object 
hemo_MEantiquewhite2_GOdata <- new("topGOdata", description = "gene enrichment",
                                   ontology = "MF",
                                   allGenes = hemo_MEantiquewhite2_factor,
                                   nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEblack_GOdata <- new("topGOdata", description = "gene enrichment",
                           ontology = "MF",
                           allGenes = hemo_MEblack_factor,
                           nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEblue_GOdata <- new("topGOdata", description = "gene enrichment",
                          ontology = "MF",
                          allGenes = hemo_MEblue_factor,
                          nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEcyan_GOdata <- new("topGOdata", description = "gene enrichment",
                          ontology = "MF",
                          allGenes = hemo_MEcyan_factor,
                          nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEdarkgreen_GOdata <- new("topGOdata", description ="gene enrichment",
                               ontology = "MF",
                               allGenes = hemo_MEdarkgreen_factor,
                               nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEdarkorange_GOdata <- new("topGOdata", description ="gene enrichment",
                                ontology = "MF",
                                allGenes = hemo_MEdarkorange_factor,
                                nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEdarkred_GOdata <- new("topGOdata", description ="gene enrichment",
                             ontology = "MF",
                             allGenes = hemo_MEdarkred_factor,
                             nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEdarkseagreen4_GOdata <- new("topGOdata", description ="gene enrichment",
                                   ontology = "MF",
                                   allGenes = hemo_MEdarkseagreen4_factor,
                                   nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEdarkslateblue_GOdata <- new("topGOdata", description ="gene enrichment",
                                   ontology = "MF",
                                   allGenes = hemo_MEdarkslateblue_factor,
                                   nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MElightcyan_GOdata <- new("topGOdata", description ="gene enrichment",
                               ontology = "MF",
                               allGenes = hemo_MElightcyan_factor,
                               nodeSize = 5, annot = annFUN.gene2GO, gene2GO =GO_universe_rna_found_geneID2GO_mapping )
hemo_MElightpink3_GOdata <- new("topGOdata", description ="gene enrichment",
                                ontology = "MF",
                                allGenes = hemo_MElightpink3_factor,
                                nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEmediumpurple3_GOdata <- new("topGOdata", description ="gene enrichment",
                                   ontology = "MF",
                                   allGenes = hemo_MEmediumpurple3_factor,
                                   nodeSize = 5, annot = annFUN.gene2GO, gene2GO =GO_universe_rna_found_geneID2GO_mapping )
hemo_MEnavajowhite2_GOdata <- new("topGOdata", description= "gene enrichment",
                                  ontology = "MF",
                                  allGenes = hemo_MEnavajowhite2_factor,
                                  nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEorangered4_GOdata <- new("topGOdata", description ="gene enrichment",
                                ontology = "MF",
                                allGenes = hemo_MEorangered4_factor,
                                nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEpaleturquoise_GOdata <- new("topGOdata", description ="gene enrichment",
                                   ontology = "MF",
                                   allGenes = hemo_MEpaleturquoise_factor,
                                   nodeSize = 5, annot = annFUN.gene2GO, gene2GO =GO_universe_rna_found_geneID2GO_mapping)
hemo_MEplum_GOdata <- new("topGOdata", description ="gene enrichment",
                                 ontology = "MF",
                                 allGenes = hemo_MEplum_factor,
                                 nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEred_GOdata <- new("topGOdata", description = "gene enrichment",
                                ontology = "MF",
                                allGenes = hemo_MEred_factor,
                                nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEyellow_GOdata <- new("topGOdata", description ="gene enrichment",
                            ontology = "MF",
                            allGenes = hemo_MEyellow_factor,
                            nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)

#nodeSize=used to prune the GO hierarchy from the terms which have less than 1 annotated genes
#annFUN.gene2GO = this function is used when the annotations are provided as a gene-to-GOs mapping.

### Perform Encrichment tests 
hemo_MEantiquewhite2_GOdata_Fisher_Weight <- runTest(hemo_MEantiquewhite2_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEblack_GOdata_Fisher_Weight <- runTest(hemo_MEblack_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEblue_GOdata_Fisher_Weight <- runTest(hemo_MEblue_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEcyan_GOdata_Fisher_Weight <- runTest(hemo_MEcyan_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEdarkgreen_GOdata_Fisher_Weight <- runTest(hemo_MEdarkgreen_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEdarkorange_GOdata_Fisher_Weight <- runTest(hemo_MEdarkorange_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEdarkred_GOdata_Fisher_Weight <- runTest(hemo_MEdarkred_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEdarkseagreen4_GOdata_Fisher_Weight <- runTest(hemo_MEdarkseagreen4_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEdarkslateblue_GOdata_Fisher_Weight <- runTest(hemo_MEdarkslateblue_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MElightcyan_GOdata_Fisher_Weight <- runTest(hemo_MElightcyan_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MElightpink3_GOdata_Fisher_Weight <- runTest(hemo_MElightpink3_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEmediumpurple3_GOdata_Fisher_Weight <- runTest(hemo_MEmediumpurple3_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEnavajowhite2_GOdata_Fisher_Weight <- runTest(hemo_MEnavajowhite2_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEorangered4_GOdata_Fisher_Weight <- runTest(hemo_MEorangered4_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEpaleturquoise_GOdata_Fisher_Weight <- runTest(hemo_MEpaleturquoise_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEplum_GOdata_Fisher_Weight <- runTest(hemo_MEplum_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEred_GOdata_Fisher_Weight <- runTest(hemo_MEred_GOdata, algorithm = "weight01", statistic = "fisher")
hemo_MEyellow_GOdata_Fisher_Weight <- runTest(hemo_MEyellow_GOdata, algorithm = "weight01", statistic = "fisher")

## Analyze enrichment test results 
# see how many results we get where weight01 gives a P-value <= 0.05
hemo_MEantiquewhite2_GOdata_summary <- summary(attributes(hemo_MEantiquewhite2_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEblack_GOdata_summary <- summary(attributes(hemo_MEblack_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEblue_GOdata_summary <- summary(attributes(hemo_MEblue_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEcyan_GOdata_summary <- summary(attributes(hemo_MEcyan_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEdarkgreen_GOdata_summary <- summary(attributes(hemo_MEdarkgreen_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEdarkorange_GOdata_summary <- summary(attributes(hemo_MEdarkorange_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEdarkred_GOdata_summary <- summary(attributes(hemo_MEdarkred_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEdarkseagreen4_GOdata_summary <- summary(attributes(hemo_MEdarkseagreen4_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEdarkslateblue_GOdata_summary <- summary(attributes(hemo_MEdarkslateblue_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MElightcyan_GOdata_summary <- summary(attributes(hemo_MElightcyan_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MElightpink3_GOdata_summary <- summary(attributes(hemo_MElightpink3_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEmediumpurple3_GOdata_summary <- summary(attributes(hemo_MEmediumpurple3_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEnavajowhite2_GOdata_summary <- summary(attributes(hemo_MEnavajowhite2_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEorangered4_GOdata_summary <- summary(attributes(hemo_MEorangered4_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEpaleturquoise_GOdata_summary <- summary(attributes(hemo_MEpaleturquoise_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEplum_GOdata_summary <- summary(attributes(hemo_MEplum_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEred_GOdata_summary <- summary(attributes(hemo_MEred_GOdata_Fisher_Weight)$score <= 0.05)
hemo_MEyellow_GOdata_summary <- summary(attributes(hemo_MEyellow_GOdata_Fisher_Weight)$score <= 0.05)

#print out the top results, though only GDC_lightblue4 is sig
hemo_MEantiquewhite2_GOdata_Res <- GenTable(hemo_MEantiquewhite2_GOdata, topgoFisher = hemo_MEantiquewhite2_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEblack_GOdata_Res <- GenTable(hemo_MEblack_GOdata, topgoFisher = hemo_MEblack_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEblue_GOdata_Res <- GenTable(hemo_MEblue_GOdata, topgoFisher = hemo_MEblue_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEcyan_GOdata_Res <- GenTable(hemo_MEcyan_GOdata, topgoFisher = hemo_MEcyan_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEdarkgreen_GOdata_Res <- GenTable(hemo_MEdarkgreen_GOdata, topgoFisher = hemo_MEdarkgreen_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEdarkorange_GOdata_Res <- GenTable(hemo_MEdarkorange_GOdata, topgoFisher = hemo_MEdarkorange_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEdarkred_GOdata_Res <- GenTable(hemo_MEdarkred_GOdata, topgoFisher = hemo_MEdarkred_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEdarkseagreen4_GOdata_Res <- GenTable(hemo_MEdarkseagreen4_GOdata, topgoFisher = hemo_MEdarkseagreen4_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEdarkslateblue_GOdata_Res <- GenTable(hemo_MEdarkslateblue_GOdata, topgoFisher = hemo_MEdarkslateblue_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MElightcyan_GOdata_Res <- GenTable(hemo_MElightcyan_GOdata, topgoFisher = hemo_MElightcyan_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MElightpink3_GOdata_Res <- GenTable(hemo_MElightpink3_GOdata, topgoFisher = hemo_MElightpink3_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEmediumpurple3_GOdata_Res <- GenTable(hemo_MEmediumpurple3_GOdata, topgoFisher = hemo_MEmediumpurple3_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEnavajowhite2_GOdata_Res <- GenTable(hemo_MEnavajowhite2_GOdata, topgoFisher = hemo_MEnavajowhite2_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEorangered4_GOdata_Res <- GenTable(hemo_MEorangered4_GOdata, topgoFisher = hemo_MEorangered4_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEpaleturquoise_GOdata_Res <- GenTable(hemo_MEpaleturquoise_GOdata, topgoFisher = hemo_MEpaleturquoise_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEplum_GOdata_Res <- GenTable(hemo_MEplum_GOdata, topgoFisher = hemo_MEplum_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEred_GOdata_Res <- GenTable(hemo_MEred_GOdata, topgoFisher = hemo_MEred_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEyellow_GOdata_Res <- GenTable(hemo_MEyellow_GOdata, topgoFisher = hemo_MEyellow_GOdata_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)

## Repeat for BP
### Make topGO data object 
hemo_MEantiquewhite2_GOdata_BP <- new("topGOdata", description = "gene enrichment",
                                   ontology = "BP",
                                   allGenes = hemo_MEantiquewhite2_factor,
                                   nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEblack_GOdata_BP <- new("topGOdata", description = "gene enrichment",
                           ontology = "BP",
                           allGenes = hemo_MEblack_factor,
                           nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEblue_GOdata_BP <- new("topGOdata", description = "gene enrichment",
                          ontology = "BP",
                          allGenes = hemo_MEblue_factor,
                          nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEcyan_GOdata_BP <- new("topGOdata", description = "gene enrichment",
                          ontology = "BP",
                          allGenes = hemo_MEcyan_factor,
                          nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEdarkgreen_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                               ontology = "BP",
                               allGenes = hemo_MEdarkgreen_factor,
                               nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEdarkorange_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                                ontology = "BP",
                                allGenes = hemo_MEdarkorange_factor,
                                nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEdarkred_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                             ontology = "BP",
                             allGenes = hemo_MEdarkred_factor,
                             nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEdarkseagreen4_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                                   ontology = "BP",
                                   allGenes = hemo_MEdarkseagreen4_factor,
                                   nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEdarkslateblue_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                                   ontology = "BP",
                                   allGenes = hemo_MEdarkslateblue_factor,
                                   nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MElightcyan_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                               ontology = "BP",
                               allGenes = hemo_MElightcyan_factor,
                               nodeSize = 5, annot = annFUN.gene2GO, gene2GO =GO_universe_rna_found_geneID2GO_mapping )
hemo_MElightpink3_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                                ontology = "BP",
                                allGenes = hemo_MElightpink3_factor,
                                nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEmediumpurple3_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                                   ontology = "BP",
                                   allGenes = hemo_MEmediumpurple3_factor,
                                   nodeSize = 5, annot = annFUN.gene2GO, gene2GO =GO_universe_rna_found_geneID2GO_mapping )
hemo_MEnavajowhite2_GOdata_BP <- new("topGOdata", description= "gene enrichment",
                                  ontology = "BP",
                                  allGenes = hemo_MEnavajowhite2_factor,
                                  nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEorangered4_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                                ontology = "BP",
                                allGenes = hemo_MEorangered4_factor,
                                nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEpaleturquoise_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                                   ontology = "BP",
                                   allGenes = hemo_MEpaleturquoise_factor,
                                   nodeSize = 5, annot = annFUN.gene2GO, gene2GO =GO_universe_rna_found_geneID2GO_mapping)
hemo_MEplum_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                          ontology = "BP",
                          allGenes = hemo_MEplum_factor,
                          nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEred_GOdata_BP <- new("topGOdata", description = "gene enrichment",
                         ontology = "BP",
                         allGenes = hemo_MEred_factor,
                         nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)
hemo_MEyellow_GOdata_BP <- new("topGOdata", description ="gene enrichment",
                            ontology = "BP",
                            allGenes = hemo_MEyellow_factor,
                            nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GO_universe_rna_found_geneID2GO_mapping)

#nodeSize=used to prune the GO hierarchy from the terms which have less than 1 annotated genes
#annFUN.gene2GO = this function is used when the annotations are provided as a gene-to-GOs mapping.

### Perform Encrichment tests 
hemo_MEantiquewhite2_GOdata_BP_Fisher_Weight <- runTest(hemo_MEantiquewhite2_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEblack_GOdata_BP_Fisher_Weight <- runTest(hemo_MEblack_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEblue_GOdata_BP_Fisher_Weight <- runTest(hemo_MEblue_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEcyan_GOdata_BP_Fisher_Weight <- runTest(hemo_MEcyan_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEdarkgreen_GOdata_BP_Fisher_Weight <- runTest(hemo_MEdarkgreen_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEdarkorange_GOdata_BP_Fisher_Weight <- runTest(hemo_MEdarkorange_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEdarkred_GOdata_BP_Fisher_Weight <- runTest(hemo_MEdarkred_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEdarkseagreen4_GOdata_BP_Fisher_Weight <- runTest(hemo_MEdarkseagreen4_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEdarkslateblue_GOdata_BP_Fisher_Weight <- runTest(hemo_MEdarkslateblue_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MElightcyan_GOdata_BP_Fisher_Weight <- runTest(hemo_MElightcyan_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MElightpink3_GOdata_BP_Fisher_Weight <- runTest(hemo_MElightpink3_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEmediumpurple3_GOdata_BP_Fisher_Weight <- runTest(hemo_MEmediumpurple3_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEnavajowhite2_GOdata_BP_Fisher_Weight <- runTest(hemo_MEnavajowhite2_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEorangered4_GOdata_BP_Fisher_Weight <- runTest(hemo_MEorangered4_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEpaleturquoise_GOdata_BP_Fisher_Weight <- runTest(hemo_MEpaleturquoise_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEplum_GOdata_BP_Fisher_Weight <- runTest(hemo_MEplum_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEred_GOdata_BP_Fisher_Weight <- runTest(hemo_MEred_GOdata_BP, algorithm = "weight01", statistic = "fisher")
hemo_MEyellow_GOdata_BP_Fisher_Weight <- runTest(hemo_MEyellow_GOdata_BP, algorithm = "weight01", statistic = "fisher")

## Analyze enrichment test results 
# see how many results we get where weight01 gives a P-value <= 0.05
hemo_MEantiquewhite2_GOdata_BP_summary <- summary(attributes(hemo_MEantiquewhite2_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEblack_GOdata_BP_summary <- summary(attributes(hemo_MEblack_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEblue_GOdata_BP_summary <- summary(attributes(hemo_MEblue_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEcyan_GOdata_BP_summary <- summary(attributes(hemo_MEcyan_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEdarkgreen_GOdata_BP_summary <- summary(attributes(hemo_MEdarkgreen_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEdarkorange_GOdata_BP_summary <- summary(attributes(hemo_MEdarkorange_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEdarkred_GOdata_BP_summary <- summary(attributes(hemo_MEdarkred_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEdarkseagreen4_GOdata_BP_summary <- summary(attributes(hemo_MEdarkseagreen4_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEdarkslateblue_GOdata_BP_summary <- summary(attributes(hemo_MEdarkslateblue_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MElightcyan_GOdata_BP_summary <- summary(attributes(hemo_MElightcyan_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MElightpink3_GOdata_BP_summary <- summary(attributes(hemo_MElightpink3_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEmediumpurple3_GOdata_BP_summary <- summary(attributes(hemo_MEmediumpurple3_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEnavajowhite2_GOdata_BP_summary <- summary(attributes(hemo_MEnavajowhite2_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEorangered4_GOdata_BP_summary <- summary(attributes(hemo_MEorangered4_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEpaleturquoise_GOdata_BP_summary <- summary(attributes(hemo_MEpaleturquoise_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEplum_GOdata_BP_summary <- summary(attributes(hemo_MEplum_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEred_GOdata_BP_summary <- summary(attributes(hemo_MEred_GOdata_BP_Fisher_Weight)$score <= 0.05)
hemo_MEyellow_GOdata_BP_summary <- summary(attributes(hemo_MEyellow_GOdata_BP_Fisher_Weight)$score <= 0.05)

#print out the top results
hemo_MEantiquewhite2_GOdata_BP_Res <- GenTable(hemo_MEantiquewhite2_GOdata_BP, topgoFisher = hemo_MEantiquewhite2_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEblack_GOdata_BP_Res <- GenTable(hemo_MEblack_GOdata_BP, topgoFisher = hemo_MEblack_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEblue_GOdata_BP_Res <- GenTable(hemo_MEblue_GOdata_BP, topgoFisher = hemo_MEblue_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEcyan_GOdata_BP_Res <- GenTable(hemo_MEcyan_GOdata_BP, topgoFisher = hemo_MEcyan_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEdarkgreen_GOdata_BP_Res <- GenTable(hemo_MEdarkgreen_GOdata_BP, topgoFisher = hemo_MEdarkgreen_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEdarkorange_GOdata_BP_Res <- GenTable(hemo_MEdarkorange_GOdata_BP, topgoFisher = hemo_MEdarkorange_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEdarkred_GOdata_BP_Res <- GenTable(hemo_MEdarkred_GOdata_BP, topgoFisher = hemo_MEdarkred_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEdarkseagreen4_GOdata_BP_Res <- GenTable(hemo_MEdarkseagreen4_GOdata_BP, topgoFisher = hemo_MEdarkseagreen4_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEdarkslateblue_GOdata_BP_Res <- GenTable(hemo_MEdarkslateblue_GOdata_BP, topgoFisher = hemo_MEdarkslateblue_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MElightcyan_GOdata_BP_Res <- GenTable(hemo_MElightcyan_GOdata_BP, topgoFisher = hemo_MElightcyan_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MElightpink3_GOdata_BP_Res <- GenTable(hemo_MElightpink3_GOdata_BP, topgoFisher = hemo_MElightpink3_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEmediumpurple3_GOdata_BP_Res <- GenTable(hemo_MEmediumpurple3_GOdata_BP, topgoFisher = hemo_MEmediumpurple3_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEnavajowhite2_GOdata_BP_Res <- GenTable(hemo_MEnavajowhite2_GOdata_BP, topgoFisher = hemo_MEnavajowhite2_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEorangered4_GOdata_BP_Res <- GenTable(hemo_MEorangered4_GOdata_BP, topgoFisher = hemo_MEorangered4_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEpaleturquoise_GOdata_BP_Res <- GenTable(hemo_MEpaleturquoise_GOdata_BP, topgoFisher = hemo_MEpaleturquoise_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEplum_GOdata_BP_Res <- GenTable(hemo_MEplum_GOdata_BP, topgoFisher = hemo_MEplum_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEred_GOdata_BP_Res <- GenTable(hemo_MEred_GOdata_BP, topgoFisher = hemo_MEred_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)
hemo_MEyellow_GOdata_BP_Res <- GenTable(hemo_MEyellow_GOdata_BP, topgoFisher = hemo_MEyellow_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher",topNodes = 20)

#### Hemocyte top interesting modules GO visualization ####

### Dotplot of significantly enriched MF GO terms from hub modules
hemo_MEantiquewhite2_GOdata_Res$group <- "antiquewhite2"
hemo_MEblack_GOdata_Res$group <- "black"
hemo_MEblue_GOdata_Res$group <- "blue"
hemo_MEcyan_GOdata_Res$group <- "cyan"
hemo_MEdarkgreen_GOdata_Res$group <- "darkgreen"
hemo_MEdarkorange_GOdata_Res$group <- "darkorange"
hemo_MEdarkred_GOdata_Res$group <- "darkred"
hemo_MEdarkseagreen4_GOdata_Res$group <- "darkseagreen4"
hemo_MEdarkslateblue_GOdata_Res$group <- "darkslateblue"
hemo_MElightcyan_GOdata_Res$group <- "lightcyan"
hemo_MElightpink3_GOdata_Res$group <- "lightpink3"
hemo_MEmediumpurple3_GOdata_Res$group <- "mediumpurple3"
hemo_MEnavajowhite2_GOdata_Res$group <- "navajowhite2"
hemo_MEorangered4_GOdata_Res$group <- "orangered4"
hemo_MEpaleturquoise_GOdata_Res$group <- "paleturquoise"
hemo_MEplum_GOdata_Res$group <- "plum"
hemo_MEred_GOdata_Res$group <- "red"
hemo_MEyellow_GOdata_Res$group <- "yellow"

hemo_MEred_GOdata_Res$type <- "Pmar_sig_all_same"
hemo_MElightcyan_GOdata_Res$type <- "Pmar_sig_all_same"
hemo_MEantiquewhite2_GOdata_Res$type <- "Pmar_sig_all_same"
hemo_MEyellow_GOdata_Res$type <- "Pmar_sig_all_diff"
hemo_MEdarkseagreen4_GOdata_Res$type <- "Pmar_sig_all_diff"

hemo_MEnavajowhite2_GOdata_Res$type <- "GDC_and_ZVAD_diff"
hemo_MEplum_GOdata_Res$type <- "GDC_and_ZVAD_down"
hemo_MEcyan_GOdata_Res$type <- "GDC_and_ZVAD_down"
hemo_MEmediumpurple3_GOdata_Res$type <- "GDC_and_ZVAD_down"

hemo_MEpaleturquoise_GOdata_Res$type <- "GDC_only"
hemo_MEdarkorange_GOdata_Res$type <- "GDC_only"
hemo_MEdarkgreen_GOdata_Res$type <- "GDC_only"
hemo_MEdarkred_GOdata_Res$type <- "GDC_only"
hemo_MElightpink3_GOdata_Res$type <- "GDC_only"
hemo_MEblack_GOdata_Res$type <- "GDC_only"

hemo_MEblue_GOdata_Res$type <- "control_Pmar"
hemo_MEdarkslateblue_GOdata_Res$type <- "ZVAD_only"
hemo_MEorangered4_GOdata_Res$type <- "ZVAD_only"

Hemo_GO_all_dotplot <- rbind(hemo_MEantiquewhite2_GOdata_Res,
                             hemo_MEblack_GOdata_Res ,
                             hemo_MEblue_GOdata_Res,
                             hemo_MEcyan_GOdata_Res,
                             hemo_MEdarkgreen_GOdata_Res,
                             hemo_MEdarkorange_GOdata_Res,
                             hemo_MEdarkred_GOdata_Res,
                             hemo_MEdarkseagreen4_GOdata_Res,
                             hemo_MEdarkslateblue_GOdata_Res,
                             hemo_MElightcyan_GOdata_Res,
                             hemo_MElightpink3_GOdata_Res,
                             hemo_MEmediumpurple3_GOdata_Res,
                             hemo_MEnavajowhite2_GOdata_Res,
                             hemo_MEorangered4_GOdata_Res,
                             hemo_MEpaleturquoise_GOdata_Res,
                             hemo_MEplum_GOdata_Res,
                             hemo_MEred_GOdata_Res,
                             hemo_MEyellow_GOdata_Res) %>% filter(topgoFisher <=0.05 & Significant >5) 

Hemo_GO_all_dotplot_plot <- ggplot(Hemo_GO_all_dotplot, aes(x = group, y = Term )) +
  geom_point(aes(size = Significant, color = as.numeric(topgoFisher))) + 
  scale_size_continuous(range = c(4,10)) +
  scale_color_viridis(option = "viridis", name = "p-value", direction = -1) + 
  facet_grid(.~type, scales = "free_x") + 
  theme_minimal() +
  labs(x = "Module Name", y = "GO Term", title = "GO Enrichment of Significant Modules") + 
  theme(panel.border = element_rect(color = "black", fill = "NA"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text.x = element_text(size = 16, face = "bold"),
        title = element_text(size = 16))

ggsave(Hemo_GO_all_dotplot_plot, device = "tiff", path = "./FIGURES/", 
       filename = "Hemo_GO_all_dotplot_plot.tiff", width = 40, height = 10, limitsize = FALSE)

### Plot as heatmap to help show overlaps
Hemo_GO_all_dotplot_heatmap <- Hemo_GO_all_dotplot %>% mutate(GO_term = paste(GO.ID, Term, sep = "_")) %>% 
  dplyr::select(GO_term, Significant, group) %>% spread(group, Significant) %>%  mutate_if(is.numeric , replace_na, replace = 0) %>% 
  column_to_rownames(., var= "GO_term") 
Hemo_GO_all_dotplot_heatmap <- as.matrix(Hemo_GO_all_dotplot_heatmap)

hemo_labels =c("antiquewhite2", "black"    ,     "blue"       ,   "cyan"         , "darkgreen"   ,  "darkorange" ,   "darkred"      , "darkseagreen4",
                "darkslateblue", "lightcyan",     "lightpink3" ,   "mediumpurple3", "navajowhite2",  "orangered4" ,   "paleturquoise", "plum"         ,
                "red"          , "yellow")
# create named vector to hold column names
hemo_column_labels = structure(paste0(hemo_labels), names = paste0(colnames(Hemo_GO_all_dotplot_heatmap )))

# set custom colors
library(circlize)
col_fun = colorRamp2(c(0, 50,100,400), c("white", "red", "orange", "yellow"))

pdf("./FIGURES/hemo_GO_all_heatmap.pdf", width = 12, height = 10)
hemo_GO_all_heatmap <- ComplexHeatmap::Heatmap(Hemo_GO_all_dotplot_heatmap, border = TRUE, 
                                         column_title_side = "bottom", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                         row_title = "GO Term", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                         row_dend_width = unit(2, "cm"),
                                         row_labels = row.names(Hemo_GO_all_dotplot_heatmap),
                                         column_labels = hemo_column_labels[colnames(Hemo_GO_all_dotplot_heatmap)],
                                         # apply split by k-meams clustering to highlight groups
                                         row_km = 3, column_km = 4, row_names_gp = gpar(fontsize = 8),
                                         column_names_gp = gpar(fontsize = 10),
                                         heatmap_legend_param = list(title = "Number Significant"),
                                         col=col_fun,
                                         rect_gp = gpar(col = "black", lwd = 1))
ComplexHeatmap::draw(hemo_GO_all_heatmap, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 100), "mm")) #bottom, left, top, right paddings
dev.off()

## Heatmap with just the interesting modules identified in my powerpoint (those that also have significant apoptosis phenotype where relevant)
Hemo_GO_all_dotplot_subset <- rbind(
  hemo_MEnavajowhite2_GOdata_Res,
  hemo_MEblack_GOdata_Res ,
  hemo_MEdarkred_GOdata_Res,
  hemo_MElightpink3_GOdata_Res,
  hemo_MEdarkgreen_GOdata_Res,
  hemo_MEblue_GOdata_Res,
  hemo_MEdarkseagreen4_GOdata_Res,
  hemo_MEdarkslateblue_GOdata_Res,
  hemo_MEorangered4_GOdata_Res) %>% filter(topgoFisher <=0.05 & Significant >1) 

### Plot as heatmap to help show overlaps
Hemo_GO_all_dotplot_subset_heatmap <- Hemo_GO_all_dotplot_subset %>% mutate(GO_term = paste(GO.ID, Term, sep = "_")) %>% 
  dplyr::select(GO_term, Significant, group) %>% spread(group, Significant) %>%  mutate_if(is.numeric , replace_na, replace = 0) %>% 
  column_to_rownames(., var= "GO_term") 
Hemo_GO_all_dotplot_subset_heatmap_mat <- as.matrix(Hemo_GO_all_dotplot_subset_heatmap)

hemo_labels =c("black","blue","darkgreen","darkred","darkseagreen4", "darkslateblue", "lightpink3","navajowhite2",  "orangered4"  )
# create named vector to hold column names
hemo_column_labels = structure(paste0(hemo_labels), names = paste0(colnames(Hemo_GO_all_dotplot_subset_heatmap )))

# set custom colors
library(circlize)
col_fun = colorRamp2(c(0, 50,100,400), c("white", "red", "orange", "yellow"))

pdf("./FIGURES/hemo_GO_all_heatmap_subset.pdf", width = 12, height = 10)
hemo_GO_all_heatmap_subset <- ComplexHeatmap::Heatmap(Hemo_GO_all_dotplot_subset_heatmap_mat, border = TRUE, 
                                               column_title_side = "bottom", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                               row_title = "GO Term", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                               row_dend_width = unit(2, "cm"),
                                               row_labels = row.names(Hemo_GO_all_dotplot_subset_heatmap_mat),
                                               column_labels = hemo_column_labels[colnames(Hemo_GO_all_dotplot_subset_heatmap)],
                                               # apply split by k-meams clustering to highlight groups
                                               row_km = 3, column_km = 4, row_names_gp = gpar(fontsize = 8),
                                               column_names_gp = gpar(fontsize = 10),
                                               heatmap_legend_param = list(title = "Number Significant"),
                                               col=col_fun,
                                               rect_gp = gpar(col = "black", lwd = 1))
ComplexHeatmap::draw(hemo_GO_all_heatmap_subset, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 100), "mm")) #bottom, left, top, right paddings
dev.off()

# plot as dotplot
Hemo_GO_all_dotplot_subset_plot <- ggplot(Hemo_GO_all_dotplot_subset, aes(x = group, y = Term )) +
  geom_point(aes(size = Significant, color = as.numeric(topgoFisher))) + 
  scale_size_continuous(range = c(4,10)) +
  scale_color_viridis(option = "viridis", name = "p-value", direction = -1) + 
  facet_grid(.~type, scales = "free_x",space = "free_x") + 
  theme_minimal() +
  labs(x = "Module Name", y = "GO Term", title = "GO Enrichment of Significant Modules") + 
  theme(panel.border = element_rect(color = "black", fill = "NA"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text.x = element_text(size = 16, face = "bold"),
        title = element_text(size = 16))

ggsave(Hemo_GO_all_dotplot_subset_plot, device = "tiff", path = "./FIGURES/", 
       filename = "Hemo_GO_all_dotplot_subset_plot_BP.tiff", width = 20, height = 12, limitsize = FALSE)


### Dotplot of just the navajowhite2 module BP and MF terms

hemo_MEnavajowhite2_GOdata_BP_Res$type <- "BP"
hemo_MEnavajowhite2_GOdata_Res$type <- "MF"

hemo_MEnavajowhite2_GOdata_combined <- rbind(hemo_MEnavajowhite2_GOdata_BP_Res,hemo_MEnavajowhite2_GOdata_Res) %>% filter(topgoFisher <=0.05)

hemo_MEnavajowhite2_GOdata_combined_dotplot <- 
  ggplot(hemo_MEnavajowhite2_GOdata_combined, aes(x = group, y = Term )) +
  geom_point(aes(size = Significant, color = as.numeric(topgoFisher))) + 
  scale_size_continuous(range = c(4,10)) +
  scale_color_viridis(option = "viridis", name = "p-value", direction = -1) + 
  facet_grid(type~., scales = "free", space="free") + 
  theme_minimal() +
  labs(x = NULL, y = "GO Term", title = "GO Enrichment") + 
  theme(panel.border = element_rect(color = "black", fill = "NA"),
        axis.text.x = ggtext::element_markdown(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text.y = element_text(size = 16, face = "bold"),
        title = element_text(size = 12))
ggsave(hemo_MEnavajowhite2_GOdata_combined_dotplot, filename = "hemo_MEnavajowhite2_GOdata_combined_dotplot.tiff", 
       path = "./FIGURES/", device = "tiff", width = 6, height = 10)

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
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4 %>% distinct(transcript_id, GS) %>% dplyr::count() # 38
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3  %>% distinct(transcript_id, GS) %>% dplyr::count() # 21
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3  %>% distinct(transcript_id, GS) %>% dplyr::count() # 27
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2  %>% distinct(transcript_id, GS) %>% dplyr::count() # 46

### ASSESS HUB GENE OVERLAPS WITH DEG RESULTS ####
# assessing the modules I've highlighted as being most significant, GDC_lightblue4, ZVADnavajowhite2, ZVADpink3
View(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4)
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot
View(perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot)

# How many overlaps between module genes and interesting modules 
#GDC
perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot[perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot$transcript_id %in% unique(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4$transcript_id),]
  # 0
perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot[perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot$transcript_id %in% unique(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2$transcript_id),]
  # 1, hypothetical protein-XM_002767251.1 TFIIS N-terminal domain profile.
perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot[perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot$transcript_id %in% unique(FilterGenes_Pmar_comb_Interpro_GDC_steelblue$transcript_id),]
  #2, methylase hypothetical protein-XM_002783013.1

# Isolate these lines in the hub genes list
FilterGenes_Pmar_comb_Interpro_GDC_steelblue_DEG <- FilterGenes_Pmar_comb_Interpro_GDC_steelblue[FilterGenes_Pmar_comb_Interpro_GDC_steelblue$transcript_id %in% perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot$transcript_id,]

# ZVAD
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot[perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot$transcript_id %in% unique(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3$transcript_id),]
# 0
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot[perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot$transcript_id %in% unique(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3$transcript_id),]
# 0
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot[perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot$transcript_id %in% unique(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2$transcript_id),]
# 0

### ASSESS FULL MODULE OVERLAPS WITH DEG RESULTS ####

# In case I'm missing something!
# How many overlaps between module genes and interesting modules 
#GDC
GDC_lightblue4_all_genes_DEG_overlap <- perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot[perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot$transcript_id %in% as.character(GDC_lightblue4_all_genes),]
#0
GDC_ZVAD_blue4_all_genes_overlap <- perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot[perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot$transcript_id %in% as.character(GDC_ZVAD_blue4_all_genes),]
#transcript_id  baseMean log2FoldChange     lfcSE       pvalue         padj                                       product         condition    log10
#8  XM_002767241.1  36.74167       3.690627 1.3135200 1.488092e-04 4.445373e-02           hypothetical protein-XM_002767241.1 P_mar_GDC_vs_Pmar 1.352092
#22 XM_002776403.1 114.53786       1.706073 0.3736873 1.971983e-07 2.421815e-04 conserved hypothetical protein-XM_002776403.1 P_mar_GDC_vs_Pmar 3.615859
#37 XM_002787026.1  88.88279       2.701711 0.4627373 3.050280e-10 6.742950e-07               Uridine phosphorylase, putative P_mar_GDC_vs_Pmar 6.171150

GDC_ZVAD_blue4_all_genes_annot_Interpro <- GDC_ZVAD_blue4_all_genes_annot[GDC_ZVAD_blue4_all_genes_annot$transcript_id %in% GDC_ZVAD_blue4_all_genes_overlap$transcript_id,]

# ZVAD
perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot[perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot$transcript_id %in% as.character(ZVAD_lightpink3_all_genes),]
#transcript_id baseMean log2FoldChange   lfcSE       pvalue       padj                                                        product          condition
#7 XM_002771133.1 31.61013       3.300932 1.01682 3.364091e-05 0.01583559 transmembrane BAX inhibitor motif-containing protein, putative P_mar_ZVAD_vs_Pmar

perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot[perk_dds_deseq_res_Pmar_ZVAD_LFC_sig_annot$transcript_id %in% as.character(ZVAD_navajowhite2_all_genes),]

### Pmar GO enrichment for interesting module intramodular hub genes ####
Perk_geneNames <- names(Perk_GO_terms_found_geneID2GO_mapping)
# topGO tutorial: https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
head(Perk_geneNames)

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

FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_GDC_steelblue_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_GDC_darkorange2_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)

FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)

FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Res <- GenTable(FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata, topgoFisher = FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 25)

#### Pmar GO enrichment for interesting module all transcripts  ####
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
ZVAD_navajowhite2_all_genes <- as.factor(ZVAD_navajowhite2_all_genes[!is.na(ZVAD_navajowhite2_all_genes)])

# adding in blue4 as a module potentially response to general P. marinus treatment
GDC_ZVAD_blue4_all_genes <- names(perk_dds_rlog_matrix)[perk_full_moduleColors =="blue4"]  
GDC_ZVAD_blue4_all_genes <- as.factor (GDC_ZVAD_blue4_all_genes[!is.na(GDC_ZVAD_blue4_all_genes)])

GDC_steelblue_all_genes_factor <- factor(as.integer(Perk_geneNames %in% GDC_steelblue_all_genes))
GDC_darkorange2_all_genes_factor <- factor(as.integer(Perk_geneNames %in% GDC_darkorange2_all_genes))
GDC_lightblue4_all_genes_factor <- factor(as.integer(Perk_geneNames %in% GDC_lightblue4_all_genes))
ZVAD_lightpink3_all_genes_factor <- factor(as.integer(Perk_geneNames %in% ZVAD_lightpink3_all_genes))
ZVAD_pink3_all_genes_factor <- factor(as.integer(Perk_geneNames %in% ZVAD_pink3_all_genes))
ZVAD_navajowhite2_all_genes_factor <- factor(as.integer(Perk_geneNames %in% ZVAD_navajowhite2_all_genes))

# adding blue4
GDC_ZVAD_blue4_all_genes_factor <- factor(as.integer(Perk_geneNames %in% GDC_ZVAD_blue4_all_genes))

names(GDC_steelblue_all_genes_factor) <- Perk_geneNames
names(GDC_darkorange2_all_genes_factor) <- Perk_geneNames
names(GDC_lightblue4_all_genes_factor) <- Perk_geneNames
names(ZVAD_lightpink3_all_genes_factor) <- Perk_geneNames
names(ZVAD_pink3_all_genes_factor) <- Perk_geneNames
names(ZVAD_navajowhite2_all_genes_factor) <- Perk_geneNames

names(GDC_ZVAD_blue4_all_genes_factor) <- Perk_geneNames

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

GDC_ZVAD_blue4_GOdata <- new("topGOdata", description = "blue4 Gene Enrichment", 
                                # I want to test MF
                                ontology = "MF",
                                # define here the genes of interest
                                allGenes = GDC_ZVAD_blue4_all_genes_factor,
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

GDC_ZVAD_blue4_GOdata_Fisher_Weight <-  runTest(GDC_ZVAD_blue4_GOdata, algorithm = "weight01", statistic = "fisher")

## Analyze enrichment test results 
# see how many results we get where weight01 gives a P-value <= 0.05
GDC_steelblue_GOdata_summary <- summary(attributes(GDC_steelblue_GOdata_Fisher_Weight)$score <= 0.05) # 9 sig
GDC_darkorange2_GOdata_summary <- summary(attributes(GDC_darkorange2_GOdata_Fisher_Weight)$score <= 0.05) # 8 sig
GDC_lightblue4_GOdata_summary <- summary(attributes(GDC_lightblue4_GOdata_Fisher_Weight)$score <= 0.05) # 4 sig
ZVAD_lightpink3_GOdata_summary <- summary(attributes(ZVAD_lightpink3_GOdata_Fisher_Weight)$score <= 0.05) # 9 sig
ZVAD_pink3_GOdata_summary <- summary(attributes(ZVAD_pink3_GOdata_Fisher_Weight)$score <= 0.05) # 3 sig
ZVAD_navajowhite2_GOdata_summary <- summary(attributes(ZVAD_navajowhite2_GOdata_Fisher_Weight)$score <= 0.05) # 11

GDC_ZVAD_blue4_GOdata_summary <- summary(attributes(GDC_ZVAD_blue4_GOdata_Fisher_Weight)$score <= 0.05) # 2

#print out the top results, though only GDC_lightblue4 is sig
GDC_steelblue_GOdata_Res <- GenTable(GDC_steelblue_GOdata, topgoFisher = GDC_steelblue_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
GDC_darkorange2_GOdata_Res <- GenTable(GDC_darkorange2_GOdata, topgoFisher = GDC_darkorange2_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
GDC_lightblue4_GOdata_Res <- GenTable(GDC_lightblue4_GOdata, topgoFisher = GDC_lightblue4_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)

ZVAD_lightpink3_GOdata_Res <- GenTable(ZVAD_lightpink3_GOdata, topgoFisher = ZVAD_lightpink3_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
ZVAD_pink3_GOdata_Res <- GenTable(ZVAD_pink3_GOdata, topgoFisher = ZVAD_pink3_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
ZVAD_navajowhite2_GOdata_Res <- GenTable(ZVAD_navajowhite2_GOdata, topgoFisher = ZVAD_navajowhite2_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)

GDC_ZVAD_blue4_GOdata_Res <- GenTable(GDC_ZVAD_blue4_GOdata, topgoFisher = GDC_ZVAD_blue4_GOdata_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)

## Repeat for BP
GDC_steelblue_GOdata_BP <- new("topGOdata", description = "GDC_steelblue Gene Enrichment", 
                            # I want to test BP
                            ontology = "BP",
                            # define here the genes of interest
                            allGenes = GDC_steelblue_all_genes_factor,
                            nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

GDC_darkorange2_GOdata_BP <- new("topGOdata", description = "GDC_darkorange2 Gene Enrichment", 
                              # I want to test BP
                              ontology = "BP",
                              # define here the genes of interest
                              allGenes = GDC_darkorange2_all_genes_factor,
                              nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

GDC_lightblue4_GOdata_BP <- new("topGOdata", description = "GDC_lightblue4 Gene Enrichment", 
                             # I want to test BP
                             ontology = "BP",
                             # define here the genes of interest
                             allGenes = GDC_lightblue4_all_genes_factor,
                             nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

ZVAD_lightpink3_GOdata_BP <- new("topGOdata", description = "ZVAD_lightpink3 Gene Enrichment", 
                              # I want to test BP
                              ontology = "BP",
                              # define here the genes of interest
                              allGenes = ZVAD_lightpink3_all_genes_factor,
                              nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

ZVAD_pink3_GOdata_BP <- new("topGOdata", description = "ZVAD_pink3 Gene Enrichment", 
                         # I want to test BP
                         ontology = "BP",
                         # define here the genes of interest
                         allGenes = ZVAD_pink3_all_genes_factor,
                         nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

ZVAD_navajowhite2_GOdata_BP <- new("topGOdata", description = "ZVAD_navajowhite2 Gene Enrichment", 
                                # I want to test BP
                                ontology = "BP",
                                # define here the genes of interest
                                allGenes = ZVAD_navajowhite2_all_genes_factor,
                                nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)

GDC_ZVAD_blue4_GOdata_BP <- new("topGOdata", description = "blue4 Gene Enrichment", 
                             # I want to test MF
                             ontology = "BP",
                             # define here the genes of interest
                             allGenes = GDC_ZVAD_blue4_all_genes_factor,
                             nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = Perk_GO_terms_found_geneID2GO_mapping)


#nodeSize=used to prune the GO hierarchy from the terms which have less than 1 annotated genes
#annFUN.gene2GO = this function is used when the annotations are provided as a gene-to-GOs mapping.

### Perform Encrichment tests 
GDC_steelblue_GOdata_BP_Fisher_Weight <- runTest(GDC_steelblue_GOdata_BP, algorithm = "weight01", statistic = "fisher")
GDC_darkorange2_GOdata_BP_Fisher_Weight <- runTest(GDC_darkorange2_GOdata_BP, algorithm = "weight01", statistic = "fisher")
GDC_lightblue4_GOdata_BP_Fisher_Weight <- runTest(GDC_lightblue4_GOdata_BP, algorithm = "weight01", statistic = "fisher")
ZVAD_lightpink3_GOdata_BP_Fisher_Weight <- runTest(ZVAD_lightpink3_GOdata_BP, algorithm = "weight01", statistic = "fisher")
ZVAD_pink3_GOdata_BP_Fisher_Weight <- runTest(ZVAD_pink3_GOdata_BP, algorithm = "weight01", statistic = "fisher")
ZVAD_navajowhite2_GOdata_BP_Fisher_Weight <- runTest(ZVAD_navajowhite2_GOdata_BP, algorithm = "weight01", statistic = "fisher")

GDC_ZVAD_blue4_GOdata_BP_Fisher_Weight <- runTest(GDC_ZVAD_blue4_GOdata_BP, algorithm = "weight01", statistic = "fisher")

## Analyze enrichment test results 
# see how many results we get where weight01 gives a P-value <= 0.05
GDC_steelblue_GOdata_BP_summary <- summary(attributes(GDC_steelblue_GOdata_BP_Fisher_Weight)$score <= 0.05) 
GDC_darkorange2_GOdata_BP_summary <- summary(attributes(GDC_darkorange2_GOdata_BP_Fisher_Weight)$score <= 0.05)
GDC_lightblue4_GOdata_BP_summary <- summary(attributes(GDC_lightblue4_GOdata_BP_Fisher_Weight)$score <= 0.05) 
ZVAD_lightpink3_GOdata_BP_summary <- summary(attributes(ZVAD_lightpink3_GOdata_BP_Fisher_Weight)$score <= 0.05)
ZVAD_pink3_GOdata_BP_summary <- summary(attributes(ZVAD_pink3_GOdata_BP_Fisher_Weight)$score <= 0.05) 
ZVAD_navajowhite2_GOdata_BP_summary <- summary(attributes(ZVAD_navajowhite2_GOdata_BP_Fisher_Weight)$score <= 0.05) 

GDC_ZVAD_blue4_GOdata_BP_summary <- summary(attributes(GDC_ZVAD_blue4_GOdata_BP_Fisher_Weight)$score <= 0.05) # 2

#print out the top results, though only GDC_lightblue4 is sig
GDC_steelblue_GOdata_BP_Res <- GenTable(GDC_steelblue_GOdata_BP, topgoFisher = GDC_steelblue_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
GDC_darkorange2_GOdata_BP_Res <- GenTable(GDC_darkorange2_GOdata_BP, topgoFisher = GDC_darkorange2_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
GDC_lightblue4_GOdata_BP_Res <- GenTable(GDC_lightblue4_GOdata_BP, topgoFisher = GDC_lightblue4_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)

ZVAD_lightpink3_GOdata_BP_Res <- GenTable(ZVAD_lightpink3_GOdata_BP, topgoFisher = ZVAD_lightpink3_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
ZVAD_pink3_GOdata_BP_Res <- GenTable(ZVAD_pink3_GOdata_BP, topgoFisher = ZVAD_pink3_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)
ZVAD_navajowhite2_GOdata_BP_Res <- GenTable(ZVAD_navajowhite2_GOdata_BP, topgoFisher = ZVAD_navajowhite2_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)

GDC_ZVAD_blue4_GOdata_BP_Res <- GenTable(GDC_ZVAD_blue4_GOdata_BP, topgoFisher = GDC_ZVAD_blue4_GOdata_BP_Fisher_Weight, orderBy = "topgoFisher", topNodes = 20)

#### Pmar top interesting modules GO visualization ####

### Dotplot of significantly enriched GO terms from hub modules plotting just the MF
  # GDC_lightblue4, ZVAD navajowhite2 and ZVAD pink 3

FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_Res$group <-"pink3"
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Res$group <-"navajowhite2"
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Res$group <-"lightblue4"
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_Res$treat <-"ZVAD-fmk"
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Res$treat <-"ZVAD-fmk"
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Res$treat <-"GDC-0152"

Pmar_GO_hub_dotplot <- rbind(FilterGenes_Pmar_comb_Interpro_ZVAD_pink3_GOdata_Res,
                             FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Res,
                             FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Res) %>% filter(topgoFisher <=0.05)

# edit GO terms where the rest of the term name is cutoff 
Pmar_GO_hub_dotplot <- Pmar_GO_hub_dotplot %>% mutate(Term = case_when(
  GO.ID == "GO:0016646" ~ "oxidoreductase activity, acting on the CH-NH group of donors, NAD or NADP as acceptor",
  GO.ID == "GO:0046933" ~ "proton-transporting ATP synthase activity, rotational mechanism",
  GO.ID == "GO:0016684" ~ "oxidoreductase activity, acting on peroxide as acceptor",
  GO.ID == "GO:0009678" ~ "pyrophosphate hydrolysis-driven proton transmembrane transporter activity",
  GO.ID == "GO:0004198" ~ "calcium-dependent cysteine-type endopeptidase activity",
  TRUE ~ Term
))

Pmar_GO_hub_dotplot_plot <- ggplot(Pmar_GO_hub_dotplot, aes(x = group, y = Term )) +
  geom_point(aes(size = Significant, color = as.numeric(topgoFisher))) + 
  scale_size_continuous(range = c(4,10)) +
  scale_color_viridis(option = "viridis", name = "p-value", direction = -1) + 
  facet_grid(.~treat, scales = "free") + 
  theme_minimal() +
  labs(x = "Module Name", y = "GO Term", title = "GO Enrichment of Intramodular Hub Genes") + 
  theme(panel.border = element_rect(color = "black", fill = "NA"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text.x = element_text(size = 16, face = "bold"),
        title = element_text(size = 16))

ggsave(Pmar_GO_hub_dotplot_plot, device = "tiff", path = "./FIGURES/", 
       filename = "Pmar_GO_hub_dotplot_plot.tiff", width = 15, height = 10)

### View interesting genes associated with particular terms

# GDC_lightblue4
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Res_sig_terms <- FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Res %>% 
  filter(topgoFisher <=0.05) %>% dplyr::select(GO.ID, Term) %>% dplyr::rename(Ontology_term = GO.ID)

FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO_enriched <- left_join(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Res_sig_terms, 
                                                                       FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GO) %>% left_join(., dplyr::select(FilterGenes_Pmar_comb_Interpro_GDC_lightblue4, -Ontology_term))
# the terms in this list are not the ones that I found interesting! 
# I think this has to do with bias in GO term identification?


### GO figure from enrichment analysis of all module genes
GDC_steelblue_GOdata_Res$group <- "steelblue"
GDC_darkorange2_GOdata_Res$group <- "darkorange2"
GDC_lightblue4_GOdata_Res$group <- "lightblue4"
ZVAD_lightpink3_GOdata_Res$group <- "lightpink3"
ZVAD_pink3_GOdata_Res$group <-"pink3"
ZVAD_navajowhite2_GOdata_Res $group <-"navajowhite2"

GDC_steelblue_GOdata_Res$treat <- "GDC-0152"
GDC_darkorange2_GOdata_Res$treat <- "GDC-0152"
GDC_lightblue4_GOdata_Res$treat <- "GDC-0152"
ZVAD_lightpink3_GOdata_Res$treat <- "ZVAD-fmk"
ZVAD_pink3_GOdata_Res$treat <-"ZVAD-fmk"
ZVAD_navajowhite2_GOdata_Res $treat <-"ZVAD-fmk"

Pmar_GO_all_dotplot <- rbind(GDC_steelblue_GOdata_Res,
                             GDC_darkorange2_GOdata_Res,
                             GDC_lightblue4_GOdata_Res,
                             ZVAD_lightpink3_GOdata_Res,
                             ZVAD_pink3_GOdata_Res,
                             ZVAD_navajowhite2_GOdata_Res) %>% filter(topgoFisher <=0.05)

# edit GO terms where the rest of the term name is cutoff 

Pmar_GO_all_dotplot_plot <- ggplot(Pmar_GO_all_dotplot, aes(x = group, y = Term )) +
  geom_point(aes(size = Significant, color = as.numeric(topgoFisher))) + 
  scale_size_continuous(range = c(4,10)) +
  scale_color_viridis(option = "viridis", name = "p-value", direction = -1) + 
  facet_grid(.~treat, scales = "free") + 
  theme_minimal() +
  labs(x = "Module Name", y = "GO Term", title = "GO Enrichment of Module Genes") + 
  theme(panel.border = element_rect(color = "black", fill = "NA"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text.x = element_text(size = 16, face = "bold"),
        title = element_text(size = 16))

ggsave(Pmar_GO_all_dotplot_plot, device = "tiff", path = "./FIGURES/", 
       filename = "Pmar_GO_all_dotplot_plot.tiff", width = 15, height = 10)

#### Assess Pmar SOD ####

## Lau et al., 2018 SOD primers are best hit for 4 iron-dependent superoxide dismutase proteins
SOD_list <- c("XM_002768746.1", "XM_002768745.1", "XM_002765900.1", "XM_002765899.1")
# search in the hub gene lists

# these SOD enzymes are not hub genes 
FilterGenes_Pmar_comb_Interpro_GDC_steelblue %>% filter(transcript_id %in% SOD_list) # 0
FilterGenes_Pmar_comb_Interpro_GDC_darkorange2 %>% filter(transcript_id %in% SOD_list) # 0
FilterGenes_Pmar_comb_Interpro_GDC_lightblue4 %>% filter(transcript_id %in% SOD_list) # 0
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3 %>% filter(transcript_id %in% SOD_list) # 0
FilterGenes_Pmar_comb_Interpro_ZVAD_pink3 %>% filter(transcript_id %in% SOD_list) # 0
FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2 %>% filter(transcript_id %in% SOD_list) # 0
 
# Check in all genes in each module
GDC_steelblue_all_genes[GDC_steelblue_all_genes %in% SOD_list]
GDC_darkorange2_all_genes[GDC_darkorange2_all_genes %in% SOD_list]
GDC_lightblue4_all_genes[GDC_lightblue4_all_genes %in% SOD_list]
ZVAD_lightpink3_all_genes[ZVAD_lightpink3_all_genes %in% SOD_list]
ZVAD_pink3_all_genes[ZVAD_pink3_all_genes %in% SOD_list]
ZVAD_navajowhite2_all_genes[ZVAD_navajowhite2_all_genes %in% SOD_list]

# no SOD in any of these modules

#### QUANTIFY MODULE ASSOCIATIONS WITH CHALLENGE AND APOPTOSIS PHENOTYPE ####

## Add phenotype data for the hemocytes that have engulfed hemocytes and to the binarized data
PCA_pheno_2020_all_samplename <-PCA_pheno_2020_all %>% mutate(Sample_Name = paste(ID, Treat, sep = "_")) 
PCA_pheno_2020_all_samplename_APOP_hemo_perk <- na.omit(PCA_pheno_2020_all_samplename[,c("Sample_Name","Percent_of_this_plot_APOP_hemo_perk")])
# Format by arcsine transform
PCA_pheno_2020_all_samplename_APOP_hemo_perk$Percent_of_this_plot_APOP_hemo_perk_arcsin <- transf.arcsin(PCA_pheno_2020_all_samplename_APOP_hemo_perk$Percent_of_this_plot_APOP_hemo_perk*0.01)

# Join apoptosis phenotype 
hemo_coldata_collapse_binarize_apop_perk <- as.data.frame(hemo_coldata_collapse_binarize) %>% rownames_to_column(., var = "Sample_Name") %>% 
  left_join(.,PCA_pheno_2020_all_samplename_APOP_hemo_perk) %>% column_to_rownames(., var = "Sample_Name")

# tutorial for this section: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf

# Recalculate Module trait correlations and Pvalue using this new phenotype dataframe
hemo_full_apop_moduleTraitCor = cor(hemo_full_MEs, hemo_coldata_collapse_binarize_apop_perk, use = "p");
hemo_full_apop_moduleTraitPvalue = corPvalueStudent(hemo_full_apop_moduleTraitCor, hemo_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
hemo_full_apop_textMatrix = paste(signif(hemo_full_apop_moduleTraitCor, 2), "\n(",
                             signif(hemo_full_apop_moduleTraitPvalue, 1), ")", sep = "");
dim(hemo_full_apop_textMatrix) = dim(hemo_full_apop_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = hemo_full_apop_moduleTraitCor,
               xLabels = names(hemo_coldata_collapse_binarize_apop_perk),
               yLabels = names(hemo_full_MEs),
               ySymbols = names(hemo_full_MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = hemo_full_apop_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with apoptosis (high correlation and low P value)?
hemo_full_apop_moduleTraitCor_df <- as.data.frame(hemo_full_apop_moduleTraitCor) %>% dplyr::select(Percent_of_this_plot_APOP_hemo_perk,Percent_of_this_plot_APOP_hemo_perk_arcsin)
colnames(hemo_full_apop_moduleTraitCor_df) <- c("Percent_of_this_plot_APOP_hemo_perk.moduleTraitCor", "Percent_of_this_plot_APOP_hemo_perk_arcsin.moduleTraitCor")
hemo_full_apop_moduleTraitPvalue_df <- as.data.frame(hemo_full_apop_moduleTraitPvalue) %>% dplyr::select(Percent_of_this_plot_APOP_hemo_perk,Percent_of_this_plot_APOP_hemo_perk_arcsin)
colnames(hemo_full_apop_moduleTraitPvalue_df) <- c("Percent_of_this_plot_APOP_hemo_perk.moduleTraitPvalue", "Percent_of_this_plot_APOP_hemo_perk_arcsin.moduleTraitPvalue")

hemo_full_apop_moduleTraitCor_Pval_df <- cbind(hemo_full_apop_moduleTraitCor_df, hemo_full_apop_moduleTraitPvalue_df) 

hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk <- hemo_full_apop_moduleTraitCor_Pval_df %>% dplyr::select(contains("Percent_of_this_plot_APOP_hemo_perk."))
hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_arcsin <- hemo_full_apop_moduleTraitCor_Pval_df %>% dplyr::select(contains("APOP_hemo_perk_arcsin."))

# Significantly correlated modules
hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig <- hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk %>% filter(Percent_of_this_plot_APOP_hemo_perk.moduleTraitPvalue <= 0.05)  %>% rownames_to_column(., "mod_names") 
nrow(hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig) #14 significant modules

hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_arcsin_sig <- hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_arcsin %>% filter(Percent_of_this_plot_APOP_hemo_perk_arcsin.moduleTraitPvalue<= 0.05)  %>% rownames_to_column(., "mod_names") 
nrow(hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_arcsin_sig) #13

# compare modules
hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_compare <-  
  full_join(hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig,hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_arcsin_sig) %>% dplyr::select(!contains("Pvalue"))

# find those shared between all 
hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_compare_shared <- drop_na(hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_compare)
# 12 shared between both...the arcine transform does not make a difference! Correlation values are very similar for both, in general slightly lower correlation for arcsine 

## Plot only those modules with significant apoptosis members as before
## Heatmap of only the significantly correlated modules that have apoptosis hub genes
# Graph and color code each the strength of association (correlation) of module eigengenes and trait

# subset hemo_full_moduleTraitCor, hemo_full_moduleTraitPvalue, hemo_full_MEs for only those modules significant in either challenge
# also remove for plotting those modules that are only correlated with ZVAD
FilterGenes_comb_apop_count_mod_names_noZVAD <- FilterGenes_comb_apop_count_mod_names[!FilterGenes_comb_apop_count_mod_names %in% "MEdarkslateblue"]
hemo_full_moduleTraitCor_sig_apop_pheno <- hemo_full_apop_moduleTraitCor[rownames(hemo_full_apop_moduleTraitCor) %in% FilterGenes_comb_apop_count_mod_names_noZVAD,]
hemo_full_moduleTraitCor_sig_apop_pheno <- hemo_full_moduleTraitCor_sig_apop_pheno[,c(1,2,8)] # keep only control, GDC, and apop_arcsin
hemo_full_moduleTraitPvalue_sig_apop_pheno <- hemo_full_apop_moduleTraitPvalue[rownames(hemo_full_apop_moduleTraitPvalue) %in% FilterGenes_comb_apop_count_mod_names_noZVAD,]
hemo_full_moduleTraitPvalue_sig_apop_pheno <- hemo_full_moduleTraitPvalue_sig_apop_pheno[,c(1,2,8)]
hemo_full_MEs_sig_apop_pheno <- hemo_full_MEs[,colnames(hemo_full_MEs) %in% FilterGenes_comb_apop_count_mod_names_noZVAD]
hemo_coldata_collapse_binarize_sig_apop_pheno <- hemo_coldata_collapse_binarize_apop_perk[,c(1,2,8)]

# Will display correlations and their p-values
hemo_full_textMatrix_sig_apop_pheno = paste(signif(hemo_full_moduleTraitCor_sig_apop_pheno, 2), "\n(",
                                      signif(hemo_full_moduleTraitPvalue_sig_apop_pheno, 1), ")", sep = "");
dim(hemo_full_textMatrix_sig_apop_pheno) = dim(hemo_full_moduleTraitCor_sig_apop_pheno)

# make plot for use in multipanel figure
sizeGrWindow(10,10)
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
pdf("./FIGURES/hemo_full_moduleTraitCor_sig_apop_pheno_heatmap.pdf", width = 2.5, height = 4)
labeledHeatmap(Matrix = hemo_full_moduleTraitCor_sig_apop_pheno,
               xLabels = c("Control vs. P. marinus", "Control vs P. marinus\nand GDC-0152 ", "Arsine Transformed\nApoptosis Percentage"),
               yLabels = names(hemo_full_MEs_sig_apop_pheno),
               ySymbols = names(hemo_full_MEs_sig_apop_pheno),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = hemo_full_textMatrix_sig_apop_pheno,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.lab = 0.3,
               #zlim = c(-1,1), 
               yColorWidth = 0.2)
dev.off()
dev.off()

### Repeat for P_marinus counts 

# format phenotype data 
perk_coldata_collapse_binarize
# Join apoptosis phenotype 
perk_coldata_collapse_binarize_apop_perk <- as.data.frame(perk_coldata_collapse_binarize) %>% rownames_to_column(., var = "Sample_Name") %>% 
  left_join(.,PCA_pheno_2020_all_samplename_APOP_hemo_perk) %>% column_to_rownames(., var = "Sample_Name")

# recalculate the correlation with traits
perk_full_apop_moduleTraitCor = cor(perk_full_MEs, perk_coldata_collapse_binarize_apop_perk, use = "p");
perk_full_apop_moduleTraitPvalue = corPvalueStudent(perk_full_apop_moduleTraitCor, perk_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
perk_full_apop_textMatrix = paste(signif(perk_full_apop_moduleTraitCor, 2), "\n(",
                             signif(perk_full_apop_moduleTraitPvalue, 1), ")", sep = "");
dim(perk_full_apop_textMatrix) = dim(perk_full_apop_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = perk_full_apop_moduleTraitCor,
               xLabels = names(perk_coldata_collapse_binarize_apop_perk),
               yLabels = names(perk_full_MEs),
               ySymbols = names(perk_full_MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = perk_full_apop_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
perk_full_apop_moduleTraitCor_df <- as.data.frame(perk_full_apop_moduleTraitCor) %>% dplyr::select(Percent_of_this_plot_APOP_hemo_perk,Percent_of_this_plot_APOP_hemo_perk_arcsin)
colnames(perk_full_apop_moduleTraitCor_df) <- c( "Percent_of_this_plot_APOP_hemo_perk.moduleTraitCor","Percent_of_this_plot_APOP_hemo_perk_arcsin.moduleTraitCor")
perk_full_apop_moduleTraitPvalue_df <- as.data.frame(perk_full_apop_moduleTraitPvalue) %>% dplyr::select(Percent_of_this_plot_APOP_hemo_perk,Percent_of_this_plot_APOP_hemo_perk_arcsin)
colnames(perk_full_apop_moduleTraitPvalue_df) <- c( "Percent_of_this_plot_APOP_hemo_perk.moduleTraitPvalue","Percent_of_this_plot_APOP_hemo_perk_arcsin.moduleTraitPvalue")

perk_full_apop_moduleTraitCor_Pval_df <- cbind(perk_full_apop_moduleTraitCor_df, perk_full_apop_moduleTraitPvalue_df) 

perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk <- perk_full_apop_moduleTraitCor_Pval_df %>% dplyr::select(contains("Percent_of_this_plot_APOP_hemo_perk"))
perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_arcsin <- perk_full_apop_moduleTraitCor_Pval_df %>% dplyr::select(contains("Percent_of_this_plot_APOP_hemo_perk_arcsin"))

# Significantly correlated modules

perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig <- perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk %>% filter(Percent_of_this_plot_APOP_hemo_perk.moduleTraitPvalue<= 0.05)  %>% rownames_to_column(., "mod_names") 
nrow(perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig) #15 significant modules

perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_arcsin_sig <- perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_arcsin %>% filter(Percent_of_this_plot_APOP_hemo_perk_arcsin.moduleTraitPvalue<= 0.05)  %>% rownames_to_column(., "mod_names") 
nrow(perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_arcsin_sig) #14 significant modules

# compare modules between these and the other modules
perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_compare <- 
  full_join(perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig,perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_arcsin_sig) %>%
  # also join with original perkinsus dataframe
  full_join(.,perk_full_moduleTraitCor_Pval_df_Pmar_sig_compare) %>%
  dplyr::select(!contains("Pvalue"))

# all mod_names significant in either challenge
perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_compare_mod_names <- perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_compare$mod_names

# find those shared between all 
perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_compare_shared <- drop_na(perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_compare)
# 15 shared, again all are almost exactly the same correlation

## Heatmap of only the significantly correlated modules 
# Graph and color code each the strength of association (correlation) of module eigengenes and trait

# subset perk_full_moduleTraitCor, perk_full_moduleTraitPvalue, perk_full_MEs for only those modules significant in either challenge
perk_full_apop_moduleTraitCor_sig <- perk_full_apop_moduleTraitCor[rownames(perk_full_apop_moduleTraitCor) %in% perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_compare_mod_names,]
perk_full_apop_moduleTraitCor_sig <- perk_full_apop_moduleTraitCor_sig[,c(-3)]
perk_full_apop_moduleTraitPvalue_sig <- perk_full_apop_moduleTraitPvalue[rownames(perk_full_apop_moduleTraitPvalue) %in% perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_compare_mod_names,]
perk_full_apop_moduleTraitPvalue_sig <- perk_full_apop_moduleTraitPvalue_sig[,c(-3)]
perk_full_apop_MEs_sig <- perk_full_MEs[,colnames(perk_full_MEs) %in% perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_compare_mod_names]
perk_coldata_collapse_binarize_apop_perk_sig <- perk_coldata_collapse_binarize_apop_perk[,c(-3)]

# Will display correlations and their p-values
perk_full_apop_textMatrix_sig = paste(signif(perk_full_apop_moduleTraitCor_sig, 2), "\n(",
                                 signif(perk_full_apop_moduleTraitPvalue_sig, 1), ")", sep = "");
dim(perk_full_apop_textMatrix_sig) = dim(perk_full_apop_moduleTraitCor_sig)

# make plot
sizeGrWindow(10,6)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = perk_full_apop_moduleTraitCor_sig,
               xLabels = c("P.mar. and GDC-0152", "P. mar. and ZVAD-fmk", "Apoptosis_Percentage", "Apoptosis Percentage Arcsine"),
               yLabels = names(perk_full_apop_MEs_sig),
               ySymbols = names(perk_full_apop_MEs_sig),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = perk_full_apop_textMatrix_sig,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.7,
               #zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("P. marinus Module-Challenge Relationships"))
# have to save manually...weird!

### Gene relationship to trait and important modules: Gene Significance and Module Membership ####

## calculate gene trait significance for each treatment
# Apoptosis  percentage 
hemo_full_apop = as.data.frame(hemo_coldata_collapse_binarize_apop_perk$Percent_of_this_plot_APOP_hemo_perk);
names(hemo_full_apop) = "Hemo_apop"
hemo_full_geneTraitSignificance_apop = as.data.frame(cor(hemo_dds_rlog_matrix,hemo_full_apop, use = "p"))
hemo_full_GSPvalue_apop = as.data.frame(corPvalueStudent(as.matrix(hemo_full_geneTraitSignificance_apop), hemo_full_nSamples))

names(hemo_full_geneTraitSignificance_apop) = paste("GS.", names(hemo_full_apop), sep="")
names(hemo_full_GSPvalue_apop) = paste("p.GS.", names(hemo_full_apop), sep="")

hemo_full_apop_moduleTraitPvalue_df

# Hemo APOP percentage
hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_list <-  as.character(unlist(hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig$mod_names))
hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_list <- str_remove(hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_list, "ME")

GS_MM_plot_hemo_apop <- function(list) {
  hemo_full_module = list 
  hemo_full_column = match(hemo_full_module, hemo_full_modNames)
  hemo_full_moduleGenes = hemo_full_moduleColors==hemo_full_module
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(hemo_full_geneModuleMembership [hemo_full_moduleGenes, hemo_full_column]),
                     abs(hemo_full_geneTraitSignificance_apop[hemo_full_moduleGenes, 1]),
                     xlab = paste("Module Membership in", hemo_full_module, "module"),
                     ylab = "Gene significance for challenge",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = hemo_full_module)
  quartz.save(paste("./FIGURES/hemo_Pmar_apop_perc",list, sep ="_"), type = "png", device = dev.cur(), dpi = 100)
}
lapply(hemo_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_list,  GS_MM_plot_hemo_apop)

## modules with a good correlation between GS and module membership (looking at those with higher than 0.4 correlation)
# high correlation of GS and MM illustrates that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait
  # >0.6
    # navajowhite2 = 0.63
    # lightpink3 = 0.61
    # tan4 = 0.64
  # 0.4 -0.6
    # royalblue = 0.43
    # palevioletred1 = 0.54
    # darkgreen = 0.44
    # mediumpurple1 = 0.56
    # darkred = 0.53

## Create plot for just navajowhite2 for use in paper
hemo_full_module = "navajowhite2" 
hemo_full_column = match(hemo_full_module, hemo_full_modNames)
hemo_full_moduleGenes = hemo_full_moduleColors==hemo_full_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
pdf("./FIGURES/hemo_navajowhite2_apop_GS_MM_plot.pdf",width=5, height = 5)
verboseScatterplot(abs(hemo_full_geneModuleMembership [hemo_full_moduleGenes, hemo_full_column]),
                   abs(hemo_full_geneTraitSignificance_apop[hemo_full_moduleGenes, 1]),
                   xlab = paste("Module Membership in", hemo_full_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module Membership vs. Gene Significance\n"),
                   cex.main = 1.0, cex.lab = 1.0, cex.axis = 1.0, col = "black")
dev.off()
dev.off()

## calculate gene trait significance for each treatment
# Pmar vs Apop percentage 
perk_full_apop = as.data.frame(perk_coldata_collapse_binarize_apop_perk$Percent_of_this_plot_APOP_hemo_perk);
names(perk_full_apop) = "Pmar_apop"
perk_full_geneTraitSignificance_apop = as.data.frame(cor(perk_dds_rlog_matrix,perk_full_apop, use = "p"))
perk_full_GSPvalue_apop = as.data.frame(corPvalueStudent(as.matrix(perk_full_geneTraitSignificance_apop), perk_full_nSamples))

names(perk_full_geneTraitSignificance_apop) = paste("GS.", names(perk_full_apop), sep="")
names(perk_full_GSPvalue_apop) = paste("p.GS.", names(perk_full_apop), sep="")

## Parasite intramodular analysis - perform for each treatment 
# P. mar APOP percentage
perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_list <-  as.character(unlist(perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig$mod_names))
perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_list <- str_remove(perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_list, "ME")

GS_MM_plot_perk_apop <- function(list) {
  perk_full_module = list 
  perk_full_column = match(perk_full_module, perk_full_modNames)
  perk_full_moduleGenes = perk_full_moduleColors==perk_full_module
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(perk_full_geneModuleMembership [perk_full_moduleGenes, perk_full_column]),
                     abs(perk_full_geneTraitSignificance_apop[perk_full_moduleGenes, 1]),
                     xlab = paste("Module Membership in", perk_full_module, "module"),
                     ylab = "Gene significance for challenge",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = perk_full_module)
  quartz.save(paste("./FIGURES/perk_Pmar_apop_perc",list, sep ="_"), type = "png", device = dev.cur(), dpi = 100)
}
lapply(perk_full_apop_moduleTraitCor_Pval_df_APOP_hemo_perk_sig_list,  GS_MM_plot_perk_apop)

## modules with a good correlation between GS and module membership (looking at those with higher than 0.4 correlation)
# high correlation of GS and MM illustrates that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait
  # pink = 0.56
  # deeppink = 0.66
  # darkseagreen4 = 0.66
  # pink3 = 0.61
  # lightblue3 = 0.44
  # darkolivegreen = 0.77
  # black = 0.61
  # darkorange = 0.76
  # purple = 0.47
  # lightblue4 = 0.4
  # grey = 0.45

### FIND INTRAMODULAR HUB GENES USING APOPTOSIS PHENOTYPE ####

# Which modules should be looked at?
# For hemocytes: modules that are significant for treatment, have >5 apoptosis transcripts, and have high (>0.4) correlation between GS and MM
# For hemocytes: same as above, except no filtering of modules for apoptosis related transcripts
# filter genes for each module of interest that have correlation between 

# annotate intramodular hub genes in each module
matrix = hemo_dds_rlog_matrix
moduleColors = hemo_full_moduleColors
lookup =   C_vir_rtracklayer_transcripts
datKME =  hemo_datKME

lookup_annot_intramodular <- function(list) {
  module_name = str_remove(list, "MM.")
  mod_list <- names(matrix)[moduleColors == module_name]
  FilterGenes = abs(GS1)> .6 & abs(datKME[list])>.8
  FilterGenes_annot <- data.frame("ID" = dimnames(data.frame(matrix))[[2]][FilterGenes]) %>% left_join(., lookup) %>% mutate(mod_names = list) %>% filter(ID %in% mod_list)
}

## Hemocytes first
# control_Pmar_GDC 

GS1 = hemo_full_geneTraitSignificance_apop
hemo_control_Pmar_GDC_navajowhite2 <- c("MM.navajowhite2")
names(hemo_control_Pmar_GDC_navajowhite2) <- c("MM.navajowhite2")

FilterGenes_control_Pmar_GDC_navajowhite2 <- lapply(hemo_control_Pmar_GDC_navajowhite2, lookup_annot_intramodular)
FilterGenes_control_Pmar_GDC_navajowhite2_df <- do.call(rbind,FilterGenes_control_Pmar_GDC_navajowhite2) %>% mutate(group = "control_Pmar_GDC")

# trait correlation for each module
hemo_full_apop_moduleTraitPvalue_df
# join with gene trait significance 
hemo_full_geneTraitSignificance_apop_join <- hemo_full_geneTraitSignificance_apop %>% rownames_to_column(.,var = "ID")
FilterGenes_control_Pmar_GDC_navajowhite2_apop_sig <- left_join(FilterGenes_control_Pmar_GDC_navajowhite2_df,hemo_full_geneTraitSignificance_apop_join) %>% dplyr::select(ID,product,gene,GS.Hemo_apop)

# find any apoptosis related 

FilterGenes_control_Pmar_GDC_navajowhite2_apop_sig_apop_transcript <- left_join(FilterGenes_control_Pmar_GDC_navajowhite2_apop_sig, C_vir_rtracklayer_apop_product_final[,c("ID","product","Name")])

#### ANNOTATE ALL GENES IN MOST INTERESTING MODULES ####

# Hemocyte list: "MEnavajowhite2", "MEblue", "MEyellow","MEdarkslateblue", "MEorangered4"
hemo_MEnavajowhite2_annot <- as.data.frame(hemo_MEnavajowhite2) %>% dplyr::rename("ID" = "hemo_MEnavajowhite2") %>% left_join(., C_vir_rtracklayer_transcripts)
hemo_MEblue_annot <- as.data.frame(hemo_MEblue) %>% dplyr::rename("ID" = "hemo_MEblue") %>% left_join(., C_vir_rtracklayer_transcripts)
hemo_MEyellow_annot <- as.data.frame(hemo_MEyellow) %>% dplyr::rename("ID" = "hemo_MEyellow") %>% left_join(., C_vir_rtracklayer_transcripts)
hemo_MEdarkslateblue_annot <- as.data.frame(hemo_MEdarkslateblue) %>% dplyr::rename("ID" = "hemo_MEdarkslateblue") %>% left_join(., C_vir_rtracklayer_transcripts)
hemo_MEorangered4_annot <- as.data.frame(hemo_MEorangered4) %>% dplyr::rename("ID" = "hemo_MEorangered4") %>% left_join(., C_vir_rtracklayer_transcripts)

# Annotate apoptosis from all 
hemo_MEnavajowhite2_annot_apop <- as.data.frame(hemo_MEnavajowhite2) %>% dplyr::rename("ID" = "hemo_MEnavajowhite2") %>% left_join(., C_vir_rtracklayer_apop_product_final) %>% filter(!is.na(product))
  # 8 total apop genes including caspase 2 and caspase 8, TNFRSF5, AP-1
hemo_MEblue_annot_apop <- as.data.frame(hemo_MEblue) %>% dplyr::rename("ID" = "hemo_MEblue") %>% left_join(., C_vir_rtracklayer_apop_product_final) %>% filter(!is.na(product))
  # 41 total apoptosis genes including many NF-kB pathway
hemo_MEyellow_annot_apop <- as.data.frame(hemo_MEyellow) %>% dplyr::rename("ID" = "hemo_MEyellow") %>% left_join(., C_vir_rtracklayer_apop_product_final) %>% filter(!is.na(product))
  # 51 total apoptosis - many IAPs, ER stress response, some TLR and TNFR
hemo_MEdarkslateblue_annot_apop <- as.data.frame(hemo_MEdarkslateblue) %>% dplyr::rename("ID" = "hemo_MEdarkslateblue") %>% left_join(., C_vir_rtracklayer_apop_product_final) %>% filter(!is.na(product))
  # 8 total apoptosis, TRAF2, cytochrome c, GADD45a, cathepsin B could suggest
hemo_MEorangered4_annot_apop <- as.data.frame(hemo_MEorangered4) %>% dplyr::rename("ID" = "hemo_MEorangered4") %>% left_join(., C_vir_rtracklayer_apop_product_final) %>% filter(!is.na(product))
  # also only 8 apoptosis transcripts, caspase7, TLR2, TLR4, TLR13, IFI27, calpain-7

## Annotate all perkinsus interesting perkinsus
GDC_lightblue4_all_genes_annot <- as.data.frame(GDC_lightblue4_all_genes) %>% dplyr::rename("transcript_id" = "GDC_lightblue4_all_genes") %>% left_join(., Perk_Interpro_GO_terms_XP)
ZVAD_lightpink3_all_genes_annot <- as.data.frame(ZVAD_lightpink3_all_genes) %>% dplyr::rename("transcript_id" = "ZVAD_lightpink3_all_genes") %>% left_join(., Perk_Interpro_GO_terms_XP)
ZVAD_navajowhite2_all_genes_annot <- as.data.frame(ZVAD_navajowhite2_all_genes) %>% dplyr::rename("transcript_id" = "ZVAD_navajowhite2_all_genes") %>% left_join(., Perk_Interpro_GO_terms_XP)

# add in blue4
GDC_ZVAD_blue4_all_genes_annot <- as.data.frame(GDC_ZVAD_blue4_all_genes) %>% dplyr::rename("transcript_id" = "GDC_ZVAD_blue4_all_genes") %>% left_join(., Perk_Interpro_GO_terms_XP)


#### EXPORT COMPILED DATA TO SPREADSHEETS ####

# GOAL:
# Export two data frames with all treatments, each for the hemocyte experiment and for the perkinsus experiment
# one table will have all the significant module hub genes, the GS and MM, and whether it is an apoptosis hub gene
# one table will have all the significantly enriched GO terms from the full module

## Hemocyte experiment
# df for all hub genes
FilterGenes_comb
# df for all apoptosis hub genes
FilterGenes_comb_apop

# add apoptosis label for the apoptosis DEGs
FilterGenes_comb_apop <- FilterGenes_comb_apop %>% mutate(apop = "apoptosis")

# join DFs
hemo_mod_list <- c("MEnavajowhite2", "MEblue", "MEyellow","MEdarkslateblue", "MEorangered4")
FilterGenes_comb_apop_labeled <- left_join(FilterGenes_comb, FilterGenes_comb_apop[,c("apop","transcript_id")]) %>% 
  dplyr::select(ID,gene,product,transcript_id,mod_names,group,moduleTraitCor,moduleTraitPvalue,GS, apop) %>%
    # filter out only the module of interest
    filter(mod_names %in% hemo_mod_list) %>% 
  distinct(ID, mod_names, group, .keep_all = TRUE)

# Compile hemocyte experiment GO enrichment for all important modules (run with all genes, not just the hub genes), for both BP and MF
# using only my pruned list of most important modules: MEnavajowhite2, MEblue, MEyellow,MEdarkslateblue, MEorangered4

Hemo_GO_export_subset <- rbind(
  hemo_MEnavajowhite2_GOdata_Res,
  hemo_MEblue_GOdata_Res,
  hemo_MEyellow_GOdata_Res,
  hemo_MEdarkslateblue_GOdata_Res,
  hemo_MEorangered4_GOdata_Res) %>% filter(topgoFisher >= 0.05) %>% mutate(GO_level = "MF") %>% dplyr::select(-type)

hemo_MEnavajowhite2_GOdata_BP_Res$group <- "navajowhite2"
hemo_MEblue_GOdata_BP_Res$group <- "blue"
hemo_MEyellow_GOdata_BP_Res$group <- "yellow"
hemo_MEdarkslateblue_GOdata_BP_Res$group <- "darkslateblue"
hemo_MEorangered4_GOdata_BP_Res$group <- "orangered4"

Hemo_GO_BP_export_subset <- rbind(
  hemo_MEnavajowhite2_GOdata_BP_Res,
  hemo_MEblue_GOdata_BP_Res,
  hemo_MEyellow_GOdata_BP_Res,
  hemo_MEdarkslateblue_GOdata_BP_Res,
  hemo_MEorangered4_GOdata_BP_Res) %>% filter(topgoFisher >= 0.05) %>% mutate(GO_level = "BP")

# combine GO data
Hemo_GO_export_subset_all <- rbind(Hemo_GO_export_subset, Hemo_GO_BP_export_subset)

## Perkinsus experiment
# all hub genes without the interproscan
FilterGenes_Pmar_comb
# all hub genes with interproscan
FilterGenes_Pmar_comb_Interpro_slim <- FilterGenes_Pmar_comb_Interpro %>% filter(signature_desc !="consensus disorder prediction") %>% filter(Dbxref != "character(0)") %>% 
  filter(!is.na(Name)) %>% filter(signature_desc != "character(0)") %>%
  distinct(transcript_id, Dbxref, .keep_all = TRUE)

# Export GO enrichment data - FOR all genes not just the hub genes 
# focusing on the following modules: lightblue4, lightpink3, navajowhite2,blue4
# export for hub genes
FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata_Res$group <- "lightpink3"
Pmar_GO_hub_export <- rbind(FilterGenes_Pmar_comb_Interpro_ZVAD_lightpink3_GOdata_Res,
                            FilterGenes_Pmar_comb_Interpro_ZVAD_navajowhite2_GOdata_Res[,-8],
                            FilterGenes_Pmar_comb_Interpro_GDC_lightblue4_GOdata_Res[,-8]) %>% filter(topgoFisher <=0.05) %>% mutate(type = "MF")
GDC_lightblue4_GOdata_BP_Res$group <- "lightblue4"
ZVAD_lightpink3_GOdata_BP_Res$group <- "lightpink3"
ZVAD_navajowhite2_GOdata_BP_Res$group <- "navajowhite2"

# export the GO enrichment results for the full module
GDC_ZVAD_blue4_GOdata_Res$group <- "blue4"
GDC_ZVAD_blue4_GOdata_Res$treat <- "GDC_ZVAD"
Pmar_GO_all_export_MF <- rbind(GDC_lightblue4_GOdata_Res,
                               GDC_ZVAD_blue4_GOdata_Res,
                               ZVAD_lightpink3_GOdata_Res,
                               ZVAD_navajowhite2_GOdata_Res) %>% filter(topgoFisher <=0.05) %>% mutate(type = "MF")

GDC_lightblue4_GOdata_BP_Res$group <- "lightblue4"
GDC_ZVAD_blue4_GOdata_BP_Res$group <- "blue4"
ZVAD_lightpink3_GOdata_BP_Res$group <- "lightpink3"
ZVAD_navajowhite2_GOdata_BP_Res$group <- "navajowhite2"
Pmar_GO_all_export_BP <- rbind(GDC_lightblue4_GOdata_BP_Res,
                               GDC_ZVAD_blue4_GOdata_BP_Res,
                               ZVAD_lightpink3_GOdata_BP_Res,
                               ZVAD_navajowhite2_GOdata_BP_Res)  %>% filter(topgoFisher <=0.05) %>% mutate(type = "BP")

# combine all 
Pmar_GO_export_all <- rbind(Pmar_GO_all_export_MF[,-8], Pmar_GO_all_export_BP)

# Files to export
write.table(FilterGenes_comb_apop_labeled, file = "FilterGenes_comb_apop_labeled_HUB_GENES.txt",sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(Hemo_GO_export_subset_all, file = "Hemo_GO_export_subset_all.txt",sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(FilterGenes_Pmar_comb_Interpro_slim, file = "FilterGenes_Pmar_comb_Interpro_slim_HUB_GENES.txt",sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(Pmar_GO_export_all, file = "Pmar_GO_export_all.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

#### Select Pmar hub genes for further analysis based on overlap with GO enrichment terms ####

# In order to narrow down my list of interesting GO terms for Perkinsus, I will interrogate enriched GO terms and which enzymes of interest are in them
# specifically things like endopeptidases, carboxypeptidases, methylases, and any redox related enzymes 
# these will be my lists of interesting parts of the cytoscape network to interrogate
# Focusing on my modules of interest
# GDC related:
#   Lightblue4
# ZVAD related:
#   Navajowhite2
#   lightpink3
# GDC and ZVAD
Pmar_GO_hub_export_all_GO_list_lightblue4 <- Pmar_GO_export_all %>% filter(group == "lightblue4") %>% dplyr::select(GO.ID)
Pmar_GO_hub_export_all_GO_list_lightblue4 <- as.character(unlist(Pmar_GO_hub_export_all_GO_list_lightblue4))

Pmar_GO_hub_export_all_GO_list_blue4 <- Pmar_GO_export_all %>% filter(group == "blue4") %>% dplyr::select(GO.ID)
Pmar_GO_hub_export_all_GO_list_blue4 <- as.character(unlist(Pmar_GO_hub_export_all_GO_list_blue4))

Pmar_GO_hub_export_all_GO_list_lightpink3 <- Pmar_GO_export_all %>% filter(group == "lightpink3") %>% dplyr::select(GO.ID)
Pmar_GO_hub_export_all_GO_list_lightpink3 <- as.character(unlist(Pmar_GO_hub_export_all_GO_list_lightpink3))

Pmar_GO_hub_export_all_GO_list_navajowhite2 <- Pmar_GO_export_all %>% filter(group == "navajowhite2") %>% dplyr::select(GO.ID)
Pmar_GO_hub_export_all_GO_list_navajowhite2 <- as.character(unlist(Pmar_GO_hub_export_all_GO_list_navajowhite2))

# get lists of hub genes that overlap with enriched GO terms in either BP or MF
Pmar_GO_hub_export_all_Interpro_hub_genes_lightblue4_list <- FilterGenes_Pmar_comb_Interpro_slim %>% filter(mod_names == "MElightblue4") %>%
  filter(grepl(paste(Pmar_GO_hub_export_all_GO_list_lightblue4, collapse = "|"), Ontology_term))
    #1 XM_002768739.1  iron-sulfur cluster assembly protein, putative XM_002768739.1      0.9089457      0.0006854917 MElightblue4 0.7920236 XP_002768785.1    CDD
    #2 XM_002766308.1                       cytochrome P450, putative XM_002766308.1      0.9089457      0.0006854917 MElightblue4 0.7677165 XP_002766354.1   Pfam
    #3 XM_002766308.1                       cytochrome P450, putative XM_002766308.1      0.9089457      0.0006854917 MElightblue4 0.7677165 XP_002766354.1 PRINTS
    #4 XM_002784868.1 DNA replication licensing factor MCM2, putative XM_002784868.1      0.9089457      0.0006854917 MElightblue4 0.7757555 XP_002784914.1   Pfam

# join back in the full Interproscan terms to look at all the domains for these proteins
Pmar_GO_hub_export_all_Interpro_hub_genes_lightblue4_all_domain <- Pmar_GO_hub_export_all_Interpro_hub_genes_lightblue4_list %>%  
  left_join(., Perk_Interpro_GO_terms_XP[,c("transcript_id","Dbxref","signature_desc")], by = "transcript_id")

# repeat for blue4
Pmar_GO_hub_export_all_Interpro_hub_genes_blue4_list <- FilterGenes_Pmar_comb_Interpro_slim %>% filter(mod_names == "MEblue4") %>%
  filter(grepl(paste(Pmar_GO_hub_export_all_GO_list_blue4, collapse = "|"), Ontology_term))
#Name                        product  transcript_id moduleTraitCor moduleTraitPvalue mod_names        GS     protein_id source Ontology_term       Dbxref
#1 XM_002772524.1 CAAX prenyl protease, putative XM_002772524.1      0.8944067       0.001134529   MEblue4 0.6877459 XP_002772570.1   Pfam  "GO:0004.... "InterPr....
#signature_desc

# join back in the full Interproscan terms to look at all the domains for these proteins
Pmar_GO_hub_export_all_Interpro_hub_genes_blue4_all_domain <- Pmar_GO_hub_export_all_Interpro_hub_genes_blue4_list %>%  
  left_join(., Perk_Interpro_GO_terms_XP[,c("transcript_id","Dbxref","signature_desc")], by = "transcript_id")

# repeat for lightpink3
Pmar_GO_hub_export_all_Interpro_hub_genes_lightpink3_list <- FilterGenes_Pmar_comb_Interpro_slim %>% filter(mod_names == "MElightpink3") %>%
  filter(grepl(paste(Pmar_GO_hub_export_all_GO_list_lightpink3, collapse = "|"), Ontology_term))
#Name                             product  transcript_id moduleTraitCor moduleTraitPvalue    mod_names        GS     protein_id          source Ontology_term       Dbxref signature_desc
#1 XM_002776978.1                hypothetical protein XM_002776978.1      0.8357802       0.005012701 MElightpink3 0.7085296 XP_002777024.1          PRINTS  "GO:0004.... "InterPr....   Pepsin (....
#2 XM_002776978.1                hypothetical protein XM_002776978.1      0.8357802       0.005012701 MElightpink3 0.7085296 XP_002777024.1 ProSitePatterns  "GO:0004.... "InterPr....   Eukaryot....
#3 XM_002786856.1 proteasome subunit alpha1, putative XM_002786856.1      0.8357802       0.005012701 MElightpink3 0.7399695 XP_002786902.1            Pfam  "GO:0004.... "InterPr....   Proteaso....
#4 XM_002785313.1         aspartyl protease, putative XM_002785313.1      0.8357802       0.005012701 MElightpink3 0.7359564 XP_002785359.1          PRINTS  "GO:0004.... "InterPr....   Pepsin (....
#5 XM_002777539.1   succinate dehydrogenase, putative XM_002777539.1      0.8357802       0.005012701 MElightpink3 0.8321954 XP_002777585.1         TIGRFAM  "GO:0000.... "InterPr....   flavo_cy....
#6 XM_002765664.1      conserved hypothetical protein XM_002765664.1      0.8357802       0.005012701 MElightpink3 0.8670879 XP_002765710.1            Pfam  "GO:0004.... "InterPr....   Aminopep....

# join back in the full Interproscan terms to look at all the domains for these proteins
Pmar_GO_hub_export_all_Interpro_hub_genes_lightpink3_all_domain <- Pmar_GO_hub_export_all_Interpro_hub_genes_lightpink3_list %>%  
  left_join(., Perk_Interpro_GO_terms_XP[,c("transcript_id","Dbxref","signature_desc")], by = "transcript_id")
 
# repeat for navajowhite2
Pmar_GO_hub_export_all_Interpro_hub_genes_navajowhite2_list <- FilterGenes_Pmar_comb_Interpro_slim %>% filter(mod_names == "MEnavajowhite2") %>%
  filter(grepl(paste(Pmar_GO_hub_export_all_GO_list_navajowhite2, collapse = "|"), Ontology_term))
  #1  XM_002780795.1                     cytochrome p450, putative XM_002780795.1      0.9145422      0.0005521049 MEnavajowhite2 0.8081779 XP_002780841.1          PRINTS  "GO:0005.... "InterPr....
  #2  XM_002780795.1                     cytochrome p450, putative XM_002780795.1      0.9145422      0.0005521049 MEnavajowhite2 0.8081779 XP_002780841.1            Pfam  "GO:0005.... "InterPr....
  #3  XM_002788329.1                          hypothetical protein XM_002788329.1      0.9145422      0.0005521049 MEnavajowhite2 0.8636239 XP_002788375.1            Pfam  "GO:0022.... "InterPr....
  #4  XM_002784845.1  succinyl-coa synthetase beta chain, putative XM_002784845.1      0.9145422      0.0005521049 MEnavajowhite2 0.6890149 XP_002784891.1            Pfam  "GO:0003824" "InterPr....
  #5  XM_002765483.1      trehalose-6-phosphate synthase, putative XM_002765483.1      0.9145422      0.0005521049 MEnavajowhite2 0.9398572 XP_002765529.1         TIGRFAM  "GO:0003.... "InterPr....
  #6  XM_002765483.1      trehalose-6-phosphate synthase, putative XM_002765483.1      0.9145422      0.0005521049 MEnavajowhite2 0.9398572 XP_002765529.1            Pfam  "GO:0003.... "InterPr....
  #7  XM_002773250.1                         SEC61-gamma, putative XM_002773250.1      0.9145422      0.0005521049 MEnavajowhite2 0.8010266 XP_002773296.1           Hamap  "GO:0006.... "InterPr....
  #8  XM_002773011.1                           Myoglobin, putative XM_002773011.1      0.9145422      0.0005521049 MEnavajowhite2 0.9740561 XP_002773057.1            Pfam  "GO:0019825" "InterPr....
  #9  XM_002771060.1 transcription elongation factor SII, putative XM_002771060.1      0.9145422      0.0005521049 MEnavajowhite2 0.9560542 XP_002771106.1 ProSitePatterns  "GO:0003.... "InterPr....
  #10 XM_002772820.1               acetolactate synthase, putative XM_002772820.1      0.9145422      0.0005521049 MEnavajowhite2 0.8030358 XP_002772866.1            Pfam  "GO:0003.... "InterPr....
  #11 XM_002773752.1                          hypothetical protein XM_002773752.1      0.9145422      0.0005521049 MEnavajowhite2 0.8418256 XP_002773798.1 ProSiteProfiles  "GO:0003.... "InterPr....

# join back in the full Interproscan terms to look at all the domains for these proteins
Pmar_GO_hub_export_all_Interpro_hub_genes_navajowhite2_all_domain <- Pmar_GO_hub_export_all_Interpro_hub_genes_navajowhite2_list %>%  
  left_join(., Perk_Interpro_GO_terms_XP[,c("transcript_id","Dbxref","signature_desc")], by = "transcript_id")

#### EXPORT MODULE HUB GENE LISTS AND INFO FOR VIEW IN CYTOSCAPE ####

## Export alternate name lists
C_vir_rtracklayer_transcripts_cytoscape_alt_names <- unique(C_vir_rtracklayer_transcripts[,c("ID","product")])
Perkinsus_rtracklayer_transcripts_cytoscape_alt_names <- unique(Perkinsus_rtracklayer_transcripts[,c("transcript_id","product")])

write.table(C_vir_rtracklayer_transcripts_cytoscape_alt_names, file = "./Cytoscape_files/C_vir_rtracklayer_transcripts_cytoscape_alt_names.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(Perkinsus_rtracklayer_transcripts_cytoscape_alt_names, file = "./Cytoscape_files/Perkinsus_rtracklayer_transcripts_cytoscape_alt_names.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Hemocyte hub gene rnaID lists
#"MEnavajowhite2", "MEblue", "MEyellow","MEdarkslateblue", "MEorangered4"

hemo_hub_gene_ids_MEnavajowhite2 <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEnavajowhite2") %>% dplyr::select(ID)
hemo_hub_gene_ids_MEblue <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEblue") %>% dplyr::select(ID)
hemo_hub_gene_ids_MEyellow <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEyellow") %>% dplyr::select(ID)
hemo_hub_gene_ids_MEdarkslateblue <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEdarkslateblue") %>% dplyr::select(ID)
hemo_hub_gene_ids_MEorangered4 <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEorangered4") %>% dplyr::select(ID)

write.table(hemo_hub_gene_ids_MEnavajowhite2, file = "./Cytoscape_files/hemo_hub_gene_ids_MEnavajowhite2.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_hub_gene_ids_MEblue, file = "./Cytoscape_files/hemo_hub_gene_ids_MEblue.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_hub_gene_ids_MEyellow, file = "./Cytoscape_files/hemo_hub_gene_ids_MEyellow.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_hub_gene_ids_MEdarkslateblue, file = "./Cytoscape_files/hemo_hub_gene_ids_MEdarkslateblue.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hemo_hub_gene_ids_MEorangered4, file = "./Cytoscape_files/hemo_hub_gene_ids_MEorangered4.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# export as labeled list indicating hub genes and apop genes
hemo_hub_gene_ids_MEnavajowhite2_apop_hub <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEnavajowhite2") %>% mutate(hub = "yes") %>%
  dplyr::select(ID, apop, hub) %>% mutate(apop = case_when(is.na(apop) ~ "non_apop", TRUE ~ "apop"))
hemo_hub_gene_ids_MEblue_apop_hub <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEblue") %>% mutate(hub = "yes") %>%
  dplyr::select(ID, apop, hub) %>% mutate(apop = case_when(is.na(apop) ~ "non_apop", TRUE ~ "apop"))
hemo_hub_gene_ids_MEyellow_apop_hub <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEyellow") %>% mutate(hub = "yes") %>%
  dplyr::select(ID, apop, hub) %>% mutate(apop = case_when(is.na(apop) ~ "non_apop", TRUE ~ "apop"))
hemo_hub_gene_ids_MEdarkslateblue_apop_hub <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEdarkslateblue")%>% mutate(hub = "yes") %>%
  dplyr::select(ID, apop, hub) %>% mutate(apop = case_when(is.na(apop) ~ "non_apop", TRUE ~ "apop"))
hemo_hub_gene_ids_MEorangered4_apop_hub <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEorangered4") %>% mutate(hub = "yes") %>%
  dplyr::select(ID, apop, hub) %>% mutate(apop = case_when(is.na(apop) ~ "non_apop", TRUE ~ "apop"))

write.table(hemo_hub_gene_ids_MEnavajowhite2_apop_hub, file = "./Cytoscape_files/hemo_hub_gene_ids_MEnavajowhite2_apop_hub.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(hemo_hub_gene_ids_MEblue_apop_hub, file = "./Cytoscape_files/hemo_hub_gene_ids_MEblue_apop_hub.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(hemo_hub_gene_ids_MEyellow_apop_hub, file = "./Cytoscape_files/hemo_hub_gene_ids_MEyellow_apop_hub.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(hemo_hub_gene_ids_MEdarkslateblue_apop_hub, file = "./Cytoscape_files/hemo_hub_gene_ids_MEdarkslateblue_apop_hub.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(hemo_hub_gene_ids_MEorangered4_apop_hub, file = "./Cytoscape_files/hemo_hub_gene_ids_MEorangered4_apop_hub.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Create dummy list that also includes oxidoreducase terms for the navajowhite2 module
hemo_hub_gene_ids_MEnavajowhite2_apop_hub_oxidoreductase <- FilterGenes_comb_apop_labeled %>% filter(mod_names == "MEnavajowhite2") %>% mutate(hub = "yes") %>%
  dplyr::select(ID, apop, hub) %>% mutate(apop = case_when(is.na(apop) ~ "non_apop", TRUE ~ "apop"))
hemo_hub_gene_ids_MEnavajowhite2_apop_hub_oxidoreductase <- hemo_hub_gene_ids_MEnavajowhite2_apop_hub_oxidoreductase %>%
  mutate(apop = case_when(
    ID == "rna44333" ~ "apop", 
    ID == "rna31750" ~ "apop",
    ID == "rna6693"  ~ "apop",
    TRUE ~ apop))

# only rna6693 was changed because the others aren't hub genes...add these to the list
navajowhite2_oxido_cytochrome <- data.frame("ID" = c("rna44333", "rna31750", "rna24080", "rna46044"), "apop" = c("apop","apop","apop","apop"), "hub" = c("no","no","no","no"))

# also add the non hub gene apoptosis transcripts cdc42 homolog and caspase 8 
hemo_MEnavajowhite2_annot
  # "rna24080", "rna46044"

# join together
hemo_hub_gene_ids_MEnavajowhite2_apop_hub_oxidoreductase <- rbind(hemo_hub_gene_ids_MEnavajowhite2_apop_hub_oxidoreductase, navajowhite2_oxido_cytochrome)

write.table(hemo_hub_gene_ids_MEnavajowhite2_apop_hub_oxidoreductase, file = "./Cytoscape_files/hemo_hub_gene_ids_MEnavajowhite2_apop_hub_oxidoreductase.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


## Repeat for Perkinsus: modules of interest lightblue4, lightpink3, navajowhite2,
perk_hub_gene_ids_MElightblue4 <- FilterGenes_Pmar_comb %>% filter(mod_names == "MElightblue4") %>% dplyr::select(Name)
perk_hub_gene_ids_MElightpink3 <- FilterGenes_Pmar_comb %>% filter(mod_names == "MElightpink3") %>% dplyr::select(Name)
perk_hub_gene_ids_MEnavajowhite2 <- FilterGenes_Pmar_comb %>% filter(mod_names == "MEnavajowhite2") %>% dplyr::select(Name)
perk_hub_gene_ids_MEblue4 <- FilterGenes_Pmar_comb %>% filter(mod_names == "MEblue4") %>% dplyr::select(Name) %>% distinct(Name)

write.table(perk_hub_gene_ids_MElightblue4, file = "./Cytoscape_files/perk_hub_gene_ids_MElightblue4.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(perk_hub_gene_ids_MElightpink3, file = "./Cytoscape_files/perk_hub_gene_ids_MElightpink3.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(perk_hub_gene_ids_MEnavajowhite2, file = "./Cytoscape_files/perk_hub_gene_ids_MEnavajowhite2.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(perk_hub_gene_ids_MEblue4 , file = "./Cytoscape_files/perk_hub_gene_ids_MEblue4.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# export with labels for hub vs enzymes of interest to use in cytoscape
# join with the lists for the enzymes of interest from "Select Pmar hub genes based on overlap with GO enrichment terms"
# just like I subset the total hemocyte networks based on apoptosis, here I will subset the 
perk_hub_gene_ids_MElightblue4_hub <- FilterGenes_Pmar_comb %>% filter(mod_names == "MElightblue4") %>% mutate(hub = "yes") %>%
  dplyr::select(Name,hub)
perk_hub_gene_ids_MElightpink3_hub <- FilterGenes_Pmar_comb %>% filter(mod_names == "MElightpink3") %>% mutate(hub = "yes") %>%
  dplyr::select(Name, hub)
perk_hub_gene_ids_MEnavajowhite2_hub <- FilterGenes_Pmar_comb %>% filter(mod_names == "MEnavajowhite2") %>% mutate(hub = "yes") %>%
  dplyr::select(Name,hub) 
perk_hub_gene_ids_MEblue4_hub <- FilterGenes_Pmar_comb %>% filter(mod_names == "MEblue4") %>% mutate(hub = "yes") %>%
  dplyr::select(Name,hub) 

write.table(perk_hub_gene_ids_MElightblue4_hub, file = "./Cytoscape_files/perk_hub_gene_ids_MElightblue4_hub.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(perk_hub_gene_ids_MElightpink3_hub, file = "./Cytoscape_files/perk_hub_gene_ids_MElightpink3_hub.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(perk_hub_gene_ids_MEnavajowhite2_hub, file = "./Cytoscape_files/perk_hub_gene_ids_MEnavajowhite2_hub.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(perk_hub_gene_ids_MEblue4_hub, file = "./Cytoscape_files/perk_hub_gene_ids_MEblue4_hub.txt",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## EXPORT GENE SIGNIFICANCE FOR APOPTOSIS PERCENTAGE FOR EACH GENE
hemo_full_geneTraitSignificance_apop_cytoscape <- hemo_full_geneTraitSignificance_apop %>% rownames_to_column(.,var = "ID")
hemo_full_GSPvalue_apop_cytoscape <- hemo_full_GSPvalue_apop %>% rownames_to_column(.,var = "ID")
perk_full_geneTraitSignificance_apop_cytoscape <- perk_full_geneTraitSignificance_apop %>% rownames_to_column(.,var = "ID")
perk_full_GSPvalue_apop_cytoscape <- perk_full_GSPvalue_apop %>% rownames_to_column(.,var = "ID")

hemo_full_geneTraitSignificance_GSPvalue <- full_join(hemo_full_geneTraitSignificance_apop_cytoscape,hemo_full_GSPvalue_apop_cytoscape )
perk_full_geneTraitSignificance_GSPvalue <- full_join(perk_full_geneTraitSignificance_apop_cytoscape,perk_full_GSPvalue_apop_cytoscape )

write.table(hemo_full_geneTraitSignificance_GSPvalue,  file = "./Cytoscape_files/hemo_full_geneTraitSignificance_GSPvalue.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(perk_full_geneTraitSignificance_GSPvalue , file = "./Cytoscape_files/perk_full_geneTraitSignificance_GSPvalue.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#### DETERMINE EDGE WEIGHT CUT OFF TO VIEW MOST INTERESTING CYTOSCAPE GENE CONNECTIONS ####

# What is the top 5% of edge weights for each file?
CytoscapeInput_edges_perk_fulllightblue4 <- read.table(file="./Cytoscape_files/CytoscapeInput-edges-perk_fulllightblue4.txt", header = TRUE, col.names = c("fromNode","toNode","weight", "direction", "fromAltName", "toAltName"))
CytoscapeInput_edges_perk_fullblue4 <- read.table(file="./Cytoscape_files/CytoscapeInput-edges-perk_fullblue4.txt", header = TRUE, col.names = c("fromNode","toNode","weight", "direction", "fromAltName", "toAltName"))

CytoscapeInput_edges_hemo_fullyellow <- read.table(file="./Cytoscape_files/CytoscapeInput-edges-hemo_fullyellow.txt", header = TRUE,col.names = c("fromNode","toNode","weight", "direction", "fromAltName", "toAltName"))
CytoscapeInput_edges_hemo_fullnavajowhite2 <- read.table(file="./Cytoscape_files/CytoscapeInput-edges-hemo_fullnavajowhite2.txt", header = TRUE,col.names = c("fromNode","toNode","weight", "direction", "fromAltName", "toAltName"))
CytoscapeInput_edges_hemo_fullblue <- read.table(file="./Cytoscape_files/CytoscapeInput-edges-hemo_fullblue.txt", header = TRUE, col.names = c("fromNode","toNode","weight", "direction", "fromAltName", "toAltName"))

CytoscapeInput_edges_perk_fulllightblue4$weight <- as.numeric(CytoscapeInput_edges_perk_fulllightblue4$weight)
CytoscapeInput_edges_perk_fullblue4$weight <- as.numeric(CytoscapeInput_edges_perk_fullblue4$weight)

CytoscapeInput_edges_hemo_fullyellow$weight <- as.numeric(CytoscapeInput_edges_hemo_fullyellow$weight)
CytoscapeInput_edges_hemo_fullnavajowhite2$weight <- as.numeric(CytoscapeInput_edges_hemo_fullnavajowhite2$weight)
CytoscapeInput_edges_hemo_fullblue$weight <- as.numeric(CytoscapeInput_edges_hemo_fullblue$weight)

quantile(CytoscapeInput_edges_perk_fulllightblue4$weight, 0.80) # 0.1320015
quantile(CytoscapeInput_edges_perk_fullblue4$weight, 0.80) #  0.1121706

quantile(CytoscapeInput_edges_hemo_fullyellow$weight, 0.90) #
quantile(CytoscapeInput_edges_hemo_fullnavajowhite2$weight, .9990) #  0.09
quantile(CytoscapeInput_edges_hemo_fullblue$weight, 0.90) #  0.3489623 


#### ANALYSIS OF CYTOSCAPE RESULTS FOR INTERESTING MODULES ####

# WITH APOP GENE SIGNIFICANCE P VALUE FILTERED TO 0.05 AND EDGE WEIHT FILTERED TO 80 PERCENTILE
# we get the following Perkinsus nodes in the lightblue4 module
perk_lightblue3_cytoscape_P05_EW80 <- read.csv("/Users/erinroberts/Documents/PhD_Research/DERMO_EXP_18_19/2020_Hemocyte_experiment/2020_Dermo_Inhibitors_main_exp/ANALYSIS_FILES/WGCNA/Cytoscape_files/CytoscapeInput-edges-perk_fulllightblue4.txt default node.csv", header= TRUE)
# filter for selected 
perk_lightblue3_cytoscape_P05_EW80_selected <- perk_lightblue3_cytoscape_P05_EW80 %>% filter(selected == "true") %>% dplyr::rename("transcript_id" = "name") %>% 
  left_join(., GDC_lightblue4_all_genes_annot)

# any overlap with GDC DEG?
perk_lightblue3_cytoscape_P05_EW80_selected$transcript_id %in% perk_dds_deseq_res_Pmar_GDC_LFC_sig_annot$transcript_id # no Overlap with DEGs!
XM_002770984.1

# any overlap at all with enriched GO terms?
perk_lightblue3_cytoscape_P05_EW80_selected %>%
  filter(grepl(paste(Pmar_GO_hub_export_all_GO_list_lightblue4, collapse = "|"), Ontology_term))
  # just the iron-sulfur cluster assembly protein iron-sulfur cluster assembly protein, putative 

## What are the most significant for apoptosis navajowhite2 genes?
hemo_MEnavajowhite2_annot <- left_join(hemo_MEnavajowhite2_annot, hemo_full_geneTraitSignificance_GSPvalue)

## What are the oxidoreductase genes in the navajowhite2 modules 
hemo_MEnavajowhite2_GOdata_Res
  # oxidoreductase activity terms: GO:0016705,GO:0016628

# Interproscan of all genes to get GO terms
hemo_MEnavajowhite2_annot_Interpro <- GO_universe_rna_found[GO_universe_rna_found$transcript_id %in% hemo_MEnavajowhite2_annot$ID,]

## look up in the list of hub genes for navajowhite2
# GO:0016705
hemo_MEnavajowhite2_annot_Interpro %>% filter(grepl("GO:0016705", Ontology_term))
  # two transcripts: rna44333, rna31750
hemo_MEnavajowhite2_annot %>% filter(ID == "rna44333" | ID == "rna31750") %>% dplyr::select(ID, product,gene)
 # ID                   product         gene
 # 1 rna44333 cytochrome P450 4F22-like LOC111099622
 # 2 rna31750 cytochrome P450 2D27-like LOC111137052

#GO:0016628
hemo_MEnavajowhite2_annot_Interpro %>% filter(grepl("GO:0016628", Ontology_term)) %>% View()
  # 1 transcript: rna6693
hemo_MEnavajowhite2_annot %>% filter(ID == "rna6693") %>% dplyr::select(ID, product,gene)
#ID                                                    product         gene
#1 rna6693 7-dehydrocholesterol reductase-like, transcript variant X2 LOC111121358
# Contains the following domains  
    # Sterol reductase family signature 1. = IPR018083



## Are either of these oxidoreductase genes also hub genes for apoptosis phenotype?

FilterGenes_control_Pmar_GDC_navajowhite2_df %>% filter(ID == "rna44333" | ID == "rna31750" | ID == "rna6693")
  # rna6693 is also a hub gene and is connected to the apoptosis genes 

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

#### MULTIPANEL FIGURE OF WGCNA RAW DATA FOR PAPER SUPPLEMENTARY FIGURE ####

