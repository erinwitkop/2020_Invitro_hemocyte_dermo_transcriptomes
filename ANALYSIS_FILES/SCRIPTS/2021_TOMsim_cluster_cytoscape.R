## R script to run TOMsimilarityFromExpr on the bluewaves cluster and export cytoscape networks
## Erin Michele Roberts, PhD Candidate

# Load libraries 
library(WGCNA)

options(stringsAsFactors = FALSE) # run every time
allowWGCNAThreads()
cor <- WGCNA::cor

#### LOAD SAVED DATA ####

# Read in data
hemo_dds_rlog_matrix <- read.table(file="/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/WGCNA/hemo_dds_rlog_matrix.table")
perk_dds_rlog_matrix <- read.table(file="/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/WGCNA/perk_dds_rlog_matrix.table")

# Load module colors 
load(file = "/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/WGCNA/hemo_full_moduleColors.RData")
load(file = "/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/WGCNA/perk_full_moduleColors.RData")

# Calculate TOM
hemo_full_TOM = TOMsimilarityFromExpr(hemo_dds_rlog_matrix, power = 7, TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 

save(hemo_full_TOM, file = "/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/WGCNA/hemo_full_TOM.RData")

perk_full_TOM = TOMsimilarityFromExpr(perk_dds_rlog_matrix, power = 7, TOMType = "signed", networkType= "signed hybrid", corType ="bicor" )

save(perk_full_TOM, file="/data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/WGCNA/perk_full_TOM.RData" )

###### CAN RUN BELOW EITHER INTERACTIVE OR AS A SCRIPT ####

## EXPORTING EACH AS SEPARATE NETWORKS BECAUSE THE INDIVIDUAL NETWORKS ARE TOO LARGE 
# $ pwd 
# /data/marine_diseases_lab/erin/2020_Hemolymph_Dermo_Transcriptome_Project/WGCNA/
# $ interactive 
# $ module load R/3.6.0-intel-2019a
# $ R

library(WGCNA)

# Read in the annotation files 
load(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/C_gig_C_vir_annotations.RData")

# Check they loaded
# nrow(C_vir_rtracklayer)

## EXPORT HEMOCYTE MODULES ##
  # EXPORT: 
      # GDC navajowhite2, 
      # Pmar vs control = blue, yellow
      # ZVAD = darkslateblue, orangered4

## EXPORT navajowhite2
hemo_full_modules = "navajowhite2"
# Select module probes
hemo_full_probes = colnames(hemo_dds_rlog_matrix)
hemo_full_inModule = is.finite(match(hemo_full_moduleColors, hemo_full_modules))
hemo_full_modProbes = hemo_full_probes[hemo_full_inModule]
hemo_full_modGenes = C_vir_rtracklayer$ID[match(hemo_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
hemo_full_modTOM = hemo_full_TOM[hemo_full_inModule, hemo_full_inModule]
dimnames(hemo_full_modTOM) = list(hemo_full_modProbes, hemo_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
hemo_full_cyt = exportNetworkToCytoscape(hemo_full_modTOM,
                                         edgeFile = paste("CytoscapeInput-edges-hemo_full", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
                                         nodeFile = paste("CytoscapeInput-nodes-hemo_full", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
                                         weighted = TRUE,
                                         threshold = 0.00, # using 0 threshold so no genes are subset out 
                                         nodeNames = hemo_full_modProbes,
                                         altNodeNames = hemo_full_modGenes,
                                         nodeAttr = hemo_full_moduleColors[hemo_full_inModule])
# The command writes results to file automatically

## EXPORT blue
hemo_full_modules = "blue"
# Select module probes
hemo_full_probes = colnames(hemo_dds_rlog_matrix)
hemo_full_inModule = is.finite(match(hemo_full_moduleColors, hemo_full_modules))
hemo_full_modProbes = hemo_full_probes[hemo_full_inModule]
hemo_full_modGenes = C_vir_rtracklayer$ID[match(hemo_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
hemo_full_modTOM = hemo_full_TOM[hemo_full_inModule, hemo_full_inModule]
dimnames(hemo_full_modTOM) = list(hemo_full_modProbes, hemo_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
hemo_full_cyt = exportNetworkToCytoscape(hemo_full_modTOM,
                                         edgeFile = paste("CytoscapeInput-edges-hemo_full", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
                                         nodeFile = paste("CytoscapeInput-nodes-hemo_full", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
                                         weighted = TRUE,
                                         threshold = 0.00, # using 0 threshold so no genes are subset out 
                                         nodeNames = hemo_full_modProbes,
                                         altNodeNames = hemo_full_modGenes,
                                         nodeAttr = hemo_full_moduleColors[hemo_full_inModule])
## EXPORT yellow
hemo_full_modules = "yellow"
# Select module probes
hemo_full_probes = colnames(hemo_dds_rlog_matrix)
hemo_full_inModule = is.finite(match(hemo_full_moduleColors, hemo_full_modules))
hemo_full_modProbes = hemo_full_probes[hemo_full_inModule]
hemo_full_modGenes = C_vir_rtracklayer$ID[match(hemo_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
hemo_full_modTOM = hemo_full_TOM[hemo_full_inModule, hemo_full_inModule]
dimnames(hemo_full_modTOM) = list(hemo_full_modProbes, hemo_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
hemo_full_cyt = exportNetworkToCytoscape(hemo_full_modTOM,
                                         edgeFile = paste("CytoscapeInput-edges-hemo_full", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
                                         nodeFile = paste("CytoscapeInput-nodes-hemo_full", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
                                         weighted = TRUE,
                                         threshold = 0.00, # using 0 threshold so no genes are subset out 
                                         nodeNames = hemo_full_modProbes,
                                         altNodeNames = hemo_full_modGenes,
                                         nodeAttr = hemo_full_moduleColors[hemo_full_inModule])

## EXPORT darkslateblue
hemo_full_modules = "darkslateblue"
# Select module probes
hemo_full_probes = colnames(hemo_dds_rlog_matrix)
hemo_full_inModule = is.finite(match(hemo_full_moduleColors, hemo_full_modules))
hemo_full_modProbes = hemo_full_probes[hemo_full_inModule]
hemo_full_modGenes = C_vir_rtracklayer$ID[match(hemo_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
hemo_full_modTOM = hemo_full_TOM[hemo_full_inModule, hemo_full_inModule]
dimnames(hemo_full_modTOM) = list(hemo_full_modProbes, hemo_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
hemo_full_cyt = exportNetworkToCytoscape(hemo_full_modTOM,
                                         edgeFile = paste("CytoscapeInput-edges-hemo_full", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
                                         nodeFile = paste("CytoscapeInput-nodes-hemo_full", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
                                         weighted = TRUE,
                                         threshold = 0.00, # using 0 threshold so no genes are subset out 
                                         nodeNames = hemo_full_modProbes,
                                         altNodeNames = hemo_full_modGenes,
                                         nodeAttr = hemo_full_moduleColors[hemo_full_inModule])

## EXPORT orangered4
hemo_full_modules = "orangered4"
# Select module probes
hemo_full_probes = colnames(hemo_dds_rlog_matrix)
hemo_full_inModule = is.finite(match(hemo_full_moduleColors, hemo_full_modules))
hemo_full_modProbes = hemo_full_probes[hemo_full_inModule]
hemo_full_modGenes = C_vir_rtracklayer$ID[match(hemo_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
hemo_full_modTOM = hemo_full_TOM[hemo_full_inModule, hemo_full_inModule]
dimnames(hemo_full_modTOM) = list(hemo_full_modProbes, hemo_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
hemo_full_cyt = exportNetworkToCytoscape(hemo_full_modTOM,
                                         edgeFile = paste("CytoscapeInput-edges-hemo_full", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
                                         nodeFile = paste("CytoscapeInput-nodes-hemo_full", paste(hemo_full_modules, collapse="-"), ".txt", sep=""),
                                         weighted = TRUE,
                                         threshold = 0.00, # using 0 threshold so no genes are subset out 
                                         nodeNames = hemo_full_modProbes,
                                         altNodeNames = hemo_full_modGenes,
                                         nodeAttr = hemo_full_moduleColors[hemo_full_inModule])

## EXPORT PERKINSUS MODULES ###
# Export Perkinsus GDC:lightblue4
# Export Perkinsus ZVAD : lightpink3, navajowhite2,

## EXPORT LIGHTBLUE4
perk_full_modules = "lightblue4"
# Select module probes
perk_full_probes = colnames(perk_dds_rlog_matrix)
perk_full_inModule = is.finite(match(perk_full_moduleColors, perk_full_modules))
perk_full_modProbes = perk_full_probes[perk_full_inModule]
perk_full_modGenes = C_vir_rtracklayer$ID[match(perk_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
perk_full_modTOM = perk_full_TOM[perk_full_inModule, perk_full_inModule]
dimnames(perk_full_modTOM) = list(perk_full_modProbes, perk_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
perk_full_cyt = exportNetworkToCytoscape(perk_full_modTOM,
                                              edgeFile = paste("CytoscapeInput-edges-perk_full", paste(perk_full_modules, collapse="-"), ".txt", sep=""),
                                              nodeFile = paste("CytoscapeInput-nodes-perk_full", paste(perk_full_modules, collapse="-"), ".txt", sep=""),
                                              weighted = TRUE,
                                              threshold = 0.00, # using 0 threshold so no genes are subset out 
                                              nodeNames = perk_full_modProbes,
                                              altNodeNames = perk_full_modGenes,
                                              nodeAttr = perk_full_moduleColors[perk_full_inModule])
# The command writes results to file automatically

## EXPORT lightpink3
perk_full_modules = "lightpink3"
# Select module probes
perk_full_probes = colnames(perk_dds_rlog_matrix)
perk_full_inModule = is.finite(match(perk_full_moduleColors, perk_full_modules))
perk_full_modProbes = perk_full_probes[perk_full_inModule]
perk_full_modGenes = C_vir_rtracklayer$ID[match(perk_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
perk_full_modTOM = perk_full_TOM[perk_full_inModule, perk_full_inModule]
dimnames(perk_full_modTOM) = list(perk_full_modProbes, perk_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
perk_full_cyt = exportNetworkToCytoscape(perk_full_modTOM,
                                         edgeFile = paste("CytoscapeInput-edges-perk_full", paste(perk_full_modules, collapse="-"), ".txt", sep=""),
                                         nodeFile = paste("CytoscapeInput-nodes-perk_full", paste(perk_full_modules, collapse="-"), ".txt", sep=""),
                                         weighted = TRUE,
                                         threshold = 0.00, # using 0 threshold so no genes are subset out 
                                         nodeNames = perk_full_modProbes,
                                         altNodeNames = perk_full_modGenes,
                                         nodeAttr = perk_full_moduleColors[perk_full_inModule])

## EXPORT navajowhite2
perk_full_modules = "navajowhite2"
# Select module probes
perk_full_probes = colnames(perk_dds_rlog_matrix)
perk_full_inModule = is.finite(match(perk_full_moduleColors, perk_full_modules))
perk_full_modProbes = perk_full_probes[perk_full_inModule]
perk_full_modGenes = C_vir_rtracklayer$ID[match(perk_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
perk_full_modTOM = perk_full_TOM[perk_full_inModule, perk_full_inModule]
dimnames(perk_full_modTOM) = list(perk_full_modProbes, perk_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
perk_full_cyt = exportNetworkToCytoscape(perk_full_modTOM,
                                         edgeFile = paste("CytoscapeInput-edges-perk_full", paste(perk_full_modules, collapse="-"), ".txt", sep=""),
                                         nodeFile = paste("CytoscapeInput-nodes-perk_full", paste(perk_full_modules, collapse="-"), ".txt", sep=""),
                                         weighted = TRUE,
                                         threshold = 0.00, # using 0 threshold so no genes are subset out 
                                         nodeNames = perk_full_modProbes,
                                         altNodeNames = perk_full_modGenes,
                                         nodeAttr = perk_full_moduleColors[perk_full_inModule])



sessionInfo()
#sessionInfo()
#R version 3.6.0 (2019-04-26)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: CentOS release 6.5 (Final)

#Matrix products: default
#BLAS:   /net/clusterhn.cluster.com/opt/software/R/3.6.0-intel-2019a/lib64/R/lib/libR.so
#LAPACK: /net/clusterhn.cluster.com/opt/software/R/3.6.0-intel-2019a/lib64/R/modules/lapack.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#[3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#
#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#  [1] WGCNA_1.69            fastcluster_1.1.25    dynamicTreeCut_1.63-1
#
#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.1            lattice_0.20-38       GO.db_3.10.0         
#[4] assertthat_0.2.1      digest_0.6.19         foreach_1.4.4        
#[7] R6_2.4.0              plyr_1.8.4            backports_1.1.4      
#[10] acepack_1.4.1         stats4_3.6.0          RSQLite_2.1.1        
#[13] ggplot2_3.1.1         pillar_1.4.1          rlang_0.3.4          
#[16] lazyeval_0.2.2        rstudioapi_0.10       data.table_1.12.2    
#[19] blob_1.1.1            S4Vectors_0.24.4      rpart_4.1-15         
#[22] Matrix_1.2-17         preprocessCore_1.48.0 checkmate_1.9.3      
#[25] splines_3.6.0         stringr_1.4.0         foreign_0.8-71       
#[28] htmlwidgets_1.3       bit_1.1-14            munsell_0.5.0        
#[31] compiler_3.6.0        xfun_0.7              pkgconfig_2.0.2      
#[34] BiocGenerics_0.32.0   base64enc_0.1-3       htmltools_0.3.6      
#[37] nnet_7.3-12           tidyselect_0.2.5      tibble_2.1.3         
#[40] gridExtra_2.3         htmlTable_1.13.1      Hmisc_4.2-0          
#[43] IRanges_2.20.2        codetools_0.2-16      matrixStats_0.54.0   
#[46] crayon_1.3.4          dplyr_0.8.1           grid_3.6.0           
#[49] gtable_0.3.0          DBI_1.0.0             magrittr_1.5         
#[52] scales_1.0.0          stringi_1.4.3         impute_1.60.0        
#[55] doParallel_1.0.14     latticeExtra_0.6-28   Formula_1.2-3        
#[58] RColorBrewer_1.1-2    iterators_1.0.10      tools_3.6.0          
#[61] bit64_0.9-7           Biobase_2.46.0        glue_1.3.1           
#[64] purrr_0.3.2           parallel_3.6.0        survival_2.44-1.1    
#[67] AnnotationDbi_1.48.0  colorspace_1.4-1      cluster_2.0.9        
#[70] memoise_1.1.0         knitr_1.23  #