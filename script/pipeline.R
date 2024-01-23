################################################################################
############### DBscan cluster + FOCALs not in CLUSTERS ########################
##############################################################################
################################################################################

library(tidyverse)
library(data.table)
library(readxl)

library(dbscan)
library(factoextra)
library(GenomicRanges)
library(proxy)

source("script/format_segment.R")
source("script/DBscan_clustering.R")

#data preparation
df <- fread("D:/analisi_in_corso/clonal_evolution/Clonal_Evolution/workfiles/FIRST_RemasterCNA_correction/5.geneCN_CALLS_samples_FINAL.txt") %>% as.data.frame()

focal <- fread("D:/analisi_in_corso/clonal_evolution/Clonal_Evolution/workfiles/FIRST_RemasterCNA_correction/focal_loci_hg19.txt") 


focal <- format_segments(focal)
format_segments(df)

outpath <- "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT Clonal Evolution - Documenti/DBscan_approach/final_plots/"

#------------------------------ clustering and HR focal analysis ---------------------------


clust_res <- DBscan_clustering (df, 
                        focal, 
                        pt_code="PAZ_[0-9]+", 
                        phase_code=c("_E", "_R"), 
                        NOISE=0.04, 
                        EPS_pair=0.03, 
                        MIN_pair=30,
                        wt_lim= c(1.8, 2.2),
                        plot=T, 
                        outpath=outpath)


#write_tsv(clust_res$cluster_centers, "workfiles/FIRST_RemasterCNA_correction/6.DBscan_long_results.txt")
#write_tsv(clust_res$focals_not_in_clusters, "workfiles/FIRST_RemasterCNA_correction/6.focals_not_in_clusters_results.txt")


# buttleship

states <- c("hT",
            "T",
            "AT",
            "A",
            "sA",
            "wt",
            "sD",
            "D",
            "HD",
            "H") %>% rev

limits_low <- c(-Inf,
                0.2,
                0.8,
                1.2,
                1.8,
                2.2,
                2.8,
                3.2,
                3.8,
                4.2)

limits_high <- c(0.2,
                 0.8,
                 1.2,
                 1.8,
                 2.2,
                 2.8,
                 3.2,
                 3.8,
                 4.2,
                 Inf)


all_states <- data.frame(states, limits_low, limits_high)

#######################################################################
################################ CLUSTER ##############################
#######################################################################

#import pts data and assign clusters 

pts_clusters <- fread("workfiles/FIRST_RemasterCNA_correction/6.DBscan_long_results.txt")

pts <- pts_clusters$pt_name %>% unique 

i=1
j=2
res <- data.frame()
for(i in 1:length(pts)) {
  
  print(pts[i])
  pt <- pts_clusters %>% filter(pt_name==pts[i])
  pt$relapse_label <- NA
  pt$diagnosis_label <- NA
  
  for( j in 1:nrow(pt)) {
    
    cluster_diag_CN <- pt$diagnosis_noise[j]
    pt$diagnosis_label[j] <- all_states %>% filter(cluster_diag_CN > limits_low & cluster_diag_CN < limits_high) %>% .$states
    
    cluster_relapse_CN <- pt$relapse_noise[j]
    pt$relapse_label[j] <- all_states %>% filter(cluster_relapse_CN > limits_low & cluster_relapse_CN < limits_high) %>% .$states
    
  }
  
  res <- rbind(res, pt)
}

write_tsv(res, "workfiles/FIRST_RemasterCNA_correction/7.DBscan_res_buttleship.txt")


########################################################
###################### FOCAL ###########################
########################################################


#import data

df <- fread("workfiles/FIRST_RemasterCNA_correction/6.focals_not_in_clusters_results.txt")

pts <- df$pt_name %>% unique

i=1
pts[66]
j=1
res <- data.frame()
for(i in 1:length(pts)) {
  
  print(i)
  
  #idx_D <- (i*2 - 1) + 6
  #idx_R <- (i*2) + 6
  #
  #D_name <- names(df)[idx_D]
  #R_name <- names(df)[idx_R]
  
  
  dat_pt <- df %>% filter(pt_name == pts[i]) #%>%  select(Chromosome=seqnames, Gene=gene_names, diagnosis = idx_D, relapse=idx_R)
  
  #dat_pt$Chromosome <- as.character(dat_pt$Chromosome)
  dat_pt$relapse_label <- NA
  dat_pt$diagnosis_label <- NA
  
  for( j in 1:nrow(dat_pt)) {
    
    cluster_diag_CN <- dat_pt$diagnosis[j]
    dat_pt$diagnosis_label[j] <- all_states %>% filter(cluster_diag_CN >= limits_low & cluster_diag_CN < limits_high) %>% .$states
    
    cluster_relapse_CN <- dat_pt$relapse[j]
    dat_pt$relapse_label[j] <- all_states %>% filter(cluster_relapse_CN >= limits_low & cluster_relapse_CN < limits_high) %>% .$states
    
    #dat_pt$pt_idx <- pts[i]
  }
  
  res <- rbind(res, dat_pt)
}

write_tsv(res, "workfiles/FIRST_RemasterCNA_correction/7.Focal_lesion_buttleship.txt")

