
library(tidyverse)
library(data.table)
library(GenomicRanges)

library(dbscan)
library(factoextra)

library(proxy)


DBscan_clustering <- function(df, 
                              focal, 
                              pt_code="PAZ_[0-9]+", 
                              phase_code=c("_E", "_R"), 
                              NOISE=0.04, 
                              EPS_pair=0.03, 
                              MIN_pair=30,
                              wt_lim= c(1.8, 2.2),
                              plot=T, 
                              outpath=".") {
  
  set.seed(123)
  output <- list()
  # create output folder and files
  ALL_CLUSTERS_CENTERS_COMPLETE_list <- list()
  
  ALL_FOCALS_ALONE_COMPLETE_list <- list()
  
  #dir.create(outpath,showWarnings = F)
  source("script/format_segment.R")
  
  focal <- format_segments(focal)
  focalGR <- makeGRangesFromDataFrame(focal, keep.extra.columns = T)
  
  
  
  i=2
  
  pts <- unique(str_extract(colnames(df), pattern = pt_code))
  pts <- pts[!is.na(pts)]
  
  df <- format_segments(df)
  
  for(i in 1:133) {
    
    print(pts[i])
    sample_D = paste0(pts[i],phase_code[1])
    sample_R = paste0(pts[i],phase_code[2])
    
    idx_D <- which(grepl(colnames(df), pattern=sample_D))
    idx_R <- which(grepl(colnames(df), pattern=sample_R))
    
    dat_pt <- df %>% select(chr, start, end, gene_names, diagnosis=all_of(idx_D), relapse=all_of(idx_R))
    
    dat_pt$chr <- as.character(dat_pt$chr)
    
    
    noise_x <- rnorm(1:nrow(df), mean = 0, sd = NOISE)
    noise_y <- rnorm(1:nrow(df), mean = 0, sd = NOISE)
    
    plot(noise_x %>% sort)
    abline(h=0.07)
    abline(h=-0.07)
    
    dat_pt$diagnosis_noise <- dat_pt$diagnosis + noise_x
    dat_pt$relapse_noise <- dat_pt$relapse + noise_y
    
    if (plot) {
      dat_pt_wt <- dat_pt %>% filter( ((diagnosis_noise> wt_lim[1] & diagnosis_noise<wt_lim[2]) & (relapse_noise>wt_lim[1] & relapse_noise<wt_lim[2])) )
    }
    
    dat_pt <- dat_pt %>% filter( !((diagnosis_noise>wt_lim[1] & diagnosis_noise<wt_lim[2]) & (relapse_noise>wt_lim[1] & relapse_noise<wt_lim[2])) )
    
    
    #==================== DBscan analysis ======================
    
    df_dbscan <- dat_pt %>% select(diagnosis_noise, relapse_noise)
    
    #______________ dbscan diagnosis and relapse _______________ 
    
    dbs <- dbscan(df_dbscan, eps = EPS_pair, minPts = MIN_pair, borderPoints = F)
    
    dat_pt$DBscan_cluster <- dbs$cluster %>% as.character() 
    
    #dat_pt$Nb_cluster <- dbs$cluster %>% as.character() 
    
    #compute cluster centers
    cluster_centers <- dat_pt %>% filter(DBscan_cluster != 0) %>%
      group_by(DBscan_cluster) %>%
      summarise(diagnosis_noise=mean(diagnosis_noise), 
                relapse_noise=mean(relapse_noise),
                n=n(),
                chrs_list=paste(
                  unique(chr), 
                  collapse = " "),
                chrs_genes=paste(
                  table(chr), 
                  collapse = " ")
      )
    
    # long-df where each chromosome is assigned a cluster
    cluster_centers_complete2<- dat_pt %>%
      group_by(DBscan_cluster) %>%
      summarise(diagnosis_noise=mean(diagnosis_noise),
                relapse_noise=mean(relapse_noise),
                n=n(),
                chr= unique(chr),
                n_genes = table(as.numeric(chr)))
    
    
    cluster_centers_complete2$pt_idx <- i
    cluster_centers_complete2$pt_name <- pts[i]
    
    ALL_CLUSTERS_CENTERS_COMPLETE_list[[i]] <- cluster_centers_complete2
    
    
    #=============== match with focal regions =================
    # looking for some HR gene out of clusters
    
    dat_ptGR <- dat_pt %>% makeGRangesFromDataFrame(keep.extra.columns = T) #focalGR
    
    focal_calls <- plyranges::join_overlap_intersect(focalGR, dat_ptGR) %>% as.data.frame()
    
    
    focal_calls_0 <- focal_calls %>% filter(DBscan_cluster==0)
    
    data_frame1 = focal_calls_0 %>% select(diagnosis_noise, relapse_noise)
    data_frame2 = cluster_centers %>% select(diagnosis_noise, relapse_noise)
    
    # creating a distance matrix
    distmat <- dist(data_frame1, data_frame2, method="euclidean")
    
    
    focal_calls_0$min_dist_clust_center <- apply(distmat, 1, min) 
    
    focal_calls_0_genes <- focal_calls_0$hgnc_symbol
    
    focal_calls_0_dist <- focal_calls_0 %>% filter(min_dist_clust_center> 0.3)
    
    try({
      focal_calls_0_dist$pt_idx <- i
      focal_calls_0_dist$pt_name <- pts[i]
    })
    
    
    focal_calls_0_dist_genes <- focal_calls_0_dist$hgnc_symbol
    
    
    ALL_FOCALS_ALONE_COMPLETE_list[[i]] <-  focal_calls_0_dist
    
    #paired phases plot
    if (plot== T) {
      
      
      if (any(dat_pt$DBscan_cluster !=0)) {
        
        cat(paste0("pt_indx= ",i, " plotted with success\n"), file = paste0(outpath,"ERRORS.log"), append = TRUE)
        
        png(paste0(outpath, pts[i],"_plot.png"), width = 10, height = 9, units = "in", res = 300)
        print(
          
          dat_pt %>% filter(DBscan_cluster != 0) %>%
            ggplot(aes(diagnosis_noise, relapse_noise)) +
            geom_abline(slope = 1, linetype=2)+
            geom_point(data = dat_pt %>% filter(DBscan_cluster == 0 & diagnosis_noise<=4.5 & relapse_noise<= 4.5), colour="black", alpha=0.1, shape=1 )+
            geom_point(data = dat_pt_wt, colour="black", alpha=0.05, shape=1 )+
            geom_point(aes(colour=DBscan_cluster), alpha=0.2) +
            ggforce::geom_mark_hull(aes(group=DBscan_cluster), expand = 0.008) +
            geom_point(data=cluster_centers, colour="red", size=2)+
            geom_vline(xintercept = c(seq(0.2,5,1), seq(0.8,4,1)), linetype=3)+
            geom_hline(yintercept = c(seq(0.2,5,1), seq(0.8,4,1)), linetype=3)+
            geom_rect(aes(xmin = 1.8, xmax = 2.2, ymin = 1.8, ymax = 2.2), fill = "white", alpha=0, color = "black") +
            coord_fixed()+
            ggtitle(paste0(pts[i], " pt idx:",i," Diagnosis: ", sample_D, " Relapse: ", sample_R),
                    subtitle = paste0("noise=",NOISE, "   EPS=",EPS_pair, "   MIN=", MIN_pair))+
            geom_point(data = dat_pt %>% filter(gene_names %in% focal_calls_0_genes), colour="blue", size=2)+
            geom_point(data = dat_pt %>% filter(gene_names %in% focal_calls_0_dist_genes), colour="green", size=2)+
            ggrepel::geom_label_repel(data = dat_pt %>% filter(gene_names %in% focal_calls_0_dist_genes), aes(label=gene_names))+
            xlim(0,4.5)+ylim(0,4.5)+
            ylab("CN Relapse")+ xlab("CN Diagnosis") + theme_light()
          
        )
        
        dev.off()
        
      } else { cat(paste0("pt_indx= ",i, " has no output plot - 0 cluster\n"), file = paste0(outpath,"ERRORS.log"), append = TRUE) }
      
    }
    
  }
  
  
  ALL_FOCALS_ALONE_COMPLETE <- Reduce(rbind, ALL_FOCALS_ALONE_COMPLETE_list)
  ALL_CLUSTERS_CENTERS_COMPLETE <- Reduce(rbind, ALL_CLUSTERS_CENTERS_COMPLETE_list)
  
  output[["focals_not_in_clusters"]] <- ALL_FOCALS_ALONE_COMPLETE
  output[["cluster_centers"]] <- ALL_CLUSTERS_CENTERS_COMPLETE
  
  output
}