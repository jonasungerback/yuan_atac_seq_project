#PCA trine

FL_vs_adult_for_PCA <- read.delim("180220_FL_vs_adult_tc_for_PCA_std_8_filtered.txt", stringsAsFactors = F)
FL_vs_adult_for_PCA_unfiltered <- read.delim("180220_FL_vs_adult_tc_for_PCA.txt",stringsAsFactors = F)

for_pca<-FL_vs_adult_for_PCA_unfiltered[,c(20:27)]
for_pca<-for_pca+0.01
for_pca<-log(for_pca)

#Individually normalize each experiment
for_pca_rep1 <-for_pca[,c(1,3,5,7)]
for_pca_rep2 <-for_pca[,c(2,4,6,8)]

gene.means_1 <- apply(for_pca_rep1,1,mean)
gene.means_2 <- apply(for_pca_rep2,1,mean)

dfl.base.c1 <- sweep(for_pca_rep1, 1, gene.means_1, "-")
dfl.base.c2 <- sweep(for_pca_rep2, 1, gene.means_2, "-")

#Merge experiment
dfl.base.c <- cbind(dfl.base.c1,dfl.base.c2)

#PCA
pca.base.c <- prcomp(t(dfl.base.c), scale=FALSE, center=FALSE, retx=TRUE, cor=TRUE)
print(summary(pca.base.c))
eigenvectors <- as.data.frame(pca.base.c$rotation)
loadings.base <- t(as.matrix(dfl.base.c)) %*% as.matrix(eigenvectors)
loadings.base <- as.data.frame(loadings.base)

#PCA1: 62.8% of the variation
#PCA2: 24.1% of the cariation.


#Plot

B<-matrix(c("BM_HSC_1", "BM_LMPP_1", "FL_HSC_1", "FL_LMPP_1", "BM_HSC_2", "BM_LMPP_2", "FL_HSC_2", "FL_LMPP_2"),nrow=8,ncol=1)
loadings.base<-cbind(loadings.base,B)
library(ggplot2)
ggplot(data = loadings.base)+geom_point(aes(PC1,PC2, colour=factor(B)), size=4)+geom_text(aes(PC1,PC2), label=rownames(loadings.base))+theme(legend.position = "bottom",legend.key=element_rect(colour="NA", fill = "NA"), legend.title=element_blank())+scale_x_continuous(limits = c(-220,160),breaks =  seq(-220,160, by = 20)) + scale_y_continuous(limits = c(-160,120),breaks =  seq(-160,120, by = 20))+theme_bw()


ggsave(filename = "FL_vs_adult_for_PCA.pdf",dpi = 600, units = c("cm"), height = 10, width = 15)



#Create a function of the code above

atac_pca_trine <- function(for_pca, sd_cutoff) {
  
  for_pca$standard_d <- apply(for_pca,1,sd)
  for_pca <- subset(for_pca, for_pca$standard_d >= sd_cutoff)
  for_pca <- for_pca[,-c(max(ncol(for_pca)))]
  
  for_pca<-for_pca+0.01
  for_pca<-log(for_pca)
  
  #Individually normalize each experiment
  for_pca_rep1 <-for_pca[,c(1,3,5,7)]
  for_pca_rep2 <-for_pca[,c(2,4,6,8)]
  
  gene.means_1 <- apply(for_pca_rep1,1,mean)
  gene.means_2 <- apply(for_pca_rep2,1,mean)
  
  dfl.base.c1 <- sweep(for_pca_rep1, 1, gene.means_1, "-")
  dfl.base.c2 <- sweep(for_pca_rep2, 1, gene.means_2, "-")
  
  #Merge experiment
  dfl.base.c <- cbind(dfl.base.c1,dfl.base.c2)
  
  #PCA
  pca.base.c <- prcomp(t(dfl.base.c), scale=FALSE, center=FALSE, retx=TRUE, cor=TRUE)
  pca_var <- summary(pca.base.c)
  print(summary(pca.base.c))
  eigenvectors <- as.data.frame(pca.base.c$rotation)
  loadings.base <- t(as.matrix(dfl.base.c)) %*% as.matrix(eigenvectors)
  loadings.base <- as.data.frame(loadings.base)
  
  #PCA1: 62.8% of the variation
  #PCA2: 24.1% of the cariation.
  
  
  #Plot
  
  B<-matrix(c("BM_HSC_1", "BM_LMPP_1", "FL_HSC_1", "FL_LMPP_1", "BM_HSC_2", "BM_LMPP_2", "FL_HSC_2", "FL_LMPP_2"),nrow=8,ncol=1)
  loadings.base<-cbind(loadings.base,B)
  library(ggplot2)
  gg_pca_plot <- ggplot(data = loadings.base)+geom_point(aes(PC1,PC2, colour=factor(B)), size=4)+geom_text(aes(PC1,PC2), label=rownames(loadings.base))+theme(legend.position = "bottom",legend.key=element_rect(colour="NA", fill = "NA"), legend.title=element_blank())+scale_x_continuous(limits = c(-220,160),breaks =  seq(-220,160, by = 20)) + scale_y_continuous(limits = c(-160,120),breaks =  seq(-160,120, by = 20))+theme_bw()
  
  plot_and_var <-list(gg_pca_plot = gg_pca_plot, pca_var = pca_var)
  return(plot_and_var)
  
}

base_plot <- atac_pca_trine(for_pca = FL_vs_adult_for_PCA_unfiltered[,c(20:27)], sd_cutoff = 8)
base_plot$gg_pca_plot

# Remove the peaks that are associated with MkP-signature genes and rerun:

source("/Users/Jonas/Documents/Arbetet/Biotoolkits/scripts/R/match_list_function.R")

megpk_genes <- read.delim("/Users/Jonas/Documents/Arbetet/KAW-project_Lund/ATAC-seq/Trine_ATAC_project/Diff_peaks_vs_gene_expression/gsea_analyses/MkP_geneset.grp", header = F)
colnames(megpk_genes)[1]<-"genes"

FL_vs_adult_for_PCA_unfiltered <- match_lists(FL_vs_adult_for_PCA_unfiltered,FL_vs_adult_for_PCA_unfiltered$Gene.Name, megpk_genes$genes,1,name_match_column = "megpk")


FL_vs_adult_for_PCA_unfiltered_no_MpK <- subset(FL_vs_adult_for_PCA_unfiltered, FL_vs_adult_for_PCA_unfiltered$megpk != 1)

no_MpK_peaks_plot <- atac_pca_trine(for_pca = FL_vs_adult_for_PCA_unfiltered_no_MpK[,c(20:27)], sd_cutoff = 0)

no_MpK_peaks_plot$gg_pca_plot

#Didn't do anything to remove those genes.
#What about removing peaks with NFE2, GATA1 and CTCFL motifs?

#In terminal:
#  gimme scan 180220_FL_vs_adult_tc_for_PCA.bed -p ../Motifs_in_diff_peaks/GM.5.0.NFE2.GATA1.CTCFL.txt -g mm10 -f 0.05 -t  > 180220_FL_vs_adult_tc_for_PCA.NFE2.GATA.CTCFL_motifs_f_0.05.txt
#  gimme scan 180220_FL_vs_adult_tc_for_PCA.bed -p ../Motifs_in_diff_peaks/GM.5.0.NFE2.GATA1.CTCFL.txt -g mm10 -f 0.05 -b  > 180220_FL_vs_adult_tc_for_PCA.NFE2.GATA.CTCFL_motifs_f_0.05.bed

  
  
FL_vs_adult_for_PCA_unfiltered_motifs <- read.delim("180220_FL_vs_adult_tc_for_PCA.NFE2.GATA.CTCFL_motifs_f_0.05.txt", stringsAsFactors = F)

FL_vs_adult_for_PCA_no_NFE2_motif <- subset(FL_vs_adult_for_PCA_unfiltered_motifs,FL_vs_adult_for_PCA_unfiltered_motifs$GM.5.0.bZIP.NFE2.0005 != 1)
FL_vs_adult_for_PCA_no_GATA1_motif <- subset(FL_vs_adult_for_PCA_unfiltered_motifs,FL_vs_adult_for_PCA_unfiltered_motifs$GM.5.0.GATA.0001 != 1)
FL_vs_adult_for_PCA_no_CTCFL_motif <- subset(FL_vs_adult_for_PCA_unfiltered_motifs,FL_vs_adult_for_PCA_unfiltered_motifs$GM.5.0.C2H2_ZF.0023 != 1)

FL_vs_adult_for_PCA_no_NFE2_motif <- as.data.frame( FL_vs_adult_for_PCA_no_NFE2_motif[,2] )
colnames(FL_vs_adult_for_PCA_no_NFE2_motif)[1]<-"PeakID"
FL_vs_adult_for_PCA_no_NFE2_motif <- merge(FL_vs_adult_for_PCA_no_NFE2_motif,FL_vs_adult_for_PCA_unfiltered, by = "PeakID", all.x = T)

FL_vs_adult_for_PCA_no_GATA1_motif <- as.data.frame( FL_vs_adult_for_PCA_no_GATA1_motif[,2] )
colnames(FL_vs_adult_for_PCA_no_GATA1_motif)[1]<-"PeakID"
FL_vs_adult_for_PCA_no_GATA1_motif <- merge(FL_vs_adult_for_PCA_no_GATA1_motif,FL_vs_adult_for_PCA_unfiltered, by = "PeakID", all.x = T)

FL_vs_adult_for_PCA_no_CTCFL_motif <- as.data.frame( FL_vs_adult_for_PCA_no_CTCFL_motif[,2] )
colnames(FL_vs_adult_for_PCA_no_CTCFL_motif)[1]<-"PeakID"
FL_vs_adult_for_PCA_no_CTCFL_motif <- merge(FL_vs_adult_for_PCA_no_CTCFL_motif,FL_vs_adult_for_PCA_unfiltered, by = "PeakID", all.x = T)  

base_plot <- atac_pca_trine(for_pca = FL_vs_adult_for_PCA_unfiltered[,c(20:27)], sd_cutoff = 0)
x11()
base_plot$gg_pca_plot


#Run PCA on this
no_NFE2_peaks_plot <- atac_pca_trine(for_pca = FL_vs_adult_for_PCA_no_NFE2_motif[,c(20:27)], sd_cutoff = 0)
x11()
no_NFE2_peaks_plot$gg_pca_plot

no_GATA1_peaks_plot <- atac_pca_trine(for_pca = FL_vs_adult_for_PCA_no_GATA1_motif[,c(20:27)], sd_cutoff = 0)
x11()
no_GATA1_peaks_plot$gg_pca_plot

no_CTCFL_peaks_plot <- atac_pca_trine(for_pca = FL_vs_adult_for_PCA_no_CTCFL_motif[,c(20:27)], sd_cutoff = 0)
x11()
no_CTCFL_peaks_plot$gg_pca_plot

ggsave(plot = base_plot$gg_pca_plot, filename = "ATAC_pca_base_all_peaks.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)
ggsave(plot = no_MpK_peaks_plot$gg_pca_plot, filename = "ATAC_pca_no_MpK_peaks.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)
ggsave(plot = no_NFE2_peaks_plot$gg_pca_plot, filename = "ATAC_pca_no_NFE2_peaks.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)
ggsave(plot = no_GATA1_peaks_plot$gg_pca_plot, filename = "ATAC_pca_no_GATA1_peaks.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)
ggsave(plot = no_CTCFL_peaks_plot$gg_pca_plot, filename = "ATAC_pca_no_CTCFL_peaks.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)


#Remove FLvsBM_HSC_diff peaks
FLvsABM_HSC_closing <- read.delim("/Users/Jonas/Documents/Arbetet/KAW-project_Lund/ATAC-seq/Trine_ATAC_project/FL_vs_adult.Down_FL_HSC_vs_BM_HSC.txt", header = T)
FLvsABM_HSC_opening <- read.delim("/Users/Jonas/Documents/Arbetet/KAW-project_Lund/ATAC-seq/Trine_ATAC_project/FL_vs_adult.Up_FL_HSC_vs_BM_HSC.txt", header = T)
colnames(FLvsABM_HSC_opening)[1] <-"PeakID"
colnames(FLvsABM_HSC_closing)[1] <-"PeakID"

FLvsABM_HSC_opening_closing <- rbind(FLvsABM_HSC_opening,FLvsABM_HSC_closing)
FL_vs_adult_for_PCA_unfiltered <- match_lists(FL_vs_adult_for_PCA_unfiltered,FL_vs_adult_for_PCA_unfiltered$PeakID, FLvsABM_HSC_opening_closing$PeakID,1,name_match_column = "BMvsFL_HSC_diff_peaks")

FL_vs_adult_for_PCA_unfiltered_no_HSC_diff <- subset(FL_vs_adult_for_PCA_unfiltered,FL_vs_adult_for_PCA_unfiltered$BMvsFL_HSC_diff_peaks != 1)
no_HSC_diff_peaks_plot <- atac_pca_trine(for_pca = FL_vs_adult_for_PCA_unfiltered_no_HSC_diff[,c(20:27)], sd_cutoff = 0)
x11()
no_HSC_diff_peaks_plot$gg_pca_plot
ggsave(plot = no_HSC_diff_peaks_plot$gg_pca_plot, filename = "ATAC_pca_no_HSC_diff_peaks.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)

#boxplot standard deviation for MpK-cats
library(ggplot2)
ggplot(data = FL_vs_adult_for_PCA_unfiltered)+geom_boxplot(aes(y = log10(Standar_dev), x = factor(megpk),fill=factor(megpk)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme_bw()
ggplot(data = FL_vs_adult_for_PCA_unfiltered)+geom_boxplot(aes(y = log10(Standar_dev), x = factor(BMvsFL_HSC_diff_peaks),fill=factor(BMvsFL_HSC_diff_peaks)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme_bw()
log10(mean(FL_vs_adult_for_PCA_unfiltered$Standar_dev[FL_vs_adult_for_PCA_unfiltered$megpk == 1]))

FL_vs_adult_for_PCA_unfiltered <- match_lists(FL_vs_adult_for_PCA_unfiltered,FL_vs_adult_for_PCA_unfiltered$PeakID, FL_vs_adult_for_PCA_no_NFE2_motif$PeakID,1,name_match_column = "NFE2")
ggplot(data = FL_vs_adult_for_PCA_unfiltered)+geom_boxplot(aes(y = log10(Standar_dev), x = factor(NFE2),fill=factor(NFE2)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme_bw()

FL_vs_adult_for_PCA_unfiltered <- match_lists(FL_vs_adult_for_PCA_unfiltered,FL_vs_adult_for_PCA_unfiltered$PeakID, FL_vs_adult_for_PCA_no_GATA1_motif$PeakID,1,name_match_column = "GATA1")
ggplot(data = FL_vs_adult_for_PCA_unfiltered)+geom_boxplot(aes(y = log10(Standar_dev), x = factor(GATA1),fill=factor(GATA1)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme_bw()

FL_vs_adult_for_PCA_unfiltered <- match_lists(FL_vs_adult_for_PCA_unfiltered,FL_vs_adult_for_PCA_unfiltered$PeakID, FL_vs_adult_for_PCA_no_CTCFL_motif$PeakID,1,name_match_column = "CTCFL")
ggplot(data = FL_vs_adult_for_PCA_unfiltered)+geom_boxplot(aes(y = log10(Standar_dev), x = factor(CTCFL),fill=factor(CTCFL)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme_bw()



