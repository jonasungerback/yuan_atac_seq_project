
##################### Creates peak files from Clust clusters
# FLvsABM_HSC

clust_FLvsABm_atac <- read.delim("Clusters_Objects.tsv", header = T)
clust_FLvsABm_atac <- clust_FLvsABm_atac[c(2:max(nrow(clust_FLvsABm_atac))),]
#Read in peak file

FLvsABM_closing <- read.delim("/Users/Jonas/Documents/Arbetet/KAW-project_Lund/ATAC-seq/Trine_ATAC_project/FL_vs_adult.Down_FL_LMPP_vs_BM_LMPP.txt", header = T)
FLvsABM_opening <- read.delim("/Users/Jonas/Documents/Arbetet/KAW-project_Lund/ATAC-seq/Trine_ATAC_project/FL_vs_adult.Up_FL_LMPP_vs_BM_LMPP.txt", header = T)
colnames(FLvsABM_opening)[1]<-"PeakID"
colnames(FLvsABM_closing)[1]<-"PeakID"
FLvsABM_opening_closing <- rbind(FLvsABM_closing,FLvsABM_opening)

C0 <- as.data.frame(clust_FLvsABm_atac[1:999, 1])
colnames(C0)[1] <- "PeakID"
C0$cluster <- "C0" 
C0 <- merge(FLvsABM_opening_closing,C0, by = "PeakID", all.y = T)

C1 <- cbind( as.data.frame(clust_FLvsABm_atac[1:453, 2]))
colnames(C1)[1] <- "PeakID"
C1$cluster <- "C1"
C1 <- merge(FLvsABM_opening_closing,C1, by = "PeakID", all.y = T)


C2 <- cbind( as.data.frame(clust_FLvsABm_atac[1:64, 3]))
colnames(C2)[1] <- "PeakID"
C2$cluster <- "C2"
C2 <- merge(FLvsABM_opening_closing,C2, by = "PeakID", all.y = T)

C3 <- cbind( as.data.frame(clust_FLvsABm_atac[1:157, 4]))
colnames(C3)[1] <- "PeakID"
C3$cluster <- "C3"
C3 <- merge(FLvsABM_opening_closing,C3, by = "PeakID", all.y = T)

C4 <- as.data.frame(clust_FLvsABm_atac[1:123, 5])
colnames(C4)[1] <- "PeakID"
C4$cluster <- "C4"
C4 <- merge(FLvsABM_opening_closing,C4, by = "PeakID", all.y = T)


C5 <- as.data.frame(clust_FLvsABm_atac[1:130, 6])
colnames(C5)[1] <- "PeakID"
C5$cluster <- "C5"
C5 <- merge(FLvsABM_opening_closing,C5, by = "PeakID", all.y = T)


C6 <- as.data.frame(clust_FLvsABm_atac[1:201, 7])
colnames(C6)[1] <- "PeakID"
C6$cluster <- "C6"
C6 <- merge(FLvsABM_opening_closing,C6, by = "PeakID", all.y = T)

C7 <- as.data.frame(clust_FLvsABm_atac[1:113, 8])
colnames(C7)[1] <- "PeakID"
C7$cluster <- "C7"
C7 <- merge(FLvsABM_opening_closing,C7, by = "PeakID", all.y = T)


C8 <- as.data.frame(clust_FLvsABm_atac[1:715, 9])
colnames(C8)[1] <- "PeakID"
C8$cluster <- "C8"
C8 <- merge(FLvsABM_opening_closing,C8, by = "PeakID", all.y = T)


C9 <- as.data.frame(clust_FLvsABm_atac[1:134, 10])
colnames(C9)[1] <- "PeakID"
C9$cluster <- "C9"
C9 <- merge(FLvsABM_opening_closing,C9, by = "PeakID", all.y = T)

C10 <- as.data.frame(clust_FLvsABm_atac[1:608, 11])
colnames(C10)[1] <- "PeakID"
C10$cluster <- "C10"
C10 <- merge(FLvsABM_opening_closing,C10, by = "PeakID", all.y = T)


C11 <- as.data.frame(clust_FLvsABm_atac[1:694, 12])
colnames(C11)[1] <- "PeakID"
C11$cluster <- "C11"
C11 <- merge(FLvsABM_opening_closing,C11, by = "PeakID", all.y = T)

C12 <- as.data.frame(clust_FLvsABm_atac[1:99, 13])
colnames(C12)[1] <- "PeakID"
C12$cluster <- "C12"
C12 <- merge(FLvsABM_opening_closing,C12, by = "PeakID", all.y = T)


C13 <- as.data.frame(clust_FLvsABm_atac[1:110, 14])
colnames(C13)[1] <- "PeakID"
C13$cluster <- "C13"
C13 <- merge(FLvsABM_opening_closing,C13, by = "PeakID", all.y = T)

C14 <- as.data.frame(clust_FLvsABm_atac[1:41, 15])
colnames(C14)[1] <- "PeakID"
C14$cluster <- "C14"
C14 <- merge(FLvsABM_opening_closing,C14, by = "PeakID", all.y = T)

C15 <- as.data.frame(clust_FLvsABm_atac[1:202, 16])
colnames(C15)[1] <- "PeakID"
C15$cluster <- "C15"
C15 <- merge(FLvsABM_opening_closing,C15, by = "PeakID", all.y = T)

C16 <- as.data.frame(clust_FLvsABm_atac[1:125, 17])
colnames(C16)[1] <- "PeakID"
C16$cluster <- "C16"
C16 <- merge(FLvsABM_opening_closing,C16, by = "PeakID", all.y = T)


C17 <- as.data.frame(clust_FLvsABm_atac[1:189, 18])
colnames(C17)[1] <- "PeakID"
C17$cluster <- "C17"
C17 <- merge(FLvsABM_opening_closing,C17, by = "PeakID", all.y = T)

C18 <- as.data.frame(clust_FLvsABm_atac[1:176, 19])
colnames(C18)[1] <- "PeakID"
C18$cluster <- "C18"
C18 <- merge(FLvsABM_opening_closing,C18, by = "PeakID", all.y = T)


#Create one final columns
clust_grps_FLvsABM_opening_closing <-rbind(C0,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18)
write.table(clust_grps_FLvsABM_opening_closing, "clust_grps_FLvsABM_LMPP_opening_closing_t_5_n_101_4.txt", sep = "\t", row.names = F)

source("/Users/Jonas/Documents/Arbetet/Biotoolkits/scripts/R/match_list_function.R")

megpk_genes <- read.delim("/Users/Jonas/Documents/Arbetet/KAW-project_Lund/ATAC-seq/Trine_ATAC_project/Diff_peaks_vs_gene_expression/gsea_analyses/MkP_geneset.grp", header = F)
colnames(megpk_genes)[1]<-"genes"

clust_grps_FLvsABM_opening_closing <- match_lists(clust_grps_FLvsABM_opening_closing,clust_grps_FLvsABM_opening_closing$Gene.Name, megpk_genes$genes,1,name_match_column = "megpk")

clust_grps_FLvsABM_opening_closing_cluster <- clust_grps_FLvsABM_opening_closing[,31:32]

#Function to calculate cat genes per cluster
collapse_clusters <- function(df, clust_col, inp){
  
  df_counts <- aggregate(inp ~ clust_col, data = df, FUN=sum)
  df_lengths <- aggregate(inp ~ clust_col, data = df, FUN=length)
  df_new <- cbind(df_counts, df_lengths[,2])
  df_new$freq <- df_new[,2]/df_new[,3]
  return(df_new)
  
}


test <- collapse_clusters(clust_grps_FLvsABM_opening_closing_cluster,clust_col = clust_grps_FLvsABM_opening_closing_cluster[,1],inp = clust_grps_FLvsABM_opening_closing_cluster[,2])


library(ggplot2)
ggplot(test, aes(x=clust_col, y=freq, fill = clust_col)) + geom_bar(stat="identity")

ggplot(test, aes(x=clust_col, y=inp, fill = clust_col)) + geom_bar(stat="identity")


