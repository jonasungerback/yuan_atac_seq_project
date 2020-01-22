#General code for histogram plot from a Homer file with only the Coverage columns left.
#Read in thed data.

#Plotting Coverage over NFE2 sites: (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42110)
setwd("~/Desktop/histogram_plots_ATAC_coverage_NFE2_peaks")
MKD5_unique_peaks_ATAC_coverage <-read.delim("NFE2_MKD5.RBC2D.MKD5_unique.for_hist.txt", stringsAsFactors = F)
colnames(MKD5_unique_peaks_ATAC_coverage)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")

#Colors
#BM_HSC: blue
#BM_LMPP: purple
#FL_HSC: red
#FL_LMPP: gold

library(ggplot2)
ggplot(data = MKD5_unique_peaks_ATAC_coverage)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(8,5.6,5.2,4.8), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "MKD5_unique_peaks_ATAC_coverage.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)


RBC2D_unique_peaks_ATAC_coverage <-read.delim("NFE2_MKD5.RBC2D.RBC2D_unique.for_hist.txt", stringsAsFactors = F)
colnames(RBC2D_unique_peaks_ATAC_coverage)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = RBC2D_unique_peaks_ATAC_coverage)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(4,3.2,2.8,2.4), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "RBC2D_unique_peaks_ATAC_coverage.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)

MKD5_RBC2D_overlapping_peaks_ATAC_coverage <-read.delim("NFE2_MKD5.RBC2D.overlapping.for_hist.txt", stringsAsFactors = F)
colnames(MKD5_RBC2D_overlapping_peaks_ATAC_coverage)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = MKD5_RBC2D_overlapping_peaks_ATAC_coverage)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(4,3.2,2.8,2.4), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "MKD5_RBC2D_overlapping_peaks_ATAC_coverage.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)

#Size 5000
RBC2D_unique_peaks_ATAC_coverage_size_5000 <-read.delim("NFE2_MKD5.RBC2D.RBC2D_unique.for_hist_size_5000.txt", stringsAsFactors = F)
colnames(RBC2D_unique_peaks_ATAC_coverage_size_5000)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = RBC2D_unique_peaks_ATAC_coverage_size_5000)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(4,3.2,2.8,2.4), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "RBC2D_unique_peaks_ATAC_coverage_size_5000.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)

MKD5_unique_peaks_ATAC_coverage_size_5000 <-read.delim("NFE2_MKD5.RBC2D.MKD5_unique.for_hist_size_5000.txt", stringsAsFactors = F)
colnames(MKD5_unique_peaks_ATAC_coverage_size_5000)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = MKD5_unique_peaks_ATAC_coverage_size_5000)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(8,5.6,5.2,4.8), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "MKD5_unique_peaks_ATAC_coverage_size_5000.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)


#Plotting ATAC-coverage over GATA1 peaks (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51337)

Gata1_meg_unique_peaks_ATAC_coverage <-read.delim("Gata1_GSM923586_PSU_ChipSeq_Megakaryo_GATA1_peaks.for_hist.txt", stringsAsFactors = F)
colnames(Gata1_meg_unique_peaks_ATAC_coverage)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = Gata1_meg_unique_peaks_ATAC_coverage)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(4,3.2,2.8,2.4), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "GSE51337_Gata1_MegEry.mega_unique_ATAC_coverage.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)


Gata1_ery_unique_peaks_ATAC_coverage <-read.delim("Gata1_GSM923575_PSU_ChipSeq_Erythrobl_GATA1_peaks.for_hist.txt", stringsAsFactors = F)
colnames(Gata1_ery_unique_peaks_ATAC_coverage)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = Gata1_ery_unique_peaks_ATAC_coverage)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(2.5,2,1.5,1), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "GSE51337_Gata1_MegEry.ery_unique_ATAC_coverage.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)

Gata1_ery_mega_shared_peaks_ATAC_coverage <-read.delim("GSE51337_Gata1_MegEry.shared.for_hist.txt", stringsAsFactors = F)
colnames(Gata1_ery_mega_shared_peaks_ATAC_coverage)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = Gata1_ery_mega_shared_peaks_ATAC_coverage)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(2.5,2,1.5,1), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "Gata1_ery_mega_shared_peaks_ATAC_coverage.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)

Gata1_ery_unique_peaks_ATAC_coverage_10000 <-read.delim("GSE51337_Gata1_MegEry.ery_unique.for_hist_10000.txt", stringsAsFactors = F)
colnames(Gata1_ery_unique_peaks_ATAC_coverage_10000)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = Gata1_ery_unique_peaks_ATAC_coverage_10000)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(2.5,2,1.5,1), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "GSE51337_Gata1_MegEry.ery_unique_ATAC_coverage.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)

Gata1_meg_unique_peaks_ATAC_coverage_1_bin <-read.delim("GSE51337_Gata1_MegEry.mega_unique.for_hist_1_bin.txt", stringsAsFactors = F)
colnames(Gata1_meg_unique_peaks_ATAC_coverage_1_bin)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data =Gata1_meg_unique_peaks_ATAC_coverage_1_bin)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(4,3.2,2.8,2.4), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "GSE51337_Gata1_MegEry.mega_unique_ATAC_coverage.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)


Gata1_ery_unique_peaks_ATAC_coverage_1_bin <-read.delim("GSE51337_Gata1_MegEry.ery_unique.for_hist_1_bin.txt", stringsAsFactors = F)
colnames(Gata1_ery_unique_peaks_ATAC_coverage_1_bin)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = Gata1_ery_unique_peaks_ATAC_coverage_1_bin)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(2.5,2,1.5,1), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "GSE51337_Gata1_MegEry.ery_unique_ATAC_coverage.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)

#Size 5000

Gata1_meg_unique_peaks_ATAC_coverage_size_5000 <-read.delim("Gata1_GSM923586_PSU_ChipSeq_Megakaryo_GATA1_peaks.5000_for_hist.txt", stringsAsFactors = F)
colnames(Gata1_meg_unique_peaks_ATAC_coverage_size_5000)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = Gata1_meg_unique_peaks_ATAC_coverage_size_5000)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(4,3.2,2.8,2.4), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "Gata1_meg_unique_peaks_ATAC_coverage_size_5000.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)


Gata1_ery_unique_peaks_ATAC_coverage_size_5000 <-read.delim("Gata1_GSM923575_PSU_ChipSeq_Erythrobl_GATA1_peaks.5000_for_hist.txt", stringsAsFactors = F)
colnames(Gata1_ery_unique_peaks_ATAC_coverage_size_5000)<-c("Distance_from_center","BM_HSC","BM_HSC_tag+","BM_HSC_tag-","BM_LMPP","BM_LMPP_tag+","BM_LMPP_tag-", "FL_HSC","FL_HSC_tag+","FL_HSC_tag-","FL_LMPP","FL_LMPP_tag+","FL_LMPP_tag-")
ggplot(data = Gata1_ery_unique_peaks_ATAC_coverage_size_5000)+geom_line(aes(Distance_from_center,BM_HSC, group = 1), color = "blue", size = 1)+geom_line(aes(Distance_from_center,BM_LMPP, group = 1), color = "purple", size = 1)+geom_line(aes(Distance_from_center,FL_HSC, group = 1), color = "red",size = 1)+geom_line(aes(Distance_from_center,FL_LMPP, group = 1), color = "gold", size = 1)+theme_bw()+annotate("text", x = c(500,500,500,500) , y = c(2.5,2,1.5,1), label = c("BM HSC","BM LMPP","FL HSC","FL LMPP"), size = 5, color = c("blue","purple","red","gold"))+ylab("Coverage")+xlab("Distance from peak center")
ggsave(filename = "Gata1_ery_unique_peaks_ATAC_coverage_size_5000.pdf", dpi = 300, width = 8, height = 6, useDingbats=FALSE)
