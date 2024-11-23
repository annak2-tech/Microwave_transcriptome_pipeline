master=read.table("C:/...Mock04h_vs_Treated04h_allgenes_deseq2_final_CPM.csv", header=TRUE, row.names=1, sep=",")

library("dplyr")
master = tibble::rownames_to_column(master, "rn")
#Na omit
master=na.omit(master)
row.names(master)= master[,16] 
#renaming ensembl gene ID header
names(master)[1]="ID"

sorted_order_master=order(master[,7], decreasing=FALSE)
master=master[sorted_order_master,]
#Add a column flagging significance to the master table
master$sig=as.factor(master$padj<0.05 & abs(master$log2FoldChange) > 1.0)
#removing symbol column as now titles
master=master[,-16]

#Adding a column of -log10p
log_10_p_4H=-log10(master$pvalue)
master$mlog10p=log_10_p_4H

EM_4H = master[8:15]
#make a scaled expression matrix
#scale function works by columns not row so need to transpose before scaling 
# also include a second transpose to flip back
#transpose returns a matrix not a data frame so cast back to data frame
EM.s_4H=na.omit(data.frame(t(scale(t(EM_4H)))))

#Make a new column in master table for mean expression
master$means=rowMeans(master[,8:15])

log_10_means_4H=log10(master$means)
master$log10means=log_10_means_4H
subset(master, padj<0.05)
master_sig_4H=subset(master, padj<0.05 & abs(log2FoldChange)>1)
# producing a list of names of the genes significantly up or down regulated
sig_genes_4H=row.names(master_sig_4H)

#Making an expression table of significant genes only
em_symbols_sig_4H=EM_4H[sig_genes_4H,] 

em_scaled_sig_4H=EM.s_4H[sig_genes_4H,]

master_sig_up_4H = subset(master_sig_4H, log2FoldChange>0)
master_sig_down_4H = subset(master_sig_4H, log2FoldChange<0)

master_sig_up_top5_4H = master_sig_up_4H[1:5,]
master_sig_down_top5_4H = master_sig_down_4H[1:5,]



#Saving tables 
write.table(master, file="C:/...master_table4h.csv", sep=",")
write.table(EM_4H, file="C:/...EM_4H.csv", sep=",")
write.table(EM.s_4H, file="C:/...EM.s_4H.csv", sep=",")
write.table(master_sig_4H, file="C:/...master_sig_4h.csv", sep=",")
write.table(master_sig_up_4H, file="C:/...master_sig_up_4h.csv", sep=",")
write.table(master_sig_down_4H, file="C:/...master_sig_down_4h.csv", sep=",")


#loading libraries
library(ggplot2)
library(ggrepel)
library(reshape2)
library(amap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(devEMF)
library(dplyr)

#Creating a custom theme
my_theme= theme(
  panel.grid = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=1),
  plot.background = element_blank(), 
  legend.background = element_rect(fill="transparent", colour=NA),
  legend.key = element_rect(fill="transparent", colour=NA),
  plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, family="sans", face="bold"),
  title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, family="sans", face="bold"),
  axis.text.y = element_text(size = 11, margin = margin(r = 5),hjust=1,vjust=0.5, family="sans", face="bold",colour="black"),
  axis.text.x = element_text(size = 11, margin = margin(t = 5),hjust=0.5,vjust=1, family="sans", face="bold",colour="black"), 
  axis.title.y = element_text(size = 12, margin = margin(r = 10),angle = 90,hjust=0.5,vjust=0.5, family="sans", face="bold"),
  axis.title.x = element_text(size = 12, margin = margin(t = 10),hjust=0.5,vjust=1, family="sans", face="bold"),
  legend.text=element_text(size=12, family="sans", face="bold"),
  legend.title=element_blank(), 
  legend.key.size=unit(1,"line"),
  plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"),
  strip.text.x = element_text(size = 12, family="sans", face="bold", vjust=1),
  panel.spacing = unit(1, "lines")
)

#Making an MA plot
ggp = ggplot(master, aes(x=log10means, y=log2FoldChange))+
  geom_point(colour="black")+
  geom_point(data=master_sig_up_4H, colour = "red")+
  geom_point(data=master_sig_down_4H, colour="blue")+
  labs(title="MA plot 4 hours", x="Log10(Mean Expression)", y="Log2 Fold Change")+
  theme_bw()+
  geom_hline(yintercept=-1, linetype="dashed", color="grey", linewidth=0.5)+
  geom_hline(yintercept=1,linetype="dashed", color="grey", linewidth=0.5)+
  geom_text_repel(data=master_sig_up_top5_4H, aes(label=row.names(master_sig_up_top5_4H)), colour="red")+
  geom_text_repel(data=master_sig_down_top5_4H, aes(label=row.names(master_sig_down_top5_4H)), colour="blue")
ggp

png("C:/...MA plot 4H.png", height=400, width=600)
print(ggp)
dev.off()

#for direction need to set this in master
levels(master$sig)
master$sig = factor(master$sig, levels=c("TRUE", "FALSE"))

#since sig only groups by yes or no, need to add an additional column for up and down
# If we pre-set these groups then when we add to ggplot function with colour set to direction, 
#will automatically be coloured based on which group they are in
#If just wanting directions of no significance, up and down shown
master_no_sig = subset(master, sig==FALSE)
master_no_sig$direction = "a"
master_sig_up_4H$direction = "b"
master_sig_down_4H$direction = "c"

master = rbind(master_no_sig, master_sig_up_4H, master_sig_down_4H)

#To include another direction that shows genes where log2Foldchange <-1 or >1 but no sig use the following code
master_no_sig_2 = subset(master, sig==FALSE & log2FoldChange >= -1 & log2FoldChange <= 1)
master_no_sig_3 = subset(master, sig==FALSE & log2FoldChange <= -1)
master_no_sig_4 = subset(master, sig==FALSE & log2FoldChange >=1)
master_no_sig_3=master_no_sig_3[,-22]
master_no_sig_4=master_no_sig_4[,-22]
master_no_sig_2$direction = "a"
master_no_sig_3$direction = "d"
master_no_sig_4$direction = "d"
master_sig_up_4H$direction = "b"
master_sig_down_4H$direction = "c"

master_2 = rbind(master_no_sig_2, master_no_sig_3, master_no_sig_4, master_sig_up_4H, master_sig_down_4H)
master = master_2


#Using my_theme on volcano plot- note no brackets when loading my_theme
#Whereas when using theme_bw() brackets are needed


ggp = ggplot(master, aes(x=log2FoldChange, y=mlog10p, colour=direction)) +geom_point() +
  scale_colour_manual(values=c("black", "red", "blue", "grey"), labels = c("NS", "Significantly up", "Significantly down", "Log2 Fold Change"))+
  labs(title="Volcano plot 4H", x="Log2 Fold Change", y="-Log10P")+
  my_theme+
  geom_vline(xintercept=-1, linetype="dashed", color="grey", linewidth=0.5)+
  geom_vline(xintercept=1,linetype="dashed", color="grey", linewidth=0.5)+
  geom_hline(yintercept=3.75,linetype="dashed", color="grey", linewidth=0.5)+
  xlim(c(-5, 10))+
  ylim(c(0,40))+
  geom_text_repel(data=master_sig_up_top5_4H, aes(label=row.names(master_sig_up_top5_4H), colour="b"), show.legend = FALSE)+
  geom_text_repel(data=master_sig_down_top5_4H, aes(label=row.names(master_sig_down_top5_4H), colour="c"), show.legend = FALSE)

ggp

png("C:/...Volcano plot 4H_shortery.png", height=400, width=600)
print(ggp)
dev.off()

#Heatmaps and clustering
#For clustering need to load the library amap
#install.packages("amap")
library(amap)


# producing a list of names of the genes significantly up or down regulated
sig_genes_4H=row.names(master_sig_4H)
sig_genes_4H_up = row.names(master_sig_up_4H)

#Making an expression table of significant genes only
em_symbols_sig_4H=em_symbols_sig_4H[sig_genes_4H,]

em_scaled_sig_4H=EM.s_4H[sig_genes_4H,]

em_symbols_sig_4H_up = em_symbols_sig_4H[sig_genes_4H_up,]
em_scaled_sig_4H_up=EM.s_4H[sig_genes_4H_up,]


#Running pheatmap
#install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(pheatmap)
pheatmap(em_scaled_sig_4H, cluster_cols=F, show_rownames=F, clustering_method="complete", color = inferno(10))

#To annotate pheatmap
colnames(em_scaled_sig_4H)
em_scaled_colnames=em_scaled_sig_4H
ann_col = data.frame(Condition=c("Mock", "Mock", "Mock", "Mock", "4H TR", "4H TR", "4H TR", "4H TR"))
rownames(ann_col) = colnames(em_scaled_colnames)
ann_row=NA
#pheatmap with annotations of sample type
#pheatmap codes, clustering, row_names
#colour=inferno(10)
#cluster_cols=F
#clustering_method="complete"
#show_colnames=F
#annotation_col=ann_col
#annotation_row=ann_row
pheatmap(em_scaled_colnames, show_colnames = F, annotation_col=ann_col, annotation_row=ann_row, color = inferno(10))
#pheatmap without annotations of sample type
pheatmap(em_scaled_sig_4H, show_colnames=F, color=inferno(10))

#pheatmap without gene names included
pheatmap(em_scaled_sig_4H, cluster_cols=F, clustering_method="complete", color=inferno(10), show_colnames=F, show_rownames=F, annotation_col=ann_col, annotation_row=ann_row)

#Making box plots of top genes in interesting pathways from metascape GSEA
#Loading results from metascape
mscape=read.table("C:/...mscape_4h.csv", header=TRUE, sep=",")
ss = read.table("C:/...ss.csv", header = TRUE, sep=",")
#Get a list of the significant genes in the most enriched gene set

enriched_gene_set_1=as.character(mscape[2,5])
candidate_genes=unlist(strsplit(enriched_gene_set_1, ","))

candidate_genes = c(candidate_genes)
gene_data = EM.s_4H[candidate_genes,]
gene_data = data.frame(t(gene_data))
gene_data$sample_group = ss$SAMPLE_GROUP
levels(gene_data$sample_group)
gene_data$sample_group= factor(gene_data$sample_group, levels=c("Mock", "TR_4H"))
#melt this to form a table with fewer columns
gene_data.m=melt(gene_data, id.vars="sample_group")

ggp=ggplot(gene_data.m, aes(x=variable, y=value, fill=sample_group)) + geom_boxplot() +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(title="Protein folding", x="Gene name", y="Scaled expression")
ggp

em_scaled_pf=EM.s_4H[candidate_genes,]
pheatmap(em_scaled_pf, cluster_cols=F, clustering_method="complete", color=inferno(10), show_colnames=F, annotation_col=ann_col, annotation_row=ann_row, main="Protein folding")

#repeating this for another upregulated pathway from metascape ('Cytokine signalling')
enriched_gene_set_1=as.character(mscape[24,5])
candidate_genes=unlist(strsplit(enriched_gene_set_1, ","))

candidate_genes = c(candidate_genes)
gene_data = EM.s_4H[candidate_genes,]
gene_data = data.frame(t(gene_data))
gene_data$sample_group = ss$SAMPLE_GROUP
levels(gene_data$sample_group)
gene_data$sample_group= factor(gene_data$sample_group, levels=c("Mock", "TR_4H"))
#melt this to form a table with fewer columns
gene_data.m=melt(gene_data, id.vars="sample_group")

ggp=ggplot(gene_data.m, aes(x=variable, y=value, fill=sample_group)) + geom_boxplot() +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(title="Cytokine signalling", x="Gene name", y="Scaled expression")
ggp

em_scaled_cyto=EM.s_4H[candidate_genes,]
pheatmap(em_scaled_cyto, cluster_cols=F, clustering_method="complete", color=inferno(10), show_colnames=F, annotation_col=ann_col, annotation_row=ann_row, main="Cytokine signalling")

#Repeating this for another upregulated pathway (NGF-stimulated transcription redundant cluster)
enriched_gene_set_1=as.character(mscape[97,5])
candidate_genes=unlist(strsplit(enriched_gene_set_1, ","))

candidate_genes = c(candidate_genes)
gene_data = EM.s_4H[candidate_genes,]
gene_data = data.frame(t(gene_data))
gene_data$sample_group = ss$SAMPLE_GROUP
levels(gene_data$sample_group)
gene_data$sample_group= factor(gene_data$sample_group, levels=c("Mock", "TR_4H"))
#melt this to form a table with fewer columns
gene_data.m=melt(gene_data, id.vars="sample_group")

ggp=ggplot(gene_data.m, aes(x=variable, y=value, fill=sample_group)) + geom_boxplot() +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(title="NGF-stimulated transcription", x="Gene name", y="Scaled expression")
ggp

em_scaled_NGF=EM.s_4H[candidate_genes,]
pheatmap(em_scaled_NGF, cluster_cols=F, clustering_method="complete", color=inferno(10), show_colnames=F, annotation_col=ann_col, annotation_row=ann_row, main="Cytokine signalling")



#Writing pathway em tables to file

write.table(em_scaled_NGF, file="C:/...em_NGF.csv", sep=",")
write.table(em_scaled_cyto, file="C:/...em_cyto.csv", sep=",")
write.table(em_scaled_pf, file="C:/...em_pf.csv", sep=",")


#Combined tables in excel and loading back in combined file
em_heatmap = read.table("C:/.../em_heat.csv", sep=",", header=TRUE, row.names=1)
pheatmap(em_heatmap, cluster_cols=F, cluster_rows=F, color=inferno(10), show_colnames=F, annotation_col=ann_col, annotation_row=ann_row)
