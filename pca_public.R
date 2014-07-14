#Statistical test-correlation test and Kruskalwallys

pca_study_public <- read.delim("published_rpkms.txt",sep=" ")

head(pca_study_public)
dim(pca_study_public)

#Filter not used,Zeros count values can be specific to tissues 
filter <- apply(pca_study_public[,2:11], 1, function(row) all(row !=0 ))
pca_study_public_nonzero <- pca_study_public[filter,]

pca <- prcomp(t(pca_study_public[,]))
#pca <- svd(pca_study_public[,])

var.percent <- ((pca$sdev)^2)/sum(pca$sdev^2) *100
#var.percent <- (pca$d^2 / sum(pca$d^2) ) *100

barplot(var.percent[1:5], xlab="PC", ylab="Percent Variance",names.arg=1:length(var.percent[1:5]), las=1,ylim=c(0,max(var.percent[1:5])+10), col="gray")


x <- pca$x
plot(pca)
summary(pca)
screeplot(pca,type=c("lines"))

sampleinfo <- read.table("sample_info_published.txt",header=TRUE)

patientsample_NRaw<-cbind(as.character(sampleinfo$Study_labels),
                          as.numeric(sampleinfo$NumberRaw),
                          as.character(sampleinfo$Tissue),
                          as.character(sampleinfo$Preparation),
                          as.character(sampleinfo$readlength),
                          as.numeric(sampleinfo$Numbermapped),
                          as.character(sampleinfo$Readtype))

samplenames <- patientsample_NRaw[,-1]
rownames(samplenames) <- patientsample_NRaw[,1]
pca_comps<-x[,0:7]


Nraw_pca.comp_na <- as.matrix(merge(samplenames,pca_comps,by="row.names",all=TRUE))
Nraw_pca.comp<-na.omit(data.frame(Nraw_pca.comp_na))


#NumberRaw
plot(as.numeric(Nraw_pca.comp[,8]),Nraw_pca.comp[,2],xlab="NRaw",ylab="PCA_C1")

pval_Nraw_pc1<-cor.test(as.numeric(Nraw_pca.comp[,8]),as.numeric(Nraw_pca.comp[,2]),method = "spearman")$p.value
pval_Nraw_pc2<-cor.test(as.numeric(Nraw_pca.comp[,9]),as.numeric(Nraw_pca.comp[,2]),method = "spearman")$p.value
pval_Nraw_pc3<-cor.test(as.numeric(Nraw_pca.comp[,10]),as.numeric(Nraw_pca.comp[,2]),method = "spearman")$p.value
pval_Nraw_pc4<-cor.test(as.numeric(Nraw_pca.comp[,11]),as.numeric(Nraw_pca.comp[,2]),method = "spearman")$p.value
pval_Nraw_pc5<-cor.test(as.numeric(Nraw_pca.comp[,12]),as.numeric(Nraw_pca.comp[,2]),method = "spearman")$p.value
pval_Nraw_pc6<-cor.test(as.numeric(Nraw_pca.comp[,13]),as.numeric(Nraw_pca.comp[,2]),method = "spearman")$p.value
pval_Nraw_pc7<-cor.test(as.numeric(Nraw_pca.comp[,14]),as.numeric(Nraw_pca.comp[,2]),method = "spearman")$p.value

cat(sprintf("Number_of_rawreads~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"PCA4=\"%f\"PCA5=\"%f\"PCA6=\"%f\"\n", pval_Nraw_pc1,pval_Nraw_pc2,pval_Nraw_pc3,pval_Nraw_pc4,pval_Nraw_pc5,pval_Nraw_pc6,pval_Nraw_pc7))
#print(pval_agecorr_pc1)


#NumberMapped
plot(as.numeric(Nraw_pca.comp[,8]),Nraw_pca.comp[,6],xlab="NRaw",ylab="PCA_C1")

pval_Nmap_pc1<-cor.test(as.numeric(Nraw_pca.comp[,8]),as.numeric(Nraw_pca.comp[,6]),method = "spearman")$p.value
pval_Nmap_pc2<-cor.test(as.numeric(Nraw_pca.comp[,9]),as.numeric(Nraw_pca.comp[,6]),method = "spearman")$p.value
pval_Nmap_pc3<-cor.test(as.numeric(Nraw_pca.comp[,10]),as.numeric(Nraw_pca.comp[,6]),method = "spearman")$p.value
pval_Nmap_pc4<-cor.test(as.numeric(Nraw_pca.comp[,11]),as.numeric(Nraw_pca.comp[,6]),method = "spearman")$p.value
pval_Nmap_pc5<-cor.test(as.numeric(Nraw_pca.comp[,12]),as.numeric(Nraw_pca.comp[,6]),method = "spearman")$p.value
pval_Nmap_pc6<-cor.test(as.numeric(Nraw_pca.comp[,13]),as.numeric(Nraw_pca.comp[,6]),method = "spearman")$p.value
pval_Nmap_pc7<-cor.test(as.numeric(Nraw_pca.comp[,14]),as.numeric(Nraw_pca.comp[,6]),method = "spearman")$p.value

cat(sprintf("Number_of_mapped_reads~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"PCA4=\"%f\"PCA5=\"%f\"PCA6=\"%f\"\n", pval_Nmap_pc1,pval_Nmap_pc2,pval_Nmap_pc3,pval_Nmap_pc4,pval_Nmap_pc5,pval_Nmap_pc6,pval_Nmap_pc7))

#Tissue
pval_tissue_pc1<-kruskal.test(as.numeric(Nraw_pca.comp[,8]) ~ as.factor(Nraw_pca.comp[,3]))$p.value
pval_tissue_pc2<-kruskal.test(as.numeric(Nraw_pca.comp[,9]) ~ as.factor(Nraw_pca.comp[,3]))$p.value
pval_tissue_pc3<-kruskal.test(as.numeric(Nraw_pca.comp[,10]) ~ as.factor(Nraw_pca.comp[,3]))$p.value
pval_tissue_pc4<-kruskal.test(as.numeric(Nraw_pca.comp[,11]) ~ as.factor(Nraw_pca.comp[,3]))$p.value
pval_tissue_pc5<-kruskal.test(as.numeric(Nraw_pca.comp[,12]) ~ as.factor(Nraw_pca.comp[,3]))$p.value
pval_tissue_pc6<-kruskal.test(as.numeric(Nraw_pca.comp[,13]) ~ as.factor(Nraw_pca.comp[,3]))$p.value
pval_tissue_pc7<-kruskal.test(as.numeric(Nraw_pca.comp[,14]) ~ as.factor(Nraw_pca.comp[,3]))$p.value

cat(sprintf("Tissues~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"PCA4=\"%f\"PCA5=\"%f\"PCA6=\"%f\"\n", pval_tissue_pc1,pval_tissue_pc2,pval_tissue_pc3,pval_tissue_pc4,pval_tissue_pc5,pval_tissue_pc6,pval_tissue_pc7))

#Library_preparation
pval_tissue_pc1<-kruskal.test(as.numeric(Nraw_pca.comp[,8]) ~ as.factor(Nraw_pca.comp[,4]))$p.value
pval_tissue_pc2<-kruskal.test(as.numeric(Nraw_pca.comp[,9]) ~ as.factor(Nraw_pca.comp[,4]))$p.value
pval_tissue_pc3<-kruskal.test(as.numeric(Nraw_pca.comp[,10]) ~ as.factor(Nraw_pca.comp[,4]))$p.value
pval_tissue_pc4<-kruskal.test(as.numeric(Nraw_pca.comp[,11]) ~ as.factor(Nraw_pca.comp[,4]))$p.value
pval_tissue_pc5<-kruskal.test(as.numeric(Nraw_pca.comp[,12]) ~ as.factor(Nraw_pca.comp[,4]))$p.value
pval_tissue_pc6<-kruskal.test(as.numeric(Nraw_pca.comp[,13]) ~ as.factor(Nraw_pca.comp[,4]))$p.value
pval_tissue_pc7<-kruskal.test(as.numeric(Nraw_pca.comp[,14]) ~ as.factor(Nraw_pca.comp[,4]))$p.value

cat(sprintf("Library_prep~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"PCA4=\"%f\"PCA5=\"%f\"PCA6=\"%f\"\n", pval_tissue_pc1,pval_tissue_pc2,pval_tissue_pc3,pval_tissue_pc4,pval_tissue_pc5,pval_tissue_pc6,pval_tissue_pc7))

#pval_interferoncorr_pc1<-kruskal.test(interferon_pca.comp[complete.cases(interferon_pca.comp),][,2]  ~ interferon_pca.comp[complete.cases(interferon_pca.comp),][,9])$p.value
#sequence_type-paired/single
pval_seqtype_pc1<-kruskal.test(as.numeric(Nraw_pca.comp[,8]) ~ as.factor(Nraw_pca.comp[,7]))$p.value
pval_seqtype_pc2<-kruskal.test(as.numeric(Nraw_pca.comp[,9]) ~ as.factor(Nraw_pca.comp[,7]))$p.value
pval_seqtype_pc3<-kruskal.test(as.numeric(Nraw_pca.comp[,10]) ~ as.factor(Nraw_pca.comp[,7]))$p.value
pval_seqtype_pc4<-kruskal.test(as.numeric(Nraw_pca.comp[,11]) ~ as.factor(Nraw_pca.comp[,7]))$p.value
pval_seqtype_pc5<-kruskal.test(as.numeric(Nraw_pca.comp[,12]) ~ as.factor(Nraw_pca.comp[,7]))$p.value
pval_seqtype_pc6<-kruskal.test(as.numeric(Nraw_pca.comp[,13]) ~ as.factor(Nraw_pca.comp[,7]))$p.value
pval_seqtype_pc7<-kruskal.test(as.numeric(Nraw_pca.comp[,14]) ~ as.factor(Nraw_pca.comp[,7]))$p.value

cat(sprintf("Sequence_type~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"PCA4=\"%f\"PCA5=\"%f\"PCA6=\"%f\"\n", pval_seqtype_pc1,pval_seqtype_pc2,pval_seqtype_pc3,pval_seqtype_pc4,pval_seqtype_pc5,pval_seqtype_pc6,pval_seqtype_pc7))

#readlength- 1x32/1x35/1x76/2x100/2x76 
pval_readlength_pc1<-kruskal.test(as.numeric(Nraw_pca.comp[,8]) ~ as.factor(Nraw_pca.comp[,5]))$p.value
pval_readlength_pc2<-kruskal.test(as.numeric(Nraw_pca.comp[,9]) ~ as.factor(Nraw_pca.comp[,5]))$p.value
pval_readlength_pc3<-kruskal.test(as.numeric(Nraw_pca.comp[,10]) ~ as.factor(Nraw_pca.comp[,5]))$p.value
pval_readlength_pc4<-kruskal.test(as.numeric(Nraw_pca.comp[,11]) ~ as.factor(Nraw_pca.comp[,5]))$p.value
pval_readlength_pc5<-kruskal.test(as.numeric(Nraw_pca.comp[,12]) ~ as.factor(Nraw_pca.comp[,5]))$p.value
pval_readlength_pc6<-kruskal.test(as.numeric(Nraw_pca.comp[,13]) ~ as.factor(Nraw_pca.comp[,5]))$p.value
pval_readlength_pc7<-kruskal.test(as.numeric(Nraw_pca.comp[,14]) ~ as.factor(Nraw_pca.comp[,5]))$p.value

cat(sprintf("readlength~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"PCA4=\"%f\"PCA5=\"%f\"PCA6=\"%f\"\n", pval_readlength_pc1,pval_readlength_pc2,pval_readlength_pc3,pval_readlength_pc4,pval_readlength_pc5,pval_readlength_pc6,pval_readlength_pc7))

