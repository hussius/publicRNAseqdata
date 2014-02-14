# Cluster reported F/RPKM values from public data sets.

# We will join the data sets on ENSEMBL ID:s, losing a lot of data in the process; but joining on gene symbols or something else would lead to an even worse loss. 

library(pheatmap)
library(org.Hs.eg.db)
library(NMF)
library(data.table)

### 
# Download data.
###

### "HPA" (Human Protein Atlas)
temp <- tempfile()
download.file(url="http://www.proteinatlas.org/download/rna.csv.zip",destfile=temp)
hpa <- read.csv(unz(temp, "rna.csv"))
unlink(temp)

hpa.heart <- hpa[hpa$Sample=="heart muscle", c("Gene", "Value")]
hpa.brain <- hpa[hpa$Sample=="cerebral cortex", c("Gene", "Value")]
hpa.kidney <- hpa[hpa$Sample=="kidney", c("Gene", "Value")]

hpa.fpkms <- merge(hpa.heart, hpa.brain, by="Gene")
hpa.fpkms <- merge(hpa.fpkms, hpa.kidney, by="Gene")
colnames(hpa.fpkms) <- c("ENSG_ID", "HPA_heart", "HPA_brain", "HPA_kidney")

length(hpa.fpkms[,1])
length(unique(hpa.fpkms[,1]))

write.table(hpa.fpkms,file="hpa_fpkms.txt",quote=F,sep="\t")

# If the table has already been saved:
# hpa.fpkms <- read.delim("hpa_fpkms.txt")

### "Altiso" Alternative isoform regulation in human tissue transcriptomes
temp <- tempfile()
download.file(url="http://genes.mit.edu/burgelab/Supplementary/wang_sandberg08/hg18.ensGene.CEs.rpkm.txt",destfile=temp)
altiso <- read.delim(temp, sep="\t")
unlink(temp)

# No kidney sample here, just use heart + brain
altiso.fpkms <- altiso[,c("X.Gene","heart","brain")]
colnames(altiso.fpkms) <- c("ENSG_ID", "AltIso_heart", "AltIso_brain")

length(altiso.fpkms[,1])
# [1] 23115
length(unique(altiso.fpkms[,1]))
# [1] 23115

write.table(altiso.fpkms,file="altiso_fpkms.txt",quote=F,sep="\t")
# altiso.fpkms <- read.delim("altiso_fpkms.txt")

### "GTEx" (Genotype-Tissue Expression)
# Note that this is a big download: 337.8 Mb (as of 2014-02-04)
temp <- tempfile()
download.file(url="http://www.broadinstitute.org/gtex/rest/file/download?portalFileId=119363&forDownload=true",destfile=temp)
header_lines <- readLines(temp, n=2)
gtex <- read.delim(temp, skip=2, sep="\t")
unlink(temp)

download.file(url="http://www.broadinstitute.org/gtex/rest/file/download?portalFileId=119273&forDownload=true",destfile="GTEx_description.txt")

metadata <- read.delim("GTEx_description.txt", sep="\t")

## GTEX: select a random heart sample 
random.heart <- sample(which(metadata$SMTS=="Heart"), size=1)
random.heart.samplename <- metadata[random.heart, "SAMPID"]
gtex.heart.fpkm <- as.numeric(gtex[,random.heart])

## GTEx: select a random brain sample
random.brain <- sample(which(metadata$SMTS=="Brain"), size=1)
random.brain.samplename <- metadata[random.heart, "SAMPID"]
gtex.brain.fpkm <- as.numeric(gtex[,random.brain])

## GTEx: select a random kidney sample
random.kidney <- sample(which(metadata$SMTS=="Kidney"), size=1)
random.kidney.samplename <- metadata[random.kidney, "SAMPID"]
gtex.kidney.fpkm <- as.numeric(gtex[,random.kidney])

## GTEx: Get gene IDs on same format as the other data sets
gtex.names <- gtex[,"Name"]
temp_list <- strsplit(as.character(gtex.names), split="\\.")
gtex.names.nodot <- unlist(temp_list)[2*(1:length(gtex.names))-1]

gtex.fpkms <- data.frame(ENSG_ID=gtex.names.nodot, GTEx_heart=gtex.heart.fpkm, GTEx_brain=gtex.brain.fpkm,GTEx_kidney=gtex.kidney.fpkm)

length(gtex.fpkms[,1])
# [1] 52576
length(unique(gtex.fpkms[,1]))
# [1] 52576

write.table(gtex.fpkms,file="gtex_fpkms.txt",quote=F,sep="\t")
# gtex.fpkms <- read.delim("gtex_fpkms.txt")

###
# Merge HPA, AltIso, GTEx into one object f. (RNA-seq Atlas will be used later)
###
fpkms <- merge(hpa.fpkms, altiso.fpkms, by="ENSG_ID")
fpkms <- merge(fpkms, gtex.fpkms, by="ENSG_ID")
gene_id <- fpkms[,1]
f <- fpkms[,2:ncol(fpkms)]
rownames(f) <- gene_id

write.table(f, file="fpkms.txt", quote=F)

# Heatmap from Pearson correlation
pheatmap(cor(f))
# Heatmap from Spearman correlation
pheatmap(cor(f, method="spearman")) # Not so bad actually
# PCA 
p <- prcomp(t(f))
plot(p$x[,2],p$x[,3],pch=".")
text(p$x[,2],p$x[,3],labels=rownames(p$x),cex=0.75)

# NMF
# need to remove all-zero rows (should maybe be done anyway, before PCA?)
f.nz <- f[-which(rowSums(f)==0),]
n <- nmf(f.nz, rank=3)
# Metagene expression profiles / mixture coefficient matrix / projections
h <- coef(n)
plot(h[1,],h[2,],col=rep(c(1,2,3),2),pch=20)
textxy(h[1,],h[2,],labs=colnames(h))

### "Atlas" (RNA-seq Atlas)
temp <- tempfile()
download.file(url="http://medicalgenomics.org/rna_seq_atlas/download?download_revision1=1",destfile=temp)
atlas <- read.delim(temp, sep="\t")
unlink(temp)

atlas.fpkms <- atlas[,c("ensembl_gene_id","heart","hypothalamus","kidney")]
colnames(atlas.fpkms) <- c("ENSG_ID","Atlas_heart","Atlas_brain","Atlas_kidney")
write.table(atlas.fpkms,file="atlas_fpkms.txt",quote=F,sep="\t")
# atlas.fpkms <- read.delim("atlas_fpkms.txt")

### Try to include Atlas also. 
length(unique(atlas$entrez_gene_id)) # 19,413
length(unique(atlas$ensembl_gene_id)) # 5,772
length(unique(atlas$hgnc_symbol)) # 21,399

# Approach (1). Merge on ENSEMBL genes as given in RNA-seq Atlas.
# Note that there are repeated ENSG ID:s in RNA-seq Atlas
# (as opposed to the other data sets)

data.dt <- data.table(atlas.fpkms)
setkey(data.dt, ENSG_ID)
temp <- data.dt[, lapply(.SD, sum), by=ENSG_ID]
collapsed <- as.data.frame(temp)
atlas.fpkms.summed <- collapsed[,2:ncol(collapsed)] 
rownames(atlas.fpkms.summed) <- collapsed[,1]

atlas.fpkms.summed <- atlas.fpkms.summed[2:nrow(atlas.fpkms.summed),]

fpkms <- merge(fpkms, atlas.fpkms.summed, by.x="ENSG_ID", by.y=0)
gene_id <- fpkms[,1]
f <- fpkms[,2:ncol(fpkms)]
rownames(f) <- gene_id

# Only 4506 genes left!
# Heatmap from Pearson correlation
pheatmap(cor(f))
# Heatmap from Spearman correlation
pheatmap(cor(f, method="spearman")) # Not so bad actually
# PCA 
p <- prcomp(t(f))
plot(p$x[,1],p$x[,2],pch=".")
text(p$x[,1],p$x[,2],labels=rownames(p$x),cex=0.75)

# log scale
library(limma)
log.fpkm <- voom(f)
pheatmap(cor(log.fpkm$E))
pheatmap(cor(log.fpkm$E, method="spearman"))

colors <- c(1,2,3,1,2,1,2,3,1,2,3)

p <- prcomp(t(log.fpkm$E))

par(mfrow=c(4,4))
for (i in 1:6){
	for(j in 1:6){
		if (i<j){ 
		plot(p$x[,i],p$x[,j],pch=20,col=colors,xlab=paste("PC", i),ylab=paste("PC", j))
		}
	}
}

# Alternatives:

# table(table(atlas$ensembl_gene_id))
# table(table(atlas$entrez_gene_id))
# table(table(atlas$hgnc_symbol))
# table(table(atlas$transcript)) # <<- almost unique

# Try to map (RefSeq) transcripts, which are almost a unique ID,
# to Entrez genes and from there to ENSEMBL
m <- org.Hs.egREFSEQ2EG
mapped_genes <- mappedkeys(m)
entrez.for.refseq <- as.list(m[mapped_genes])
remapped.entrez <- entrez.for.refseq[atlas$transcript]
rem.entr <- as.character(remapped.entrez)
names(rem.entr) <- names(remapped.entrez)
# Results in just 4121 ... maybe not a good idea

# Try to map Entrez symbols to ENSEMBL to recover more ENSG IDs than already present in the table. Does it matter which direction?
# Take the names from "fpkms" and convert to Entrez, or vice versa?
m <- org.Hs.egENSEMBL
mapped_genes <- mappedkeys(m)
ensg.for.entrez <- as.list(m[mapped_genes])
remapped.ensg <- ensg.for.entrez[as.character(atlas$entrez_gene_id)]

atlas.fpkms$remapped_ensg <- as.character(remapped.ensg)

# And add expression values
data.dt <- data.table(atlas.fpkms[,2:ncol(atlas.fpkms)])
setkey(data.dt, remapped_ensg)
temp <- data.dt[, lapply(.SD, sum), by=remapped_ensg]
collapsed <- as.data.frame(temp)
atlas.fpkms.summed <- collapsed[,2:ncol(collapsed)] 
rownames(atlas.fpkms.summed) <- collapsed[,1]

# Many more recovered!!
fpkms <- merge(fpkms, atlas.fpkms.summed, by.x="ENSG_ID", by.y=0)
gene_id <- fpkms[,1]
f <- fpkms[,2:ncol(fpkms)]
rownames(f) <- gene_id

# 13537 common IDs, cool.

# Heatmap from Pearson correlation
pheatmap(cor(f))
# Heatmap from Spearman correlation
pheatmap(cor(f, method="spearman")) # Not so bad actually
# PCA 
p <- prcomp(t(f))
plot(p$x[,1],p$x[,2],pch=".")
text(p$x[,1],p$x[,2],labels=rownames(p$x),cex=0.75)

# log scale
library(limma)
log.fpkm <- voom(f)
pheatmap(cor(log.fpkm$E))
pheatmap(cor(log.fpkm$E, method="spearman"))

colors <- c(1,2,3,1,2,1,2,3,1,2,3)

p <- prcomp(t(log.fpkm$E)) # PC1 vs PC5 pretty good here!
p <- prcomp(t(f)) # No outstanding combo

par(mfrow=c(4,4))
for (i in 1:6){
	for(j in 1:6){
		if (i<j){ 
		plot(p$x[,i],p$x[,j],pch=20,col=colors,xlab=paste("PC", i),ylab=paste("PC", j))
		}
	}
}


