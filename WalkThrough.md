Walkthrough with all data analyses that comes with the manuscript
Comparing published FPKM/RPKM values and reprocessed FPKM values.
========================================================

Prepare by loading packages etc.

```r
library(pheatmap)
library(reshape)
library(gplots)
```

```
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
library(ops)
```

```
## 
## Attaching package: 'ops'
## 
## The following object is masked from 'package:stats':
## 
##     filter
```

```r
library(calibrate)
```

```
## Loading required package: MASS
```

```r
library(biomaRt)
library(sva)
```

```
## Loading required package: corpcor
## Loading required package: mgcv
## Loading required package: nlme
## This is mgcv 1.8-3. For overview type 'help("mgcv-package")'.
```

Data from four different public sources are downloaded and brain, heart and kidney samples are extracted. 

Start with Human Protein Atlas (HPA):


```r
#temp <- tempfile()
#download.file(url="http://www.proteinatlas.org/download/rna.csv.zip",destfile=temp)
#hpa <- read.csv(unz(temp, "rna.csv"))
#unlink(temp)

#hpa.heart <- hpa[hpa$Sample=="heart muscle", c("Gene", "Value")]
#hpa.brain <- hpa[hpa$Sample=="cerebral cortex", c("Gene", "Value")]
#hpa.kidney <- hpa[hpa$Sample=="kidney", c("Gene", "Value")]

#hpa.fpkms <- merge(hpa.heart, hpa.brain, by="Gene")
#hpa.fpkms <- merge(hpa.fpkms, hpa.kidney, by="Gene")
#colnames(hpa.fpkms) <- c("ENSG_ID", "HPA_heart", "HPA_brain", "HPA_kidney")
```

Check if the identifiers are unique and write table to file.


```r
#length(hpa.fpkms[,1])
#length(unique(hpa.fpkms[,1]))

#write.table(hpa.fpkms,file="hpa_fpkms.txt",quote=F,sep="\t")
```

Next dataset is from the article "Alternative isoform regulation in human tissue transcriptomes.." by Wang et.al
AltIso:


```r
#temp <- tempfile()
#download.file(url="http://genes.mit.edu/burgelab/Supplementary/wang_sandberg08/hg18.ensGene.CEs.rpkm.txt",destfile=temp)
#altiso <- read.delim(temp, sep="\t")
#unlink(temp)
```

There is no kidney sample here, so just use heart and brain


```r
#altiso.fpkms <- altiso[,c("X.Gene","heart","brain")]
#colnames(altiso.fpkms) <- c("ENSG_ID", "AltIso_heart", "AltIso_brain")
```

Check uniqueness of IDs.


```r
#length(altiso.fpkms[,1])
#length(unique(altiso.fpkms[,1]))

#write.table(altiso.fpkms,file="altiso_fpkms.txt",quote=F,sep="\t")
```

Next dataset is derived from "GTEx": Genotype-Tissue Expression

This is a big download: 337.8 Mb (as of 2014-02-04)
We also add some code to randomly select one sample from each tissue type; there are many biological replicates in this data set.


```r
#temp <- tempfile()
#download.file(url="http://www.gtexportal.org/home/rest/file/download?portalFileId=175729&forDownload=true",destfile=temp)
#header_lines <- readLines(temp, n=2)
#gtex <- read.delim(temp, skip=2, sep="\t")
#unlink(temp)

#write.table(gtex, file="gtex_all.txt",   quote=F, sep="\t")

#download.file(url="http://www.gtexportal.org/home/rest/file/download?portalFileId=175707&forDownload=true",destfile="GTEx_description.txt")

#metadata <- read.delim("GTEx_description.txt", sep="\t")
```

The metadata table seems to contain entries that are not in the RPKM table.


```r
#samp.id <- gsub('-','.',metadata$SAMPID)
#eligible.samples <- which(samp.id %in% colnames(gtex))
#metadata <- metadata[eligible.samples,]
```

Select random heart, kidney and brain samples.


```r
#random.heart <- sample(which(metadata$SMTS=="Heart"), size=1)
#random.heart.samplename <- gsub('-','.',metadata[random.heart, "SAMPID"])
#gtex.heart.fpkm <- as.numeric(gtex[,random.heart.samplename])

#random.brain <- sample(which(metadata$SMTS=="Brain"), size=1)
#random.brain.samplename <- gsub('-','.',metadata[random.brain, "SAMPID"])
#gtex.brain.fpkm <- as.numeric(gtex[,random.brain.samplename])

#random.kidney <- sample(which(metadata$SMTS=="Kidney"), size=1)
#random.kidney.samplename <- gsub('-','.',metadata[random.kidney, "SAMPID"])
#gtex.kidney.fpkm <- as.numeric(gtex[,random.kidney.samplename])
```

Get gene IDs on same format as the other data sets by removing the part after the dot; check ID uniqueness and write to file.


```r
#gtex.names <- gtex[,"Name"]
#temp_list <- strsplit(as.character(gtex.names), split="\\.")
#gtex.names.nodot <- unlist(temp_list)[2*(1:length(gtex.names))-1]

#gtex.fpkms <- data.frame(ENSG_ID=gtex.names.nodot, GTEx_heart=gtex.heart.fpkm, GTEx_brain=gtex.brain.fpkm,GTEx_kidney=gtex.kidney.fpkm)

#length(gtex.fpkms[,1])
#length(unique(gtex.fpkms[,1]))

#write.table(gtex.fpkms,file="gtex_fpkms.txt",quote=F,sep="\t")
```

*RNA-seq Atlas*


```r
#temp <- tempfile()
#download.file(url="http://medicalgenomics.org/rna_seq_atlas/download?download_revision1=1",destfile=temp)
#atlas <- read.delim(temp, sep="\t")
#unlink(temp)

#atlas.fpkms <- atlas[,c("ensembl_gene_id","heart","hypothalamus","kidney")]
#colnames(atlas.fpkms) <- c("ENSG_ID","Atlas_heart","Atlas_brain","Atlas_kidney")
#write.table(atlas.fpkms,file="atlas_fpkms.txt",quote=F,sep="\t")
```

Combining F/RPKM values from public data sets
---------------------------------------------

We will join the data sets on ENSEMBL ID:s, losing a lot of data in the process - but joining on gene symbols or something else would lead to an even worse loss. 


```r
library(org.Hs.eg.db) # for transferring gene identifiers
```

```
## Loading required package: AnnotationDbi
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Loading required package: GenomeInfoDb
## 
## Attaching package: 'AnnotationDbi'
## 
## The following object is masked from 'package:MASS':
## 
##     select
## 
## Loading required package: DBI
```

```r
library(data.table) # for collapsing transcript RPKMs
library(pheatmap) # for nicer visualization
library(edgeR) # for TMM normalization
```

```
## Loading required package: limma
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
#hpa.fpkms <- read.delim("hpa_fpkms.txt")
#altiso.fpkms <- read.delim("altiso_fpkms.txt")
#gtex.fpkms <- read.delim("gtex_fpkms.txt")
#atlas.fpkms <- read.delim("atlas_fpkms.txt")
```

The RNA-seq Atlas data set uses many different identifiers, while the other all use ENSG as the primary identifier

Approach 1: Merge on ENSEMBL genes (ENSG) as given in RNA-seq Atlas. Note that there are repeated ENSG ID:s in RNA-seq Atlas, as opposed to the other data sets, so we need to do something about that. In this case, we just sum the transcripts that belong to each ENSG gene. We use data.table for this.


```r
#data.dt <- data.table(atlas.fpkms)
#setkey(data.dt, ENSG_ID)
#temp <- data.dt[, lapply(.SD, sum), by=ENSG_ID]
#collapsed <- as.data.frame(temp)
#atlas.fpkms.summed <- collapsed[,2:ncol(collapsed)] 
#rownames(atlas.fpkms.summed) <- collapsed[,1]

#atlas.fpkms.summed <- atlas.fpkms.summed[2:nrow(atlas.fpkms.summed),]
```

Finally, combine all the data sets into a data frame.


```r
#fpkms <- merge(hpa.fpkms, altiso.fpkms, by="ENSG_ID")
#fpkms <- merge(fpkms, gtex.fpkms, by="ENSG_ID")
#fpkms <- merge(fpkms, atlas.fpkms.summed, by.x="ENSG_ID", by.y=0)
#gene_id <- fpkms[,1]
#f <- fpkms[,2:ncol(fpkms)]
#rownames(f) <- gene_id
```

Check how many ENSG IDs we have left.


```r
#dim(f)
```

Approach 2: Try to map Entrez symbols to ENSEMBL to recover more ENSG IDs than already present in the table. 


```r
#m <- org.Hs.egENSEMBL
#mapped_genes <- mappedkeys(m)
#ensg.for.entrez <- as.list(m[mapped_genes])
#remapped.ensg <- ensg.for.entrez[as.character(atlas$entrez_gene_id)]

#atlas.fpkms$remapped_ensg <- as.character(remapped.ensg)

# And add expression values
#data.dt <- data.table(atlas.fpkms[,2:ncol(atlas.fpkms)])
#setkey(data.dt, remapped_ensg)
#temp <- data.dt[, lapply(.SD, sum), by=remapped_ensg]
#collapsed <- as.data.frame(temp)
#atlas.fpkms.summed <- collapsed[,2:ncol(collapsed)] 
#rownames(atlas.fpkms.summed) <- collapsed[,1]
```

Combine data sets again


```r
#fpkms <- merge(hpa.fpkms, altiso.fpkms, by="ENSG_ID")
#fpkms <- merge(fpkms, gtex.fpkms, by="ENSG_ID")
#fpkms <- merge(fpkms, atlas.fpkms.summed, by.x="ENSG_ID", by.y=0)
#gene_id <- fpkms[,1]
#f <- fpkms[,2:ncol(fpkms)]
#rownames(f) <- gene_id
#write.table(f, file = 'published_rpkms.txt', quote=F)

#instead of downloading everytime:
```
Download the data from local file:


```r
published <- read.delim("published_rpkms.txt",sep=" ")

sampleinfo_published <- read.table("sample_info_published.txt",header=TRUE)
```

The published FPKM values are first filtered by removing all lines where FPKM is less or equal to 0.01 in all samples:


```r
published.nozero <- published[-which(rowSums(published[,])<=0.01),]
```

Heatmap of Spearman correlations between published expression profiles (# genes = 13,323), Figure 1b:


```r
pheatmap(cor(published.nozero, method="spearman")) 
```

![plot of chunk :published heatmap spearman](figure/:published heatmap spearman.png) 

Alternatively, one could use Pearson correlation (not shown in paper):


```r
pheatmap(cor(published.nozero))
```

![plot of chunk :published heatmap pearson](figure/:published heatmap pearson.png) 

PCA analysis of published FPKM values, Figure 1c:


```r
colors <- c("indianred", "dodgerblue", "forestgreen",
            "indianred", "dodgerblue",
            "indianred", "dodgerblue", "forestgreen", 
            "indianred", "dodgerblue", "forestgreen")
  
p <- prcomp(t(published.nozero))

shapes <- c(rep(15,3),rep(16,2),rep(17,3),rep(8,3))

plot(p$x[,1],p$x[,2],pch=shapes,cex=1.5,col=colors,xlab=paste("PC1 58% of variance"),ylab=paste("PC2 13% of variance"),main="Published FPKM/RPKM values \n n=13,323")
legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20)
legend("top",legend=c("HPA","AltIso","GTEx","Atlas"),col="black",pch=c(15,16,17,8),ncol=2)
```

![plot of chunk :published PCA](figure/:published PCA.png) 

We can plot all pairwise combinations of principal components 1 to 5. (not shown in paper):


```r
par(mfrow=c(4,4))
for (i in 1:6){
  for(j in 1:6){
  	if (i<j){ 
		plot(p$x[,i],p$x[,j],pch=20,col=colors,xlab=paste("PC",i),ylab=paste("PC",j),main="Published FPKM values \n n=13,323")
		}
	}
}
```

![plot of chunk :published pairwise PCA](figure/:published pairwise PCA.png) 

Look a bit closer at PCs 1-3 in prcomp (not shown in paper):


```r
      #PC1:
      load.pc1 <- p$rotation[,1][order(p$rotation[,1])]
      extreme.pc1 <- c(tail(load.pc1), head(load.pc1))

      extreme.pc1.ensg <- names(extreme.pc1)
      ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl") #select the ensembl database
      extreme.pc1.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc1.ensg,
                           mart=ensembl)
      
      q <- extreme.pc1.symbols[,2]
      names(q) <- extreme.pc1.symbols[,1]

      fpkm.pc1 <- cbind(q[extreme.pc1.ensg],published.nozero[extreme.pc1.ensg,])
      names(fpkm.pc1)[names(fpkm.pc1) == 'q[extreme.pc1.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc1, names.arg=q[extreme.pc1.ensg],las=2,main="Genes w highest absolute loadings in PC1 (raw RPKM)")
```

![plot of chunk :published PC 1-3](figure/:published PC 1-31.png) 

```r
      print(fpkm.pc1)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000136872       ALDOB       0.4       2.4     1824.1         0.00
## ENSG00000123560        PLP1      11.9    1295.2        0.7         2.56
## ENSG00000118785        SPP1       8.3     368.2     1251.1         7.74
## ENSG00000120885         CLU     336.9    3545.9      631.6        58.72
## ENSG00000087086         FTL     436.5     487.2     1523.1       463.86
## ENSG00000131095        GFAP       3.0    1372.0        1.4         1.70
## ENSG00000111245        MYL2    5291.5       0.2        2.6     16780.36
## ENSG00000118194       TNNT2    2693.8       7.0       20.8      8074.39
## ENSG00000198125          MB    3937.3       0.9        2.3      7065.54
## ENSG00000175084         DES    3403.6       2.4       10.2      3160.05
## ENSG00000159251       ACTC1    2914.4       1.6        0.4      3457.43
## ENSG00000092054        MYH7    1802.4       2.1        0.8      4117.17
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000136872         0.49  1.630e-01     0.2298   1927.8312       1.452
## ENSG00000123560      1544.93  2.432e+00  1164.7164      0.4163       5.086
## ENSG00000118785       309.26  5.712e-01   317.7870   2806.6497      16.293
## ENSG00000120885      1682.03  2.138e+01   307.3469    365.1711      29.357
## ENSG00000087086       407.35  5.738e+02  2099.7219   5310.5093      16.930
## ENSG00000131095       834.13  1.671e+00  3633.8174      0.8538       0.549
## ENSG00000111245         1.03  1.577e+04    10.2006      2.3143    1064.500
## ENSG00000118194         9.37  2.028e+03     1.6755     11.8667   10930.800
## ENSG00000198125         1.35  2.706e+03     2.6915      0.6978    1415.554
## ENSG00000175084         2.65  5.524e+03    13.1029      4.9042    1675.890
## ENSG00000159251         0.42  5.320e+03     3.6530      1.0412     310.068
## ENSG00000092054         0.83  4.714e+03     6.4870      0.6887    2137.510
##                 Atlas_brain Atlas_kidney
## ENSG00000136872       0.125      788.910
## ENSG00000123560     779.755        0.513
## ENSG00000118785     121.874      878.535
## ENSG00000120885     394.746       47.363
## ENSG00000087086       7.700       40.892
## ENSG00000131095     691.908        0.135
## ENSG00000111245       0.000        0.358
## ENSG00000118194       8.713       43.633
## ENSG00000198125       0.032        0.000
## ENSG00000175084       4.500        6.252
## ENSG00000159251       0.421        0.311
## ENSG00000092054       0.635        0.354
```

```r
      #PC2:
      load.pc2 <- p$rotation[,2][order(p$rotation[,2])]
      extreme.pc2 <- c(tail(load.pc2), head(load.pc2))
      
      extreme.pc2.ensg <- names(extreme.pc2)
      extreme.pc2.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc2.ensg,
                           mart=ensembl)
      
      q <- extreme.pc2.symbols[,2]
      names(q) <- extreme.pc2.symbols[,1]

      fpkm.pc2 <- cbind(q[extreme.pc2.ensg],published.nozero[extreme.pc2.ensg,])
      names(fpkm.pc2)[names(fpkm.pc2) == 'q[extreme.pc2.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc2, names.arg=q[extreme.pc2.ensg],las=2,main="Genes w highest absolute loadings in PC2 (raw RPKM)")
```

![plot of chunk :published PC 1-3](figure/:published PC 1-32.png) 

```r
      print(fpkm.pc2)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000125971     DYNLRB1     207.8     164.7      170.0       1533.0
## ENSG00000160808        MYL3    1213.4       3.2       24.9       3909.8
## ENSG00000087086         FTL     436.5     487.2     1523.1        463.9
## ENSG00000075624        ACTB     423.6     786.5      521.8        377.9
## ENSG00000111245        MYL2    5291.5       0.2        2.6      16780.4
## ENSG00000167996        FTH1     363.8     366.9      498.9       3563.3
## ENSG00000118194       TNNT2    2693.8       7.0       20.8       8074.4
## ENSG00000140416        TPM1    4937.5      61.9      128.6       4379.1
## ENSG00000148677      ANKRD1    2451.5       0.2        0.3       1088.9
## ENSG00000175206        NPPA    6693.0       8.1        0.1        193.7
## ENSG00000188257     PLA2G2A      51.2       0.8        2.3        210.7
## ENSG00000106631        MYL7    4606.5       0.2        0.1        896.0
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000125971      2293.04     21.458    42.5226     32.5899       20.87
## ENSG00000160808         0.84   3996.971     5.8288     19.8493       80.67
## ENSG00000087086       407.35    573.762  2099.7219   5310.5093       16.93
## ENSG00000075624      3277.02    283.548   782.8235    653.1090       67.06
## ENSG00000111245         1.03  15765.179    10.2006      2.3143     1064.50
## ENSG00000167996      9086.14    337.276  1151.6777   1828.0887       43.50
## ENSG00000118194         9.37   2027.960     1.6755     11.8667    10930.80
## ENSG00000140416        68.39   1132.401    17.7402     31.1861     4228.46
## ENSG00000148677         0.10   1471.472     2.0740      0.5602     3863.61
## ENSG00000175206         1.74    137.694     0.8982      2.7192      311.40
## ENSG00000188257         0.00      2.314     0.4622      0.4419     2351.22
## ENSG00000106631         0.00    214.278     1.5850      0.5706      115.29
##                 Atlas_brain Atlas_kidney
## ENSG00000125971      42.315       23.413
## ENSG00000160808       0.204        2.634
## ENSG00000087086       7.700       40.892
## ENSG00000075624     140.560       69.191
## ENSG00000111245       0.000        0.358
## ENSG00000167996      53.459       66.455
## ENSG00000118194       8.713       43.633
## ENSG00000140416     105.073      189.505
## ENSG00000148677       0.403        0.263
## ENSG00000175206       0.445        0.131
## ENSG00000188257       3.057       24.728
## ENSG00000106631       0.000        0.000
```

```r
      #PC3:
      load.pc3 <- p$rotation[,3][order(p$rotation[,3])]
      extreme.pc3 <- c(tail(load.pc3), head(load.pc3))
      
      extreme.pc3.ensg <- names(extreme.pc3)
      extreme.pc3.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc3.ensg,
                           mart=ensembl)
      
      q <- extreme.pc3.symbols[,2]
      names(q) <- extreme.pc3.symbols[,1]

      fpkm.pc3 <- cbind(q[extreme.pc3.ensg],published.nozero[extreme.pc3.ensg,])
      names(fpkm.pc3)[names(fpkm.pc3) == 'q[extreme.pc3.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc3, names.arg=q[extreme.pc3.ensg],las=2,main="Genes w highest absolute loadings in PC3 (raw RPKM)")
```

![plot of chunk :published PC 1-3](figure/:published PC 1-33.png) 

```r
      print(fpkm.pc3)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000111245        MYL2    5291.5       0.2        2.6      16780.4
## ENSG00000101608      MYL12A    2033.5      32.2      214.9        529.0
## ENSG00000175206        NPPA    6693.0       8.1        0.1        193.7
## ENSG00000175084         DES    3403.6       2.4       10.2       3160.1
## ENSG00000159251       ACTC1    2914.4       1.6        0.4       3457.4
## ENSG00000129991       TNNI3    2157.5       0.5        0.1       2698.4
## ENSG00000167996        FTH1     363.8     366.9      498.9       3563.3
## ENSG00000118194       TNNT2    2693.8       7.0       20.8       8074.4
## ENSG00000125971     DYNLRB1     207.8     164.7      170.0       1533.0
## ENSG00000075624        ACTB     423.6     786.5      521.8        377.9
## ENSG00000071082       RPL31     632.0     348.9      561.3       2239.9
## ENSG00000106211       HSPB1     521.6      37.7      196.6       3971.4
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000111245         1.03   15765.18    10.2006      2.3143     1064.50
## ENSG00000101608        28.70    1968.66    54.9994    141.2811      300.41
## ENSG00000175206         1.74     137.69     0.8982      2.7192      311.40
## ENSG00000175084         2.65    5524.18    13.1029      4.9042     1675.89
## ENSG00000159251         0.42    5319.64     3.6530      1.0412      310.07
## ENSG00000129991         0.47    5751.82     4.0958      0.6091      999.56
## ENSG00000167996      9086.14     337.28  1151.6777   1828.0887       43.50
## ENSG00000118194         9.37    2027.96     1.6755     11.8667    10930.80
## ENSG00000125971      2293.04      21.46    42.5226     32.5899       20.87
## ENSG00000075624      3277.02     283.55   782.8235    653.1090       67.06
## ENSG00000071082      1678.90      38.42    99.8055     51.5868      378.92
## ENSG00000106211       643.46     333.52    86.1187    201.8283       61.19
##                 Atlas_brain Atlas_kidney
## ENSG00000111245       0.000        0.358
## ENSG00000101608      15.476       24.862
## ENSG00000175206       0.445        0.131
## ENSG00000175084       4.500        6.252
## ENSG00000159251       0.421        0.311
## ENSG00000129991       0.675        0.132
## ENSG00000167996      53.459       66.455
## ENSG00000118194       8.713       43.633
## ENSG00000125971      42.315       23.413
## ENSG00000075624     140.560       69.191
## ENSG00000071082     484.379      259.057
## ENSG00000106211       6.808       17.175
```
             
Try Anova on a "melted" expression matrix with some metadata, Figure 1d:


```r
m <- melt(published.nozero)
```

```
## Using  as id variables
```

```r
colnames(m) <- c("sample_ID","RPKM")

meta <- sampleinfo_published[,c("Study","Tissue","Preparation","NumberRaw","Numbermapped","Readtype")]
rownames(meta) <- colnames(published.nozero)
tissue <- rep(meta$Tissue, each=nrow(published.nozero))
study <- rep(meta$Study, each=nrow(published.nozero))
prep <- rep(meta$Preparation, each=nrow(published.nozero))
layout <- rep(meta$Readtype, each=nrow(published.nozero))
raw <- rep(meta$NumberRaw, each=nrow(published.nozero))
mapped <- rep(meta$Numbermapped, each=nrow(published.nozero))
data <- data.frame(m, tissue=tissue, study=study, prep=prep, layout=layout,nraw=raw,nmapped=mapped)
fit <- lm(RPKM ~ layout + prep + nraw + study + tissue, data=data)
a <- anova(fit)
maxval = 100

barplot(100*a$"Sum Sq"[1:5]/sum(a$"Sum Sq"[1:5]),names.arg=rownames(a[1:5,]),main="Anova, published FPKM/RPKM values",ylim=c(0,maxval))
```

![plot of chunk :published anova](figure/:published anova.png) 

Try log2 transformation of the published FPKM values:


```r
pseudo <- 1
published.log <- log2(published.nozero + pseudo)
```

Heatmap of Spearman correlations between published expression profiles with log2 values, Figure 2a:


```r
pheatmap(cor(published.log),method="spearman")
```

![plot of chunk :published log heatmap spearman](figure/:published log heatmap spearman.png) 

PCA analysis of log2 published FPKM values, Figure 2b:


```r
colors <- c("indianred", "dodgerblue", "forestgreen",
            "indianred", "dodgerblue",
            "indianred", "dodgerblue", "forestgreen", 
            "indianred", "dodgerblue", "forestgreen")

p.log <- prcomp(t(published.log))

shapes <- c(rep(15,3),rep(16,2),rep(17,3),rep(8,3))

plot(p.log$x[,1],p.log$x[,2],pch=shapes,col=colors,xlab=paste("PC1 31% of variance"),ylab=paste("PC2 27% of variance"),main="log2 Published FPKM/RPKM values \n n=13,323")
legend("bottomright",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("top",legend=c("HPA","AltIso","GTEx","Atlas"),col="black",pch=c(15,16,17,8),ncol=2)
```

![plot of chunk :published log PCA 1&2](figure/:published log PCA 1&2.png) 

Figure 2c:


```r
plot(p.log$x[,2],p.log$x[,3],pch=shapes,col=colors,xlab=paste("PC2 27% of variance"),ylab=paste("PC3 19% of variance"),main="log2 Published FPKM values \n n=13,323")
legend("bottomright",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("topright",legend=c("HPA","AltIso","GTEx","Atlas"),col="black",pch=c(15,16,17,8),ncol=2)
```

![plot of chunk :published log PCA 2&3](figure/:published log PCA 2&3.png) 

Cross-validattion by leaving out one of the studies, in this case AltIso, Figure 2d: 


```r
p.loo <- published.log[,-c(4,5)]
colors.loo <- colors[-c(4,5)]
p.loo.log <- prcomp(t(p.loo))

p.add <- published.log[,c(4,5)]
projection <- t(p.add) %*% p.loo.log$rotation

p.original.plus.new <- rbind(p.loo.log$x, projection)
col.original.plus.new <- c(colors.loo, colors[c(4,5)])
plot(p.original.plus.new[,2],p.original.plus.new[,3],pch=c(rep(15,3),rep(17,3),rep(8,3),rep(22,nrow(projection))),col=col.original.plus.new,xlab="PC2",ylab="PC3",main="log2 Published FPKM/RPKM values; AltIso projected onto existing PCs \n n=13,323",xlim=c(-150,100))

legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("top",legend=c("HPA","GTEx","Atlas","AltIso"),col="black",pch=c(15,17,8,22),ncol=2)
```

![plot of chunk :leave-one-out-pca](figure/:leave-one-out-pca.png) 

We can plot all pairwise combinations of principal components 1 to 5. (not shown in paper):


```r
par(mfrow=c(4,4))
for (i in 1:6){
  for(j in 1:6){
    if (i<j){ 
		plot(p.log$x[,i],p.log$x[,j],pch=shapes,col=colors,xlab=paste("PC",i),ylab=paste("PC",j),main="log2 Published FPKM values \n n=13323")
		}
	}
}
```

![plot of chunk :published log PCA pairwise](figure/:published log PCA pairwise.png) 

Look a bit closer at PCs 1-3 in prcomp (not shown in paper):


```r
     #PC1:
     load.pc1 <- p.log$rotation[,1][order(p.log$rotation[,1])]
     extreme.pc1 <- c(tail(load.pc1), head(load.pc1))

     extreme.pc1.ensg <- names(extreme.pc1)
     extreme.pc1.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc1.ensg,
                           mart=ensembl)
     
      q <- extreme.pc1.symbols[,2]
      names(q) <- extreme.pc1.symbols[,1]

      fpkm.log.pc1 <- cbind(q[extreme.pc1.ensg],published.nozero[extreme.pc1.ensg,])
      names(fpkm.log.pc1)[names(fpkm.log.pc1) == 'q[extreme.pc1.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc1, names.arg=q[extreme.pc1.ensg],las=2,main="Genes w highest absolute loadings in PC1 (log2 RPKM)")
```

![plot of chunk :published log PC 1-3](figure/:published log PC 1-31.png) 

```r
      print(fpkm.log.pc1)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000171560         FGA       0.0       0.0        0.7         0.00
## ENSG00000148677      ANKRD1    2451.5       0.2        0.3      1088.94
## ENSG00000135218        CD36     499.4       2.6        8.0       132.46
## ENSG00000057593          F7       0.2       0.3        0.3         0.00
## ENSG00000118194       TNNT2    2693.8       7.0       20.8      8074.39
## ENSG00000188257     PLA2G2A      51.2       0.8        2.3       210.68
## ENSG00000105372       RPS19     419.2     268.3      453.6       118.30
## ENSG00000063177       RPL18     455.9     237.8      458.0       205.67
## ENSG00000100097      LGALS1     498.7     113.2       77.2       389.82
## ENSG00000105640      RPL18A      45.8      31.3       68.1       505.42
## ENSG00000105583     WDR83OS     132.7      52.6       99.0        23.17
## ENSG00000087086         FTL     436.5     487.2     1523.1       463.86
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000171560         0.00  0.000e+00     0.0000      6.4728      63.271
## ENSG00000148677         0.10  1.471e+03     2.0740      0.5602    3863.610
## ENSG00000135218         0.08  1.551e+02     0.5914      0.3047     894.983
## ENSG00000057593         0.30  6.089e-03     0.1806      0.4763      26.652
## ENSG00000118194         9.37  2.028e+03     1.6755     11.8667   10930.800
## ENSG00000188257         0.00  2.314e+00     0.4622      0.4419    2351.218
## ENSG00000105372       140.77  1.181e+02   324.9812    306.6752       3.360
## ENSG00000063177       284.33  1.004e+02   159.0170    134.4763       3.984
## ENSG00000100097       161.11  1.838e+02   216.8468     58.2166       2.138
## ENSG00000105640       422.37  4.883e+01    55.6920     60.6585       0.000
## ENSG00000105583        26.35  6.529e+01    62.2844     80.6909       0.411
## ENSG00000087086       407.35  5.738e+02  2099.7219   5310.5093      16.930
##                 Atlas_brain Atlas_kidney
## ENSG00000171560       0.000       25.045
## ENSG00000148677       0.403        0.263
## ENSG00000135218       8.763       17.456
## ENSG00000057593      15.712       19.155
## ENSG00000118194       8.713       43.633
## ENSG00000188257       3.057       24.728
## ENSG00000105372       2.945        3.123
## ENSG00000063177       5.143        2.603
## ENSG00000100097       1.154        0.332
## ENSG00000105640       0.000        0.000
## ENSG00000105583       0.352        0.452
## ENSG00000087086       7.700       40.892
```

```r
      #PC2:
      load.pc2 <- p.log$rotation[,2][order(p.log$rotation[,2])]
      extreme.pc2 <- c(tail(load.pc2), head(load.pc2))
      
      extreme.pc2.ensg <- names(extreme.pc2)
      extreme.pc2.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc2.ensg,
                           mart=ensembl)

      q <- extreme.pc2.symbols[,2]
      names(q) <- extreme.pc2.symbols[,1]
      
      fpkm.log.pc2 <- cbind(q[extreme.pc2.ensg],published.nozero[extreme.pc2.ensg,])
      names(fpkm.log.pc2)[names(fpkm.log.pc2) == 'q[extreme.pc2.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc2, names.arg=q[extreme.pc2.ensg],las=2,main="Genes w highest absolute loadings in PC2 (log2 RPKM)")
```

![plot of chunk :published log PC 1-3](figure/:published log PC 1-32.png) 

```r
      print(fpkm.log.pc2)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000168314        MOBP       0.0     138.0        0.1         0.00
## ENSG00000104833      TUBB4A       0.3     235.1        1.5         2.94
## ENSG00000104435       STMN2       0.1     235.9        0.2         0.00
## ENSG00000132639      SNAP25       1.2     802.9        1.7         0.15
## ENSG00000123560        PLP1      11.9    1295.2        0.7         2.56
## ENSG00000131095        GFAP       3.0    1372.0        1.4         1.70
## ENSG00000111245        MYL2    5291.5       0.2        2.6     16780.36
## ENSG00000114854       TNNC1    3087.9       0.3       14.0      4100.03
## ENSG00000160808        MYL3    1213.4       3.2       24.9      3909.81
## ENSG00000198125          MB    3937.3       0.9        2.3      7065.54
## ENSG00000104879         CKM    1172.0       0.1        2.6      3321.37
## ENSG00000159251       ACTC1    2914.4       1.6        0.4      3457.43
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000168314       100.74  6.864e-03    486.669     0.02784       0.176
## ENSG00000104833      2038.44  5.334e-01    387.557     1.91307       0.408
## ENSG00000104435       176.88  1.790e-02    217.035     0.15558       0.105
## ENSG00000132639       535.77  3.758e-02    117.303     0.65681       0.346
## ENSG00000123560      1544.93  2.432e+00   1164.716     0.41628       5.086
## ENSG00000131095       834.13  1.671e+00   3633.817     0.85376       0.549
## ENSG00000111245         1.03  1.577e+04     10.201     2.31427    1064.500
## ENSG00000114854         0.69  2.618e+03      3.013     5.23761     316.135
## ENSG00000160808         0.84  3.997e+03      5.829    19.84928      80.668
## ENSG00000198125         1.35  2.706e+03      2.692     0.69784    1415.554
## ENSG00000104879         0.26  2.775e+03      4.498     1.39116      57.889
## ENSG00000159251         0.42  5.320e+03      3.653     1.04124     310.068
##                 Atlas_brain Atlas_kidney
## ENSG00000168314     153.127        0.139
## ENSG00000104833      82.655        0.389
## ENSG00000104435     139.464        0.000
## ENSG00000132639     284.680        0.180
## ENSG00000123560     779.755        0.513
## ENSG00000131095     691.908        0.135
## ENSG00000111245       0.000        0.358
## ENSG00000114854       0.763        1.014
## ENSG00000160808       0.204        2.634
## ENSG00000198125       0.032        0.000
## ENSG00000104879       0.000        0.060
## ENSG00000159251       0.421        0.311
```

```r
      #PC3:
      load.pc3 <- p.log$rotation[,3][order(p$rotation[,3])]
      extreme.pc3 <- c(tail(load.pc3), head(load.pc3))
      
      extreme.pc3.ensg <- names(extreme.pc3)
      extreme.pc3.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc3.ensg,
                           mart=ensembl)
      
      q <- extreme.pc3.symbols[,2]
      names(q) <- extreme.pc3.symbols[,1]
      
      fpkm.log.pc3 <- cbind(q[extreme.pc3.ensg],published.nozero[extreme.pc3.ensg,])
      names(fpkm.log.pc3)[names(fpkm.log.pc3) == 'q[extreme.pc3.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc3, names.arg=q[extreme.pc3.ensg],las=2,main="Genes w highest absolute loadings in PC3 (log2 RPKM)")
```

![plot of chunk :published log PC 1-3](figure/:published log PC 1-33.png) 

```r
      print(fpkm.log.pc3)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000111245        MYL2    5291.5       0.2        2.6      16780.4
## ENSG00000101608      MYL12A    2033.5      32.2      214.9        529.0
## ENSG00000175206        NPPA    6693.0       8.1        0.1        193.7
## ENSG00000175084         DES    3403.6       2.4       10.2       3160.1
## ENSG00000159251       ACTC1    2914.4       1.6        0.4       3457.4
## ENSG00000129991       TNNI3    2157.5       0.5        0.1       2698.4
## ENSG00000167996        FTH1     363.8     366.9      498.9       3563.3
## ENSG00000118194       TNNT2    2693.8       7.0       20.8       8074.4
## ENSG00000125971     DYNLRB1     207.8     164.7      170.0       1533.0
## ENSG00000075624        ACTB     423.6     786.5      521.8        377.9
## ENSG00000071082       RPL31     632.0     348.9      561.3       2239.9
## ENSG00000106211       HSPB1     521.6      37.7      196.6       3971.4
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000111245         1.03   15765.18    10.2006      2.3143     1064.50
## ENSG00000101608        28.70    1968.66    54.9994    141.2811      300.41
## ENSG00000175206         1.74     137.69     0.8982      2.7192      311.40
## ENSG00000175084         2.65    5524.18    13.1029      4.9042     1675.89
## ENSG00000159251         0.42    5319.64     3.6530      1.0412      310.07
## ENSG00000129991         0.47    5751.82     4.0958      0.6091      999.56
## ENSG00000167996      9086.14     337.28  1151.6777   1828.0887       43.50
## ENSG00000118194         9.37    2027.96     1.6755     11.8667    10930.80
## ENSG00000125971      2293.04      21.46    42.5226     32.5899       20.87
## ENSG00000075624      3277.02     283.55   782.8235    653.1090       67.06
## ENSG00000071082      1678.90      38.42    99.8055     51.5868      378.92
## ENSG00000106211       643.46     333.52    86.1187    201.8283       61.19
##                 Atlas_brain Atlas_kidney
## ENSG00000111245       0.000        0.358
## ENSG00000101608      15.476       24.862
## ENSG00000175206       0.445        0.131
## ENSG00000175084       4.500        6.252
## ENSG00000159251       0.421        0.311
## ENSG00000129991       0.675        0.132
## ENSG00000167996      53.459       66.455
## ENSG00000118194       8.713       43.633
## ENSG00000125971      42.315       23.413
## ENSG00000075624     140.560       69.191
## ENSG00000071082     484.379      259.057
## ENSG00000106211       6.808       17.175
```

To further validate the above results, indicating that tissue specificity appears mainly in PC 2 and 3, we will extract the 500 genes with highest loadings in each component and plot the corresponding published FPKM values in a heatmap (not shown in paper):


```r
     load.pc1 <- abs(p.log$rotation[,1])[order(abs(p.log$rotation[,1]),decreasing=TRUE)]
     top.pc1 <- names(load.pc1[1:500]) 

     load.pc2 <- abs(p.log$rotation[,2])[order(abs(p.log$rotation[,2]),decreasing=TRUE)]
     top.pc2 <- names(load.pc2[1:500])

     load.pc3 <- abs(p.log$rotation[,3])[order(abs(p.log$rotation[,3]),decreasing=TRUE)]
     top.pc3 <- names(load.pc3[1:500])

     pheatmap(cor(published[top.pc1,]),method="spearman")
```

![plot of chunk :published log top 500 loadings heatmap](figure/:published log top 500 loadings heatmap1.png) 

```r
     pheatmap(cor(published[top.pc2,]),method="spearman")
```

![plot of chunk :published log top 500 loadings heatmap](figure/:published log top 500 loadings heatmap2.png) 

```r
     pheatmap(cor(published[top.pc3,]),method="spearman")
```

![plot of chunk :published log top 500 loadings heatmap](figure/:published log top 500 loadings heatmap3.png) 

Try Anova on a "melted" expression matrix with logged values and some metadata, Figure 2e:


```r
n <- melt(published.log)
```

```
## Using  as id variables
```

```r
colnames(n) <- c("sample_ID","RPKM")
meta <- sampleinfo_published[,c("Study","Tissue","Preparation","NumberRaw","Numbermapped","Readtype")]
rownames(meta) <- colnames(published.log)
tissue <- rep(meta$Tissue, each=nrow(published.log))
study <- rep(meta$Study, each=nrow(published.log))
prep <- rep(meta$Preparation, each=nrow(published.log))
layout <- rep(meta$Readtype, each=nrow(published.log))
raw <- rep(meta$NumberRaw, each=nrow(published.log))
mapped <- rep(meta$Numbermapped, each=nrow(published.log))
data <- data.frame(n, tissue=tissue, study=study, prep=prep, layout=layout,nraw=raw,nmapped=mapped)
fit <- lm(RPKM ~ layout + prep + nraw + study + tissue, data=data)
b <- anova(fit)
maxval = 100

barplot(100*b$"Sum Sq"[1:5]/sum(b$"Sum Sq"[1:5]),names.arg=rownames(a[1:5,]),main="Anova, log2 published FPKM/RPKM values",ylim=c(0,maxval))
```

![plot of chunk :published log anova](figure/:published log anova.png) 

Combat analysis is performed on log2 values (n=13,323):


```r
meta <- data.frame(study=c(rep("HPA",3),rep("AltIso",2),rep("GTex",3),rep("Atlas",3)),tissue=c("Heart","Brain","Kidney","Heart","Brain","Heart","Brain","Kidney","Heart","Brain","Kidney"))

batch <- meta$study
design <- model.matrix(~1,data=meta)
combat <- ComBat(dat=published.log,batch=batch,mod=design,numCovs=NULL,par.prior=TRUE)
```

```
## Found 4 batches
## Found 0  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

Heatmap of Spearman correlations between published expression profiles after combat run (# genes = 13,323), Figure 3a:


```r
pheatmap(cor(combat, method="spearman")) 
```

![plot of chunk :published log combat heatmap spearman](figure/:published log combat heatmap spearman.png) 

PCA analysis of published FPKM values after combat run, Figure 3b:


```r
colors <- c("indianred", "dodgerblue", "forestgreen",
            "indianred", "dodgerblue",
            "indianred", "dodgerblue", "forestgreen", 
            "indianred", "dodgerblue", "forestgreen")

p.combat <- prcomp(t(combat))

shapes <- c(rep(15,3),rep(16,2),rep(17,3),rep(8,3))

plot(p.combat$x[,1],p.combat$x[,2],pch=shapes,col=colors,xlab=paste("PC1 54% of variance"),ylab=paste("PC2 38% of variance"),main="Published FPKM values \n COMBAT \n n=13,323")
legend("bottomright",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("topright",legend=c("HPA","AltIso","GTEx","Atlas"),col="black",pch=c(15,16,17,8),ncol=2)
```

![plot of chunk :published log combat PCA](figure/:published log combat PCA.png) 

We can plot all pairwise combinations of principal components 1 to 5. (not shown in paper):


```r
par(mfrow=c(4,4))
for (i in 1:6){
  for(j in 1:6){
    if (i<j){ 
		plot(p.combat$x[,i],p.combat$x[,j],pch=shapes,col=colors,xlab=paste("PC",i),ylab=paste("PC",j),main="Published FPKM values \n COMBAT \ n=13323")
		}
	}
}
```

![plot of chunk :published log combat PCA pairwise](figure/:published log combat PCA pairwise.png) 

Look a bit closer at PCs 1-3 in prcomp (not shown in paper):


```r
     #PC1:
     load.pc1 <- p.combat$rotation[,1][order(p.combat$rotation[,1])]
     extreme.pc1 <- c(tail(load.pc1), head(load.pc1))

     extreme.pc1.ensg <- names(extreme.pc1)
     extreme.pc1.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc1.ensg,
                           mart=ensembl)

      q <- extreme.pc1.symbols[,2]
      names(q) <- extreme.pc1.symbols[,1]

      fpkm.combat.pc1 <- cbind(q[extreme.pc1.ensg],published.nozero[extreme.pc1.ensg,])
      names(fpkm.combat.pc1)[names(fpkm.combat.pc1) == 'q[extreme.pc1.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc1, names.arg=q[extreme.pc1.ensg],las=2,main="Genes w highest absolute loadings in PC1 (ComBat log2 RPKM)")
```

![plot of chunk :published log combat PC 1-3](figure/:published log combat PC 1-31.png) 

```r
      print(fpkm.combat.pc1)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000104833      TUBB4A       0.3     235.1        1.5         2.94
## ENSG00000168314        MOBP       0.0     138.0        0.1         0.00
## ENSG00000104435       STMN2       0.1     235.9        0.2         0.00
## ENSG00000132639      SNAP25       1.2     802.9        1.7         0.15
## ENSG00000123560        PLP1      11.9    1295.2        0.7         2.56
## ENSG00000131095        GFAP       3.0    1372.0        1.4         1.70
## ENSG00000111245        MYL2    5291.5       0.2        2.6     16780.36
## ENSG00000198125          MB    3937.3       0.9        2.3      7065.54
## ENSG00000114854       TNNC1    3087.9       0.3       14.0      4100.03
## ENSG00000148677      ANKRD1    2451.5       0.2        0.3      1088.94
## ENSG00000118194       TNNT2    2693.8       7.0       20.8      8074.39
## ENSG00000129991       TNNI3    2157.5       0.5        0.1      2698.36
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000104833      2038.44  5.334e-01    387.557     1.91307       0.408
## ENSG00000168314       100.74  6.864e-03    486.669     0.02784       0.176
## ENSG00000104435       176.88  1.790e-02    217.035     0.15558       0.105
## ENSG00000132639       535.77  3.758e-02    117.303     0.65681       0.346
## ENSG00000123560      1544.93  2.432e+00   1164.716     0.41628       5.086
## ENSG00000131095       834.13  1.671e+00   3633.817     0.85376       0.549
## ENSG00000111245         1.03  1.577e+04     10.201     2.31427    1064.500
## ENSG00000198125         1.35  2.706e+03      2.692     0.69784    1415.554
## ENSG00000114854         0.69  2.618e+03      3.013     5.23761     316.135
## ENSG00000148677         0.10  1.471e+03      2.074     0.56016    3863.610
## ENSG00000118194         9.37  2.028e+03      1.675    11.86674   10930.800
## ENSG00000129991         0.47  5.752e+03      4.096     0.60909     999.559
##                 Atlas_brain Atlas_kidney
## ENSG00000104833      82.655        0.389
## ENSG00000168314     153.127        0.139
## ENSG00000104435     139.464        0.000
## ENSG00000132639     284.680        0.180
## ENSG00000123560     779.755        0.513
## ENSG00000131095     691.908        0.135
## ENSG00000111245       0.000        0.358
## ENSG00000198125       0.032        0.000
## ENSG00000114854       0.763        1.014
## ENSG00000148677       0.403        0.263
## ENSG00000118194       8.713       43.633
## ENSG00000129991       0.675        0.132
```

```r
      #PC2:
      load.pc2 <- p.combat$rotation[,2][order(p.combat$rotation[,2])]
      extreme.pc2 <- c(tail(load.pc2), head(load.pc2))

      extreme.pc2.ensg <- names(extreme.pc2)
      extreme.pc2.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc2.ensg,
                           mart=ensembl)
      q <- extreme.pc2.symbols[,2]
      names(q) <- extreme.pc2.symbols[,1]

      fpkm.combat.pc2 <- cbind(q[extreme.pc2.ensg],published.nozero[extreme.pc2.ensg,])
      names(fpkm.combat.pc2)[names(fpkm.combat.pc2) == 'q[extreme.pc2.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc2, names.arg=q[extreme.pc2.ensg],las=2,main="Genes w highest absolute loadings in PC2 (ComBat log2 RPKM)")
```

![plot of chunk :published log combat PC 1-3](figure/:published log combat PC 1-32.png) 

```r
      print(fpkm.combat.pc2)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000095932    C19orf77       0.1       2.2      424.6         0.00
## ENSG00000145692        BHMT       0.0       1.0      538.5         0.00
## ENSG00000124253        PCK1       0.5       0.7      489.6         0.93
## ENSG00000164825       DEFB1       0.1       0.6      662.2         0.84
## ENSG00000169344        UMOD       0.1       1.5     1425.7         0.00
## ENSG00000136872       ALDOB       0.4       2.4     1824.1         0.00
## ENSG00000129991       TNNI3    2157.5       0.5        0.1      2698.36
## ENSG00000092054        MYH7    1802.4       2.1        0.8      4117.17
## ENSG00000111245        MYL2    5291.5       0.2        2.6     16780.36
## ENSG00000148677      ANKRD1    2451.5       0.2        0.3      1088.94
## ENSG00000198125          MB    3937.3       0.9        2.3      7065.54
## ENSG00000159251       ACTC1    2914.4       1.6        0.4      3457.43
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000095932         3.06  0.000e+00    0.75338    237.5330       0.000
## ENSG00000145692         0.26  1.949e-02    0.04664    210.6384       0.052
## ENSG00000124253         0.00  1.512e-02    1.43248   1107.3146       0.923
## ENSG00000164825         0.22  1.320e-01    0.37882    726.4011       0.000
## ENSG00000169344         0.00  0.000e+00    0.00000    285.9097       0.000
## ENSG00000136872         0.49  1.630e-01    0.22985   1927.8312       1.452
## ENSG00000129991         0.47  5.752e+03    4.09582      0.6091     999.559
## ENSG00000092054         0.83  4.714e+03    6.48698      0.6887    2137.510
## ENSG00000111245         1.03  1.577e+04   10.20060      2.3143    1064.500
## ENSG00000148677         0.10  1.471e+03    2.07400      0.5602    3863.610
## ENSG00000198125         1.35  2.706e+03    2.69152      0.6978    1415.554
## ENSG00000159251         0.42  5.320e+03    3.65304      1.0412     310.068
##                 Atlas_brain Atlas_kidney
## ENSG00000095932       0.042      109.653
## ENSG00000145692       0.106      116.957
## ENSG00000124253       0.099      127.729
## ENSG00000164825       0.082      179.618
## ENSG00000169344       0.000      371.344
## ENSG00000136872       0.125      788.910
## ENSG00000129991       0.675        0.132
## ENSG00000092054       0.635        0.354
## ENSG00000111245       0.000        0.358
## ENSG00000148677       0.403        0.263
## ENSG00000198125       0.032        0.000
## ENSG00000159251       0.421        0.311
```

```r
      #PC3:
      load.pc3 <- p.combat$rotation[,3][order(p.combat$rotation[,3])]
      extreme.pc3 <- c(tail(load.pc3), head(load.pc3))

      extreme.pc3.ensg <- names(extreme.pc3)
      extreme.pc3.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc3.ensg,
                           mart=ensembl)

      q <- extreme.pc3.symbols[,2]
      names(q) <- extreme.pc3.symbols[,1]

      fpkm.combat.pc3 <- cbind(q[extreme.pc3.ensg],published.nozero[extreme.pc3.ensg,])
      names(fpkm.combat.pc3)[names(fpkm.combat.pc3) == 'q[extreme.pc3.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc3, names.arg=q[extreme.pc3.ensg],las=2,main="Genes w highest absolute loadings in PC3 (ComBat log2 RPKM)")
```

![plot of chunk :published log combat PC 1-3](figure/:published log combat PC 1-33.png) 

```r
      print(fpkm.combat.pc3)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000163631         ALB       0.8       0.8      327.6         0.00
## ENSG00000173432        SAA1       0.6       1.3        1.2         5.20
## ENSG00000161610        HCRT       0.0       0.0        0.1         0.00
## ENSG00000140575      IQGAP1      22.5      15.9       40.6         6.29
## ENSG00000188257     PLA2G2A      51.2       0.8        2.3       210.68
## ENSG00000183395        PMCH       0.2       0.0        0.0         0.13
## ENSG00000129824      RPS4Y1      40.8      53.7       60.9        66.40
## ENSG00000170345         FOS     147.3     116.7      288.4         6.23
## ENSG00000104888     SLC17A7       1.1     263.4        0.5         0.84
## ENSG00000103316        CRYM     124.5      58.4       43.4        30.53
## ENSG00000070808      CAMK2A       8.8     134.7        1.2         1.63
## ENSG00000170579      DLGAP1       0.4      57.5        0.8         0.00
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000163631         0.00  6.669e-02     0.6302    12.09975      37.193
## ENSG00000173432         0.00  2.748e-02     0.5260     9.07668     110.718
## ENSG00000161610         0.00  0.000e+00   165.0033     0.04222       0.000
## ENSG00000140575         4.86  1.251e+01   331.6092    14.63416      28.961
## ENSG00000188257         0.00  2.314e+00     0.4622     0.44191    2351.218
## ENSG00000183395         0.14  0.000e+00  2045.6382     0.00000       0.000
## ENSG00000129824        45.76  8.935e+01     0.6486    86.30981      30.176
## ENSG00000170345        20.94  2.698e+02     8.9939    29.33131      15.149
## ENSG00000104888       311.50  1.172e+00     0.5592     0.32666       0.044
## ENSG00000103316        38.56  6.342e+01     4.4874    48.14844       7.586
## ENSG00000070808       322.44  1.065e+00     3.2227     0.08715       1.640
## ENSG00000170579        56.90  4.434e-03     1.1074     0.42644       0.241
##                 Atlas_brain Atlas_kidney
## ENSG00000163631       0.152       24.347
## ENSG00000173432       0.000       35.264
## ENSG00000161610       0.197        0.000
## ENSG00000140575       8.308       23.509
## ENSG00000188257       3.057       24.728
## ENSG00000183395       0.050        0.000
## ENSG00000129824      42.832       14.097
## ENSG00000170345      41.809       34.855
## ENSG00000104888      11.004        0.019
## ENSG00000103316      20.478        5.674
## ENSG00000070808     131.757        0.144
## ENSG00000170579      53.706        0.905
```

Revisit Anova with combated values, Figure 3c:


```r
o <- melt(combat)
```

```
## Using  as id variables
```

```r
colnames(o) <- c("sample_ID","combat")
meta <- sampleinfo_published[,c("Study","Tissue","Preparation","NumberRaw","Numbermapped","Readtype")]
rownames(meta) <- colnames(combat)
tissue <- rep(meta$Tissue, each=nrow(combat))
study <- rep(meta$Study, each=nrow(combat))
prep <- rep(meta$Preparation, each=nrow(combat))
layout <- rep(meta$Readtype, each=nrow(combat))
raw <- rep(meta$NumberRaw, each=nrow(combat))
mapped <- rep(meta$Numbermapped, each=nrow(combat))
data <- data.frame(o, tissue=tissue, study=study, prep=prep, layout=layout,nraw=raw,nmapped=mapped)
fit <- lm(combat ~ layout + prep + nraw + study + tissue, data=data)
c <- anova(fit)
maxval = 100

barplot(100*c$"Sum Sq"[1:5]/sum(c$"Sum Sq"[1:5]),names.arg=rownames(c[1:5,]),main="Anova Combat",ylim=c(0,maxval))
```

![plot of chunk :published log combat anova](figure/:published log combat anova.png) 

Comparing FPKMs for FASTQ files reprocessed with TopHat and Cufflinks:


```r
cufflinks <- read.delim("fpkm_table_tophat.txt")

sampleinfo_cufflinks <- read.delim("sample_info_reprocessed.txt")
```

First, we will restrict the data set to only include protein coding genes using the ensembl based R package biomaRt:


```r
gene_ids <- as.vector(cufflinks[,1])

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl") #select the ensembl database

gene_type <- getBM(attributes=c("ensembl_gene_id", "gene_biotype"), 
                   filters = "ensembl_gene_id",
                   values=gene_ids,
                   mart=ensembl)

pc <- subset(gene_type[,1],gene_type[,2]=="protein_coding")

cufflinks_pc <- cufflinks[match(pc,cufflinks[,1]),]
```

Let's remove all lines where FPKM is close to zero in all samples before we proceed with this version of the data set:


```r
cufflinks_pc_nozero <- cufflinks_pc[-which(rowSums(cufflinks_pc[,3:16])<=0.01),]
```

Heatmap of Spearman correlations between reprocessed expression profiles (# genes = 19,475), Figure 4a:


```r
pheatmap(cor(cufflinks_pc_nozero[,3:16], method="spearman")) 
```

![plot of chunk :cufflinks heatmap spearman](figure/:cufflinks heatmap spearman.png) 

Let's look at a few PCA plots, Figure 4b:


```r
cufflinks_fpkms <- cufflinks_pc_nozero[,3:16]
rownames(cufflinks_fpkms) <- cufflinks_pc_nozero[,1]

p.cufflinks <- prcomp(t(cufflinks_fpkms))

colors <- c("dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred")          

shapes_cufflinks <- c(rep(11,3),rep(8,3),rep(17,3),rep(15,3),rep(16,2))

plot(p.cufflinks$x[,1],p.cufflinks$x[,2],pch=shapes_cufflinks,col=colors,xlab=paste("PC1 87% of variance"),ylab=paste("PC2 7.7% of variance"),main="Reprocessed FPKM values \n n=19,475")
legend("bottomright",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("top",legend=c("EoGE","Atlas","BodyMap","HPA","AltIso"),col="black",pch=c(11,8,17,15,16),ncol=2)
```

![plot of chunk :cufflinks PCA](figure/:cufflinks PCA.png) 

We can plot all pairwise combinations of principal components 1 to 5 (not shown in paper):


```r
colors <- c("dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred")

par(mfrow=c(4,4))
for (i in 1:6){
  for(j in 1:6){
    if (i<j){ 
      plot(p.cufflinks$x[,i],p.cufflinks$x[,j],pch=shapes_cufflinks,col=colors,xlab=paste("PC",i),ylab=paste("PC",j),main="Cufflinks FPKM values \n n=19475")
		}
	}
}
```

![plot of chunk :cufflinks PCA pairwise](figure/:cufflinks PCA pairwise.png) 

Look at PCA loadings for PC1-3 (not shown in paper):


```r
      #PC1:
      load.pc1 <- p.cufflinks$rotation[,1][order(p.cufflinks$rotation[,1])]
      extreme.pc1 <- c(tail(load.pc1), head(load.pc1))
      
      extreme.pc1.ensg <- names(extreme.pc1)
      extreme.pc1.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc1.ensg,
                           mart=ensembl)
      
      q <- extreme.pc1.symbols[,2]
      names(q) <- extreme.pc1.symbols[,1]
      
      fpkm.cuff.pc1 <- cbind(q[extreme.pc1.ensg],cufflinks_fpkms[extreme.pc1.ensg,])
      names(fpkm.cuff.pc1)[names(fpkm.cuff.pc1) == 'q[extreme.pc1.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc1, names.arg=q[extreme.pc1.ensg],las=2,main="Genes w highest absolute loadings in PC1 (raw Cufflinks FPKM)")
```

![plot of chunk :cufflinks PC 1-3](figure/:cufflinks PC 1-31.png) 

```r
      print(fpkm.cuff.pc1)
```

```
##                 Gene Symbol EoGE_brain EoGE_heart EoGE_kidney Atlas_brain
## ENSG00000198888      MT-ND1  9.783e+03  1.297e+04     12927.4    6609.070
## ENSG00000198804      MT-CO1  9.525e+03  3.874e+04      3772.8    4331.870
## ENSG00000198938      MT-CO3  1.038e+04  2.799e+04     13574.3    2072.480
## ENSG00000198899     MT-ATP6  1.454e+04  4.370e+04     27038.3    8168.860
## ENSG00000198840      MT-ND3  1.115e+04  1.127e+04     16370.8    4281.000
## ENSG00000228253     MT-ATP8  3.594e+04  1.255e+04     29657.8    5519.300
## ENSG00000060138        YBX3  9.072e+00  2.253e+02       243.2    2266.430
## ENSG00000211445        GPX3  3.092e+01  5.421e+02     18989.3      11.297
## ENSG00000166598     HSP90B1  1.021e+02  6.230e+01       197.9    1934.760
## ENSG00000130203        APOE  8.039e+02  6.579e+01      7062.5      97.184
## ENSG00000136872       ALDOB  2.326e-01  2.628e-01      3923.9       0.142
## ENSG00000109971       HSPA8  6.278e+02  4.847e+02       508.5    2205.090
##                 Atlas_heart Atlas_kidney BodyMap_brain BodyMap_heart
## ENSG00000198888   11043.300       9043.7     1.264e+04     2.880e+04
## ENSG00000198804   12564.700      12148.2     1.416e+04     4.935e+04
## ENSG00000198938    3861.350       4634.5     1.246e+04     3.890e+04
## ENSG00000198899   16937.600      16023.4     2.621e+04     4.914e+04
## ENSG00000198840   11015.100      11389.9     2.173e+04     3.456e+04
## ENSG00000228253   12604.100      11548.2     2.719e+05     5.177e+04
## ENSG00000060138   29731.300       2030.8     4.952e+01     3.201e+02
## ENSG00000211445     392.659       2265.2     2.872e+01     3.661e+02
## ENSG00000166598    5957.270       2593.6     1.054e+02     9.520e+01
## ENSG00000130203       9.908        131.7     2.369e+02     7.501e+00
## ENSG00000136872       1.949        839.9     4.263e-01     2.485e-01
## ENSG00000109971    2103.800       1277.0     4.202e+02     2.989e+02
##                 BodyMap_kidney HPA_brain HPA_heart HPA_kidney AltIso_brain
## ENSG00000198888        12167.3 3.654e+03 2.356e+04    6011.14    2.617e+04
## ENSG00000198804        30835.8 1.124e+04 6.494e+04   18088.50    1.694e+04
## ENSG00000198938        29411.4 7.538e+03 5.671e+04   15282.40    2.023e+04
## ENSG00000198899        33007.1 8.507e+03 5.368e+04   16901.40    3.770e+04
## ENSG00000198840        32774.3 5.552e+03 4.185e+04   12978.90    6.827e+04
## ENSG00000228253       172497.0 1.553e+04 1.430e+05   34584.40    1.938e+05
## ENSG00000060138          123.3 1.763e+01 2.088e+02      41.42    2.459e+01
## ENSG00000211445         3920.1 8.354e+01 7.295e+02    8203.31    2.393e+01
## ENSG00000166598          415.5 1.441e+02 1.008e+02     192.41    6.632e+01
## ENSG00000130203          318.5 7.025e+02 1.755e+01    1933.39    2.607e+02
## ENSG00000136872          492.2 2.343e-01 2.951e-01    3730.37    2.186e-01
## ENSG00000109971          605.5 6.475e+02 2.972e+02     716.93    2.934e+02
##                 AltIso_heart
## ENSG00000198888    3.015e+04
## ENSG00000198804    3.748e+04
## ENSG00000198938    3.172e+04
## ENSG00000198899    6.117e+04
## ENSG00000198840    1.231e+05
## ENSG00000228253    3.279e+05
## ENSG00000060138    3.083e+02
## ENSG00000211445    4.083e+02
## ENSG00000166598    6.668e+01
## ENSG00000130203    2.136e+01
## ENSG00000136872    1.616e-01
## ENSG00000109971    1.233e+02
```

```r
      #PC2:
      load.pc2 <- p.cufflinks$rotation[,2][order(p.cufflinks$rotation[,2])]
      extreme.pc2 <- c(tail(load.pc2), head(load.pc2))
      
      extreme.pc2.ensg <- names(extreme.pc2)
      extreme.pc2.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc2.ensg,
                           mart=ensembl)
      
      q <- extreme.pc2.symbols[,2]
      names(q) <- extreme.pc2.symbols[,1]
      
      fpkm.cuff.pc2 <- cbind(q[extreme.pc2.ensg],cufflinks_fpkms[extreme.pc2.ensg,])
      names(fpkm.cuff.pc2)[names(fpkm.cuff.pc2) == 'q[extreme.pc2.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc2, names.arg=q[extreme.pc2.ensg],las=2,main="Genes w highest absolute loadings in PC2 (raw Cufflinks FPKM)")
```

![plot of chunk :cufflinks PC 1-3](figure/:cufflinks PC 1-32.png) 

```r
      print(fpkm.cuff.pc2)
```

```
##                 Gene Symbol EoGE_brain EoGE_heart EoGE_kidney Atlas_brain
## ENSG00000123560        PLP1   2829.730      8.547       0.000       474.8
## ENSG00000131095        GFAP   1848.190      9.741       5.943       934.2
## ENSG00000060138        YBX3      9.072    225.311     243.229      2266.4
## ENSG00000120885         CLU   5495.280    275.906    2845.960      3538.5
## ENSG00000197971         MBP   6099.710      9.180      28.511      2304.7
## ENSG00000228253     MT-ATP8  35940.500  12550.600   29657.800      5519.3
## ENSG00000198804      MT-CO1   9524.530  38741.000    3772.820      4331.9
## ENSG00000198899     MT-ATP6  14543.900  43699.600   27038.300      8168.9
## ENSG00000198938      MT-CO3  10377.900  27986.400   13574.300      2072.5
## ENSG00000212907     MT-ND4L   9785.530  61755.400   17136.400      3655.7
## ENSG00000198840      MT-ND3  11154.600  11274.100   16370.800      4281.0
## ENSG00000198886      MT-ND4   9043.370  31881.500   17645.400     10485.3
##                 Atlas_heart Atlas_kidney BodyMap_brain BodyMap_heart
## ENSG00000123560   2.276e+00    3.284e-01       1993.67         3.180
## ENSG00000131095   9.644e-01    5.573e+00       2412.65         1.405
## ENSG00000060138   2.973e+04    2.031e+03         49.52       320.109
## ENSG00000120885   1.706e+02    3.301e+02       2412.68        94.740
## ENSG00000197971   1.240e+01    2.300e+01       3361.86         6.559
## ENSG00000228253   1.260e+04    1.155e+04     271870.00     51767.400
## ENSG00000198804   1.256e+04    1.215e+04      14160.20     49349.400
## ENSG00000198899   1.694e+04    1.602e+04      26213.00     49140.100
## ENSG00000198938   3.861e+03    4.634e+03      12464.40     38902.200
## ENSG00000212907   8.938e+03    1.239e+04      13244.60     22257.200
## ENSG00000198840   1.102e+04    1.139e+04      21727.90     34560.500
## ENSG00000198886   2.365e+04    2.227e+04      13637.50     30344.000
##                 BodyMap_kidney HPA_brain HPA_heart HPA_kidney AltIso_brain
## ENSG00000123560      2.661e-01   1789.87 1.369e+01  2.369e-01      2404.65
## ENSG00000131095      4.526e-01   2752.72 2.589e+00  5.729e-01      1472.64
## ENSG00000060138      1.233e+02     17.63 2.088e+02  4.142e+01        24.59
## ENSG00000120885      5.582e+02   5895.80 2.144e+02  7.291e+02      6333.52
## ENSG00000197971      1.286e+01   2097.86 1.046e+01  1.988e+01     15257.10
## ENSG00000228253      1.725e+05  15531.10 1.430e+05  3.458e+04    193776.00
## ENSG00000198804      3.084e+04  11235.40 6.494e+04  1.809e+04     16936.60
## ENSG00000198899      3.301e+04   8507.14 5.368e+04  1.690e+04     37697.30
## ENSG00000198938      2.941e+04   7538.32 5.671e+04  1.528e+04     20230.10
## ENSG00000212907      1.949e+04   7248.37 3.878e+04  1.292e+04     21268.10
## ENSG00000198840      3.277e+04   5552.03 4.185e+04  1.298e+04     68273.50
## ENSG00000198886      2.339e+04   7150.81 3.850e+04  1.045e+04     21114.90
##                 AltIso_heart
## ENSG00000123560         3.64
## ENSG00000131095        11.56
## ENSG00000060138       308.35
## ENSG00000120885       206.26
## ENSG00000197971        11.95
## ENSG00000228253    327884.00
## ENSG00000198804     37479.40
## ENSG00000198899     61171.10
## ENSG00000198938     31723.40
## ENSG00000212907     29519.50
## ENSG00000198840    123074.00
## ENSG00000198886     25640.50
```

```r
      #PC3:
      load.pc3 <- p.cufflinks$rotation[,3][order(p.cufflinks$rotation[,3])]
      extreme.pc3 <- c(tail(load.pc3), head(load.pc3))
      
      extreme.pc3.ensg <- names(extreme.pc3)
      extreme.pc3.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc3.ensg,
                           mart=ensembl)
      
      q <- extreme.pc3.symbols[,2]
      names(q) <- extreme.pc3.symbols[,1]
      
      fpkm.cuff.pc3 <- cbind(q[extreme.pc3.ensg],cufflinks_fpkms[extreme.pc3.ensg,])
      names(fpkm.cuff.pc3)[names(fpkm.cuff.pc3) == 'q[extreme.pc3.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc3, names.arg=q[extreme.pc3.ensg],las=2,main="Genes w highest absolute loadings in PC3 (raw Cufflinks FPKM)")
```

![plot of chunk :cufflinks PC 1-3](figure/:cufflinks PC 1-33.png) 

```r
      print(fpkm.cuff.pc3)
```

```
##                 Gene Symbol EoGE_brain EoGE_heart EoGE_kidney Atlas_brain
## ENSG00000060138        YBX3  9.072e+00      225.3   2.432e+02   2.266e+03
## ENSG00000111245        MYL2  0.000e+00     9422.8   1.395e-01   3.665e-01
## ENSG00000118194       TNNT2  3.716e+01     6316.9   1.101e+02   1.452e+01
## ENSG00000198888      MT-ND1  9.783e+03    12970.4   1.293e+04   6.609e+03
## ENSG00000198695      MT-ND6  3.821e+03     5483.3   6.604e+03   1.162e+04
## ENSG00000198840      MT-ND3  1.115e+04    11274.1   1.637e+04   4.281e+03
## ENSG00000198804      MT-CO1  9.525e+03    38741.0   3.773e+03   4.332e+03
## ENSG00000198938      MT-CO3  1.038e+04    27986.4   1.357e+04   2.072e+03
## ENSG00000228253     MT-ATP8  3.594e+04    12550.6   2.966e+04   5.519e+03
## ENSG00000212907     MT-ND4L  9.786e+03    61755.4   1.714e+04   3.656e+03
## ENSG00000175206        NPPA  2.841e+00      451.7   3.932e-01   6.998e-01
## ENSG00000106631        MYL7  6.932e-01     4082.4   2.477e-01   0.000e+00
##                 Atlas_heart Atlas_kidney BodyMap_brain BodyMap_heart
## ENSG00000060138     29731.3    2.031e+03     4.952e+01         320.1
## ENSG00000111245      2115.7    2.477e+00     5.078e-01        6847.0
## ENSG00000118194     18168.9    6.913e+01     1.142e+01        2903.3
## ENSG00000198888     11043.3    9.044e+03     1.264e+04       28802.2
## ENSG00000198695     17693.0    2.026e+04     8.165e+03       12493.1
## ENSG00000198840     11015.1    1.139e+04     2.173e+04       34560.5
## ENSG00000198804     12564.7    1.215e+04     1.416e+04       49349.4
## ENSG00000198938      3861.3    4.634e+03     1.246e+04       38902.2
## ENSG00000228253     12604.1    1.155e+04     2.719e+05       51767.4
## ENSG00000212907      8938.0    1.239e+04     1.324e+04       22257.2
## ENSG00000175206       425.6    1.603e-01     9.122e-01         262.6
## ENSG00000106631       317.6    0.000e+00     1.877e-01         462.3
##                 BodyMap_kidney HPA_brain HPA_heart HPA_kidney AltIso_brain
## ENSG00000060138      1.233e+02    17.635 2.088e+02      41.42    2.459e+01
## ENSG00000111245      0.000e+00     0.000 5.943e+00       0.00    1.765e+00
## ENSG00000118194      2.374e+01     6.302 3.838e+03      19.30    2.617e+01
## ENSG00000198888      1.217e+04  3653.900 2.356e+04    6011.14    2.617e+04
## ENSG00000198695      1.559e+04  4505.170 2.130e+04    5803.01    1.413e+04
## ENSG00000198840      3.277e+04  5552.030 4.185e+04   12978.90    6.827e+04
## ENSG00000198804      3.084e+04 11235.400 6.494e+04   18088.50    1.694e+04
## ENSG00000198938      2.941e+04  7538.320 5.671e+04   15282.40    2.023e+04
## ENSG00000228253      1.725e+05 15531.100 1.430e+05   34584.40    1.938e+05
## ENSG00000212907      1.949e+04  7248.370 3.878e+04   12920.00    2.127e+04
## ENSG00000175206      3.065e-01    11.069 2.553e+04       0.00    3.203e+00
## ENSG00000106631      0.000e+00     0.000 1.001e+04       0.00    0.000e+00
##                 AltIso_heart
## ENSG00000060138        308.3
## ENSG00000111245      15922.7
## ENSG00000118194      12212.1
## ENSG00000198888      30152.7
## ENSG00000198695      28571.2
## ENSG00000198840     123074.0
## ENSG00000198804      37479.4
## ENSG00000198938      31723.4
## ENSG00000228253     327884.0
## ENSG00000212907      29519.5
## ENSG00000175206        737.4
## ENSG00000106631        851.3
```
Mitochondrially encoded genes have relatively high expresseion levels, FPKM values of several thousands.

Anova analysis of different batch factors, Figure 4c:


```r
p <- melt(cufflinks_fpkms)
```

```
## Using  as id variables
```

```r
colnames(p) <- c("sample_ID","Cuff_FPKM")
meta <- sampleinfo_cufflinks[,c("Study","Tissue","Preparation","NumberRaw","Numbermapped","Readtype")]
rownames(meta) <- colnames(cufflinks_fpkms)
tissue <- rep(meta$Tissue, each=nrow(cufflinks_fpkms))
study <- rep(meta$Study, each=nrow(cufflinks_fpkms))
prep <- rep(meta$Preparation, each=nrow(cufflinks_fpkms))
layout <- rep(meta$Readtype, each=nrow(cufflinks_fpkms))
nraw <- rep(meta$NumberRaw, each=nrow(cufflinks_fpkms))
data <- data.frame(p, tissue=tissue, study=study, prep=prep, layout=layout, nraw=nraw)
fit <- lm(Cuff_FPKM ~ layout + prep + nraw + study + tissue, data=data)
d <- anova(fit)
maxval = 100

barplot(100*d$"Sum Sq"[1:5]/sum(d$"Sum Sq"[1:5]),names.arg=rownames(d[1:5,]),main="Anova, Cufflinks FPKM",ylim=c(0,maxval))
```

![plot of chunk :cufflinks anova](figure/:cufflinks anova.png) 

Try log2 transformation of the reprocessed FPKM values:


```r
pseudo <- 1
cufflinks_log <- log2(cufflinks_fpkms + pseudo)
```

Heatmap of Spearman correlations between log2 reprocessed cufflinks FPKM values, Figure 4d:


```r
pheatmap(cor(cufflinks_log) ,method="spearman")
```

![plot of chunk :cufflinks log heatmap spearman](figure/:cufflinks log heatmap spearman.png) 

PCA analysis of log2 reprocessed cufflinks FPKM values, Figure 4e:


```r
colors <- c("dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred")          

p.log.cufflinks <- prcomp(t(cufflinks_log))

shapes_cufflinks <- c(rep(11,3),rep(8,3),rep(17,3),rep(15,3),rep(16,2))
plot(p.log.cufflinks$x[,1],p.log.cufflinks$x[,2],pch=shapes_cufflinks,col=colors,xlab=paste("PC1 33% of variance"),ylab=paste("PC2 25% of variance"),main="log2 reprocessed cufflinks FPKM values \n n=19,475")
legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("top",legend=c("EoGE","Atlas","BodyMap","HPA","AltIso"),col="black",pch=c(11,8,17,15,16),ncol=2)
```

![plot of chunk :cufflinks log PCA 1&2](figure/:cufflinks log PCA 1&2.png) 

Figure 4f:


```r
plot(p.log.cufflinks$x[,2],p.log.cufflinks$x[,3],pch=shapes_cufflinks,col=colors,xlab=paste("PC2 25% of variance"),ylab=paste("PC3 22% of variance"),main="log2 reprocessed cufflinks FPKM values \n n=19,475")
legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("top",legend=c("EoGE","Atlas","BodyMap","HPA","AltIso"),col="black",pch=c(11,8,17,15,16),ncol=2)
```

![plot of chunk :cufflinks log PCA 2&3](figure/:cufflinks log PCA 2&3.png) 

We can plot all pairwise combinations of principal components 1 to 5. (not shown in paper):


```r
par(mfrow=c(4,4))
for (i in 1:6){
  for(j in 1:6){
    if (i<j){ 
  	plot(p.log.cufflinks$x[,i],p.log.cufflinks$x[,j],pch=shapes_cufflinks,col=colors,xlab=paste("PC",i),ylab=paste("PC",j),main="log2 reprocessed FPKM values \n n=19475")
		}
	}
}
```

![plot of chunk :cufflinks log PCA pairwise](figure/:cufflinks log PCA pairwise.png) 

Cross-validattion by leaving out one of the studies, in this case AltIso, Figure 4g: 


```r
colors <- c("dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred", "forestgreen",
            "dodgerblue", "indianred") 

p.loo <- cufflinks_log[,-c(13,14)]
colors.loo <- colors[-c(13,14)]
p.loo.log <- prcomp(t(p.loo))

p.add <- cufflinks_log[,c(13,14)]
projection <- t(p.add) %*% p.loo.log$rotation

p.original.plus.new <- rbind(p.loo.log$x, projection)
col.original.plus.new <- c(colors.loo, colors[c(13,14)])
plot(p.original.plus.new[,2],p.original.plus.new[,3],pch=c(shapes_cufflinks[1:12],rep(22,nrow(projection))),col=col.original.plus.new,xlab="PC2",ylab="PC3",main="log2 Cufflinks FPKM values; AltIso projected onto existing PCs \n n=19,475,",xlim=c(-150,100))

legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("topleft",legend=c("HPA","GTEx","Atlas","AltIso"),col="black",pch=c(11,8,17,15,22),ncol=2)
```

![plot of chunk :leave-one-out-pca cufflinks](figure/:leave-one-out-pca cufflinks.png) 

Look a bit closer at PCs 1-3 in prcomp for the logged FPKM values from cufflinks (not shown in paper):


```r
      #PC1
      load.pc1 <- p.log.cufflinks$rotation[,1][order(p.log.cufflinks$rotation[,1])]
      extreme.pc1 <- c(tail(load.pc1), head(load.pc1))
      
      extreme.pc1.ensg <- names(extreme.pc1)
      extreme.pc1.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc1.ensg,
                           mart=ensembl)
      
      q <- extreme.pc1.symbols[,2]
      names(q) <- extreme.pc1.symbols[,1]
      
      fpkm.logcuff.pc1 <- cbind(q[extreme.pc1.ensg],cufflinks_fpkms[extreme.pc1.ensg,])
      names(fpkm.logcuff.pc1)[names(fpkm.logcuff.pc1) == 'q[extreme.pc1.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc1, names.arg=q[extreme.pc1.ensg],las=2,main="Genes w highest absolute loadings in PC1 (log2 Cufflinks FPKM)")
```

![plot of chunk :cufflinks log PC 1-3](figure/:cufflinks log PC 1-31.png) 

```r
      print(fpkm.combat.pc1)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000104833      TUBB4A       0.3     235.1        1.5         2.94
## ENSG00000168314        MOBP       0.0     138.0        0.1         0.00
## ENSG00000104435       STMN2       0.1     235.9        0.2         0.00
## ENSG00000132639      SNAP25       1.2     802.9        1.7         0.15
## ENSG00000123560        PLP1      11.9    1295.2        0.7         2.56
## ENSG00000131095        GFAP       3.0    1372.0        1.4         1.70
## ENSG00000111245        MYL2    5291.5       0.2        2.6     16780.36
## ENSG00000198125          MB    3937.3       0.9        2.3      7065.54
## ENSG00000114854       TNNC1    3087.9       0.3       14.0      4100.03
## ENSG00000148677      ANKRD1    2451.5       0.2        0.3      1088.94
## ENSG00000118194       TNNT2    2693.8       7.0       20.8      8074.39
## ENSG00000129991       TNNI3    2157.5       0.5        0.1      2698.36
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000104833      2038.44  5.334e-01    387.557     1.91307       0.408
## ENSG00000168314       100.74  6.864e-03    486.669     0.02784       0.176
## ENSG00000104435       176.88  1.790e-02    217.035     0.15558       0.105
## ENSG00000132639       535.77  3.758e-02    117.303     0.65681       0.346
## ENSG00000123560      1544.93  2.432e+00   1164.716     0.41628       5.086
## ENSG00000131095       834.13  1.671e+00   3633.817     0.85376       0.549
## ENSG00000111245         1.03  1.577e+04     10.201     2.31427    1064.500
## ENSG00000198125         1.35  2.706e+03      2.692     0.69784    1415.554
## ENSG00000114854         0.69  2.618e+03      3.013     5.23761     316.135
## ENSG00000148677         0.10  1.471e+03      2.074     0.56016    3863.610
## ENSG00000118194         9.37  2.028e+03      1.675    11.86674   10930.800
## ENSG00000129991         0.47  5.752e+03      4.096     0.60909     999.559
##                 Atlas_brain Atlas_kidney
## ENSG00000104833      82.655        0.389
## ENSG00000168314     153.127        0.139
## ENSG00000104435     139.464        0.000
## ENSG00000132639     284.680        0.180
## ENSG00000123560     779.755        0.513
## ENSG00000131095     691.908        0.135
## ENSG00000111245       0.000        0.358
## ENSG00000198125       0.032        0.000
## ENSG00000114854       0.763        1.014
## ENSG00000148677       0.403        0.263
## ENSG00000118194       8.713       43.633
## ENSG00000129991       0.675        0.132
```

```r
      #PC2
      load.pc2 <- p.log.cufflinks$rotation[,2][order(p.log.cufflinks$rotation[,2])]
      extreme.pc2 <- c(tail(load.pc2), head(load.pc2))
      
      extreme.pc2.ensg <- names(extreme.pc2)
      extreme.pc2.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc2.ensg,
                           mart=ensembl)
      
      q <- extreme.pc2.symbols[,2]
      names(q) <- extreme.pc2.symbols[,1]
      
      fpkm.logcuff.pc2 <- cbind(q[extreme.pc2.ensg],cufflinks_fpkms[extreme.pc2.ensg,])
      names(fpkm.logcuff.pc2)[names(fpkm.logcuff.pc2) == 'q[extreme.pc2.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc2, names.arg=q[extreme.pc2.ensg],las=2,main="Genes w highest absolute loadings in PC2 (log2 Cufflinks FPKM)")
```

![plot of chunk :cufflinks log PC 1-3](figure/:cufflinks log PC 1-32.png) 

```r
      print(fpkm.combat.pc2)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000095932    C19orf77       0.1       2.2      424.6         0.00
## ENSG00000145692        BHMT       0.0       1.0      538.5         0.00
## ENSG00000124253        PCK1       0.5       0.7      489.6         0.93
## ENSG00000164825       DEFB1       0.1       0.6      662.2         0.84
## ENSG00000169344        UMOD       0.1       1.5     1425.7         0.00
## ENSG00000136872       ALDOB       0.4       2.4     1824.1         0.00
## ENSG00000129991       TNNI3    2157.5       0.5        0.1      2698.36
## ENSG00000092054        MYH7    1802.4       2.1        0.8      4117.17
## ENSG00000111245        MYL2    5291.5       0.2        2.6     16780.36
## ENSG00000148677      ANKRD1    2451.5       0.2        0.3      1088.94
## ENSG00000198125          MB    3937.3       0.9        2.3      7065.54
## ENSG00000159251       ACTC1    2914.4       1.6        0.4      3457.43
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000095932         3.06  0.000e+00    0.75338    237.5330       0.000
## ENSG00000145692         0.26  1.949e-02    0.04664    210.6384       0.052
## ENSG00000124253         0.00  1.512e-02    1.43248   1107.3146       0.923
## ENSG00000164825         0.22  1.320e-01    0.37882    726.4011       0.000
## ENSG00000169344         0.00  0.000e+00    0.00000    285.9097       0.000
## ENSG00000136872         0.49  1.630e-01    0.22985   1927.8312       1.452
## ENSG00000129991         0.47  5.752e+03    4.09582      0.6091     999.559
## ENSG00000092054         0.83  4.714e+03    6.48698      0.6887    2137.510
## ENSG00000111245         1.03  1.577e+04   10.20060      2.3143    1064.500
## ENSG00000148677         0.10  1.471e+03    2.07400      0.5602    3863.610
## ENSG00000198125         1.35  2.706e+03    2.69152      0.6978    1415.554
## ENSG00000159251         0.42  5.320e+03    3.65304      1.0412     310.068
##                 Atlas_brain Atlas_kidney
## ENSG00000095932       0.042      109.653
## ENSG00000145692       0.106      116.957
## ENSG00000124253       0.099      127.729
## ENSG00000164825       0.082      179.618
## ENSG00000169344       0.000      371.344
## ENSG00000136872       0.125      788.910
## ENSG00000129991       0.675        0.132
## ENSG00000092054       0.635        0.354
## ENSG00000111245       0.000        0.358
## ENSG00000148677       0.403        0.263
## ENSG00000198125       0.032        0.000
## ENSG00000159251       0.421        0.311
```

```r
      #PC3
      load.pc3 <- p.log.cufflinks$rotation[,3][order(p.log.cufflinks$rotation[,3])]
      extreme.pc3 <- c(tail(load.pc3), head(load.pc3))
      
      extreme.pc3.ensg <- names(extreme.pc3)
      extreme.pc3.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc3.ensg,
                           mart=ensembl)

      q <- extreme.pc3.symbols[,2]
      names(q) <- extreme.pc3.symbols[,1]
  
      fpkm.logcuff.pc3 <- cbind(q[extreme.pc3.ensg],cufflinks_fpkms[extreme.pc3.ensg,])
      names(fpkm.logcuff.pc3)[names(fpkm.logcuff.pc3) == 'q[extreme.pc3.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc3, names.arg=q[extreme.pc3.ensg],las=2,main="Genes w highest absolute loadings in PC3 (log2 Cufflinks FPKM)")
```

![plot of chunk :cufflinks log PC 1-3](figure/:cufflinks log PC 1-33.png) 

```r
      print(fpkm.combat.pc3)
```

```
##                 Gene Symbol HPA_heart HPA_brain HPA_kidney AltIso_heart
## ENSG00000163631         ALB       0.8       0.8      327.6         0.00
## ENSG00000173432        SAA1       0.6       1.3        1.2         5.20
## ENSG00000161610        HCRT       0.0       0.0        0.1         0.00
## ENSG00000140575      IQGAP1      22.5      15.9       40.6         6.29
## ENSG00000188257     PLA2G2A      51.2       0.8        2.3       210.68
## ENSG00000183395        PMCH       0.2       0.0        0.0         0.13
## ENSG00000129824      RPS4Y1      40.8      53.7       60.9        66.40
## ENSG00000170345         FOS     147.3     116.7      288.4         6.23
## ENSG00000104888     SLC17A7       1.1     263.4        0.5         0.84
## ENSG00000103316        CRYM     124.5      58.4       43.4        30.53
## ENSG00000070808      CAMK2A       8.8     134.7        1.2         1.63
## ENSG00000170579      DLGAP1       0.4      57.5        0.8         0.00
##                 AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart
## ENSG00000163631         0.00  6.669e-02     0.6302    12.09975      37.193
## ENSG00000173432         0.00  2.748e-02     0.5260     9.07668     110.718
## ENSG00000161610         0.00  0.000e+00   165.0033     0.04222       0.000
## ENSG00000140575         4.86  1.251e+01   331.6092    14.63416      28.961
## ENSG00000188257         0.00  2.314e+00     0.4622     0.44191    2351.218
## ENSG00000183395         0.14  0.000e+00  2045.6382     0.00000       0.000
## ENSG00000129824        45.76  8.935e+01     0.6486    86.30981      30.176
## ENSG00000170345        20.94  2.698e+02     8.9939    29.33131      15.149
## ENSG00000104888       311.50  1.172e+00     0.5592     0.32666       0.044
## ENSG00000103316        38.56  6.342e+01     4.4874    48.14844       7.586
## ENSG00000070808       322.44  1.065e+00     3.2227     0.08715       1.640
## ENSG00000170579        56.90  4.434e-03     1.1074     0.42644       0.241
##                 Atlas_brain Atlas_kidney
## ENSG00000163631       0.152       24.347
## ENSG00000173432       0.000       35.264
## ENSG00000161610       0.197        0.000
## ENSG00000140575       8.308       23.509
## ENSG00000188257       3.057       24.728
## ENSG00000183395       0.050        0.000
## ENSG00000129824      42.832       14.097
## ENSG00000170345      41.809       34.855
## ENSG00000104888      11.004        0.019
## ENSG00000103316      20.478        5.674
## ENSG00000070808     131.757        0.144
## ENSG00000170579      53.706        0.905
```

Seems to yield a heart vs. brain separation

To further validate the above results, indicating that tissue specificity appears mainly in PC 3, we will extract the 500 genes with highest loadings in each component and plot the corresponding cufflinks FPKM values in a heatmap (not shown in paper):


```r
     cufflinks_values <- cufflinks_pc_nozero[,3:16]
     rownames(cufflinks_values) <- cufflinks_pc_nozero[,1]

     load.pc1 <- abs(p.log.cufflinks$rotation[,1])[order(abs(p.log.cufflinks$rotation[,1]),decreasing=TRUE)]
     top.pc1 <- names(load.pc1[1:500])
     load.pc2 <- abs(p.log.cufflinks$rotation[,2])[order(abs(p.log.cufflinks$rotation[,2]),decreasing=TRUE)]
     top.pc2 <- names(load.pc2[1:500])
     load.pc3 <- abs(p.log.cufflinks$rotation[,3])[order(abs(p.log.cufflinks$rotation[,3]),decreasing=TRUE)]
     top.pc3 <- names(load.pc3[1:500])
     
     pheatmap(cor(cufflinks_values[top.pc1,]),method="spearman")
```

![plot of chunk :cufflinks log top 500 loadings heatmap](figure/:cufflinks log top 500 loadings heatmap1.png) 

```r
     pheatmap(cor(cufflinks_values[top.pc2,]),method="spearman")
```

![plot of chunk :cufflinks log top 500 loadings heatmap](figure/:cufflinks log top 500 loadings heatmap2.png) 

```r
     pheatmap(cor(cufflinks_values[top.pc3,]),method="spearman")
```

![plot of chunk :cufflinks log top 500 loadings heatmap](figure/:cufflinks log top 500 loadings heatmap3.png) 
    
Try Anova on a "melted" expression matrix with logged cufflinks values and some metadata (not shown in paper):


```r
q <- melt(cufflinks_log[,])
```

```
## Using  as id variables
```

```r
colnames(q) <- c("sample_ID","logFPKM")
meta <- sampleinfo_cufflinks[,c("Study","Tissue","Preparation","NumberRaw","Numbermapped","Readtype")]
rownames(meta) <- colnames(cufflinks_log)
tissue <- rep(meta$Tissue, each=nrow(cufflinks_log))
study <- rep(meta$Study, each=nrow(cufflinks_log))
prep <- rep(meta$Preparation, each=nrow(cufflinks_log))
layout <- rep(meta$Readtype, each=nrow(cufflinks_log))
data <- data.frame(q, tissue=tissue, study=study, prep=prep, layout=layout)
fit <- lm(logFPKM ~ layout + prep + nraw + study + tissue, data=data)
e <- anova(fit)
maxval = 100

barplot(100*e$"Sum Sq"[1:5]/sum(e$"Sum Sq"[1:5]),names.arg=rownames(e[1:5,]),main="Anova, Cufflinks log2 FPKM",ylim=c(0,maxval))
```

![plot of chunk :cufflinks log anova](figure/:cufflinks log anova.png) 

Combat analysis for removal of batch effects (n=19,475):


```r
meta <- data.frame(study=c(rep("EoGE",3),rep("Atlas",3),rep("BodyMap",3),rep("HPA",3),rep("AltIso",2)),tissue=c("Brain","Heart","Kidney","Brain","Heart","Kidney","Brain","Heart","Kidney","Brain","Heart","Kidney","Brain","Heart"),prep=c(rep("poly-A",3),rep("rRNA-depl",3),rep("poly-A",8)),layout=c(rep("PE",3),rep("SE",3),rep("PE",6),rep("SE",2)))

batch <- meta$study
design <- model.matrix(~1,data=meta)
combat.cufflinks <- ComBat(dat=cufflinks_log,batch=batch,mod=design,numCovs=NULL,par.prior=TRUE)
```

```
## Found 5 batches
## Found 0  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

```r
rownames(combat.cufflinks) <- rownames(cufflinks_log)
```

Heatmap of Spearman correlations between reprocessed cufflinks FPKM values after ComBat run, Figure 4h:


```r
pheatmap(cor(combat.cufflinks),method="spearman")
```

![plot of chunk :cufflinks log combat heatmap spearman](figure/:cufflinks log combat heatmap spearman.png) 

PCA analysis on reprocessed cufflinks FPKM values after ComBat run, Figure 4i:


```r
p.combat.cufflinks <- prcomp(t(combat.cufflinks))

shapes_cufflinks <- c(rep(11,3),rep(8,3),rep(17,3),rep(15,3),rep(16,2))


plot(p.combat.cufflinks$x[,1],p.combat.cufflinks$x[,2],pch=shapes_cufflinks,col=colors,xlab=paste("PC1 54% of variance"),ylab=paste("PC2 37% of variance"),main="Cufflinks FPKM values \n COMBAT \n n=19,475")
legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("top",legend=c("EoGE","Atlas","BodyMap","HPA","AltIso"),col="black",pch=c(11,8,17,15,16),ncol=2)
```

![plot of chunk :cufflinks log combat PCA](figure/:cufflinks log combat PCA.png) 

We can plot all pairwise combinations of principal components 1 to 5. (not shown in paper):


```r
par(mfrow=c(4,4))
for (i in 1:6){
  for(j in 1:6){
    if (i<j){ 
  	plot(p.combat.cufflinks$x[,i],p.combat.cufflinks$x[,j],pch=shapes_cufflinks,col=colors,xlab=paste("PC",i),ylab=paste("PC",j),main="Cufflinks FPKM values \n COMBAT \ n=19,475")
		}
	}
}
```

![plot of chunk :cufflinks log combat PCA pairwise](figure/:cufflinks log combat PCA pairwise.png) 

Look a bit closer at PCs 1-3 in prcomp (not shown in paper):


```r
      #PC1:
      load.pc1 <- p.combat.cufflinks$rotation[,1][order(p.combat.cufflinks$rotation[,1])]
      extreme.pc1 <- c(tail(load.pc1), head(load.pc1))
      
      extreme.pc1.ensg <- names(extreme.pc1)
      extreme.pc1.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc1.ensg,
                           mart=ensembl)
      
      q <- extreme.pc1.symbols[,2]
      names(q) <- extreme.pc1.symbols[,1]

      fpkm.combatcuff.pc1 <- cbind(q[extreme.pc1.ensg],cufflinks_fpkms[extreme.pc1.ensg,])
      names(fpkm.combatcuff.pc1)[names(fpkm.combatcuff.pc1) == 'q[extreme.pc1.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc1, names.arg=q[extreme.pc1.ensg],las=2,main="Genes w highest absolute loadings in PC1 (COMBAT Cufflinks FPKM)")
```

![plot of chunk :cufflinks log combat PC 1-3](figure/:cufflinks log combat PC 1-31.png) 

```r
      print(fpkm.combatcuff.pc1)
```

```
##                 Gene Symbol EoGE_brain EoGE_heart EoGE_kidney Atlas_brain
## ENSG00000106631        MYL7     0.6932  4.082e+03      0.2477   0.000e+00
## ENSG00000111245        MYL2     0.0000  9.423e+03      0.1395   3.665e-01
## ENSG00000175084         DES     1.7880  1.305e+04      4.6120   5.345e+00
## ENSG00000159251       ACTC1     0.7814  8.115e+03      2.4721   4.595e-01
## ENSG00000114854       TNNC1     0.9556  5.476e+03     37.9980   2.029e+00
## ENSG00000198125          MB     0.4224  6.147e+03      1.0122   8.953e-02
## ENSG00000132639      SNAP25  1281.4400  3.197e-01      0.4139   1.057e+03
## ENSG00000123560        PLP1  2829.7300  8.547e+00      0.0000   4.748e+02
## ENSG00000131095        GFAP  1848.1900  9.741e+00      5.9433   9.342e+02
## ENSG00000197971         MBP  6099.7100  9.180e+00     28.5108   2.305e+03
## ENSG00000104435       STMN2   370.6450  1.170e-01      0.0000   2.420e+02
## ENSG00000125462     C1orf61   267.1770  0.000e+00      0.0000   2.315e+01
##                 Atlas_heart Atlas_kidney BodyMap_brain BodyMap_heart
## ENSG00000106631    317.5820       0.0000        0.1877     4.623e+02
## ENSG00000111245   2115.7100       2.4769        0.5078     6.847e+03
## ENSG00000175084   1877.9900       6.6131        7.8118     3.316e+03
## ENSG00000159251    306.7440       0.3626        0.4032     2.516e+03
## ENSG00000114854    760.2540       2.1117        0.5150     1.921e+03
## ENSG00000198125    645.0390       0.0000        0.2452     4.562e+03
## ENSG00000132639      1.7793       0.1911      989.4640     1.474e-01
## ENSG00000123560      2.2760       0.3284     1993.6700     3.180e+00
## ENSG00000131095      0.9644       5.5727     2412.6500     1.405e+00
## ENSG00000197971     12.3986      23.0044     3361.8600     6.559e+00
## ENSG00000104435      0.2339       0.0000      180.8680     5.872e-03
## ENSG00000125462      0.2242       0.4376      184.8080     5.941e-03
##                 BodyMap_kidney HPA_brain HPA_heart HPA_kidney AltIso_brain
## ENSG00000106631         0.0000    0.0000 1.001e+04    0.00000    0.000e+00
## ENSG00000111245         0.0000    0.0000 5.943e+00    0.00000    1.765e+00
## ENSG00000175084        28.8400    0.3434 2.853e+03   15.32580    1.244e+01
## ENSG00000159251         4.9243    0.2034 4.205e+03    0.45761    1.214e+00
## ENSG00000114854         6.2579    0.0000 1.663e+03    8.55966    5.488e-01
## ENSG00000198125         0.2971    0.7570 6.001e+03    0.61315    8.099e-01
## ENSG00000132639         0.6910  383.4000 2.514e+00    0.92989    5.958e+02
## ENSG00000123560         0.2661 1789.8700 1.369e+01    0.23695    2.405e+03
## ENSG00000131095         0.4526 2752.7200 2.589e+00    0.57293    1.473e+03
## ENSG00000197971        12.8627 2097.8600 1.046e+01   19.88150    1.526e+04
## ENSG00000104435         0.0000   91.1251 3.292e-02    0.03517    1.676e+02
## ENSG00000125462         0.0310 1182.5300 0.000e+00    0.00000    2.999e+02
##                 AltIso_heart
## ENSG00000106631    8.513e+02
## ENSG00000111245    1.592e+04
## ENSG00000175084    1.120e+04
## ENSG00000159251    1.726e+03
## ENSG00000114854    4.131e+03
## ENSG00000198125    7.280e+03
## ENSG00000132639    9.816e-02
## ENSG00000123560    3.640e+00
## ENSG00000131095    1.156e+01
## ENSG00000197971    1.195e+01
## ENSG00000104435    2.997e-02
## ENSG00000125462    1.383e-01
```

```r
      #PC2:
      load.pc2 <- p.combat.cufflinks$rotation[,2][order(p.combat.cufflinks$rotation[,2])]
      extreme.pc2 <- c(tail(load.pc2), head(load.pc2))
      
      extreme.pc2.ensg <- names(extreme.pc2)
      extreme.pc2.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc2.ensg,
                           mart=ensembl)
      
      q <- extreme.pc2.symbols[,2]
      names(q) <- extreme.pc2.symbols[,1]
      
      fpkm.combatcuff.pc2 <- cbind(q[extreme.pc2.ensg],cufflinks_fpkms[extreme.pc2.ensg,])
      names(fpkm.combatcuff.pc2)[names(fpkm.combatcuff.pc2) == 'q[extreme.pc2.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc2, names.arg=q[extreme.pc2.ensg],las=2,main="Genes w highest absolute loadings in PC2 (COMBAT Cufflinks FPKM)")
```

![plot of chunk :cufflinks log combat PC 1-3](figure/:cufflinks log combat PC 1-32.png) 

```r
      print(fpkm.combatcuff.pc2)
```

```
##                 Gene Symbol EoGE_brain EoGE_heart EoGE_kidney Atlas_brain
## ENSG00000095932    C19orf77    6.95853     0.0000   1019.6500     0.06484
## ENSG00000164825       DEFB1    0.00000     0.0000   1469.9200     0.14473
## ENSG00000162366    PDZK1IP1    0.00000     0.3384   4259.2200     0.00000
## ENSG00000137731       FXYD2    0.08533     1.0573   2744.2100     0.00000
## ENSG00000136872       ALDOB    0.23265     0.2628   3923.8700     0.14197
## ENSG00000169344        UMOD    0.00000     0.0000   1380.6500     0.00000
## ENSG00000198125          MB    0.42237  6146.6400      1.0122     0.08953
## ENSG00000092054        MYH7    1.11922  7826.2300      0.1244     0.95501
## ENSG00000129991       TNNI3    1.06130  3832.3200      0.0000     1.88834
## ENSG00000106631        MYL7    0.69323  4082.3600      0.2477     0.00000
## ENSG00000175206        NPPA    2.84121   451.6720      0.3932     0.69976
## ENSG00000148677      ANKRD1    0.00000   400.9800      0.0000     0.47166
##                 Atlas_heart Atlas_kidney BodyMap_brain BodyMap_heart
## ENSG00000095932   0.000e+00     237.9050       5.01882     0.000e+00
## ENSG00000164825   0.000e+00     300.8200       0.44295     9.557e-01
## ENSG00000162366   1.766e-01      34.1620       0.22068     6.883e-02
## ENSG00000137731   3.148e-02      83.4440       0.07880     2.428e-01
## ENSG00000136872   1.949e+00     839.9170       0.42626     2.485e-01
## ENSG00000169344   0.000e+00     701.1950       0.00000     0.000e+00
## ENSG00000198125   6.450e+02       0.0000       0.24523     4.562e+03
## ENSG00000092054   3.061e+03       0.4458       0.45322     2.772e+03
## ENSG00000129991   3.112e+03       0.4260       0.42347     1.481e+03
## ENSG00000106631   3.176e+02       0.0000       0.18768     4.623e+02
## ENSG00000175206   4.256e+02       0.1603       0.91225     2.626e+02
## ENSG00000148677   4.320e+03       0.4388       0.04642     7.954e+02
##                 BodyMap_kidney HPA_brain HPA_heart HPA_kidney AltIso_brain
## ENSG00000095932      1.781e+02   0.95308 0.000e+00  6.261e+02      5.21370
## ENSG00000164825      2.895e+02   0.00000 0.000e+00  9.020e+02      0.25451
## ENSG00000162366      5.258e+02   0.00000 2.178e-01  1.305e+03      0.00000
## ENSG00000137731      1.244e+03   0.24916 7.389e-01  1.503e+03      0.01764
## ENSG00000136872      4.922e+02   0.23425 2.951e-01  3.730e+03      0.21865
## ENSG00000169344      1.357e+03   0.00000 0.000e+00  2.096e+03      0.00000
## ENSG00000198125      2.971e-01   0.75697 6.001e+03  6.131e-01      0.80990
## ENSG00000092054      5.451e-02   4.28417 3.286e+02  7.335e-02      0.64258
## ENSG00000129991      0.000e+00   0.13866 1.104e+03  0.000e+00      1.06883
## ENSG00000106631      0.000e+00   0.00000 1.001e+04  0.000e+00      0.00000
## ENSG00000175206      3.065e-01  11.06860 2.553e+04  0.000e+00      3.20302
## ENSG00000148677      2.788e-01   0.04946 1.858e+03  7.979e-02      0.31433
##                 AltIso_heart
## ENSG00000095932       0.5796
## ENSG00000164825       1.0971
## ENSG00000162366       0.2821
## ENSG00000137731       0.1487
## ENSG00000136872       0.1616
## ENSG00000169344       0.0000
## ENSG00000198125    7280.0900
## ENSG00000092054    2699.2100
## ENSG00000129991    3405.6700
## ENSG00000106631     851.3080
## ENSG00000175206     737.3900
## ENSG00000148677     793.0480
```

```r
      #PC3:      
      load.pc3 <- p.combat.cufflinks$rotation[,3][order(p.combat.cufflinks$rotation[,3])]
      extreme.pc3 <- c(tail(load.pc3), head(load.pc3))
      
      extreme.pc3.ensg <- names(extreme.pc3)
      extreme.pc3.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                           filters = "ensembl_gene_id",
                           values=extreme.pc3.ensg,
                           mart=ensembl)
      
      q <- extreme.pc3.symbols[,2]
      names(q) <- extreme.pc3.symbols[,1]

      fpkm.combatcuff.pc3 <- cbind(q[extreme.pc3.ensg],cufflinks_fpkms[extreme.pc3.ensg,])
      names(fpkm.combatcuff.pc3)[names(fpkm.combatcuff.pc3) == 'q[extreme.pc3.ensg]'] <- 'Gene Symbol'
      
      barplot(extreme.pc3, names.arg=q[extreme.pc3.ensg],las=2,main="Genes w highest absolute loadings in PC3 (COMBAT Cufflinks FPKM)")
```

![plot of chunk :cufflinks log combat PC 1-3](figure/:cufflinks log combat PC 1-33.png) 

```r
      print(fpkm.combatcuff.pc3)
```

```
##                 Gene Symbol EoGE_brain EoGE_heart EoGE_kidney Atlas_brain
## ENSG00000198838        RYR3     16.211    0.01536      0.5697     68.2328
## ENSG00000110195       FOLR1      0.157    0.55841    243.4410     10.6218
## ENSG00000223865    HLA-DPB1     15.340   74.86100    132.6340     18.6981
## ENSG00000163220      S100A9      0.629   14.06570     65.3565     67.6230
## ENSG00000268396                  0.000    0.00000      0.0000    303.2290
## ENSG00000118271         TTR      1.781    0.62967      6.3687   4138.2800
## ENSG00000266956                  0.000   40.50640     22.4281      0.0000
## ENSG00000257017          HP      2.282   23.74320     13.0780      3.5407
## ENSG00000175445         LPL      2.607  361.14400     11.9540      4.4369
## ENSG00000005513        SOX8     53.742    1.20108      0.1120      0.7289
## ENSG00000125676       THOC2     10.123    3.51788     12.8518      0.0000
## ENSG00000269027                  0.000    0.00000      0.0000      0.0000
##                 Atlas_heart Atlas_kidney BodyMap_brain BodyMap_heart
## ENSG00000198838   1.305e-01       0.7971        4.1883       0.03303
## ENSG00000110195   1.154e-01      29.8561        8.4481       0.22343
## ENSG00000223865   1.341e+01       9.9471       36.3018      17.66610
## ENSG00000163220   1.949e+02       7.5356       23.3908      44.53170
## ENSG00000268396   0.000e+00       0.0000        0.0000       0.00000
## ENSG00000118271   4.042e+00       2.9596     1200.0000       0.07613
## ENSG00000266956   1.758e+02     119.3440        0.1838       2.88450
## ENSG00000257017   2.590e+03     921.4000        0.8874      12.62420
## ENSG00000175445   5.376e+02      14.7952        6.9406     145.16000
## ENSG00000005513   1.509e-02       0.0000       14.7852       0.16574
## ENSG00000125676   4.852e+01      52.4633       14.2816       8.78908
## ENSG00000269027   0.000e+00     471.7360        0.0000       0.00000
##                 BodyMap_kidney HPA_brain HPA_heart HPA_kidney AltIso_brain
## ENSG00000198838         0.2603   1.09129   0.46409    0.39848      38.0721
## ENSG00000110195        44.1110   0.35560   0.09313  130.67000       0.4714
## ENSG00000223865        41.0103   3.99889  49.81610   53.94210      18.3877
## ENSG00000163220        42.3277   1.54815  17.65110    4.91085      14.2365
## ENSG00000268396         0.0000   0.00000   0.00000    0.00000       0.0000
## ENSG00000118271         1.3065   0.00000   0.36730    2.10049      20.5621
## ENSG00000266956        83.2357   0.01193   0.55508    0.11698       0.1413
## ENSG00000257017        26.6850   0.17245   5.04010    0.00000       0.7159
## ENSG00000175445         5.4742  75.32840  18.72190   10.35480       2.4762
## ENSG00000005513         0.1370 119.15600   0.62632    0.08259      26.2390
## ENSG00000125676        19.1089  24.73970  11.04990   17.04860       8.8476
## ENSG00000269027         0.0000   0.00000   0.00000    0.00000       0.0000
##                 AltIso_heart
## ENSG00000198838       0.3282
## ENSG00000110195       0.4818
## ENSG00000223865      49.6954
## ENSG00000163220      74.1273
## ENSG00000268396       0.0000
## ENSG00000118271       0.3213
## ENSG00000266956     153.1870
## ENSG00000257017      16.9893
## ENSG00000175445      88.1243
## ENSG00000005513       0.1488
## ENSG00000125676       6.3458
## ENSG00000269027       0.0000
```

Revisit ANOVA analysis on reprocessed Cufflinks FPKM values after ComBat run, Figure 4j:


```r
q <- melt(combat.cufflinks)
```

```
## Using  as id variables
```

```r
colnames(q) <- c("sample_ID","combat")
meta <- sampleinfo_cufflinks[,c("Study","Tissue","Preparation","NumberRaw","Numbermapped","Readtype")]
rownames(meta) <- colnames(combat.cufflinks)
tissue <- rep(meta$Tissue, each=nrow(combat.cufflinks))
study <- rep(meta$Study, each=nrow(combat.cufflinks))
prep <- rep(meta$Preparation, each=nrow(combat.cufflinks))
layout <- rep(meta$Readtype, each=nrow(combat.cufflinks))
raw <- rep(meta$NumberRaw, each=nrow(combat.cufflinks))
mapped <- rep(meta$Numbermapped, each=nrow(combat.cufflinks))
data <- data.frame(q, tissue=tissue, study=study, prep=prep, layout=layout,nraw=raw,nmapped=mapped)
fit <- lm(combat ~ prep + nraw + layout + study + tissue, data=data)
f <- anova(fit)
maxval = 100
barplot(100*f$"Sum Sq"[1:5]/sum(f$"Sum Sq"[1:5]),names.arg=rownames(f[1:5,]),main="Anova Cufflinks Combat",ylim=c(0,maxval))
```

![plot of chunk :cufflinks log combat anova](figure/:cufflinks log combat anova.png) 
