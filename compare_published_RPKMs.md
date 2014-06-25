Downloading the F/RPKM data
---------------------------

Prepare by defining functions etc.

```r
library(edgeR)
```

```
## Loading required package: limma
```

```r
library(pheatmap)
library(gplots)
```

```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
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

normalize.voom <- function(counts) {
    require(limma)
    return(voom(counts)$E)
}

cpm.tmm <- function(counts, groups = NA) {
    require(edgeR)
    if (is.na(groups)) {
        d <- DGEList(counts = counts)
    } else {
        d <- DGEList(counts = counts, group = groups)
    }
    d <- calcNormFactors(d, method = "TMM")
    return(cpm(d, normalized.lib.sizes = TRUE))
}

log2.cpm <- function(x) {
    x.pseud <- x + 0.5
    libsizes = colSums(x)
    libsize.pseud <- libsizes + 1
    new.mtx <- log2(1e+06 * mapply("/", as.data.frame(x.pseud), libsize.pseud))
    rownames(new.mtx) <- rownames(x)
    return(as.matrix(new.mtx))
}

do.SVD = function(m, comp.1 = 1, comp.2 = 2) {
    # returns eig.cell
    s <- svd(m)
    ev <- s$d^2/sum(s$d^2)
    return(s$u[, c(comp.1, comp.2)])
}

project.SVD <- function(m, eig.cell) {
    return(t(m) %*% eig.cell)
}

plot.SVD <- function(m, comp.1 = 1, comp.2 = 2, groups = rep("blue", ncol(m)), 
    title = "") {
    eig <- do.SVD(m, comp.1, comp.2)
    proj <- project.SVD(m, eig)
    xminv <- min(proj[, 1])  # - .2 * abs(min(proj[,1]))
    xmaxv <- max(proj[, 1])  # + .2 * abs(max(proj[,1]))
    yminv <- min(proj[, 2])  # - .2 * abs(min(proj[,2]))
    ymaxv <- max(proj[, 2])  # + .2 * abs(max(proj[,2]))
    plot(proj, pch = 20, col = "white", xlim = c(xminv, xmaxv), ylim = c(yminv, 
        ymaxv), xaxt = "n", yaxt = "n", xlab = "PC1", ylab = "PC2", main = title)
    
    points(proj, col = as.character(groups), pch = 20)  # , #pch=c(rep(15,3),rep(17,3),rep(19,3),rep(18,3),rep(20,2)), cex=2)
    textxy(proj[, 1], proj[, 2], labs = colnames(m))
}

loadings.SVD <- function(m, comp = 1, gene.ids = rownames(m)) {
    s <- svd(m)
    l <- s$u[, comp]
    names(l) <- gene.ids
    l.s <- l[order(l)]
    return(l.s)
}

plot.loadings.SVD <- function(m, comp = 1, cutoff = 0.1, gene.ids = rownames(m)) {
    l <- loadings.SVD(m, comp, gene.ids)
    barplot(l[abs(l) > cutoff], las = 2, main = paste("PC", comp, "cutoff", 
        cutoff), cex.names = 0.6)
}

plotPC <- function(matrix, a, b, desc, colors) {
    eig <- do.SVD(matrix, a, b)
    proj <- project.SVD(matrix, eig)
    xminv <- min(proj[, 1]) - 0.2 * abs(min(proj[, 1]))
    xmaxv <- max(proj[, 1]) + 0.2 * abs(max(proj[, 1]))
    yminv <- min(proj[, 2]) - 0.2 * abs(min(proj[, 2]))
    ymaxv <- max(proj[, 2]) + 0.2 * abs(max(proj[, 2]))
    plot(proj, pch = 20, xlim = c(xminv, xmaxv), ylim = c(yminv, ymaxv), xaxt = "n", 
        yaxt = "n", xlab = paste0("PC", a), ylab = paste("PC", b), col = colors, 
        main = desc)
    textxy(proj[, 1], proj[, 2], labs = rownames(proj))
}
```


Here, we download data from various public sources and extract the brain, heart and kidney samples.

"HPA": Human Protein Atlas

```r
# temp <- tempfile()
# download.file(url='http://www.proteinatlas.org/download/rna.csv.zip',destfile=temp)
# hpa <- read.csv(unz(temp, 'rna.csv')) unlink(temp)

# hpa.heart <- hpa[hpa$Sample=='heart muscle', c('Gene', 'Value')] hpa.brain
# <- hpa[hpa$Sample=='cerebral cortex', c('Gene', 'Value')] hpa.kidney <-
# hpa[hpa$Sample=='kidney', c('Gene', 'Value')]

# hpa.fpkms <- merge(hpa.heart, hpa.brain, by='Gene') hpa.fpkms <-
# merge(hpa.fpkms, hpa.kidney, by='Gene') colnames(hpa.fpkms) <-
# c('ENSG_ID', 'HPA_heart', 'HPA_brain', 'HPA_kidney')
```


Check if the identifiers are unique and write table to file.

```r
# length(hpa.fpkms[,1]) length(unique(hpa.fpkms[,1]))

# write.table(hpa.fpkms,file='hpa_fpkms.txt',quote=F,sep='\t')
```


"Altiso": Alternative isoform regulation in human tissue transcriptomes

```r
# temp <- tempfile()
# download.file(url='http://genes.mit.edu/burgelab/Supplementary/wang_sandberg08/hg18.ensGene.CEs.rpkm.txt',destfile=temp)
# altiso <- read.delim(temp, sep='\t') unlink(temp)
```


There is no kidney sample here, so just use heart + brain


```r
# altiso.fpkms <- altiso[,c('X.Gene','heart','brain')]
# colnames(altiso.fpkms) <- c('ENSG_ID', 'AltIso_heart', 'AltIso_brain')
```


Check uniqueness of IDs.


```r
# length(altiso.fpkms[,1]) length(unique(altiso.fpkms[,1]))

# write.table(altiso.fpkms,file='altiso_fpkms.txt',quote=F,sep='\t')
```


"GTEx": Genotype-Tissue Expression

This is a big download: 337.8 Mb (as of 2014-02-04)
We also add some code to randomly select one sample from each tissue type; there are many biological replicates in this data set.


```r
# temp <- tempfile()
# download.file(url='http://www.broadinstitute.org/gtex/rest/file/download?portalFileId=119363&forDownload=true',destfile=temp)
# header_lines <- readLines(temp, n=2) gtex <- read.delim(temp, skip=2,
# sep='\t') unlink(temp)

# write.table(gtex, file='gtex_all.txt', quote=F, sep='\t')

# download.file(url='http://www.broadinstitute.org/gtex/rest/file/download?portalFileId=119273&forDownload=true',destfile='GTEx_description.txt')

# metadata <- read.delim('GTEx_description.txt', sep='\t')
```


The metadata table seems to contain entries that are not in the RPKM table.


```r
# samp.id <- gsub('-','.',metadata$SAMPID) eligible.samples <- which(samp.id
# %in% colnames(gtex)) metadata <- metadata[eligible.samples,]
```


Select random heart, kidney and brain samples.


```r
# random.heart <- sample(which(metadata$SMTS=='Heart'), size=1)
# random.heart.samplename <- gsub('-','.',metadata[random.heart, 'SAMPID'])
# gtex.heart.fpkm <- as.numeric(gtex[,random.heart.samplename])

# random.brain <- sample(which(metadata$SMTS=='Brain'), size=1)
# random.brain.samplename <- gsub('-','.',metadata[random.brain, 'SAMPID'])
# gtex.brain.fpkm <- as.numeric(gtex[,random.brain.samplename])

# random.kidney <- sample(which(metadata$SMTS=='Kidney'), size=1)
# random.kidney.samplename <- gsub('-','.',metadata[random.kidney,
# 'SAMPID']) gtex.kidney.fpkm <- as.numeric(gtex[,random.kidney.samplename])
```


Get gene IDs on same format as the other data sets by removing the part after the dot; check ID uniqueness and write to file.


```r
# gtex.names <- gtex[,'Name'] temp_list <-
# strsplit(as.character(gtex.names), split='\\.') gtex.names.nodot <-
# unlist(temp_list)[2*(1:length(gtex.names))-1]

# gtex.fpkms <- data.frame(ENSG_ID=gtex.names.nodot,
# GTEx_heart=gtex.heart.fpkm,
# GTEx_brain=gtex.brain.fpkm,GTEx_kidney=gtex.kidney.fpkm)

# length(gtex.fpkms[,1]) length(unique(gtex.fpkms[,1]))

# write.table(gtex.fpkms,file='gtex_fpkms.txt',quote=F,sep='\t')
```


*RNA-seq Atlas*


```r
# temp <- tempfile()
# download.file(url='http://medicalgenomics.org/rna_seq_atlas/download?download_revision1=1',destfile=temp)
# atlas <- read.delim(temp, sep='\t') unlink(temp)

# atlas.fpkms <-
# atlas[,c('ensembl_gene_id','heart','hypothalamus','kidney')]
# colnames(atlas.fpkms) <-
# c('ENSG_ID','Atlas_heart','Atlas_brain','Atlas_kidney')
# write.table(atlas.fpkms,file='atlas_fpkms.txt',quote=F,sep='\t')
```


Combining F/RPKM values from public data sets
---------------------------------------------

We will join the data sets on ENSEMBL ID:s, losing a lot of data in the process - but joining on gene symbols or something else would lead to an even worse loss. 


```r
library(org.Hs.eg.db)  # for transferring gene identifiers
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
## The following object is masked from 'package:limma':
## 
##     plotMA
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, duplicated, eval, evalq, Filter, Find, get,
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
library(data.table)  # for collapsing transcript RPKMs
library(pheatmap)  # for nicer visualization
library(edgeR)  # for TMM normalization

# hpa.fpkms <- read.delim('hpa_fpkms.txt') altiso.fpkms <-
# read.delim('altiso_fpkms.txt') gtex.fpkms <- read.delim('gtex_fpkms.txt')
# atlas.fpkms <- read.delim('atlas_fpkms.txt')
```


The RNA-seq Atlas data set uses many different identifiers, while the other all use ENSG as the primary identifier

Approach 1: Merge on ENSEMBL genes (ENSG) as given in RNA-seq Atlas. Note that there are repeated ENSG ID:s in RNA-seq Atlas, as opposed to the other data sets, so we need to do something about that. In this case, we just sum the transcripts that belong to each ENSG gene. We use data.table for this.


```r
# data.dt <- data.table(atlas.fpkms) setkey(data.dt, ENSG_ID) temp <-
# data.dt[, lapply(.SD, sum), by=ENSG_ID] collapsed <- as.data.frame(temp)
# atlas.fpkms.summed <- collapsed[,2:ncol(collapsed)]
# rownames(atlas.fpkms.summed) <- collapsed[,1]

# atlas.fpkms.summed <- atlas.fpkms.summed[2:nrow(atlas.fpkms.summed),]
```


Finally, combine all the data sets into a data frame.


```r
# fpkms <- merge(hpa.fpkms, altiso.fpkms, by='ENSG_ID') fpkms <-
# merge(fpkms, gtex.fpkms, by='ENSG_ID') fpkms <- merge(fpkms,
# atlas.fpkms.summed, by.x='ENSG_ID', by.y=0) gene_id <- fpkms[,1] f <-
# fpkms[,2:ncol(fpkms)] rownames(f) <- gene_id
```


Check how many ENSG IDs we have left.


```r
# dim(f)
```


Approach 2: Try to map Entrez symbols to ENSEMBL to recover more ENSG IDs than already present in the table. 


```r
# m <- org.Hs.egENSEMBL mapped_genes <- mappedkeys(m) ensg.for.entrez <-
# as.list(m[mapped_genes]) remapped.ensg <-
# ensg.for.entrez[as.character(atlas$entrez_gene_id)]

# atlas.fpkms$remapped_ensg <- as.character(remapped.ensg)

# And add expression values data.dt <-
# data.table(atlas.fpkms[,2:ncol(atlas.fpkms)]) setkey(data.dt,
# remapped_ensg) temp <- data.dt[, lapply(.SD, sum), by=remapped_ensg]
# collapsed <- as.data.frame(temp) atlas.fpkms.summed <-
# collapsed[,2:ncol(collapsed)] rownames(atlas.fpkms.summed) <-
# collapsed[,1]
```


Combine data sets again


```r
# fpkms <- merge(hpa.fpkms, altiso.fpkms, by='ENSG_ID') fpkms <-
# merge(fpkms, gtex.fpkms, by='ENSG_ID') fpkms <- merge(fpkms,
# atlas.fpkms.summed, by.x='ENSG_ID', by.y=0) gene_id <- fpkms[,1] f <-
# fpkms[,2:ncol(fpkms)] rownames(f) <- gene_id write.table(f, file =
# 'published_rpkms.txt', quote=F)
```


Check how many ENSG IDs we have left.


```r
f <- read.delim("published_rpkms.txt", sep = " ")
# dim(f)
```


This looks much better. Let's proceed with this version of the data set. Start by a few correlation heat maps:

**Figure 1B**

Heatmap of Spearman correlations between published expression profiles (# genes = 13,537)


```r
pheatmap(cor(f, method = "spearman"))
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 


The brain samples are in a separate cluster, whereas the heart and kidney ones are intermixed.

Alternatively, one could use Pearson correlation (not shown in paper):


```r
pheatmap(cor(f))
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 


Sometimes the linear (Pearson) correlation works better on log values.  (not shown in paper)


```r
# pseudo <- 0.125 f.log <- log2(f+pseudo)
f.log <- log2.cpm(f)
pheatmap(cor(f.log))
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 

```r
write.table(f.log, file = "published_rpkms_log2cpm.txt", quote = F)
```


What if we drop the genes that have less than FPKM 1 on average? (not shown in paper)


```r
f.nolow <- f[-which(rowMeans(f) < 1), ]
# pheatmap(cor(log2(f.nolow+pseudo)))
pheatmap(cor(log2.cpm(f.nolow)))
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 


What if we use TMM normalization?

















































