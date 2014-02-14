Downloading the F/RPKM data
===========================

Here, we download data from various public sources and extract the brain, heart and kidney samples.

"HPA": Human Protein Atlas

```r
temp <- tempfile()
download.file(url = "http://www.proteinatlas.org/download/rna.csv.zip", destfile = temp)
hpa <- read.csv(unz(temp, "rna.csv"))
unlink(temp)

hpa.heart <- hpa[hpa$Sample == "heart muscle", c("Gene", "Value")]
hpa.brain <- hpa[hpa$Sample == "cerebral cortex", c("Gene", "Value")]
hpa.kidney <- hpa[hpa$Sample == "kidney", c("Gene", "Value")]

hpa.fpkms <- merge(hpa.heart, hpa.brain, by = "Gene")
hpa.fpkms <- merge(hpa.fpkms, hpa.kidney, by = "Gene")
colnames(hpa.fpkms) <- c("ENSG_ID", "HPA_heart", "HPA_brain", "HPA_kidney")
```


Check if the identifiers are unique and write table to file.

```r
length(hpa.fpkms[, 1])
```

```
## [1] 20314
```

```r
length(unique(hpa.fpkms[, 1]))
```

```
## [1] 20314
```

```r

write.table(hpa.fpkms, file = "hpa_fpkms.txt", quote = F, sep = "\t")
```


"Altiso": Alternative isoform regulation in human tissue transcriptomes











