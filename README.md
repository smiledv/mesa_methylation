### Finding significant CpG sites located inside the promoter of genes:

```
#BiocManager::install("Rsubread")
library(Rsubread)

# we defined promoter 3kb upstream and 10kb downstream around transcription start site (STT)
GeneBody <- promoterRegions("hg19", upstream=3000, downstream=10000)

library("GenomicRanges")

gr1 <- GRanges(seqnames=GeneBody$Chr, ranges=IRanges(start=GeneBody$Start, 
                                                     end=GeneBody$End), strand=GeneBody$Strand)

# we make a data frame of our significant sites including chrom, start, end (start+1)
sig_site <- read.csv(file = "MESA-meth-seq-exam5-vs-exercm5c_sd-full-model-with-cell-contamination-adjustment-pvalue_05.csv", header = T)

sig_site <- sig_site[sig_site$qvalue < 0.05,]
sig_site <- sig_site[!is.na(sig_site$qvalue),]
sig_site <- sig_site[-(which(sig_site$chr %in% "X")),]
sig_site <- sig_site[,c("chr","pos")]
sig_site$end <- sig_site$pos+1
colnames(sig_site)[2] <- "start"
sig_site$chr <- paste("chr", sig_site$chr, sep = "")

gr2 <- GRanges(seqnames=sig_site$chr, ranges=IRanges(start=sig_site$start, 
                                                     end=sig_site$end),strand = "*")

# finding orverlaps
result <- GRanges(subsetByOverlaps(gr2, gr1))
# so we found 80 overlaps with promoters

# we need to find the genes of these promoter and we have the annotations in the csv file
#BiocManager::install("Repitools")
library(Repitools)

#changing the overlapping sites to a data frame
overlap <- annoGR2DF(result)
overlap$chr <- gsub("chr","",as.character(overlap$chr))
overlap$id <- paste(overlap$chr, overlap$start, sep = "_")
overlap <- overlap[,c("id","width")]

#extracting columns we need from the csv file
sig_site <- read.csv(file = "MESA-meth-seq-exam5-vs-exercm5c_sd-full-model-with-cell-contamination-adjustment-pvalue_05.csv", header = T)
sig_site <- sig_site[sig_site$qvalue < 0.05,]
sig_site <- sig_site[!is.na(sig_site$qvalue),]
sig_site <- sig_site[-(which(sig_site$chr %in% "X")),]
sig_site <- sig_site[,c("id","chr","pos","symbol","beta")]

#merging the sig sites in promoter with the csv file
sig_site_prom <- merge(sig_site,overlap,by=c("id"))
sig_site_prom <- sig_site_prom[,-c(6)]

write.csv(sig_site_prom, "sig_site_prom.csv", quote = F, row.names=F)
```
        
