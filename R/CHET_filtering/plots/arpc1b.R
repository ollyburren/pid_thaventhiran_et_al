library(Gviz)
library(GenomicInteractions)
library(biomaRt)
library(data.table)
library(magrittr)

## plot a given region


cand_genes <- readRDS("__REDACTED__/pid/RESULTS/pid_cand_gene_chet_JUNE_2018.RDS")

e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

## hnisz

#plot.genes <- c('ENSG00000241685','ENSG00000130429')
plot.genes <- c('ENSG00000130429')


gene <- BiomartGeneRegionTrack(
    genome="hg19", name="Genes", transcriptAnnotation="symbol",
    mart=e75.genemart,
    #collapseTranscripts="meta",
    stackHeight=0.2, filters=list(with_ox_refseq_mrna=T,ensembl_gene_id=plot.genes),
    fill="black",col="black",cex.group = 1,fontcolor.group="black",just.group="left")

chr <- cand_genes[gene=='ARPC1B',]$chr
variant <- cand_genes[gene=='ARPC1B',]$hgvc
del <- cand_genes[gene=='ARPC1B',]$sv.coords

variant.start <- gsub(".*:([0-9]+).*","\\1",variant) %>% as.numeric
del.start <- gsub(".*:([^\\-]+)\\-(.*)","\\1",del) %>% as.numeric
del.end <- gsub(".*:([^\\-]+)\\-(.*)","\\2",del) %>% as.numeric

starts <- c(variant.start,del.start)
ends <- c(variant.start+1,del.end)
varinames <- c(variant,del)



mut<-AnnotationTrack(name="Variants",start=starts,end=ends,chr=chr,group=varinames,
    col="red",fill="red", groupAnnotation = "group",cex.group=1,fontcolor.group="red",just.group = "above")
    #tracks$se <- AnnotationTrack(name="SE Hnisz",se.gr,col="orange",fill="orange")

## load in exac coverage data

## region coords

region.start <- c(ranges(gene) %>% start,ranges(mut) %>% start) %>% min
region.end <- c(ranges(gene) %>% end,ranges(mut) %>% end) %>% max



cov <- fread('zcat__REDACTED__/tmp/Panel.chr7.coverage.txt.gz')
setnames(cov,names(cov) %>% make.names)
my.cov <- cov[X.chrom==7 & between(pos,region.start,region.end),]
exac.gr <- with(my.cov,GRanges(seqnames=Rle('chr7'),ranges=IRanges(start=pos,width=1L),mean=mean))
exac <-  DataTrack(exac.gr, name = "ExAC Coverage %",type="histogram",col.histogram='blue',fill.histogram ='blue')


tracks <- list()
tracks$axis<-GenomeAxisTrack()
tracks$gene <- gene
tracks$mut <- mut
tracks$exac <- exac
feature(tracks$mut) <- varinames

pdf("~/tmp/arpc1b.pdf")
plotTracks(tracks,main=sprintf("%s",'ARPC1B'),
min.width = 1,sizes=c(0.5,1,0.5,0.5),
col.title="white",background.title="black")
dev.off()
