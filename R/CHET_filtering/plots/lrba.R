library(Gviz)
library(GenomicInteractions)
library(biomaRt)
library(data.table)
library(magrittr)
library(cowplot)
library(GenomicInteractions)

gname <- 'LRBA'
ensg <- 'ENSG00000198589'
pat.id <- 'F00642'
hic.file <- '__REDACTED__/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab'

DT <- data.table(chr=c('4','4'),start=c(151199051,151870474),end=c(151199051,152019030),type=c('SNP','DELETION'),mut.name=c('4:151199051C>A','4:151870474_152019030del'))

e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

gene <- BiomartGeneRegionTrack(
    genome="hg19", name="Genes", transcriptAnnotation="symbol",
    mart=e75.genemart,
    #collapseTranscripts="meta",
    stackHeight=0.2, filters=list(with_ox_refseq_mrna=T,ensembl_gene_id=ensg),
    fill="black",col="black",cex.group = 1,fontcolor.group="black",just.group="left")

## we can get the range from here
#extend by 500kb either side
e.region<-1e5
region.gr <- GRanges(seqnames=Rle('4'),ranges=IRanges(start=start(head(gene@range,1))-e.region,end=end(tail(gene@range,1))+e.region))
(load('__REDACTED__/pid/hnisz_SE_annotation_all_merged.RData'))
se.gr<-with(merge.results.se,GRanges(seqnames=Rle(se.gr.seqnames),ranges=IRanges(start=se.gr.start,end=se.gr.end),uid=uid))
se.gr<-reduce(subsetByOverlaps(se.gr,region.gr))
seqlevels(se.gr) <- '4'

header <- system(sprintf("head -1 %s",hic.file),intern=TRUE) %>% strsplit(.,"\t") %>% unlist
hic <- fread(sprintf("grep %s %s",ensg,hic.file))
setnames(hic,header)

bait.gr <- with(hic,GRanges(baitChr,IRanges(start=baitStart,end=baitEnd),idx=1:nrow(hic)))
oe.gr <- with(hic,GRanges(oeChr,IRanges(start=oeStart,end=oeEnd),idx=1:nrow(hic)))


mut.gr <- with(DT,GRanges(chr,IRanges(start,end)))

hics <- subsetByOverlaps(oe.gr,mut.gr)

seqlevels(bait.gr) <- 'chr4'
seqlevels(hics) <- 'chr4'
seqlevels(oe.gr) <- 'chr4'


interaction <- GenomicInteractions(bait.gr[hics$idx], hics, counts=1)
all.interaction <- GenomicInteractions(bait.gr, oe.gr, counts=1)



mut.gr <- with(DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),type=type))

tracks<-list()
tracks$axis<-GenomeAxisTrack(region.gr)
tracks$gene<-gene
tracks$mut<-AnnotationTrack(name="Variants",start=DT$start,end=DT$end,chr=DT$chr,group=DT$mut.name,
col="red",fill="red", groupAnnotation = "group",cex.group=1,fontcolor.group="red",just.group = "above")
tracks$se <- AnnotationTrack(name="SE Hnisz",se.gr,col="orange",fill="orange")
feature(tracks$mut) <- DT$mut.name
tracks$interactions <- InteractionTrack(name='Javierre et al', interaction)
displayPars(tracks$interactions)$anchor.height <- 0.1
displayPars(tracks$interactions)$col.anchors.line <- "black"
displayPars(tracks$interactions)$col.anchors.fill <- "darkgreen"
displayPars(tracks$interactions)$col.interactions <- "darkgreen"
#tracks$all.interactions <- InteractionTrack(name='All Javierre et al', all.interaction)
#displayPars(tracks$all.interactions)$anchor.height <- 0.1
#displayPars(tracks$all.interactions)$col.anchors.line <- "darkgrey"
#displayPars(tracks$all.interactions)$col.anchors.fill <- "lightgrey"
#displayPars(tracks$all.interactions)$col.interactions <- "grey"

pdf("~/tmp/LRBA.pdf")
plotTracks(tracks,main=sprintf("%s(%s)",gname,pat.id),
min.width = 1,sizes=c(1,1,2,1,2),
col.title="white",background.title="black")
dev.off()
