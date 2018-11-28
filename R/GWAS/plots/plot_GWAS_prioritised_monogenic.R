### prototype plotter for linear chicp integration

library(Gviz)
library(GenomicInteractions)
library(biomaRt)
library(data.table)
library(magrittr)


# full results with posterior

#pmi <- readRDS('__REDACTED__/pid/GWAS/ABDEF.RDS')
pmi <- fread('__REDACTED__/abdef_ai/pmi/ABDEF.pmi')

e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

## input is a set of top hits
#         Ensembl_ID        Biotype Gene_Name        Type     Score Gene_Chr
#1: ENSG00000175354 protein_coding     PTPN2 interaction 0.7474821       18
#   Gene_Start  Gene_End
#1:   12785477  12929642

## get from output of COGS
res<-fread("__REDACTED__/abdef_ai/gene_score/ABDEF.pmi_prioritised.tab")
#top.hits<-res[res$overall_ppi>0.5 & biotype=='protein_coding',]
all.hits <- res[biotype=='protein_coding',]

ensembl_archive <- 'feb2014.archive.ensembl.org'
ensembl <- useMart(host=ensembl_archive, biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
gene.details<-getBM(c("ensembl_gene_id","gene_biotype","external_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position"),
  mart=e75.genemart,filters=list(ensembl_gene_id=all.hits$ensg))
gene.details<-data.table(gene.details)
setkey(gene.details,ensembl_gene_id)
setkey(all.hits,ensg)
all.hits<-all.hits[gene.details]
all.hits<-all.hits[order(chromosome_name,start_position,overall_ppi),]
all.hits<-all.hits[,c(1:3,7,11,16:18),with=FALSE]
setnames(all.hits,c('Ensembl_ID','Biotype','Gene_Name','Type','Score','Gene_Chr','Gene_Start','Gene_End'))
all.hits<-all.hits[order(Score,Gene_Chr,Gene_Start,decreasing=TRUE),]
#write.table(all.hits,file="~/tmp/CVID.csv",sep=",",row.names=FALSE,quote=FALSE)
#write.table(all.hits[Score>0.5,],file="~/tmp/CVID_gt_0.5.csv",sep=",",row.names=FALSE,quote=FALSE)


top.hits<-all.hits[all.hits$Score>0.4 & Biotype=='protein_coding',]

## Hi-C data
hic<-fread("__REDACTED__/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab")
baits.gr<-unique(hic[,.(baitChr,baitStart,baitEnd,baitID,ensg),key='baitID']) %>% with(.,GRanges(seqnames=Rle(paste0('chr',baitChr)),ranges=IRanges(start=baitStart,end=baitEnd),id=baitID,ensg=ensg))
pirs.gr<-unique(hic[,.(oeChr,oeStart,oeEnd,oeID,ensg),key='oeID']) %>% with(.,GRanges(seqnames=Rle(paste0('chr',oeChr)),ranges=IRanges(start=oeStart,end=oeEnd),id=oeID))




## step one is to extend the region around gene by some distance

getGenes<-function(gr){
    genesToGet<-unique(subsetByOverlaps(baits.gr,gr)$ensg)
    BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=e75.genemart,
        collapseTranscripts="meta",shape = "arrow",
        stackHeight=0.2, filters=list(with_ox_refseq_mrna=T,ensembl_gene_id=genesToGet))
}

getInteraction<-function(gr){
    oeid<-subsetByOverlaps(pirs.gr,gr)$id
    bids<-baits.gr[baits.gr$ensg==gr$ensg,]$id
    fhic<-subset(hic,(oeID %in% oeid & baitID %in% bids) & biotype=='protein_coding')
    b.gr<-with(fhic,GRanges(seqnames=Rle(paste0('chr',baitChr)), IRanges(start=baitStart, end=baitEnd)))
    o.gr<-with(fhic,GRanges(seqnames=Rle(paste0('chr',oeChr)), IRanges(start=oeStart, end=oeEnd)))
    GenomicInteractions(b.gr, o.gr, counts=1)

}

getSNPS <- function(gr){
  rchr <- gsub("chr","",as.character(seqnames(gr)))
  pmi.reg <- subset(pmi,chr==rchr & between(position,start(gr),end(gr)))
  pmi.reg[,mlp:=-log10(p.val)]
  with(pmi.reg,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L),mlp=mlp,ppi=ppi))
}

## next plot -log10(p)


pdf("~/tmp/abdef_ai_tophits.pdf")
pdf("~/tmp/stim1.pdf")
ext.dist<-1e5


for(gidx in 1:nrow(top.hits)){
#gidx <- 2
message(gidx)
qh <- top.hits[gidx,]

region.gr<-with(qh,GRanges(seqnames=Rle(paste0('chr',Gene_Chr)),ranges=IRanges(start=Gene_Start-ext.dist,end=Gene_End+ext.dist),ensg=Ensembl_ID,gscore=Score,name=Gene_Name))

maintitle <- sprintf("%s:%s",as.character(seqnames(region.gr)),region.gr$name)
tracks <- list()
tracks[['axis']]<-GenomeAxisTrack(name=as.character(seqnames(region.gr)))
tracks[['Genes']] <- getGenes(region.gr)
tracks[['Interactions']] <- InteractionTrack(name="Javierre",getInteraction(region.gr))
## mlp tracks
snp.gr <- getSNPS(region.gr)
mlp.gr<-snp.gr
mcols(mlp.gr)<-mcols(mlp.gr)[['mlp']]
ppi.gr <- snp.gr
mcols(ppi.gr)<-mcols(ppi.gr)[['ppi']]

tracks[['-log10(P)']] <- DataTrack(name="-log10(P)",snp.gr)
tracks[['Posterior']] <- DataTrack(name="PPi",ppi.gr,col='red')
plotTracks(main=maintitle,tracks)
}
dev.off()
