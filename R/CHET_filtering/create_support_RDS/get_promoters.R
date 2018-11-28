library(data.table)
library(GenomicRanges)

# region around a TSS to allocate as a promoter region - This will extend 5' of TSS and could encompass coding region of a proximal but different gene
# or in the case of multiple transcripts the same gene.
PROMOTER_REGION<-500
OUT_DIR<-'__REDACTED__/pid/ANNOTATIONS'

## gencode
# wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.annotation.gtf.gz

g<-fread("__REDACTED__/pid/support/gencode.v26lift37.annotation.gtf")
## how to compute tss ?
g.transcript<-subset(g,V3=='transcript' & grepl("protein_coding",V9))
info<-do.call('rbind',lapply(strsplit(g.transcript$V9,";"),function(f){
        gid<-sub('gene_id \"(.*)\"','\\1',f[grep("gene_id",f)])
        tid<-sub(' transcript_id \"(.*)\"','\\1',f[grep("transcript_id",f)])
        gt<-sub(' gene_type \"(.*)\"','\\1',f[grep("gene_type",f)])
        tt<-sub(' transcript_type \"(.*)\"','\\1',f[grep("transcript_type",f)])
        gn<-sub(' gene_name \"(.*)\"','\\1',f[grep("gene_name",f)])
        return(cbind(gid,tid,gt,tt,gn))
}))



final.t<-cbind(g.transcript,info)
final.tpc<-subset(final.t,gt=='protein_coding' & tt=='protein_coding')
final.tpc$ensg<-sub("(ENSG[0-9]+)\\..*","\\1",final.tpc$gid)
final.tpc[final.tpc$V7=='+',V5:=V4]
final.tpc[final.tpc$V7=='-',V4:=V5]
final.tpc[,uid:=paste(V1,V4,ensg,sep=':')]
## only want distinct tss
ufinal.tpc<-unique(final.tpc,by='uid')
tss.gr<-with(ufinal.tpc,GRanges(seqnames=Rle(sub("chr","",V1,fixed=TRUE)),ranges=IRanges(start=V4,end=V4),strand=V7,ensg=ensg,name=gn))
prom.gr<-promoters(tss.gr,upstream=PROMOTER_REGION,downstream=0)
saveRDS(prom.gr,file=file.path(OUT_DIR,'promoters.RDS'))
