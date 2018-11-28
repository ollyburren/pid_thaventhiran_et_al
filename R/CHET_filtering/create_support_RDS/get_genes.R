library(data.table)
library(GenomicRanges)

OUT_DIR<-'__REDACTED__/pid/ANNOTATIONS'

g<-fread("__REDACTED__/pid/support/gencode.v26lift37.annotation.gtf")
g.gene<-subset(g,V3=='gene'& grepl("protein_coding",V9))
info<-do.call('rbind',lapply(strsplit(g.gene$V9,";"),function(f){
        gid<-sub('gene_id \"(.*)\"','\\1',f[grep("gene_id",f)])
        tid<-sub(' transcript_id \"(.*)\"','\\1',f[grep("transcript_id",f)])
        gt<-sub(' gene_type \"(.*)\"','\\1',f[grep("gene_type",f)])
        tt<-sub(' transcript_type \"(.*)\"','\\1',f[grep("transcript_type",f)])
        gn<-sub(' gene_name \"(.*)\"','\\1',f[grep("gene_name",f)])
        return(cbind(gid,tid,gt,tt,gn))
}))

final<-cbind(g.gene,info)
final.pc<-subset(final,gt=='protein_coding')
final.pc$ensg<-sub("(ENSG[0-9]+)\\..*","\\1",final.pc$gid)
final.pc[,uid:=paste(V1,V4,ensg,sep=':')]
ufinal.pc<-unique(final.pc,by='uid')
g.gr<-with(ufinal.pc,GRanges(seqnames=Rle(sub("chr","",V1,fixed=TRUE)),ranges=IRanges(start=V4,end=V5),strand=V7,ensg=ensg,name=gn))
saveRDS(g.gr,file=file.path(OUT_DIR,'genes.RDS'))
