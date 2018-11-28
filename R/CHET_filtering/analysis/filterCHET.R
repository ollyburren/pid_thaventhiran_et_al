library(data.table)
library(GenomicRanges)
library(magrittr)


DATA.DIR<-'__REDACTED__/pid/ANNOTATIONS/'

prom<-readRDS(file.path(DATA.DIR,'promoters.RDS'))

## exon regions
exon<-readRDS(file.path(DATA.DIR,'exons.RDS'))

## pchic superenhancer (for time being)
pse<-readRDS(file.path(DATA.DIR,'pchic_hnisz','se.gr.RDS'))

genes<-readRDS(file.path(DATA.DIR,'genes.RDS'))

#mevents<-lapply(manta.files,readRDS)
mevents <- list()

del.gr <- readRDS('__REDACTED__/pid/pid_paper_dels/pid_del_strict_v4_0.03.RDS')
uid <- paste(as.character(seqnames(del.gr)),del.gr$uid,sep=':')
del.gr$uid<-uid

mevents[['del']] <- del.gr
names(mevents)<-gsub("pid\\_([^_]+).*","\\1",basename(manta.files))
## drop bnd for time being
mevents<-mevents[names(mevents) != 'bnd']

computeFeatOl<-function(e.gr,f.gr){
  out<-logical(length=length(e.gr))
  ol<-unique(as.matrix(findOverlaps(e.gr,f.gr))[,1])
  out[ol]<-TRUE
  out
}

fl<-list(prom=prom,exon=exon,pse=pse)

annoMatrix<-function(e){
  tmp<-do.call('cbind',lapply(fl,function(f.gr){
    computeFeatOl(e,f.gr)
  }))
  ## get a list of all those events that overlap at least one features
  keep<-which(rowSums(tmp)!=0)
  info<-cbind(mcols(e),DataFrame(tmp))[keep,]
  e<-e[keep,]
  mcols(e)<-info
  e
}

res<-lapply(mevents,annoMatrix)

## for these want feature genes (i.e. a denormalised list of CNV's and genes and reasons)

mergeFeatures<-function(f,flist=fl){
  message("Merging features")
  mf<-mcols(f)
  rbindlist(lapply(names(fl),function(tname){
    message(sprintf("Merging %s",tname))
    p<-f[mf[[tname]],]
    m<-data.table(as.data.frame(mergeByOverlaps(p,fl[[tname]])))
    #m<-m[,.(p.seqnames,p.start,p.end,p.gt=ifelse(is.na(p.MANTA_GT),p.CANVAS_GT,p.MANTA_GT),AF_postfilt,p.uid,p.bridge_id,ensg,name)]
    m<-m[,.(p.seqnames,p.start,p.end,p.gt=ifelse(is.na(p.MANTA_GT),p.CANVAS_GT,p.MANTA_GT),AF_postfilt,p.uid,p.bridge_id,ensg,name)]
    m$ftype<-tname
    m
  }))
}

all<-lapply(res,function(x) unique(mergeFeatures(x)))

snps<-rbindlist(lapply(list.files(path='__REDACTED__/pid/VCF2/exonic_only/no_monomorphs/lof/gene_events_AUGUST_2018/RDS',full.names=TRUE),readRDS))
setnames(snps,'#CHROM','CHROM')
snps[,hgvc:=sprintf("%s:%d%s>%s",CHROM,POS,REF,ALT)]

## merge the events with related SNPs to get a really large list
pid.genes.2015<-fread("__REDACTED__/pid/all.pid.genes.csv")


if(FALSE){
  ## method for preparing
  library(biomaRt)
  library(magrittr)
  ## these were given to me by Will
  pid.2017<-scan("__REDACTED__/pid/IUIS_2017.csv",character())
  ensembl_archive="feb2014.archive.ensembl.org" # Final V37 Ensembl
  ensembl <- useMart(host=ensembl_archive, biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
  pid.genes.2017<-getBM(c("ensembl_gene_id","gene_biotype","external_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position","strand"),
    filters=c("hgnc_symbol"),
    values=list(hgnc_symbol=pid.2017), mart=ensembl
  ) %>% data.table
  pid.genes.2017 <- pid.genes.2017[gene_biotype!='LRG_gene',]
  #some are missing have asked Will about this.
  #pid.2017[!pid.2017 %in% pid.genes.2017$hgnc_symbol]
  write.table(pid.genes.2017,file="__REDACTED__/pid/all.pid.genes.IUIS_2017.tab",quote=FALSE,row.names=FALSE,sep="\t")
}
pid.genes.2017 <- fread("__REDACTED__/pid/all.pid.genes.IUIS_2017.tab")

## these are lists of anon id's that are available on reasonable request but I have ommitted for privacy reasons
poss.diag<- c('please contact the authors')
def.diag<-c('please contact the authors')
ws<-lapply(names(all),function(n){
    x<-all[[n]]
    x[,p.coords:=sprintf("%s:%d-%d",p.seqnames,p.start,p.end)]
    tmp<-unique(merge(x,snps,by.x=c('ensg','p.bridge_id'),by.y=c('Gene','ind'),allow.cartesian=TRUE))
    tmp<-subset(tmp,BIOTYPE=='protein_coding' & gt %in% c(2,3) & minOPR>0.8)
    tmp<-unique(tmp[,.(ensg,name,ftype,p.bridge_id,p.seqnames,p.uid,p.coords,p.gt,AF_postfilt,gt,hgvc,GNOMAD_AF,AF_WGS10K,Existing_variation,Consequence,CADD_PHRED)])
    tmp$pid_2015<-tmp$ensg %in% pid.genes.2015$ensembl_gene_id
    tmp$pid_2017<-tmp$ensg %in% pid.genes.2017$ensembl_gene_id
    #tmp[,sv.gt:=p.gt+1]
    tmp$sv.gt<-0
    tmp[tmp$p.gt %in% c('1/0','0/1')]$sv.gt<-2
    tmp[tmp$p.gt=='1/1']$sv.gt<-3
    tmp[tmp$p.gt=='0/0']$sv.gt<-1
    tmp$p.gt<-NULL
    tmp$sv.type<-n
    setnames(tmp,c('ensg','gene','feat.type','individual','chr','sv.uid','sv.coords','sv.AF','snp.gt','hgvc','gnomad_af','wgs10k_af','exist.var','consequence','cadd','pid_2015','pid_2017','sv.gt','sv.type'))
    setcolorder(tmp,c('sv.type','sv.uid','ensg','gene','feat.type','individual','chr','sv.gt','sv.coords','sv.AF','snp.gt','hgvc','gnomad_af','wgs10k_af','exist.var','consequence','cadd','pid_2015','pid_2017'))
    tmp$poss.diag<-tmp$individual %in% poss.diag
    tmp$def.diag<-tmp$individual %in% def.diag
    ##fix numeric columns
    for(n in c('snp.gt','gnomad_af','wgs10k_af','cadd')){
      tmp[[n]]<-as.numeric(tmp[[n]])
    }
    tmp[['chr']] <- as.numeric(as.character(tmp[['chr']]))
    tmp
})

## filters

## take a look at frameshift

getByCons<-function(DT,term){
  DT[grep(term,DT$consequence),]
}

filtByAF<-function(DT,GAF,WAF){
  GNOMAD_MAF <- ifelse(DT$gnomad_af>0.5,1-DT$gnomad_af,DT$gnomad_af)
  WGS10K_MAF <- ifelse(DT$wgs10k_af>0.5,1-DT$wgs10k_af,DT$wgs10k_af)
  tDT <- data.table(GNOMAD_MAF,WGS10K_MAF)
  idx <- which(is.na(tDT$GNOMAD_MAF) | (tDT$GNOMAD_MAF<GAF & tDT$WGS10K_MAF<WAF))
  DT[idx,]
}

filtByCADD<-function(DT,score){
  DT[as.numeric(DT$cadd)>score,]
}

filtByType<-function(DT,type){
  DT[DT$feat.type==type,]
}

filtByPID<-function(DT){
  DT[DT$pid_2015 ==TRUE | DT$pid_2017 ==TRUE ,]
}

dels <- ws[[1]]
dels[,tun:=paste(sv.uid,hgvc,sep=':')]
mergecols <- c('ensg','gene','feat.type')
cnames <- names(dels)[!names(dels) %in% mergecols ]
dels.merge <- dels[,lapply(mergecols, function(y) paste(unique(get(`y`)),sep=',',collapse=',')),by=tun]
setnames(dels.merge,c('V1','V2','V3'),c('ensg_merge','gene_merge','feat_merge'))

dels <- merge(dels,dels.merge,by.x='tun',by.y='tun')
dels.merge <- dels[,list(pid_2015_merge=sum(pid_2015)>0,pid_2017_merge=sum(pid_2017)>0),by=tun]
dels <- merge(dels,dels.merge,by.x='tun',by.y='tun')

## be careful here as when removing duplicates we might remove that a row contains a pid gene - others should be OK as SV specific


dels <- dels[!duplicated(tun),]
dels[,c('ensg','gene','feat.type','pid_2015','pid_2017'):=list(ensg_merge,gene_merge,feat_merge,pid_2015_merge,pid_2017_merge)]
dels[,c('ensg_merge','gene_merge','feat_merge','pid_2015_merge','pid_2017_merge'):=list(NULL,NULL,NULL,NULL,NULL)]

## add in Hana's annotations as to whether a validated deletion
# she sent an excel spreadsheet which I then converted to a file of confirmed (by her - Hana_confirm column) deletions

cdel<-scan("~/tmp/confirmed.txt",character())

dels[,hana_confirmed:=sv.coords %in% cdel]


ws[[1]] <- dels

all.rare<-lapply(ws,filtByAF,GAF=1e-3,WAF=1e-3)
all.hi.cadd.rare<-lapply(all.rare,filtByCADD,score=20)
all.hi.cadd.rare[[1]][,sv.start:=as.numeric(gsub(".*:([^\\-]+)\\-.*","\\1",sv.coords))]
all.hi.cadd.rare[[1]][is.na(chr),chr:=23]


by.pid<-rbindlist(lapply(all.hi.cadd.rare,filtByPID))
by.pid <- by.pid[order(chr,sv.start),]
by.pid[,sv.start:=NULL]
by.pid[,tun:=NULL]
## remove those with hom snps as unlikely to be a CHET
by.pid <- by.pid[snp.gt!=3 & sv.gt!=3,]

library(xlsx)
write.xlsx(by.pid, file="__REDACTED__/pid/RESULTS/pid_cand_gene_chet_AUGUST_2018.xlsx", sheetName="PID Candidate Gene CHET")
saveRDS(by.pid,file="__REDACTED__/pid/RESULTS/pid_cand_gene_chet_AUGUST_2018.RDS")

final.list<-rbindlist(all.hi.cadd.rare)
final.list <- final.list[order(chr,sv.start),]
final.list[,sv.start:=NULL]
final.list[,tun:=NULL]
final.list <- final.list[snp.gt!=3 & sv.gt!=3,]
write.xlsx(final.list, file="__REDACTED__/pid/RESULTS/all_CHET_AUGUST_2018.xlsx", sheetName="All CHET")

## compute summary statistics for flowchart explaining filtering

## total number of SV deletions considered
total.sv <- unique(res[[1]]$uid) %>% length
## events with LOF MAF < 10^-3 GNOMAD
gnomad.lof <- unique(all.rare[[1]]$tun) %>% length
## high CADD rare
cadd <- unique(all.hi.cadd.rare[[1]]$tun) %>% length
## pid CHET
pid <- nrow(by.pid)


unique(all[[1]]$p.uid) %>% length

## process SV homs


homs<-rbindlist(lapply(names(all),function(n){
  tmp<-subset(all[[n]],p.gt=="1/1")
  #tmp[,sv.gt:=p.gt+1]
  tmp$sv.gt<-0
  tmp[tmp$p.gt %in% c('1/0','0/1')]$sv.gt<-2
  tmp[tmp$p.gt=='1/1']$sv.gt<-3
  tmp[tmp$p.gt=='0/0']$sv.gt<-1
  tmp$p.gt<-NULL
  tmp$p.gt<-NULL
  tmp[,p.coords:=sprintf("%s:%d-%d",p.seqnames,p.start,p.end)]
  tmp$sv.type<-n
  tmp$poss.diag<-tmp$p.bridge_id %in% poss.diag
  tmp$def.diag<-tmp$p.bridge_id %in% def.diag
  tmp$pid_2015<-tmp$ensg %in% pid.genes.2015$ensembl_gene_id
  tmp$pid_2017<-tmp$ensg %in% pid.genes.2017$ensembl_gene_id
  #tmp<-tmp[,.(sv.type,p.coords,AF_postfilt,ensg,name,ftype,p.bridge_id,p.seqnames,sv.gt,pid_2015,pid_2017,poss.diag,def.diag,p.uid,p.cluster)]
  tmp<-tmp[,.(sv.type,p.coords,AF_postfilt,ensg,name,ftype,p.bridge_id,p.seqnames,sv.gt,pid_2015,pid_2017,poss.diag,def.diag,p.uid)]
  setnames(tmp,c('sv.type','sv.coords','AF','ensg','gene','feat.type','individual','chr','sv.gt','pid_2015','pid_2017','poss.diag','def.diag','sv.uid'))
  tmp$chr<-as.numeric(tmp$chr)
  tmp
}))


## merge events
mergecols <- c('ensg','gene','feat.type')
homs.merge <- homs[,lapply(mergecols, function(y) paste(unique(get(`y`)),sep=',',collapse=',')),by=sv.uid]
setnames(homs.merge,c('V1','V2','V3'),c('ensg_merge','gene_merge','feat_merge'))

homs <- merge(homs,homs.merge,by.x='sv.uid',by.y='sv.uid')
homs.merge <- homs[,list(pid_2015_merge=sum(pid_2015)>0,pid_2017_merge=sum(pid_2017)>0),by=sv.uid]
homs <- merge(homs,homs.merge,by.x='sv.uid',by.y='sv.uid')



homs <- homs[!duplicated(sv.uid),]
homs[,c('ensg','gene','feat.type','pid_2015','pid_2017'):=list(ensg_merge,gene_merge,feat_merge,pid_2015_merge,pid_2017_merge)]
homs[,c('ensg_merge','gene_merge','feat_merge','pid_2015_merge','pid_2017_merge'):=list(NULL,NULL,NULL,NULL,NULL)]
homs[,sv.start:=as.numeric(gsub(".*:([^\\-]+)\\-.*","\\1",sv.coords))]
homs <- homs[order(chr,sv.start),]
homs[,sv.start:=NULL]
## this would fit with prevalence of 50 per 100,000
homs<-homs[AF<sqrt(1e-04),]
homs[,hana_confirmed:=sv.coords %in% cdel]
homs[hana_confirmed==TRUE,]


write.xlsx(homs, file="__REDACTED__/pid/RESULTS/hom_sv_AUGUST_2018.xlsx", sheetName="all")
saveRDS(homs,file="__REDACTED__/pid/RESULTS/hom_sv_AUGUST_2018.RDS")
