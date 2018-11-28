library(Gviz)
library(GenomicInteractions)
library(biomaRt)
library(data.table)
library(magrittr)
library(cowplot)
library(GenomicInteractions)
library(simGWAS)





## load in meta analysis results

meta <- fread("__REDACTED__/analysis/pid/GWAS/lango-allen-li-ellinghaus-meta_analysis.tab")
## load in immunochip regions
ic.regions <- fread("__REDACTED__/DATA/JAVIERRE_ICHIP/support/dense.ic.regions.guessfm.bed")
setnames(ic.regions,c('chr','start','end','name'))
ic.regions[,c('start','end'):=list(start-1e5,end=end+1e5)]
ic.gr <- with(ic.regions,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end)))
meta.gr <- with(meta,GRanges(seqnames=Rle(CHR),ranges=IRanges(start=BP,width=1L)))
ol <- findOverlaps(meta.gr,ic.gr) %>% as.matrix
meta[ol[,1],region:=ol[,2]]
## get gw regions

regions <- meta[,list(min.p=min(overall.P)),by='region'][min.p<5e-8,]$region

regionID <- regions[1]

r.gr <- with(meta[region==regionID,],GRanges(seqnames=Rle(CHR),ranges=IRanges(start=BP,width=1L),mlp=-log10(overall.P),pid=paste(CHR,BP,sep=':')))
## grab index.snp <-
isnp <- r.gr[which.max(r.gr$mlp),]$pid

## need to compute LD

## stuff to compute LD
gethap <- function(cmd) {
    y=fread(cmd)
    ha <- simGWAS:::vcf2haps(as.matrix(y[,-c(1:9)]))
    rownames(ha) <- paste(gsub("chr","",y$V1),y$V2,sep=':')
    t(ha)
}
cor2 <- function (x) {
    1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}


BCF_TOOLS='/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools'
BCF_FILE='__REDACTED__/Data/reference/UK10K/chr16.bcf.gz'
bcf_maf=0.01

ld.cmd <- sprintf(
    "%s view -H __REDACTED__/Data/reference/UK10K/BCF/chr%s.bcf.gz --min-af %f:minor --max-alleles 2 --min-alleles 2 -r chr%s:%s-%s -Ov",
    BCF_TOOLS,ic.regions$chr[regionID],bcf_maf,ic.regions$chr[regionID],ic.regions$start[regionID],ic.regions$end[regionID])

h <- gethap(ld.cmd)
use <- apply(h,2,var)>0
h <- h[,use,drop=FALSE]
LD <- cor2(h)[,isnp]


ldColGrad <- c('darkblue','lightblue','green','orange','red')

range <- seq(0,1,length.out=6)
dat <- data.table(rowy=1:length(r.gr),mlp=r.gr$mlp,breaky=cut(LD[r.gr$pid]^2,breaks=range),pid=r.gr$pid)
dat$col <- ldColGrad[as.numeric(dat$breaky)]
## some of the LD info is missing remove for time being
#missing.ld <- dat[which(is.na(breaky)),]
#miss.gr <- r.gr[missing.ld$rowy,]
#r.gr <- r.gr[-missing.ld$rowy,]
#miss.dat <- dat[missing.ld$rowy,]
#dat <- dat[-missing.ld$rowy,]

dat[,breaky:=paste('group',as.numeric(dat$breaky)+1,sep='.')]
dt<-dcast(dat,rowy ~ breaky,value.var='mlp')[,-1]
## add in missing
setnames(dt,'group.NA','group.1')
#miss.snp.dat <- rep(NA,length(r.gr))
#miss.snp.dat[miss.dat$rowy] <- miss.dat$mlp
#r.gr$group.1 <- miss.snp.dat
## add in index snp
#index.snp.dat<-rep(NA,length(r.gr))
#index.snp.dat[dat[pid==isnp,]$rowy] <- dat[pid==isnp,]$mlp
#r.gr$group.7 <- index.snp.dat
#dt[miss.dat$rowy,group.1:=miss.dat$mlp]
dt[dat[pid==isnp,]$rowy,group.7:=dat[pid==isnp,]$mlp]
dt<-dt[,paste('group',1:7,sep='.'),with=FALSE] %>% as.data.frame %>% DataFrame


mcols(r.gr) <- dt

groups <- paste('group',1:7,sep='.')

library(biomaRt)

e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

gene <- BiomartGeneRegionTrack(
    genome="hg19", name="Genes", transcriptAnnotation="symbol",
    mart=e75.genemart,
    collapseTranscripts="meta",
    stackHeight=0.2, filters=list(with_ox_refseq_mrna=T,chromosome_name=ic.regions$chr[regionID],
    #start=start(query.region),end=end(query.region),biotype='protein_coding',hgnc_symbol=list('BCL6','LPP')),
    start=ic.regions$start[regionID]-1e5,end=ic.regions$end[regionID]+1e5,biotype='protein_coding'),
    fill="black",col="black",cex.group = 0.7,fontcolor.group="black",just.group="left")

pstart <- start(r.gr) %>% min
pend <- ranges(gene)[ranges(gene)$symbol=='SOCS1',] %>% end %>% max

#plot.region.gr <- GRanges(seqnames=Rle(unique(seqnames(gene)) %>% gsub("chr","",.)),ranges=IRanges(start=start(gene) %>% min,end=end(gene) %>% max))
plot.region.gr <- GRanges(seqnames=Rle(unique(seqnames(gene)) %>% gsub("chr","",.)),ranges=IRanges(start=pstart ,end=pend))
## add in atac.seq

getBed <- function(f){
  DT <- fread(f)
  if(length(names(DT))==5){
    setnames(DT,c('chr','start','end','name','score'))
  }else{
    setnames(DT,c('chr','start','end'))
  }
  DT.gr <- with(DT,GRanges(seqnames=Rle(gsub("chr","",chr)),ranges=IRanges(start=start,end=end))) %>% subsetByOverlaps(.,plot.region.gr)
  seqlevels(DT.gr) <- seqlevels(plot.region.gr)
  DT.gr
}

getReg <- function(f1,f2,maxgap=100){
  atac <- getBed(f1)
  h3k27ac <- getBed(f2)
  ol<-findOverlaps(atac,h3k27ac,maxgap=maxgap) %>% as.matrix
  union(atac[ol[,1]],h3k27ac[ol[,2]])
}

acd4 <- getReg("__REDACTED__/analysis/regulome/BED/aCD4/aCD4_ATAC.bed",
"__REDACTED__/analysis/regulome/H3K27ac/allPeaks/CD4_activated.bed")

rcd4 <- getReg("__REDACTED__/analysis/regulome/BED/rCD4/rCD4_ATAC.bed",
'__REDACTED__/analysis/regulome/H3K27ac/allPeaks/CD4_resting.bed')

nb <- getReg("__REDACTED__/analysis/regulome/ATACseq/NaiveBcell/naiveBcell_ATAC.bed",
'__REDACTED__/analysis/regulome/H3K27ac/allPeaks/CD19_NaiveBcell.bed')

mon <- getReg("__REDACTED__/analysis/regulome/ATACseq/Monocytes/Monocytes_ATAC.bed",
'__REDACTED__/analysis/regulome/H3K27ac/allPeaks/Monocytes.bed')



## add in GWAS catalogue hits

gwas.DT <- fread("__REDACTED__/pid/gwas_cat/gwas_cat_23_07_2018.tab")
setnames(gwas.DT,make.names(names(gwas.DT)))
gwas.DT[,CHR_POS:=as.numeric(CHR_POS)]
gwas.gr <- with(gwas.DT[P.VALUE<5e-8 & !is.na(CHR_POS),],GRanges(seqnames=Rle(CHR_ID),ranges=IRanges(start=as.numeric(CHR_POS),width=1L),disease=DISEASE.TRAIT))


## add in ImmunoBase
if(FALSE){
diseases <- c('AA','AS','ATD','JIA','NAR','PSC','CEL','CRO','IGE','MS','PBC','PSO','RA','SJO','SLE','T1D','UC','VIT')
for(d in diseases){
  cmd <- sprintf("curl -s https://www.immunobase.org/downloads/regions-files-archives/latest_1.11/Hs_GRCh37-%s-assoc_variantsBED | grep -v 'track' | grep -ve '^#' > __REDACTED__/pid/immunobase/%s.bed",d,d)
  message(cmd)
  system(cmd)
}
}

fnames <- list.files(path="__REDACTED__/pid/immunobase/",pattern="*.bed",full.names=TRUE)
imb.gr <- lapply(fnames,function(f){
  DT <- fread(f)[,disease:=basename(f) %>% gsub("\\.bed","",.)]
  setnames(DT,c('chr','start','pos','id','disease'))
  DT[,.(chr=gsub("^chr","",chr),pos,id,disease)]
}) %>% rbindlist %>% with(.,GRanges(seqnames=Rle(chr),ranges=IRanges(start=pos,width=1L),disease=disease)) %>% subsetByOverlaps(.,plot.region.gr)
seqlevels(imb.gr) <- seqlevels(plot.region.gr)

## pchic

hic<-fread("__REDACTED__/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab")
idx<-hic[,.(Total_CD4_NonActivated>5,Total_CD4_Activated>5,Naive_B>5,Total_B>5,Naive_CD4>5,Monocytes>5)] %>% rowSums >1
hic.f <- hic[idx,][ensg %in% ranges(gene)$gene,]
baits.gr<-hic.f[,.(baitChr,baitStart,baitEnd,baitID,ensg),key='baitID'] %>% with(.,GRanges(seqnames=Rle(baitChr),ranges=IRanges(start=baitStart,end=baitEnd),id=baitID,ensg=ensg))
pirs.gr<-hic.f[,.(oeChr,oeStart,oeEnd,oeID,ensg),key='oeID'] %>% with(.,GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,end=oeEnd),id=oeID))


if(FALSE){
  ol <- findOverlaps(pirs.gr,union(rcd4,nb) %>% union(.,mon) %>% union(.,imb.gr)) %>% as.matrix
  baits.gr <- baits.gr[ol[,1],]
  pirs.gr <- pirs.gr[ol[,1],]
}


seqlevels(baits.gr) <- seqlevels(baits.gr) %>% paste0('chr',.)
seqlevels(pirs.gr) <- seqlevels(pirs.gr) %>% paste0('chr',.)
interactions <- GenomicInteractions(baits.gr, pirs.gr, counts=1)


tracks <- list()
tracks$axis <- GenomeAxisTrack(name="Chr16")
tracks$meta <- DataTrack(name="Meta -log10(P)",r.gr,type=c("p"),groups=groups,col=c('grey',ldColGrad,'purple'),cex=1.5,legend=FALSE,col.axis="black")
#tracks$westra <- DataTrack(name="Westra",westra.37.gr,type=c("p"),cex=1.5,legend=FALSE)
tracks$gene <- gene
displayPars(tracks$gene)$collapseTranscripts <- TRUE
displayPars(tracks$gene)$shape <- "arrow"
displayPars(tracks$gene)$fontsize.group <- 16
tracks$gwas <- AnnotationTrack(name="GWAS",imb.gr,stacking="dense",col="black",fill="black",group=imb.gr$disease, just.group="above", groupAnnotation = "group",col.line="white")
#tracks$ataca <- AnnotationTrack(name="rCD4",rcd4,stacking="dense",col="black",fill="black")
tracks$cd4 <- AnnotationTrack(name="CD4",union(acd4,rcd4),stacking="dense",col="darkgreen",fill="darkgreen")
tracks$nb <- AnnotationTrack(name="B",nb,stacking="dense",col="orange",fill="orange")
tracks$mon <- AnnotationTrack(name="Mon",mon,stacking="dense",col="red",fill="red")
tracks$interaction <- InteractionTrack(name="pcHi-C",interactions,chromosome="chr16")
displayPars(tracks$interaction)$plot.outside <- FALSE
displayPars(tracks$interaction)$col.anchors.fill<-'black'
displayPars(tracks$interaction)$col.interactions<-'black'
displayPars(tracks$interaction)$col.anchors.line <- 'black'
pdf("~/tmp/clec16a.pdf",useDingbats=FALSE)
plotTracks(tracks, background.title = "white",fontcolor.title="black",col.border.title="black",cex.title = 1,from=pstart,to=pend,sizes=c(0.5,2,1,0.5,0.5,0.5,0.5,0.5))
dev.off()
