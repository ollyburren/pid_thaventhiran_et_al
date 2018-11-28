## Figure for paper

## want to layer three things
## MAF > 0.5 base GWAS
## Meta-analysis where we have it
## TACI region for rare variant GWAS

## BASE GWAS

library(data.table)
library(magrittr)
library(rtracklayer)

DT.abrare <- fread("__REDACTED__/analysis/pid/GWAS/AbDef.maf05.opr98.chr1-23.assoc.logistic.plink")
DT.abrare[,pid:=paste(CHR,BP,sep=':')]
setnames(DT.abrare,paste('rare',names(DT.abrare),sep='.'))
setkey(DT.abrare,'rare.pid')

## how many duplicated by position

duplicated(DT.abrare$rare.pid) %>% sum
#15585
## remove duplicated by position
DT.abrare <- DT.abrare[!rare.pid %in% DT.abrare[duplicated(rare.pid),]$rare.pid,]

DT.ab <- fread("__REDACTED__/analysis/pid/GWAS/AbDeficiency2.chr1-22.assoc.logistic.plink")
DT.ab[,pid:=paste(CHR,BP,sep=':')]
setkey(DT.ab,'pid')


duplicated(DT.ab$pid) %>% sum
#4631
DT.ab <- DT.ab[!pid %in% DT.ab[duplicated(pid),]$pid,]

## compare two GWAS - complete agreement
m<-DT.ab[DT.abrare,nomatch=0]

## next conduct the meta analysis

DT.li <- fread("__REDACTED__/tmp/CVID_Ichip/CVID_QCed_assoc.logistic")
DT.li[,pid:=paste(CHR,BP,sep=':')]
setnames(DT.li,paste('li',names(DT.li),sep='.'))
setkey(DT.li,'li.pid')

## are there duplicates ?
duplicated(DT.li$li.pid) %>% sum
#13
  DT.li <- DT.li[!li.pid %in% DT.li[duplicated(li.pid),]$li.pid,]

  m.DT <- DT.ab[DT.li,nomatch=0]

  ## check that A1 match
  m.DT[,allele.check:=1]
  m.DT[A1!=li.A1,allele.check:=-1]
  m.DT[,SE:=log(OR)/STAT]
  ## compute weights
  m.DT[,c('w','w.li'):=list(1/SE^2,1/li.SE^2)]
  ## compute overall SE
  m.DT[,overall.SE:=sqrt(1/(w+w.li))]
  m.DT[,overall.beta:=((log(OR)*w)+(log(li.OR) * w.li * allele.check))/(w+w.li)]
  m.DT[,overall.Z:=overall.beta/overall.SE]
  # compute pvalues
  m.DT[,overall.P:=2*pnorm(abs(overall.Z),lower.tail=FALSE)]

  META <- m.DT[,.(CHR,SNP,BP,A1,OR,P,li.SNP,li.A1,li.OR,li.P,overall.OR=exp(overall.beta),overall.P)]

## what happens if we perform the meta-analysis with rare ?

if(FALSE){
  m.DT <- DT.abrare[DT.li,nomatch=0]

  ## check that A1 match
  m.DT[,allele.check:=1]
  m.DT[rare.A1!=li.A1,allele.check:=-1]
  m.DT[,SE:=log(rare.OR)/rare.STAT]
  ## compute weights
  m.DT[,c('w','w.li'):=list(1/SE^2,1/li.SE^2)]
  ## compute overall SE
  m.DT[,overall.SE:=sqrt(1/(w+w.li))]
  m.DT[,overall.beta:=((log(rare.OR)*w)+(log(li.OR) * w.li * allele.check))/(w+w.li)]
  m.DT[,overall.Z:=overall.beta/overall.SE]
  # compute pvalues
  m.DT[,overall.P:=2*pnorm(abs(overall.Z),lower.tail=FALSE)]
}



## for plot we need p.rare,p.ab,p.meta

DT.abrare.c <- DT.abrare[,.(CHR=rare.CHR,BP=rare.BP,pid=rare.pid,P=rare.P,type='RARE')]
DT.ab.c <- DT.ab[,.(CHR,BP,pid,P,type='COMMON')]
DT.meta.c <- m.DT[,.(CHR,BP,pid,P=overall.P,type='META')]
if(FALSE)
  DT.meta.c <- m.DT[,.(CHR=rare.CHR,BP=rare.BP,pid=rare.pid,P=overall.P,type='META')]

DT.o <- rbindlist(list(DT.abrare.c,DT.ab.c,DT.meta.c))
DT.o[,mlp:=-log10(P)]

DT.o <- melt(DT.o,id.vars=c('CHR','BP','type'),measure.vars='mlp')

DT.o<-dcast(DT.o,CHR+BP~type+variable)


taci.chr <- 17
taci.start <- 16695860
taci.end <- 17022299

DT.o[,c('is_meta','is_taci','is_common'):=list(!is.na(COMMON_mlp) & !is.na(META_mlp),
  CHR==taci.chr & between(BP,taci.start,taci.end),
  !is.na(COMMON_mlp))]

DT.o[,plot.mlp:=COMMON_mlp]
## add meta
DT.o[is_meta==TRUE,plot.mlp:=META_mlp]
## add taci
DT.o[is_taci==TRUE,plot.mlp:=RARE_mlp]
DT.o<-DT.o[!is.na(plot.mlp),]
DT.o[,cat:='common']
DT.o[is_meta==TRUE,cat:='meta']
DT.o[is_taci==TRUE,cat:='taci']

## for plotting we remove nominally associated variants

#DT.o[plot.mlp>-log(0.1),]

### my code for plotting manhattan adapted from qqman !
DT.o[,c('CHR','BP'):=list(as.numeric(CHR),as.integer(BP))]
DT.o <- DT.o[order(as.numeric(CHR),as.numeric(BP)),]

## we want some spacing between chromosomes for formatting reasons

DT.o[,pos:=as.double(BP)]
DT.o[CHR>14,pos:=BP+1e8]
DT.o[,pos:=cumsum(pos)]

## add ld.blocks

r.gr <- fread("__REDACTED__/rds/hpc-work/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37.bed") %>%
with(.,GRanges(seqnames=V1,ranges=IRanges(start=V2,end=V3)))
d.gr <- with(DT.o,GRanges(seqnames=CHR,ranges=IRanges(start=BP,width=1L)))
ol <- findOverlaps(d.gr,r.gr) %>% as.matrix

DT.o[ol[,1],ld.id:=ol[,2]]

## compute ticks

ticks <- DT.o[,list(ticks=median(pos)),by=CHR]

## add gene annotations

top.hits <- DT.o[DT.o[, .I[plot.mlp == max(plot.mlp)], by=ld.id]$V1]
#top.hits <- top.hits[is_taci==TRUE | is_meta==TRUE,]
## filter so include those above a certain threshold

sig.hits<-top.hits[plot.mlp>-log10(1e-5)]
library(GenomicRanges)

sig.gr <- with(sig.hits,GRanges(seqnames=Rle(CHR),ranges=IRanges(start=BP,width=1L)))

library(biomaRt)
window.size<-1e6
ensembl_archive="feb2014.archive.ensembl.org" # Final V37 Ensembl
ensembl <- useMart(host=ensembl_archive, biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
sig.genes<-getBM(c("ensembl_gene_id","gene_biotype","external_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position","strand"),
filters=c("chromosome_name","start","end","biotype"),
values=list(chromosome_name=sig.hits$CHR,
  start=sig.hits$BP-window.size,end=sig.hits$BP+window.size,biotype='protein_coding'), mart=ensembl) %>% data.table

sig.genes.gr <- with(sig.genes,GRanges(seqnames=Rle(chromosome_name),
ranges=IRanges(start=start_position,end=end_position),gene=hgnc_symbol))
##

sig.hits[,c('nearest.gene','nearest.gene.dist'):=list(sig.genes.gr$gene[nearest(sig.gr,sig.genes.gr)],
  mcols(distanceToNearest(sig.gr,sig.genes.gr))$distance)
  ]

sig.hits[CHR==6,c('nearest.gene','nearest.gene.dist'):=list('MHC',0)]
sig.hits[,glabel:=sprintf("italic(%s)",nearest.gene)]
sig.hits <- sig.hits[sig.hits[, .I[plot.mlp == max(plot.mlp)], by=nearest.gene]$V1]

## replace these with prioritised genes
sig.hits[grep('CLEC16A',glabel),glabel:='italic(SOCS1)']
sig.hits[grep('DDC',glabel),glabel:='italic(IKZF1)']

## include genes
label.genes <- c('EOMES','MHC','CLEC16A','TNFRSF13B','PTPN2','ETS1','DDC')


##add annotations to plot object

#plot <- merge(DT.o,sig.hits[nearest.gene %in% label.genes ,.(CHR,BP,glabel,nearest.gene.dist)],by.x=c('CHR','BP'),by.y=c('CHR','BP'),all.x=TRUE)
plot <- merge(DT.o,sig.hits[ ,.(CHR,BP,glabel,nearest.gene.dist)],by.x=c('CHR','BP'),by.y=c('CHR','BP'),all.x=TRUE)


library(ggplot2)
library(cowplot)
library(ggrepel)


pp1<-ggplot(plot[plot.mlp>-log10(0.005),],aes(x=pos,y=plot.mlp,col=cat)) +
geom_point(size=0.2) +
scale_x_continuous(name="Chromosome", breaks=ticks$ticks, labels=(ticks$CHR)) +
xlab("Chr") + ylab("-log10(p)") +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(size=FALSE,label=FALSE) +
scale_size_discrete(range=c(0.5,2)) + geom_hline(yintercept=-log10(1e-5),col='blue',alpha=0.8) +
geom_hline(yintercept=-log10(5e-8),col='red',alpha=0.8) +
#geom_text_repel(aes(label=glabel),size=4,direction="y",point.padding=0.1,nudge_y=1,parse=TRUE)
geom_text_repel(aes(label=glabel),size=4,data = plot[!is.na(glabel),],
  nudge_y = 20 - plot[!is.na(glabel),]$plot.mlp,
  segment.size  = 0.2,
  direction     = "x",
  parse=TRUE,
  angle=-90,
  show.legend  = F
) + coord_cartesian(ylim=c(-log10(0.005),20)) + theme(axis.text.x = element_text(size=8),legend.position="top") +
scale_colour_manual("",values=c("dodgerblue","firebrick1","purple"),
labels=c(common="BRIDGE (MAF>5%)",meta="Meta-analysis\nLi et al.",taci="BRIDGE (MAF>0.5%)")) +
guides(colour = guide_legend(override.aes = list(size=5)))

save_plot("~/tmp/all_suggestive_abdef_man.pdf",pp1,base_width=8)


plot <- merge(DT.o,sig.hits[nearest.gene %in% label.genes ,.(CHR,BP,glabel,nearest.gene.dist)],by.x=c('CHR','BP'),by.y=c('CHR','BP'),all.x=TRUE)


pp2<-ggplot(plot[plot.mlp>-log10(0.005),],aes(x=pos,y=plot.mlp,col=cat)) +
geom_point(size=0.2) +
scale_x_continuous(name="Chromosome", breaks=ticks$ticks, labels=(ticks$CHR)) +
xlab("Chr") + ylab("-log10(p)") +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(size=FALSE,label=FALSE) +
scale_size_discrete(range=c(0.5,2)) + geom_hline(yintercept=-log10(1e-5),col='blue',alpha=0.8) +
geom_hline(yintercept=-log10(5e-8),col='red',alpha=0.8) +
#geom_text_repel(aes(label=glabel),size=4,direction="y",point.padding=0.1,nudge_y=1,parse=TRUE)
geom_text_repel(aes(label=glabel),size=4,data = plot[!is.na(glabel),],
  nudge_y = 22 - plot[!is.na(glabel),]$plot.mlp,
  segment.size  = 0.2,
  point.padding = 1,
  direction     = "x",
  parse=TRUE,
  angle=-90,
  show.legend  = F
) +  coord_cartesian(ylim=c(-log10(0.005),22)) + theme(axis.text.x = element_text(size=8),legend.position="top") +
scale_colour_manual("",values=c("dodgerblue","firebrick1","purple"),
labels=c(common="BRIDGE (MAF>5%)",meta="Meta-analysis\nLi et al.",taci="BRIDGE (MAF>0.5%)")) +
guides(colour = guide_legend(override.aes = list(size=5)))


save_plot("~/tmp/figure_abdef_man.pdf",pp2,base_width=8)
