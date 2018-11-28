## test code to obtain pvalues for enrichment if found.

library(data.table)
library(rtracklayer)
library(wgsea)
manifest.file <- '__REDACTED__/pid/support/CondGWAS.tsv'
m.DT <- fread(manifest.file)
gwas.dir <- '__REDACTED__/pid/GWAS/summary_stats/'

## JUNE 2018 LATEST GWAS MAF > 0.05% MAF
if(FALSE){
abdf <- fread("__REDACTED__/analysis/pid/GWAS/AbDef.maf05.opr98.chr1-23.assoc.logistic.plink")
maf.thresh=0.01
UK10 <- readRDS("__REDACTED__/DATA/UK10K/UK10K_0.005_MAF.RDS")
UK10 <- UK10[,.(CHROM,POS,MAF)]
UK10[CHROM=='X',CHROM:='23']
UK10 <- UK10[,CHROM:=as.numeric(CHROM)]
DT.c <- merge(abdf,UK10,by.x=c('CHR','BP'),by.y=c('CHROM','POS'))[,.(ID=SNP,CHR,BP,BETA=log(OR),SE=log(OR)/STAT,P=P,MAF)]
DT.c <- DT.c[MAF>maf.thresh,]
out <- DT.c[,.(chr=CHR,start=BP-1,end=BP,id=ID,p.val=P)]
out <- out[order(chr,start),]
options(scipen=999)
write.table(out,file="__REDACTED__/DATA/JAVIERRE_GWAS/gwas/abdef_june.bed",col.names=FALSE,row.names=FALSE,quote=FALSE)
options(scipen=0)
}


if(FALSE){
  library(parallel)
  res <- mclapply(1:nrow(m.DT),function(i){
    #fname <- sprintf("%s.RDS",m.DT$label[i])
    message(sprintf("Processing %s",m.DT$label[i]))
    cmd <- sprintf("zcat %s",file.path(gwas.dir,m.DT$filename[i]))
    DT.f <- fread(cmd)
    setnames(DT.f,c('chr','start','end','id','p.val'))
    DT.f <- DT.f[!is.na(p.val),]
    DT.c <- merge(DT.f,UK10,by.x=c('chr','end'),by.y=c('CHROM','POS'))
    DT.c <- DT.c[MAF>maf.thresh,]
    DT.maf <- DT.c[,.(id,chr,position=end,maf=MAF,p.val)]
    rm.idx<-with(DT.maf,which((between(position,25e6,45e6) & chr==6)))
    DT.maf<-DT.maf[-rm.idx,]
    data.table(trait=m.DT$label[i],study.snps=nrow(DT.f),include.snps=nrow(DT.maf))
  },mc.cores=8) %>% rbindlist
}

maf.thresh=0.01

UK10 <- readRDS("__REDACTED__/DATA/UK10K/UK10K_0.005_MAF.RDS")
UK10 <- UK10[,.(CHROM,POS,MAF)]
UK10[CHROM=='X',CHROM:='23']
UK10 <- UK10[,CHROM:=as.numeric(CHROM)]

ld.gr<-import.bed("__REDACTED__/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37_ordered.bed")

addLDBlock<-function(DT,ld.gr){
  dt.gr<-with(DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=position,width=1L),idx=1:nrow(DT)))
  ol<-as.matrix(findOverlaps(dt.gr,ld.gr))
  DT[ol[,1],ld:=ol[,2]]
}

## robustly sum logs
logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form)
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

Var.data.cc <- function(f, N, s) {
    1 / (2 * N * f * (1 - f) * s * (1 - s))
}

## compute approx bayes factors and resultant posterior probabilities
## based on the assumption of one causal variant in a region
approx.bf.p <- function(p,f, N, s,pi_i,type='CC') {
    if(type=="QUANT") {
      sd.prior <- 0.15
      V <- Var.data(f, N)
    } else {
      sd.prior <- 0.2
      V <- Var.data.cc(f, N, s)
    }
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ## tABF - to add one we create another element at the end of zero for which pi_i is 1
    tABF <- c(lABF,0)
    vpi_i<-c(rep(pi_i,length(lABF)),1)
    sBF <- logsum(tABF + log(vpi_i))
    exp(lABF+log(pi_i)-sBF)
}

maf.thresh=0.05
out.dir <- '__REDACTED__/pid/GWAS_enrich_june_0.05/'

AbDef <- readRDS(file.path(out.dir,'ABDEF.RDS'))
AbDef[,pid:=paste(chr,position,sep=':')]
AbDef[,list(ppi.sum=sum(ppi)),by=ld][ppi.sum>0.5,]$ld
AbDef<-AbDef[,.(pid,ab.ppi=ppi,ab.p=p.val,ld)]

for(i in 1:nrow(m.DT)){
  fname <- sprintf("%s.RDS",m.DT$label[i])
  if(file.exists(file.path(out.dir,fname)))
    next
  message(sprintf("Processing %s",m.DT$label[i]))
  cmd <- sprintf("zcat %s",file.path(gwas.dir,m.DT$filename[i]))
  DT.f <- fread(cmd)
  setnames(DT.f,c('chr','start','end','id','p.val'))
  DT.f <- DT.f[!is.na(p.val),]
  DT.c <- merge(DT.f,UK10,by.x=c('chr','end'),by.y=c('CHROM','POS'))
  DT.c <- DT.c[MAF>maf.thresh,]
  DT.maf <- DT.c[,.(id,chr,position=end,maf=MAF,p.val)]
  rm.idx<-with(DT.maf,which((between(position,25e6,45e6) & chr==6)))
  DT.maf<-DT.maf[-rm.idx,]
  #ret[,c('trait','cases','controls':=list(m.DT$label[i],m.DT$cases[i],m.DT$comtrols[i])]
  DT.maf<-addLDBlock(DT.maf,ld.gr)
  DT.maf <- DT.maf[!is.na(ld),]
  ## need to add MAF
  if(is.na(m.DT$controls[i])){
    message("QUANT")
    DT.maf[,p.val:=as.numeric(p.val)]
    DT.maf[,ppi:=approx.bf.p(p.val,maf,m.DT$cases[i],0,1e-4,'QUANT'),by=ld]
  }else{
    message("CC")
    N<- m.DT$cases[i] + m.DT$controls[i]
    DT.maf[,ppi:=approx.bf.p(p.val,maf,N,m.DT$cases[i]/N,1e-4),by=ld]
  }

  saveRDS(DT.maf,file=file.path(out.dir,fname))
  message(sprintf("Saved to %s",file.path(out.dir,fname)))
  summ <- DT.maf[,list(total=sum(ppi),min.p=min(p.val),chr=unique(chr)),by=ld]
  fname <- sprintf("ld_summary_%s.RDS",m.DT$label[i])
  saveRDS(summ,file=file.path(out.dir,fname))
  message(sprintf("Saved to %s",file.path(out.dir,fname)))
}


rotate<-function(v){
	l=length(v)
	vc<-l
	m<-matrix(0,nrow=l,ncol=l)
	m[,1]<-v
	for(i in 2:l){
		pm<-m[,i-1]
		m[,i]<-c(pm[l],pm[-l])
	}
	return(m)
}
n.perms <- 1e4
pmi.thresh <- 0.95
##circ permutation method
bfiles <- list.files(path=out.dir,pattern="ld_summary.*",full.names=TRUE)
all.traits <- lapply(bfiles,readRDS)
names(all.traits)<-gsub("ld_summary_([^.]+)\\..*","\\1",basename(bfiles))
## we only take some forward for assessment
AbDef <- all.traits[['ABDEF']]
all.traits <- all.traits[names(all.traits) %in% c('T2D','CAD','T1D','UC','RA','AllS','SLE','CD','AST')]
all.Wrot.raw <- mclapply(seq_along(all.traits),function(i){
  trait <- names(all.traits)[i]
  message(trait)
  DT <- all.traits[[i]]
  DT <- DT[ld %in% intersect(ld,AbDef$ld),]
  DTa <- AbDef[ld %in% intersect(ld,DT$ld),]
  DT[,region:=FALSE]

  DT[total>pmi.thresh,region:=TRUE]
  m <- merge(DT,AbDef,by.x='ld',by.y='ld')
  # keep block structure
  ## we only consider LD blocks that are in both
  W<-wilcoxon(m$total.y,m$region)
  all.perms <- lapply(split(m[,.(ld,region)],m$chr.y),function(ch){
    rotate(ch$region)
  })

  perms <- do.call('rbind',lapply(all.perms,function(M){
    M[,sample.int(nrow(M),n.perms,replace=TRUE)]
  }))

  Wstar <- apply(perms,2,function(x){
    wilcoxon(m$total.y,x==1)
  })
  message(sprintf("%s:%f,%f",trait,sum(m$region),sum(!m$region)))
  enrich<-Z.value(W=W,Wstar=Wstar,n.in=sum(m$region),n.out=sum(!m$region))
  data.table(trait=trait,Z=enrich$Z.empirical$statistic,P=enrich$Z.empirical$p.value,
    Z.t=enrich$Z.theoretical$statistic,P.t=enrich$Z.theoretical$p.value,region.count=sum(m$region))
},mc.cores=8)
## work out enrichment or not

all.Wrot <- rbindlist(all.Wrot.raw)
all.Wrot <- all.Wrot[order(all.Wrot$Z),]
all.Wrot[,xlab:=trait]
all.Wrot[trait=='AllS',xlab:='Allergy']
all.Wrot[trait=='AST',xlab:='Asthma']
all.Wrot[,fill:='Other']
all.Wrot[trait %in% c('RA','SLE','T1D','UC','CD'),fill:="Autoimmune/Autoinflammatory"]
all.Wrot[trait %in% c('AllS','AST'),fill:='Allergy/Asthma']
all.Wrot$xlab<-factor(all.Wrot$xlab,levels=all.Wrot$xlab)

sig.thresh <- -log10(0.05/length(unique(all.Wrot$trait)))

library(cowplot)
library(ggrepel)

pz<-ggplot(all.Wrot,aes(x=xlab,y=Z,fill=fill)) + geom_bar(stat="identity") +
xlab("Trait") + ylab("AbDef-Assoc. PID\nEnrichment Z-score") +
scale_fill_discrete(name = "Trait Category") +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#dev.off()

save_plot("~/tmp/gwas_enrich_barplot_rotate_z.pdf",pz,base_width=8)
