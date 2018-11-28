
library(data.table)
library(magrittr)

BASE.DIR <- '__REDACTED__/analysis/pid/GWAS'

DT.ab <- fread(file.path(BASE.DIR,'AbDef.maf05.opr98.chr1-23.assoc.logistic.plink'))
DT.ab[,pid:=paste(CHR,BP,sep=':')]
setkey(DT.ab,'pid')
duplicated(DT.ab$pid) %>% sum
#4631
DT.ab <- DT.ab[!pid %in% DT.ab[duplicated(pid),]$pid,]
## next conduct the meta analysis
DT.li <- fread(file.path(BASE.DIR,'CVID_QCed_assoc.logistic'))
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
write.table(META,file=file.path(BASE.DIR,'lango-allen-li-ellinghaus-meta_analysis_october.tab'),quote=FALSE,row.names=FALSE,sep="\t")


## create supplementary table

library(GenomicRanges)

ld.regions <- fread("__REDACTED__/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37_ordered.bed")
setnames(ld.regions,c('chr','start','end','name'))
ld.gr <- with(ld.regions,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end)))
ol.mhc <- GRanges(seqnames=Rle(6),ranges=IRanges(start=25e6,end=35e6)) %>% findOverlaps(ld.gr,.) %>% as.matrix
meta.gr <- with(META,GRanges(seqnames=Rle(CHR),ranges=IRanges(start=BP,width=1L)))
ol <- findOverlaps(meta.gr,ld.gr) %>% as.matrix
META[ol[,1],region:=ol[,2]]
## collapse MHC
META[region %in% ol.mhc,region:=1209]

## get min SNP by region
by.reg <- META[META[,.I[which.min(overall.P)],by=region]$V1][overall.P<1e-5,][order(overall.P,CHR),]

out <- by.reg[,.(chr=CHR,position_GRCh37=BP,rsid=SNP,effect_allele=A1,or_meta=signif(overall.OR,digits=2),p.value=signif(overall.P,digits=2),region)]

## use COGS to create a list of putative causal genes

ld.gr <- load_ld_regions('__REDACTED__/DATA/rCOGS/hapmap_recomb.bed')
maf.DT <- load_ref_maf('__REDACTED__/DATA/rCOGS/uk10k_0.001.tab',min.maf=0.001)
tmp.file <- tempfile()
#write.table(META[region %in% by.reg$region & !is.na(overall.P),.(chr=CHR,pos=BP,p=overall.P)],file=tmp.file,row.names=FALSE,quote=FALSE,sep="\t")
write.table(META[!is.na(overall.P),.(chr=CHR,pos=BP,p=overall.P)],file=tmp.file,row.names=FALSE,quote=FALSE,sep="\t")
gwas.DT<-load_gwas(tmp.file,maf.DT,ld.gr,n.cases=733,n.controls=9958-733)
by.reg2 <- gwas.DT[gwas.DT[,.I[which.min(p)],by=ld]$V1][p<1e-5,][order(chr,pos),]

## just include suggestive regions
gwas.filt.DT <- gwas.DT[ld %in% by.reg2$ld & chr!=6,]
feature.sets <- make_pchic('__REDACTED__/DATA/rCOGS/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab',biotype.filter='protein_coding')
all.genes <- lapply(feature.sets,function(g) unique(g$ensg)) %>% do.call('c',.) %>% unique
feature.sets[['VProm']] <- make_vprom('__REDACTED__/DATA/rCOGS/Digest_Human_HindIII.tab','__REDACTED__/DATA/rCOGS/HindIII_baits_e75.tab',all.genes)
digest.gr <- load_digest('__REDACTED__/DATA/rCOGS/Digest_Human_HindIII.tab')

res <- lapply(split(gwas.filt.DT,gwas.filt.DT$ld),function(gs){
  tmp <- compute_cogs(gs,csnps.gr,digest.gr,feature.sets)
  tmp[,ld:=unique(gs$ld)]
}) %>% rbindlist

library(biomaRt)
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','entrezgene'), filters = 'ensembl_gene_id', values =res$ensg, mart =ensembl) %>% data.table
DT <- merge(genedesc,res,by.x='ensembl_gene_id',by.y='ensg')
cand<-DT[cogs>0.05,]
cand[,cogs:=signif(cogs,digits=2)]
cand <- cand[,list(cand_genes=paste(external_gene_name,sep=',',collapse=',')),by=ld]
out <- merge(out,cand,by.x='region',by.y='ld',all.x=TRUE)[,region:=NULL]
out <- out[order(p.value),]
library(xlsx)
write.xlsx(out, file="__REDACTED__/pid/RESULTS/meta_gwas_hits.xlsx", sheetName="meta")
