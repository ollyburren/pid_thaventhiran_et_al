library(data.table)
library(magrittr)
library(GenomicRanges)
library(rCOGS)

## for any locus with a suggestive gwas hit prioritise genes
## what are the loci for meta analysis ?

meta <- fread("__REDACTED__/analysis/pid/GWAS/lango-allen-li-ellinghaus-meta_analysis.tab")
meta<-meta[SNP!='rs7785984',]
## load in immunochip regions
ic.regions <- fread("__REDACTED__/DATA/JAVIERRE_ICHIP/support/dense.ic.regions.guessfm.bed")
setnames(ic.regions,c('chr','start','end','name'))

ic.gr <- with(ic.regions,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end)))
meta.gr <- with(meta,GRanges(seqnames=Rle(CHR),ranges=IRanges(start=BP,width=1L)))
ol <- findOverlaps(meta.gr,ic.gr) %>% as.matrix
meta[ol[,1],region:=ol[,2]]

regions <- meta[,list(min.p=min(overall.P)),by='region'][min.p<1e-5,]$region
## create a custom region file

sugg.regions <- ic.regions[regions,]

write.table(sugg.regions[,.(chr,start,end)],file="__REDACTED__/pid/abdef_region_cogs/sugg_regions_october.tab",quote=FALSE,row.names=FALSE,sep="\t")

## process ichip data

ic.files <- list.files(path="__REDACTED__/DATA/JAVIERRE_ICHIP/gwas/",pattern="*.gz$",full.names=TRUE)
for(f in ic.files){
  #f <- ic.files[1]
  trait <- gsub("hg19_gwas_ic_([^4]+).*","\\1",basename(f)) %>% gsub("\\_$","",.)
  DT <- fread(sprintf("zcat %s",f))
  setnames(DT,c('chr','start','pos','info','p'))
  DT[,c('rs','cases','controls'):=tstrsplit(info,':')]
  n.cases <- DT$cases[1]
  n.controls <- DT$controls[1]
  fname <- sprintf("%s:%s:%s.tab",trait,n.cases,n.controls) %>% file.path('__REDACTED__/pid/abdef_region_cogs/gwas',.)
  write.table(DT[,.(chr,pos,p=signif(p,digits=3))],file=fname,quote=FALSE,row.names=FALSE,sep="\t")
}

## next we can get genes priorised for a given trait in a region

ld.gr <- load_ld_regions('__REDACTED__/pid/abdef_region_cogs/sugg_regions_october.tab')
maf.DT <- load_ref_maf('__REDACTED__/rCOGS/inst/extdata/uk10k_reference_maf.tab')
feature.sets <- make_pchic('__REDACTED__/rCOGS/inst/extdata/javierre_pchic_ensembl75.tab',biotype.filter='protein_coding')
all.genes <- lapply(feature.sets,function(g) unique(g$ensg)) %>% do.call('c',.) %>% unique
feature.sets[['VProm']] <- make_vprom('__REDACTED__/rCOGS/inst/extdata/hindIII_digest_grch37.tab','__REDACTED__/rCOGS/inst/extdata/javierre_design_ensembl75.tab',all.genes)
csnps.gr <- make_csnps('__REDACTED__/rCOGS/inst/extdata/coding_snps_vep_ensembl75.tab')
digest.gr <- load_digest('__REDACTED__/rCOGS/inst/extdata/hindIII_digest_grch37.tab')

## next we load in each gwas and attempt gene prioritisation
gwas_files <- list.files(path='__REDACTED__/pid/abdef_region_cogs/gwas',full.names=TRUE)
prior.genes <- lapply(gwas_files,function(f){
  message(f)
  vars <- basename(f) %>% gsub("\\.tab","",.) %>% strsplit(.,':') %>% unlist
  trait <- vars[1]
  n.cases <- vars[2] %>% as.numeric
  n.ctrls <- vars[3] %>% as.numeric
  gwas.DT <- load_gwas(f,maf.DT,ld.gr,n.cases=n.cases,n.controls=n.ctrls)
  out.DT <- compute_cogs(gwas.DT,csnps.gr,digest.gr,feature.sets)[,trait:=trait]
  out.DT
}) %>% rbindlist

filt <- prior.genes[cogs>0.5,]
filt.genes <- filt$ensg

supp.table <- prior.genes[ensg %in% filt.genes,]

## get gene information using biomaRt

library(biomaRt)

ensembl_archive="feb2014.archive.ensembl.org" # Final V37 Ensembl
ensembl <- useMart(host=ensembl_archive, biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
sig.genes<-getBM(c("ensembl_gene_id","gene_biotype","external_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position","strand"),
filters=c("ensembl_gene_id"),
values=list(ensembl_gene_id=filt$ensg), mart=ensembl) %>% data.table
sig.genes<-sig.genes[,.(ensembl_gene_id,hgnc_symbol,chromosome_name,start_position,end_position)]
M <- merge(prior.genes,sig.genes,by.x='ensg',by.y='ensembl_gene_id')
dat <- melt(M,id.vars=c('ensg','trait','hgnc_symbol','chromosome_name','start_position','end_position'),measure.vars='cogs') %>%
dcast(ensg + hgnc_symbol + chromosome_name + start_position + end_position ~ variable + trait)

## add in pLI scores

DT.pli <- fread("__REDACTED__/pid/pLI_scores.tab")

DT.out <- merge(dat,DT.pli,by.x='hgnc_symbol',by.y='gene')

## add in PID candidate genes

pid.2017<-scan("__REDACTED__/pid/IUIS_2017.csv",character())
library(xlsx)
DT.out[,iuid2017:=hgnc_symbol %in% pid.2017]
library(xlsx)
write.xlsx(DT.out, file="__REDACTED__/pid/RESULTS/abdef_region_prioritised_genes_october.xlsx", sheetName="Prioritised Genes")
## to create a supplemental table of results





## James likes this so work up into a nicer figure


library(magrittr)
library(data.table)
DT.out <- read.xlsx(file="__REDACTED__/pid/RESULTS/abdef_region_prioritised_genes_october.xlsx", sheetName="Prioritised Genes") %>% data.table
DT.out[,cogs_ibd_jostins:=NULL]
DT.out <- DT.out[,-1]
cnames <- names(DT.out)
midx <- grep("^cogs",cnames)
M <- melt(DT.out,id.vars=cnames[-midx],measure.vars=cnames[midx])
gene.factor <- unique(M[,.(hgnc_symbol,pLI)])[order(pLI,decreasing=FALSE),]$hgnc_symbol
M[,gene:=factor(hgnc_symbol,levels=gene.factor)]

## cluster genes by disease relatedness
hc.disease <- DT.out[,cnames[midx],with=FALSE] %>% as.matrix %>% t %>% dist %>% hclust
dlevels <- gsub("cogs\\_","",cnames[midx]) %>% gsub("([^\\_]+)\\_.*","\\1",.) %>% toupper
dlevels[hc.disease$order]
M[,trait:=gsub("cogs\\_","",variable) %>% gsub("([^\\_]+)\\_.*","\\1",.) %>% toupper %>% factor(.,levels=dlevels[hc.disease$order])]

M[,overall.score:=value * pLI]

## filter out traits with no COGS scores > 0.5

keep.trait <- M[,list(max(value)),by=trait][V1>0.5,]$trait

library(cowplot)
library(RColorBrewer)

myPalette <- colorRampPalette(c("blue","red"))
sc <- scale_fill_gradient2("Prioritisation\nScore",low="blue",high="red",limits=c(0, 1))
keep.genes <- M[,list(max(overall.score)),by=gene][V1>mean(M$overall.score),]$gene


plotM <- M[trait %in% keep.trait & hgnc_symbol %in% keep.genes,]

hilight.genes <- c('SOCS1','PTPN2','ETS1')

labs <- sapply(strsplit(as.character(plotM$gene), ", "), function(x) {
    parse(text = paste("italic('", x, "')", sep = ""))

})

pp<-ggplot(plotM,aes(x=trait,y=gene,fill=overall.score)) + geom_tile(col='black',size=0.1) + sc +
xlab("ImmunoChip Trait") + ylab("Gene") + scale_y_discrete(labels = labs, breaks = plotM$gene) +
theme(axis.line=element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

save_plot("~/tmp/prioritised_genes_october.pdf",pp)

suppM <- M[trait %in% keep.trait,][,c('value','variable'):=list(NULL,NULL)]
suppM <- melt(suppM,id.vars=c('hgnc_symbol' , 'ensg' , 'chromosome_name' , 'start_position' , 'end_position' , 'iuid2017' , 'gene','trait','pLI'),measured.vars='overall.score')
suppM[,value:=signif(value,digits=3)]
supp.DT <- dcast(suppM, hgnc_symbol + ensg + chromosome_name + start_position + end_position + iuid2017 + gene ~ variable + trait)
supp.DT <- supp.DT[order(chromosome_name,start_position),]
write.xlsx(supp.DT, file="__REDACTED__/pid/RESULTS/abdef_region_prioritised_genes_october_supplementary_table.xlsx", sheetName="Prioritised Genes")
