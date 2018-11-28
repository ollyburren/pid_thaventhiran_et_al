library(data.table)
library(GenomicRanges)
library(magrittr)


## load in promoter capture data
CHIC.THRESH<-5


## Dennis super enhancer retrieval and integrations with PCHIC

data.dir <- '__REDACTED__/pid'
SE.data.dir <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers'

se.files <- list.files(path=SE.data.dir,pattern='*.bed',full.names=TRUE)

## create two lists one of narrow and one of broad peaks

btl <- split(se.files,gsub(".*\\_([^\\_]+)\\_peaks.*","\\1",basename(se.files)))

## for this analysis where we want linked SE use broad peaks

bp<-btl[['broad']]

pchic<-fread(file.path(data.dir,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab"))
pchic.f<-subset(pchic,biotype=='protein_coding')

## assign each pchic to a SE type

desc<-list()
desc['Monocytes'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.Monocytes_merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Macrophages_M0'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.Macrophage_merged_1e-5_broad_peaks.SuperEnhancers.bed'
## make assumption that this is MCF stim
desc['Macrophages_M1'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.InflammatoryMacrophage_merged_1e-5_broad_peaks.SuperEnhancers.bed'
## this is IL13 stim
desc['Macrophages_M2'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.AlternativelyActivatedMacrophage_merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Neutrophils'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.Neutrophils_Merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Megakaryocytes'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.MK_merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Endothelial_precursors'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.BOEC_progenitors_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Erythroblasts'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.EB_Merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Foetal_thymus'] <- ''
desc['Naive_CD4'] <- ''
desc['Total_CD4_MF'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.CD4_CentralMemoryTcell_merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Total_CD4_Activated'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.CD4Proliferating_Merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Total_CD4_NonActivated'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.CD4NonActive_Merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Naive_CD8'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.CentralmemoryCD8positiveAlphaBetaTcell_merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Total_CD8'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.CD8positiveAlphaBetaTcell_merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Naive_B'] <- '__REDACTED__/cbrcref/SuperEnhancerCalling/FinalFiles/SuperEnhancers/Peaks.NaiveBcell_merged_1e-5_broad_peaks.SuperEnhancers.bed'
desc['Total_B'] <- ''

desc<-desc[sapply(desc,function(x) x!='')]

## get cell types and see if we can link to PCHIC

### next we can read each of these in


pchicm<-data.table::melt(pchic.f,id.vars=c('ensg','name','baitID','oeID','oeChr','oeStart','oeEnd'),measure.vars=names(pchic.f)[16:32])
pchicm<-subset(pchicm,value>CHIC.THRESH)

## annotate with whether overlaps a superenhancer
p.grl<-with(unique(pchicm,by=c('seName','variable')),GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,oeEnd),baitID=baitID,oeID=oeID,ct=variable,ensg=ensg,gname=name,pchic.score=value)) %>% split(.,.$ct)

library(rtracklayer)
pchse<-lapply(seq_along(desc),function(i){
  pchic.tissue<-names(desc)[i]
  se.fname<-desc[[i]]
  gr<-import.bed(con=se.fname)
  seqlevels(gr)<-gsub('chr','',seqlevels(gr))
  tmp<-mergeByOverlaps(gr,p.grl[[pchic.tissue]])
  df<-data.table(as.data.frame(tmp))[,c(1:4,6,7,10:13,15:20),with=FALSE]
  setnames(df,c('seChr','seStart','seEnd','seWidth','seName','seScore','pchicChr','pchicStart','pchicEnd','pchicWidth','baitID','oeID','pchicTissue','ensg','name','pchicScore'))

})

pchse<-rbindlist(pchse)

## useful for plotting at a later stage

se.gr<-with(unique(pchse,by=c('seName','pchicTissue')),GRanges(seqnames=Rle(seChr),ranges=IRanges(start=seStart,end=seEnd),uid=pchicTissue))
saveRDS(se.gr,file='__REDACTED__/pid/bp_se.RDS')


res<-data.table::dcast(pchse,pchicChr+pchicStart+pchicEnd+ensg+name+baitID+oeID~pchicTissue)
res$variable<-'match'
res$value<-6
## convert to Granges objects
out.gr<-with(res,GRanges(seqnames=Rle(pchicChr),ranges=IRanges(start=pchicStart,end=pchicEnd)))
## get meta data as DataFrame
mcol.df<-cbind(res[,.(oeID,ensg,name,baitID,variable,value)],res[,c(8:21),with=FALSE])

mcols(out.gr)<-DataFrame(as.data.frame(mcol.df))

saveRDS(out.gr,file="__REDACTED__/pid/ANNOTATIONS/pchic_blueprint/se.gr.RDS")
