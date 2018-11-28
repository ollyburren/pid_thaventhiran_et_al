library(data.table)
library(GenomicRanges)


## load in promoter capture data
CHIC.THRESH<-5
data.dir<-'__REDACTED__/pid'
mfiles<-list.files(path=file.path(data.dir,'hnisz'),pattern='*.csv',full.names = TRUE)
## select CD files first off
cd.files<-mfiles[grep('CD',basename(mfiles))]
desc<-vector(mode="character",length=length(cd.files))
desc[grep('CD34',basename(cd.files))]<-'Endothelial_precursors'
## CD4
desc[grep('CD4_Memory',basename(cd.files))]<-'Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated'
desc[grep('BI_CD4p_CD225int_CD127p_Tmem',basename(cd.files))]<-'Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated'
desc[grep('BI_CD4p_CD25-_CD45ROp_Memory',basename(cd.files))]<-'Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated'
desc[grep('CD4_Naive',basename(cd.files))]<-'Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated,Naive_CD4'
desc[grep('BI_CD4p_CD25-_CD45RAp_Naive',basename(cd.files))]<-'Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated,Naive_CD4'
## stimulated T-Cells
desc[grep('BI_CD4p_CD25-_Il17-_PMAstim_Th',basename(cd.files))]<-'Total_CD4_Activated'
desc[grep('BI_CD4p_CD25-_Il17p_PMAstim_Th17',basename(cd.files))]<-'Total_CD4_Activated'
## CD8
desc[grep('BI_CD8_Memory',basename(cd.files))]<-'Total_CD8,Naive_CD8'
desc[grep('CD8_primary',basename(cd.files))]<-'Total_CD8,Naive_CD8'
desc[grep('BI_CD8_Naive',basename(cd.files))]<-'Naive_CD8'
desc[grep('CD14',basename(cd.files))]<-'Macrophages_M0,Macrophages_M1,Macrophages_M2,Monocytes'
desc[grep('CD19',basename(cd.files))]<-'Naive_B,Total_B'
## pro thymocytes
desc[grep('CD3\\.',basename(cd.files))]<-'Naive_B,Total_B,Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated,Naive_CD4,Total_CD8,Naive_CD8'
hnisz<-data.table(f=cd.files,pchic.tissue=desc)
##remove those that don't have a mapping
hnisz<-subset(hnisz,pchic.tissue!='')

pchic.cutoff<-5
pchic<-fread(file.path(data.dir,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab"))
pchic.f<-subset(pchic,biotype=='protein_coding')
library(reshape2)
pchicm<-data.table::melt(pchic.f,id.vars=c('ensg','name','baitID','oeID','oeChr','oeStart','oeEnd'),measure.vars=names(pchic.f)[16:32])

## annotate with whether overlaps a superenhancer


p.gr<-with(unique(pchicm,by='oeID'),GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,oeEnd),oeID=oeID))

#p.gr<-with(pchicm,GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,oeEnd),id=1:nrow(pchicm)))

hnisz_ol<-lapply(cd.files,function(f){
  message(sprintf("Processing %s",f))
  DT<-fread(f,skip=1L)
	setnames(DT,c('enh_id','chr','start','end','gene','enh_rank','super','chip.den','read.den'))
	## first off retrieve SE that are found in matched tissues.
	se.gr<-with(DT[DT$super==1,],GRanges(seqnames=Rle(sub('^chr','',chr)),ranges=IRanges(start=start,end=end),id=enh_id))
  enh.gr<-with(DT[DT$super==0,],GRanges(seqnames=Rle(sub('^chr','',chr)),ranges=IRanges(start=start,end=end),id=enh_id))
  se.ol<-as.matrix(findOverlaps(p.gr,se.gr))
  enh.ol<-as.matrix(findOverlaps(p.gr,enh.gr))
  full<-integer(length=length(p.gr))
  full[unique(enh.ol[,1])]<-1
  full[unique(se.ol[,1])]<-2
  full
})


names(hnisz_ol)<-sub("\\.csv","",basename(cd.files))
all.hnisz.pchic<-do.call('cbind',hnisz_ol)
## remove those that have no overlaps
no.overlaps.idx<-which(rowSums(all.hnisz.pchic)==0)
all.hnisz<-cbind(all.hnisz.pchic,oeID=p.gr$oeID)
## get rid of rows with no overlaps
filt<-all.hnisz[-no.overlaps.idx,]


mergeCols<-function(tomerge,cname,filt){
  tmp<-integer(length=nrow(filt))
  idx<-which(apply(filt[,tomerge],1,function(x) any(x==1)))
  tmp[idx]<-1
  idx<-which(apply(filt[,tomerge],1,function(x) any(x==2)))
  tmp[idx]<-2
  filt<-filt[,!colnames(filt) %in% tomerge]
  cnames<-colnames(filt)
  filt<-cbind(filt,tmp)
  colnames(filt)<-c(cnames,cname)
  filt
}

merge.list<-list(
  CD34=c('BI_CD34_Primary_RO01480','BI_CD34_Primary_RO01536','BI_CD34_Primary_RO01549','CD34_adult'),
  CD8_Naive=c('BI_CD8_Naive_8pool','BI_CD8_Naive_7pool'),
  CD4_Naive=c('BI_CD4_Naive_Primary_7pool','BI_CD4_Naive_Primary_8pool'),
  CD4_Memory=c('BI_CD4_Memory_Primary_7pool','BI_CD4_Memory_Primary_8pool')
      )

for(n in names(merge.list)){
  message(sprintf("Processing %s",n))
  filt<-mergeCols(merge.list[[n]],n,filt)
}

filt<-data.table(filt)
setkey(filt,oeID)

pchicm.int<-subset(pchicm,value>5)
## add in enhancer/superenhancer overlaps
setkey(pchicm.int,oeID)
setnames(pchicm,'variable','tissue')
setnames(pchicm,'value','chicago.score')
final.merge<-merge(pchicm.int,filt,by.x='oeID',by.y='oeID',all.x=TRUE)
## next filter so we have four files.
m<-as.matrix(final.merge[,10:26])
#1 all interactions for any tissue regardless of overlap with hnisz
#2 all interactions for se
se.idx<-which(apply(m,1,function(x) any(x==2)))
#3 all interactions for enhancer only (i.e not convered above)
enh.idx<-which(apply(m,1,function(x) any(x==1)))
enh.idx<-setdiff(enh.idx,se.idx)
#4 all interactions for se or enhancer
se_enh.idx<-which(rowSums(m,na.rm=TRUE)>0)
df<-DataFrame(data.frame(final.merge[,!names(final.merge) %in% c('oeChr','oeStart','oeEnd'),with=FALSE]))

cGR<-function(index){
  tmp<-final.merge[index,]
  gr<-with(tmp,GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,end=oeEnd)))
  mcols(gr)<-df[index,]
  gr
}

saveRDS(cGR(se.idx),file="__REDACTED__/pid/ANNOTATIONS/pchic_hnisz/se.gr.RDS")
saveRDS(cGR(enh.idx),file="__REDACTED__/pid/ANNOTATIONS/pchic_hnisz/enh.gr.RDS")
saveRDS(cGR(se_enh.idx),file="__REDACTED__/pid/ANNOTATIONS/pchic_hnisz/se_enh.gr.RDS")
