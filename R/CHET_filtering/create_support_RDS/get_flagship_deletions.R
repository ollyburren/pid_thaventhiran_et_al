library(data.table)
library(GenomicRanges)


## get lenient CNV's from Keren

#DATA.DIR <- '/scratch/kjn29/open/SV_pipeline/flagship_13k/dels_v2/dataset/lenient_nonGEL'
FLAG_FILE <- '__REDACTED__/analysis/pid/flagship-del/flagship-deletions.txt'
OUT.DIR<-'__REDACTED__/pid/flagship_dels/'
FREQ<-0.03
MIN.WIDTH <- 10

## get list of ID's to include in analysis1

#pid.samples <- fread('__REDACTED__/pid/support/PID_MainTable_PipelineStatus_cohort20170614_forCohortPaper_29_11_2017.tsv')
#setnames(pid.samples,make.names(names(pid.samples)))
#pid.samples <- pid.samples[!is.na(Is.Affected),][['WGS.ID']]
#pid.samples <- pid.samples[PIDnoVasc_Index==1,][['WGS.ID']]
#write(pid.samples,file='__REDACTED__/pid/support/PID_895_26_06_2018.txt')
pid.samples <- scan('__REDACTED__/pid/support/PID_895_26_06_2018.txt',character())
cnv <- fread(FLAG_FILE)[,uid:=1:.N]
## compute deletion totals
## first recode genotypes 0/1 1/0 = 1 and 1/1 =2
cnv[gt_Manta %in% c('1/0','0/1'),gtm:=1]
cnv[gt_Manta == '1/1',gtm:=2]
cnv[gt_Canvas %in% c('1/0','0/1'),gtc:=1]
cnv[gt_Canvas == '1/1',gtc:=2]
cnv[,gt:=max(gtm,gtc,na.rm=TRUE),by=uid]
total_chr <- length(unique(cnv$id)) * 2
deletion<-cnv[,list(AF=sum(gt)/total_chr),by=deletion]
M <- merge(cnv,deletion,by.x='deletion',by.y='deletion')


cnv.f <- M[id %in% pid.samples,]
cnv.f[,c('START','END'):=list(min(start_Manta,start_Canvas,na.rm=TRUE),max(end_Manta,end_Canvas,na.rm=TRUE)),by=uid]
cnv.f<-cnv.f[,.(BRIDGE_ID=id,CHROM=chr,START,END,GT=gt,AF=AF,uid,cluster=deletion)]
cnv.f[,MAF:=ifelse(AF>0.5,1-AF,AF)]
gr <- with(cnv.f,GRanges(seqnames=Rle(CHROM),ranges=IRanges(start=START,end=END),GT=GT,uid=uid,bridge_id=BRIDGE_ID,MAF=MAF,cluster=cluster))
gr <- gr[width(gr)>MIN.WIDTH & gr$MAF<FREQ,]



saveRDS(gr,file=file.path(OUT.DIR,'pid_flag_0.03.RDS'))

## no longer filter on allele frequency as this dataset does not have that information
