library(data.table)
library(GenomicRanges)

## ONLY NEED TO DO ONCE FOR A GIVEN SET OF ANNOTATIONS
if(FALSE){
  OUT_DIR<-'/__REDACTED__/pid/ANNOTATIONS'
  g<-fread("/__REDACTED__/pid/support/gencode.v26lift37.annotation.gtf")
  g.exon<-subset(g,V3=='exon'& grepl("protein_coding",V9))
  info<-do.call('rbind',lapply(strsplit(g.exon$V9,";"),function(f){
          gid<-sub('gene_id \"(.*)\"','\\1',f[grep("gene_id",f)])
          tid<-sub(' transcript_id \"(.*)\"','\\1',f[grep("transcript_id",f)])
          gt<-sub(' gene_type \"(.*)\"','\\1',f[grep("gene_type",f)])
          tt<-sub(' transcript_type \"(.*)\"','\\1',f[grep("transcript_type",f)])
          gn<-sub(' gene_name \"(.*)\"','\\1',f[grep("gene_name",f)])
          return(cbind(gid,tid,gt,tt,gn))
        }))

        final<-cbind(g.exon,info)
        final.pc<-subset(final,gt=='protein_coding' & tt=='protein_coding')
        final.pc$ensg<-sub("(ENSG[0-9]+)\\..*","\\1",final.pc$gid)
        final.pc[,uid:=paste(V1,V4,ensg,sep=':')]
        ufinal.pc<-unique(final.pc,by='uid')
        e.gr<-with(ufinal.pc,GRanges(seqnames=Rle(sub("chr","",V1,fixed=TRUE)),ranges=IRanges(start=V4,end=V5),strand=V7,ensg=ensg,name=gn))
        saveRDS(e.gr,file=file.path(OUT_DIR,'exons.RDS'))

        ## to filter VCF file  create region files
        library(rtracklayer)
        exon_files<-lapply(split(e.gr,seqnames(e.gr)),function(gr){
          tmp<-reduce(gr)
          chr<-unique(seqnames(tmp))
          message(sprintf("Processing %s",chr))
          ofile<-sprintf("chr%s_gencode.v26lift37.bed",chr)
          ofile<-file.path(OUT_DIR,'exons_by_chr',ofile)
          export.bed(gr,con=ofile)
          ofile
        })
}

### CODE BELOW HERE GETS ALL SNPS OVERLAPPING EXONS FROM GENCODE
exon_files<-list.files(path='/__REDACTED__/pid/ANNOTATIONS/exons_by_chr',pattern='*.bed',full.names=TRUE)
exon_files <- split(exon_files,gsub("^chr([^\\_]+)\\_.*","\\1",basename(exon_files)))
## create list of suitable bcftool calls
bcftools_bin<-'~/bin/bcftools-1.4/bcftools'
pid_samples<-'/__REDACTED__/pid/support/PID_895_29_11_2017.txt'
cmd<-"%s view %s -R %s -S %s -Ob > %s"
# get a list of VCF files
vcf_dir<-vcf.dir<-'/scratch/WGS10K/data/release/20170614G-A/merged-vcf/no_hgmd'
fname<-'chr%s_agg3_dedup_vep.bcf'
odir<-'/__REDACTED__/pid/VCF/exonic_only/'



all.cmds<-sapply(names(exon_files),function(n){
  if(n %in% c('M','Y')){
    return(NA)
  }
  efile<-exon_files[[n]]
  vcf<-file.path(vcf_dir,sprintf(fname,n))
  ofile<-sprintf("chr%s_gencode.v26lift37.bcf",n)
  ofile<-file.path(odir,ofile)
  sprintf(cmd,bcftools_bin,vcf,efile,pid_samples,ofile)
})

all.cmds<-all.cmds[!is.na(all.cmds)]
write(all.cmds,file=file.path("/__REDACTED__/pid/VCF",'exons.txt'))

## run on the QUEUE.
# Finally fiter out those rows that are monomorphic in the PID cohort.
#INDIR=/__REDACTED__/pid/VCF/exonic_only
#OUTDIR=/__REDACTED__/pid/VCF/exonic_only/no_monomorphs
#for i in `\ls $INDIR`;do
#  ofile=$(basename $i)
#  echo "~/bin/bcftools-1.4/bcftools view $INDIR/$i -i 'AC!=0' -Ob -o $OUTDIR/$ofile"
#done

## finally on these get the LOF
#INDIR=/__REDACTED__/pid/VCF/exonic_only/no_monomorphs
#OUTDIR=/__REDACTED__/pid/VCF/exonic_only/no_monomorphs/lof/
#for i in `\ls $INDIR`;do
#  ofile=$(basename $i)
#  echo "~/bin/bcftools-1.4/bcftools view $INDIR/$i | grep -E '^#|HGMD_CLASS=DM[^?]|splice|frame|start|stop|missense' | ~/bin/bcftools-1.4/bcftools convert - -Ob -o $OUTDIR/$ofile"
#done

## note that the resultant files are not sorted and may contain duplicate lines (not sure) - they need to be sorted and indexed
#INDIR=/__REDACTED__/pid/VCF/exonic_only/no_monomorphs/lof/
#cd $INDIR
#BCFTOOLS=__REDACTED__/bin/bcftools-1.4/bcftools
#for i in `\ls *.bcf`;do
#  echo "Processing $i"
#  # dump header to a file
#  $BCFTOOLS view $i -h > tmp.h
  # now dump bcf file to vcf sort and add header
#  $BCFTOOLS view $i -H | sort -k2n >> tmp.h
#  $BCFTOOLS view tmp.h -Ob -o ./sorted/$i
#  $BCFTOOLS index ./sorted/$i
#  rm tmp.h
#done
