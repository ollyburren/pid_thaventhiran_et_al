library(data.table)
library(optparse)

## get a list of all genes that we need to process

option_list = list(
        make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file", metavar="character"),
        make_option(c("-gl", "--genelist"), type="character", default=10,
                    help="Genelist", metavar="character")
        )

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)
if (is.null(args$out)){
        print_help(opt_parser)
        stop("At least one argument must be supplied (output file)", call.=FALSE)
}

print(args)


#bcftools_bin<-'~/bin/bcftools-1.4/bcftools'
bcftools_bin <- '/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools'

## format of this file is that each variant has much data associated with it.
## these routines allow the deconstruction of the annotation data

processVEP<-function(v,vep.header){
	ret<-do.call('rbind',lapply(1:length(v),function(x){
		tmp<-strsplit(v[x],'|',fixed=TRUE)[[1]]
		m<-matrix(tmp,ncol=length(vep.header)-1,byrow=TRUE)
		cbind(m,x)
	}))
	ret<-data.table(ret)
	setnames(ret,c(head(vep.header,-1),'row'))
	ret$row<-as.numeric(ret$row)
	ret
}


createInfoDT<-function(str){
	blah<-lapply(strsplit(str,";",fixed=TRUE),function(x) {
		tmp<-strsplit(x,"=",fixed=TRUE)
		#ret<-lapply(tmp,'[[',2)
    ret<-lapply(tmp,'[',2)
		#names(ret)<-sapply(tmp,'[[',1)
    names(ret)<-sapply(tmp,'[',1)
		#ret<-data.table(data.frame(ret))
		ret
	})
	## get a list of possible headers
	pos.header<-unique(c(do.call('c',lapply(blah,names)),'GNOMAD_AF'))
	tlist<-list()
	for(n in pos.header){
		tmp<-sapply(blah,'[[',n)
		if(class(tmp)=='list'){
			tmp<-sapply(tmp,function(x){
				if(length(x)==0){
					return(NA)
				}
				return(x)
			})
			## annoyingly the odd variant is not annotated
			if(n != 'ANN')
				tmp<-as.numeric(tmp)
		}
		tlist[[n]]<-tmp
	}
	#tlist[['ANN']]<-NULL
	tlist<-data.table(data.frame(tlist,stringsAsFactors=FALSE))
	tlist$row<-1:nrow(tlist)
	tlist
}

robustDTfread<-function(cmd){
	tryCatch(as.data.frame(fread(cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE)),warning=function(w) {print(sprintf("Warning=%s fread CMD=%s",w,cmd))},error=function(e){print(sprintf("Error=%s fread CMD=%s",e,cmd));return(NA)})
}

robustDTfread<-function(cmd){
	tryCatch(as.data.frame(fread(cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE)),error=function(e){print(sprintf("Error=%s fread CMD=%s",e,cmd));return(NA)})
}
## this is prototype command for bcftools
## finish converting to BCFTools
## look at genes outside of the PID list
## reread the bioarxiv paper on reg variation
## reread BeviMed paper.

getPossDamSNPs<-function(region,individuals,...){
	ar<-list(...)
	chr<-sub("([^:]+):.*","\\1",region)
	start<-as.numeric(sub("[^:]+:([^\\-]+)\\-.*","\\1",region))
	end<-as.numeric(sub("[^:]+:[^\\-]+\\-(.*)","\\1",region))
	ind<-paste(individuals,sep=',',collapse=',')
	#args<-list(...)
	## load in the merged vcf file.
	vcf.dir<-'__REDACTED__/data/release/latest/merged-vcf/no_hgmd/gnomad/'
	data.dir<-'__REDACTED__/pid/'
	fname<-'chr%s_agg3_dedup_vep_gnomad.bcf'
	vf<-file.path(vcf.dir,sprintf(fname,chr))
	header_cmd<-sprintf("%s view -s %s --header-only %s",bcftools_bin,ind,vf)
	my.pipe<-pipe(header_cmd)
	header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
	close(my.pipe)
	cnames<-unlist(strsplit(header,"\t"))
	## should be able to do in one pass but this is so quick !
	my.pipe<-pipe(sprintf("%s | egrep -e '%s' | sed -e s'%s'",header_cmd,'^##INFO=<ID=ANN','/^.*Format: \\([^\\][^\\]*\\)\\.*/\\1/'))
	## variant annotation header is another thing to capture
	vep<-scan(my.pipe,what=character(),sep="\n",quiet=TRUE)
	close(my.pipe)
	vep<-gsub('"','',vep)
	vep.header<-strsplit(vep,'\\|')[[1]]
	if("egrep" %in% names(ar)){
		body_cmd<-sprintf("%s view %s -r %s -s %s -Ov -H | egrep -e '%s'",bcftools_bin,vf,region,ind,ar[['egrep']])
	}else{
		body_cmd<-sprintf("%s view %s -r %s -s %s -Ov -H",bcftools_bin,vf,region,ind)
	}
	message(body_cmd)
	#tmp<-as.data.frame(fread(body_cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE))
	tmp<-robustDTfread(body_cmd)
	if(any(is.na(tmp)))
		return(NA)
	colnames(tmp)<-cnames
	info<-tmp[,1:9]
	info.dt<-createInfoDT(info$INFO)
	## possibly only interested in GNOMAD_AF to start with
	obj<-list(map=info,info=info.dt)
	#get VEP stuff - note that we get multiple entries per SNP.
	vep<-processVEP(obj$info$ANN,vep.header)
	## get SNPs that have a high chance to mess up the protein
	interesting.vep<-vep[grep(ar$egrep,vep$Consequence,ignore.case=TRUE),c('SYMBOL','Gene','BIOTYPE','Existing_variation','Consequence','CADD_PHRED','SIFT','Amino_acids','PolyPhen','EXON','row'),with=FALSE]
	if(nrow(interesting.vep)==0)
		return(NA)
	## get info on GNOMAD allele freq
	#c.idx<-as.numeric(unique(subset(vep,IMPACT %in% c('HIGH','MODERATE') & CANONICAL=='YES' & BIOTYPE=='protein_coding')$row))
	## use snpMatrix to find out which ones have variation
	#csum<-col.summary(obj$gt[,c.idx])
	gaf<-data.frame(id=1:nrow(tmp),GNOMAD_AF=info.dt$GNOMAD_AF)
	merge(gaf,interesting.vep,by.x='id',by.y='row')
}

in.dir<-'__REDACTED__/pid/VCF2/exonic_only/no_monomorphs/lof/gene_events_AUGUST_2018'
fpat<-'%s_lof.bcf'
parseVCF<-function(ensg){
  message(sprintf("Processing %s",ensg))
	vf<-file.path(in.dir,sprintf(fpat,ensg))
	if(!file.exists(vf))
		return(NA)
	header_cmd<-sprintf("%s view --header-only %s",bcftools_bin,vf)
	my.pipe<-pipe(header_cmd)
	header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
	close(my.pipe)
	cnames<-unlist(strsplit(header,"\t"))
	## should be able to do in one pass but this is so quick !
	my.pipe<-pipe(sprintf("%s | egrep -e '%s' | sed -e s'%s'",header_cmd,'^##INFO=<ID=ANN','/^.*Format: \\([^\\][^\\]*\\)\\.*/\\1/'))
	## variant annotation header is another thing to capture
	vep<-scan(my.pipe,what=character(),sep="\n",quiet=TRUE)
	close(my.pipe)
	vep<-gsub('"','',vep)
	vep.header<-strsplit(vep,'\\|')[[1]]
	body_cmd<-sprintf("%s view %s -Ov -H",bcftools_bin,vf)
	#message(body_cmd)
	tmp<-robustDTfread(body_cmd)
	if(any(is.na(tmp)))
		return(NA)
	colnames(tmp)<-cnames
	info<-tmp[,1:9]
	info.dt<-createInfoDT(info$INFO)

	vep<-processVEP(info.dt$ANN,vep.header)
	all.info<-merge(info.dt,vep,by.x='row',by.y='row')[,.(AF_WGS10K,GNOMAD_AF,minOPR,SYMBOL,Gene,BIOTYPE,Existing_variation,Consequence,CADD_PHRED,Amino_acids,row)]
	## next remove those that cannot be loss of function
	all.info<-all.info[grep("HGMD_CLASS=DM[^?]|splice|frame|start|stop|missense",all.info$Consequence),]
	inds<-names(tmp)[10:length(names(tmp))]
	stub<-rbindlist(lapply(inds,function(i){
			igt<-tmp[[i]]
			igt<-sub("0\\/0.*","1",igt)
			igt<-sub("(0\\/1).*|(1\\/0).*","2",igt)
			igt<-sub("1\\/1.*","3",igt)
			igt<-sub("[0-9]\\/[0-9]","0",igt)
			data.table(ind=i,gt=igt,row=(1:length(igt)))
	}))
	M <- merge(stub,all.info,by.x='row',by.y='row',allow.cartesian=TRUE)
  data.table(cbind(info[M$row,c(1,2,4,5)],M))
}



genes<-scan(args$genelist,character())

ar<-lapply(genes,parseVCF)
## remove those that did not have any consequence or whatever
ar<-rbindlist(ar[sapply(ar,function(x) any(class(x)=='data.table'))])
saveRDS(ar,file=args$out)
