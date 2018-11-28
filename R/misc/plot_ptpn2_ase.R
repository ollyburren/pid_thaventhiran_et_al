library(cowplot)

## raw data
DT <- data.table(allele=c('min','maj','min','maj'),atype=c('gDNA','gDNA','cDNA','cDNA'),perc_col=c(0.44117647,0.55882353,0.02857143,0.97142857))

DT[,c('allele','atype'):=list(factor(allele,levels=c('min','maj'),labels=c('G(Minor) - FAM','A(Major) - VIC')),factor(atype,levels=c('gDNA','cDNA'))   )]

pp<-ggplot(DT,aes(x=allele,y=perc_col,fill=allele)) + geom_bar(stat="identity") +
facet_wrap(~atype) + xlab("PTPN2 rs2847297") + ylab("% of Colonies") + geom_hline(yintercept=0.5,lty=2,alpha=0.5) +
theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_fill_manual(name="Allele",values=c('blue','red')) +
theme(legend.position = c(0.01, 0.9))
save_plot("~/tmp/ptpn2_ase.pdf",pp,base_width=4)
