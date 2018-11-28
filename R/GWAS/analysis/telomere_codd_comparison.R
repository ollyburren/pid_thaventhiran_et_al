DT <- fread("__REDACTED__/telomere/summary/telomeres_length_all.assoc.linear")
DT$SNP<-NULL

DT[,c('CHR','BP','BETA','STAT','P'):=list(as.numeric(CHR),as.numeric(BP),as.numeric(BETA),as.numeric(STAT),as.numeric(P))]
DT[,SE:=BETA/STAT]

codd <- fread("~/tmp/CODD.tsv")

m <- merge(codd,DT,by.x=c('CHR','POS_37'),by.y=c('CHR','BP'))

tab <- m[,.(GENE,SNP,CHR,BP=POS_37,codd_A1=EFFECT_ALLELE,codd_A2=OTHER_ALLELE,codd_BETA=BETA.x,codd_SE=SE.x,codd_P=P.x,A1,BETA=BETA.y,SE=SE.y,P=P.y)]
tab[codd_A2==A1,flip:=TRUE]
tab[codd_A1==A1,flip:=FALSE]
tab[flip==TRUE,BETA:=BETA*-1]

tab <- tab[,.(SNP,chr=CHR,position=BP,effect_allele=codd_A2,codd_beta=codd_BETA,codd_P,beta=BETA,P,codd_gene=GENE)]
tab <- tab[order(codd_P),]
library(xlsx)
write.xlsx(tab, file="__REDACTED__/pid/RESULTS/telomere.xlsx", sheetName="codd_comparison")
