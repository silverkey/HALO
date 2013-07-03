library('biomaRt')

args = commandArgs(trailingOnly=T)

# For ensembl genes use the following mart:
# ENSEMBL_MART_ENSEMBL
mart = args[1]
sp1 = args[2]
sp2 = args[3]

db = paste(sp1,'gene_ensembl',sep='_')
hom = paste(sp2,'homolog_ensembl_gene',sep='_')
inter = paste(sp2,'inter_paralog_ensembl_gene',sep='_')
outfile = paste(sp1,sp2,'homologs.txt',sep='_')

ensembl = useMart(mart,host="ensembl.org",dataset=db)

tab1 = getBM(attributes=c('ensembl_gene_id',hom),mart=ensembl)
tab2 = getBM(attributes=c('ensembl_gene_id',inter),mart=ensembl,bmHeader=FALSE)
colnames(tab1) = c('sp1','sp2')
colnames(tab2) = c('sp1','sp2')
tab = rbind(tab1,tab2)

write.table(unique(tab),file=outfile,sep='\t',row.names=F,quote=F)
