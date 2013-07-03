library('biomaRt')

args = commandArgs(trailingOnly=T)

# For ensembl genes use the following mart:
# ENSEMBL_MART_ENSEMBL
mart = args[1]
sp1 = args[2]
sp2 = args[3]

db = paste(sp1,'gene_ensembl',sep='_')
hom = paste(sp2,'homolog_ensembl_gene',sep='_')
outfile = paste(sp1,sp2,'homologs.txt',sep='_')

ensembl = useMart(mart,host="ensembl.org",dataset=db)

# trubripes_inter_paralog_ensembl_gene return an error
tab = getBM(attributes=c('ensembl_gene_id',hom),mart=ensembl)

write.table(tab,file=outfile,sep='\t',row.names=F,quote=F)
