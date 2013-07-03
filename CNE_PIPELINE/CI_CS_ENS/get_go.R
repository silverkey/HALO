library('biomaRt')

args = commandArgs(trailingOnly=T)

# For ensembl genes use the following mart:
# ENSEMBL_MART_ENSEMBL
mart = args[1]
sp = args[2]

db = paste(sp,'gene_ensembl',sep='_')

outfile = paste(sp,'go.txt',sep='_')

ensembl = useMart(mart,host="ensembl.org",dataset=db)

tab = getBM(attributes=c('ensembl_gene_id','go_id'),mart=ensembl)

colnames(tab) = c('gid','goid')

write.table(unique(tab),file=outfile,sep='\t',row.names=F,quote=F)
