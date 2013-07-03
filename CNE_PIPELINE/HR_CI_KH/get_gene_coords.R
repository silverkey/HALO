library('biomaRt')

args = commandArgs(trailingOnly=T)

# For ensembl genes use the following mart:
# ENSEMBL_MART_ENSEMBL
mart = args[1]
org = args[2]

db = paste(org,'gene_ensembl',sep='_')
outfile = paste(org,'gene_coord.txt',sep='_')

ensembl = useMart(mart,host="ensembl.org",dataset=db)

tab = getBM(attributes=c('ensembl_gene_id',
                         'chromosome_name',
                         'start_position',
                         'end_position',
                         'strand',
                         'gene_biotype',
                         'description'),mart=ensembl)

write.table(tab,file=outfile,sep='\t',row.names=F,quote=F)
