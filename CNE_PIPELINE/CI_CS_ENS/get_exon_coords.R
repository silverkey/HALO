library('biomaRt')

args = commandArgs(trailingOnly=T)

# For ensembl genes use the following mart:
# ENSEMBL_MART_ENSEMBL
mart = args[1]
org = args[2]

db = paste(org,'gene_ensembl',sep='_')
outfile = paste(org,'exon_coord.txt',sep='_')

ensembl = useMart(mart,host="ensembl.org",dataset=db)

tab = getBM(attributes=c('ensembl_gene_id',
                         'ensembl_exon_id',
                         'chromosome_name',
                         'exon_chrom_start',
                         'exon_chrom_end',
                         'strand'),mart=ensembl)

write.table(tab,file=outfile,sep='\t',row.names=F,quote=F)
