library(GenomicFeatures)
library(Biostrings)

# Function read GFF
# The fname variable is used to store information from the type of the hit gff
# It is not really used at the moment to indicate the type of query because at
# moment it is a single file. In the future we might decide to change this
# by specifing more queries at the same time.
get.ranges.from.gff = function(gff.filename,fname) {
  gff = read.table(file=gff.filename,sep='\t',comment.char='#',quote='')
  colnames(gff) = c('seqid','source','type','start','end','score','strand','phase','attributes')
  ranges = GRanges(seqnames=gff$seqid,
                   ranges=IRanges(start=gff$start,end=gff$end),
                   strand=gff$strand)
  mcols(ranges) = gff[,c('source','type','score','phase','attributes')]
  if(!missing(fname)) mcols(ranges)$fname = fname
  ranges
}







hom = paste(sp2,'homolog_ensembl_gene',sep='_')





library('biomaRt')
db = 'hsapiens_gene_ensembl'
mart = useMart('ensembl',dataset=db)
hom = getBM(attributes=c('ensembl_gene_id','mmusculus_homolog_ensembl_gene'),mart=mart)
ipa = getBM(attributes=c('ensembl_gene_id','mmusculus_inter_paralog_ensembl_gene'),mart=mart)

