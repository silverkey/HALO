#!/usr/bin/R

# VERSION: 0.3

library(GenomicFeatures)
library(Biostrings)

# Target file
targetfile = '../target'
# Query GFF file
qfile = '../BLAST/TR_vs_HS_blastn_HS.gff'

# Classify the matches that overlap GFFs
# Vector based on the "type" field of the loaded gffs
genic = c('gene')
exonic = c('exon')
# String based on the last column of the "target" file
#repetitive = 'repeat'
RFAM = 'RFAM'
tRNAscan = 'tRNAscan'
#ncrna = 'CiCs.ncRNA'
#mirna = 'miRNA'

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

# Function to assign the mapping bsed on the GFF overlaps
# It takes the ID of the match to test
assign.class = function(x) {
  class = 'NA'
  types = unique(res[res$index==x,]$h.type)
  fnames = unique(res[res$index==x,]$h.fname)
#  if(sum(fnames == repetitive)) {
#    class = 'repeat'
#  }
#  if(sum(fnames == mirna)) {
#    class = 'miRNA'
#  }
  if(sum(fnames == tRNAscan)) {
    class = 'tRNA'
  }
#  if(sum(fnames == ncrna)) {
#    class = 'ncRNA'
#  }
  if(sum(fnames == RFAM)) {
    class = 'RFAM'
  }
  if(class == 'NA') {
    if(sum(types %in% genic)) {
      if(sum(types %in% exonic)) {
        class = 'exonic'
      }
      else {
        class = 'intronic'
      }
    }
  }
  as.character(class)
}

get.attribute.value = function(df,key) {
  key = paste(key,'=',sep='')
  val = lapply(as.character(df$attribute),function(x)sub(key,'',unlist(strsplit(x,';'))[grep(key,unlist(strsplit(x,';')))]))
  unlist(val)
}


# N.B.: the ID field in the attribute of the GFF has been
# created as numeric autoincrement starting by 1 in the perl
# script comparative_blast_to_gff.pl
# In such a way the row number of qranges object corresponds
# to the ID of the blast match
qranges = get.ranges.from.gff(qfile)

# Create the datafrae in which to put the results of the overlaps
res = data.frame()

target = read.table(file=targetfile,sep='\t',head=T,comment.char='',quote='')
for(i in 1:nrow(target)) {
  hfile = paste(target[i,'folder'],target[i,'file'],sep='/')
  fname = target[i,'fname']
  hranges = get.ranges.from.gff(hfile,fname)
  overlap = findOverlaps(qranges,hranges,ignore.strand=T)
  qdata = as.data.frame(qranges[queryHits(overlap)])
  colnames(qdata) = paste('q',colnames(qdata),sep='.')
  qdata$index = queryHits(overlap)
  hdata =as.data.frame(hranges[subjectHits(overlap)])
  colnames(hdata) = paste('h',colnames(hdata),sep='.')
  res = rbind(res,cbind(qdata,hdata))
}

# Classify all the matches that do not overlap against the GFFs
# as intergenic
o.index = unique(res$index)
no.overlap = as.data.frame(qranges[-o.index])
no.overlap$index = get.attribute.value(no.overlap,'ID')
no.overlap$mapping = 'intergenic'

# Attach the mapping to the matches
o.mapping = sapply(o.index,assign.class,simplify="array")
mdf = data.frame(o.index,o.mapping)
colnames(mdf) = c('index','mapping')

mdf = rbind(mdf,no.overlap[,c('index','mapping')])
gffdf = as.data.frame(qranges)
gffdf$index = get.attribute.value(gffdf,'ID')
gffdf = merge(gffdf,mdf)
gffdf$attributes = paste(gffdf$attributes,gffdf$mapping,sep=';map=')

# Write an output table with all the overlaps
write.table(res,file='GFF_OVERLAPS.xls',sep="\t",row.names=F,quote=F)
write.table(gffdf,file='GFF_MAPPINGS.xls',sep="\t",row.names=F,quote=F)
gffcol = c('seqnames','source','type','start','end','score','strand','phase','attributes')
write.table(gffdf[,gffcol],file='GFF_MAPPINGS.gff',sep="\t",row.names=F,quote=F)

# Write an output table with the classification for each conserved fragment

# intergenic: NA

# genic: gene, mRNA
  # exonic: exon, CDS, five_prime_UTR, three_prime_UTR
  # intronic: NA
