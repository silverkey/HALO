#!/usr/bin/R

# VERSION: 0.2

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

# Function to assign the mapping bsed on the GFF overlaps
# It takes the ID of the match to test
assign.class = function(x) {
  class = 'NA'
  types = res[res$index==x,]$h.type
  if(sum(types %in% genic)) {
    if(sum(types %in% exonic)) {
      class = 'exonic'
    }
    else {
      class = 'intronic'
    }
  }
  if(sum(types %in% repetitive)) {
      class = 'repeat'
  }
  as.character(class)
}

get.attribute.value = function(df,key) {
  key = paste(key,'=',sep='')
  val = lapply(as.character(df$attribute),function(x)sub(key,'',unlist(strsplit(x,';'))[grep(key,unlist(strsplit(x,';')))]))
  unlist(val)

}

# Target file
targetfile = 'target'
target = read.table(file=targetfile,sep='\t',head=T,comment.char='',quote='')

# Load the query GFF
qfile = 'HR_vs_CI_ghost_rm_megablast_CI.gff'

# N.B.: the ID field in the attribute of the GFF has been
# created as numeric autoincrement starting by 1 in the perl
# script comparative_blast_to_gff.pl
# In such a way the row number of qranges object corresponds
# to the ID of the blast match
qranges = get.ranges.from.gff(qfile)

# Create the datafrae in which to put the results of the overlaps
res = data.frame()

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

# Classify the matches that overlap GFFs
genic = c('gene','mRNA')
exonic = c('exon','CDS','five_prime_UTR','three_prime_UTR')
repetitive = c('dispersed_repeat')

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
