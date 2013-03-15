library(GenomicFeatures)
library(Biostrings)

# Function read GFF
get.ranges.from.gff = function(gff.filename,fname) {
  gff = read.table(file=gff.filename,sep='\t',comment.char='#',quote='')
  colnames(gff) = c('seqid','source','type','start','end','score','strand','phase','attributes')
  ranges = GRanges(seqnames=gff$seqid,
                   ranges=IRanges(start=gff$start,end=gff$end),
                   strand=gff$strand)
  mcols(ranges) = gff[,c('source','type','attributes')]
  if(!missing(fname)) mcols(ranges)$fname = fname
  ranges
}

# Target file
targetfile = 'target'
target = read.table(file=targetfile,sep='\t',head=T,comment.char='',quote='')

# Load the query GFF
qfile = 'HR_vs_CI_ghost_rm_blast_ws6_Ci.gff'
qranges = get.ranges.from.gff(qfile)

res = data.frame()

for(i in 1:nrow(target)) {
  hfile = paste(target[i,'folder'],target[i,'file'],sep='/')
  fname = target[i,'fname']
  hranges = get.ranges.from.gff(hfile,fname)
  overlap = findOverlaps(qranges,hranges,ignore.strand=T)
  qdata = as.data.frame(qranges[queryHits(overlap)])
  hdata =as.data.frame(hranges[subjectHits(overlap)])
  res = rbind(res,cbind(qdata,hdata))
}
