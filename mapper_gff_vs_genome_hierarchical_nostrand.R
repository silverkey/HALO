# JACK THE MAPPER!!!
#
# VERSION gff_vs_genome_hierarchical_nostrand
#
# Takes a GFF file from a blast of 2 genomes and calculates the overlap
# between the features and a given genome.
#
# The GFF should be referred to the genome of which the script will
# load the annotations to calculate the overlap.
# It consider an hyerarchi of overlap defined in the variable overlap.order.
#
# Therefore it will return only the annotations associated to the first in
# hyerarchi. I.E. if a match overlap with a CDS and with an intron the script
# will report only one of the associations as defined in the order variable.

# ----------------- #
#   RUNTIME VARS    #
# ----------------- #

# The directory and the file GFF inside it to annotate
sel.dir = '/Users/remo/ANALYSIS/HALO/ANALYSIS/BLASTRES'
gff.file = 'HR_vs_CI_rm_top_CI.gff'

# Info about the GFF to analyze
gff.filename = 'HR_vs_CI_rm_top_CI.gff'
att.list = c('ID','Target','E','O')

# The directory containing the database or in which you will create it
dbdir = '/Users/remo/ANALYSIS/HALO/ANALYSIS/BLASTRES'
organism = 'cintestinalis'
download = 'F'

# Length of promoters
s.prox = 1000
e.prox = 0

# Completely noncodings transcripts will overlap exons but not cds nor utr
# We use a hierarchy order for which every time a range is overlapping a feature than
# it is associated to the feature and cut-out from the ranges, so that it cannot overlapping
# with the next order feature.
# Order to test overlap in this script (this is defined in the overlap.order variable:
# 1) cds
# 2) utr
# 3) exon
# 4) promoter
# 5) intron

# Define the order to test the overlap:
overlap.order = c('cds','utr5','utr3','exon','promoter','intron')

# ----------------- #
#     FUNCTIONS     #
# ----------------- #
get.promoters = function(transcripts,s.prox,e.prox) {
  pos = strand(transcripts) == '+'
  promoter.start = as.vector(ifelse(pos,start(transcripts)-s.prox,end(transcripts)-e.prox))
  promoter.end = as.vector(ifelse(pos,start(transcripts)+e.prox,end(transcripts)+s.prox))
  promoter.strand = ifelse(pos,'+','-')
  promoter.chr = seqnames(transcripts)
  promoter.tx_name = unlist(mcols(transcripts[,'tx_name'])[,1])
  promoter.gene_id = unlist(mcols(transcripts[,'gene_id'])[,1])
  promoters = GRanges(seqnames=promoter.chr,
                      ranges=IRanges(start=promoter.start,end=promoter.end),
                      strand=promoter.strand,
                      tx_name=promoter.tx_name,gene_id=promoter.gene_id)
  promoters
}

# Used for transcripts
calculate.range.overlap = function(ranges,features,name,dir) {
  overlap = findOverlaps(ranges,features,ignore.strand=T)
  overlap = as.data.frame(overlap)
  id = as.character(mcols(ranges[overlap$queryHits])$ID)
  gene = as.character(mcols(features[overlap$subjectHits])$gene_id)
  transcript = as.character(mcols(features[overlap$subjectHits])$tx_name)
  strand = as.character(strand(features[overlap$subjectHits]))
  res = as.data.frame(unique(cbind(id,name,gene,transcript,strand,dir)))
  if(ncol(res)<6) no.res = ranges
  if(ncol(res)==6) no.res = ranges[-overlap$queryHits]
  list(res=res,no.res=no.res)
}

# Used for promoters, exons, introns
calculate.rangelist.overlap = function(ranges,features,name,txid.gid,dir) {
  overlap = findOverlaps(ranges,features,ignore.strand=T)
  overlap = as.data.frame(overlap)
  id = as.character(mcols(ranges[overlap$queryHits])$ID)
  gene = as.character(txid.gid[names(features[overlap$subjectHits]),2])
  transcript = as.character(txid.gid[names(features[overlap$subjectHits]),1])
  strand = as.character(txid.gid[names(features[overlap$subjectHits]),3])
  res = as.data.frame(unique(cbind(id,name,gene,transcript,strand,dir)))
  if(ncol(res)<6) no.res = ranges
  if(ncol(res)==6) no.res = ranges[-overlap$queryHits]
  list(res=res,no.res=no.res)
}

get.ranges.from.gff = function(gff.filename,att.list) {
  gff = read.table(file=gff.filename,sep='\t',comment.char='',quote='')
  colnames(gff) = c('seqid','source','type','start','end','score','strand','phase','attribute')
  att.df = c()
  for(i in 1:length(att.list)) {
    a = att.list[i]
    a = paste(a,'=',sep='')
    val = unlist(lapply(as.character(gff$attribute),function(x)sub(a,'',unlist(strsplit(x,';'))[grep(a,unlist(strsplit(x,';')))])))
    att.df = cbind(att.df,val)
  }
  att.df = as.data.frame(att.df)
  colnames(att.df) = att.list
  att.df$score = gff$score
  ranges = GRanges(seqnames=gff$seqid,
                   ranges=IRanges(start=gff$start,end=gff$end),
                   strand=gff$strand)
  mcols(ranges) = att.df
  ranges
}

# ----------------- #
#      SCRIPT       #
# ----------------- #

library("GenomicFeatures")

if(download == 'T') {
  transdb = makeTranscriptDbFromBiomart(biomart="ensembl",dataset=paste(organism,"gene_ensembl",sep='_'))
  saveDb(transdb,file=paste(organism,"sqlite",sep='.'))
} else {
  transdb = loadDb(file=paste(organism,"sqlite",sep='.'))
}

setwd(sel.dir)

# Build Features
transcripts = transcripts(transdb,columns=c("tx_id","tx_name","gene_id")) # range
cds = cdsBy(transdb,by='tx',use.names=T)                                  # rangelist
exon = exonsBy(transdb,by='tx',use.names=T)                               # rangelist
intron = intronsByTranscript(transdb,use.names=T)                         # rangelist
utr5 = fiveUTRsByTranscript(transdb,use.names=T)                          # rangelist
utr3 = threeUTRsByTranscript(transdb,use.names=T)                         # rangelist
promoter = get.promoters(transcripts,s.prox,e.prox)                       # range

# Build a comfortable table to associate transcripts->genes->strands
txid.gid = as.data.frame(cbind(
                         unlist(mcols(transcripts[,'tx_name'])[,1]),
                         unlist(mcols(transcripts[,'gene_id'])[,1]),
                         as.character(strand(transcripts))))
colnames(txid.gid) = c('tx_id','gene_id','strand')
rownames(txid.gid)=txid.gid$tx_id

# Build table for results
features.res.cn = c('ID','overlap','gene','trans','strand','dir')
features.res = matrix(ncol=length(features.res.cn),nrow=0)

# Build the ranges of the fragments you want to test the overlap with features
ranges = get.ranges.from.gff(gff.filename,att.list)
rtab = as.data.frame(ranges)

# Start the overlap analysis
for(i in 1:length(overlap.order)) {
  fname = overlap.order[i]
  features = get(fname)
  if(class(features)=='GRanges') {
    lres = calculate.range.overlap(ranges,features,fname,'NA')
    ranges = lres$no.res
    if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)
  }
  if(class(features)=='GRangesList') {
    lres = calculate.rangelist.overlap(ranges,features,fname,txid.gid,'NA')
    ranges = lres$no.res
    if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)
  }
}

colnames(features.res) = features.res.cn

intergenic = as.data.frame(lres$no.res)
if(nrow(intergenic)>=1) {
  intergenic$overlap='intergenic'
  intergenic$gene=NA
  intergenic$trans=NA
  intergenic$strand=NA
  intergenic$dir=NA
  features.res = rbind(features.res,intergenic[,features.res.cn])
}

map = unique(features.res[,c('ID','overlap')])
final = merge(rtab,map,by.x='ID',by.y='ID',sort=F)
write.table(final,file='RESULTS_FEATURE_ANNOTATED.xls',sep="\t",row.names=F,quote=F)
