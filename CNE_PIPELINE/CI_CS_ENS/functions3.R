library(GenomicFeatures)
library(Biostrings)
library(parallel)

analyze.blast = function(blast.file,ensg.file,ense.file,rfam.file,trna.file,rep.file) {

  rfam = NA
  trna = NA
  rep  = NA
  
  gff.colnames = c('seqname','source','type','start','end','score','strand','phase','attributes')

  blast = read.table(file=blast.file,sep='\t',head=F,quote='',comment.char='')
  colnames(blast) = gff.colnames
  #blast = blast[1:100,]
  
  ensg  = read.table(file=ensg.file,sep='\t',head=T,quote='',comment.char='')
  colnames(ensg) = c('gid','seqname','start','end','strand','biotype','description')

  ense  = read.table(file=ense.file,sep='\t',head=T,quote='',comment.char='')
  colnames(ense) = c('gid','eid','seqname','start','end','strand')

  if(!is.na(rfam.file)) rfam  = read.table(file=rfam.file,sep='\t',head=F,quote='',comment.char='#')
  if(!is.na(trna.file)) trna  = read.table(file=trna.file,sep='\t',head=F,quote='',comment.char='#')
  if(!is.na(rep.file)) rep = read.table(file=rep.file,sep='\t',head=F,quote='',comment.char='#')
  
  if(!is.na(rfam.file)) colnames(rfam) = gff.colnames
  if(!is.na(trna.file)) colnames(trna) = gff.colnames
  if(!is.na(rep.file)) colnames(rep) = gff.colnames
  
  lnc = c('3prime_overlapping_ncrna','lincRNA','non_coding','antisense','sense_intronic','sense_overlapping')
  prc = c('protein_coding')
  mir = c('miRNA')
  
  coding = ensg[ensg$biotype=='protein_coding',]
  
  blast.r = get.ranges.from.df(blast)  
  ensg.r  = get.ranges.from.df(ensg)
  ense.r  = get.ranges.from.df(ense)

  if(!is.na(rfam.file)) rfam.r = get.ranges.from.df(rfam)
  if(!is.na(trna.file)) trna.r = get.ranges.from.df(trna)
  if(!is.na(rep.file)) rep.r = get.ranges.from.df(rep)
  
  bl.ge = findOverlaps(blast.r,ensg.r,ignore.strand=T)
  bl.ex = findOverlaps(blast.r,ense.r,ignore.strand=T)

  if(!is.na(rfam.file)) bl.rf = findOverlaps(blast.r,rfam.r,ignore.strand=T)
  if(!is.na(trna.file)) bl.tr = findOverlaps(blast.r,trna.r,ignore.strand=T)
  if(!is.na(rep.file)) bl.re = findOverlaps(blast.r,rep.r,ignore.strand=T)
  
  omap = mclapply(1:length(blast.r),
                  map.on.ensg,ensg=ensg,ense=ense,rfam=rfam,trna=trna,rep=rep,
                              bl.ge=bl.ge,bl.ex=bl.ex,bl.rf=bl.rf,bl.tr=bl.tr,bl.re=bl.re,
                              prc=prc,mir=mir,lnc=lnc,
                  mc.cores=mc.cores)

  map = as.data.frame(matrix(unlist(omap),ncol=3,byrow=T))
  colnames(map) = c('mapping','overlap','class')
  
  fl.id = get.flanking(blast.r,coding,mc.cores)
  
  list(map=map,fl.id=fl.id)
}

get.flanking = function(r,gtab,mc.cores) {

  g = get.ranges.from.df(gtab)

  strand(r) = '+'
  strand(g) = '+'
  
  # Get the overlapping feature in chrtr
  ov = findOverlaps(r,g)
  ov = as.data.frame(ov)
  index = data.frame(index=1:length(r))
  data = merge(index,ov,by.x='index',by.y='queryHits',all.x=T)
  colnames(data) = c('index','ov')
  
  # Get the next feature in chrtr (precede is referred to the query i.e. ranges[i])
  up = precede(r,g)
  up = data.frame(index=1:length(r),up=up)
  data = merge(data,up,by.x='index',by.y='index',all.x=T)
  
  # Get the previous feature in chrtr (follow is referred to the query i.e. ranges[i])
  down = follow(r,g)
  down = data.frame(index=1:length(r),down=down)
  data = merge(data,down,by.x='index',by.y='index',all.x=T)
  
  fl.in = mclapply(unique(data$index),function(x) as.numeric(na.omit(unique(as.numeric(unlist(c(data[data$index==x,2:4])))))),mc.cores=mc.cores)
  fl.id = mclapply(fl.in,function(x) as.character(gtab[x,'gid']),mc.cores=mc.cores)
  
  fl.id
}

# Function read GFF
# The fname variable is used to store information from the type of the hit gff
# It is not really used at the moment to indicate the type of query because at
# moment it is a single file. In the future we might decide to change this
# by specifing more queries at the same time.
get.ranges.from.gff = function(gff.filename,fname) {
  gff = read.table(file=gff.filename,sep='\t',comment.char='#',quote='')
  colnames(gff) = c('seqname','source','type','start','end','score','strand','phase','attributes')
  ranges = GRanges(seqnames=gff$seqname,
                   ranges=IRanges(start=gff$start,end=gff$end),
                   strand=gff$strand)
  mcols(ranges) = gff[,c('source','type','score','phase','attributes')]
  if(!missing(fname)) mcols(ranges)$fname = fname
  ranges
}

get.ranges.from.df = function(df,fname) {
  ranges = GRanges(seqnames=df$seqname,
                   ranges=IRanges(start=df$start,end=df$end),
                   strand=df$strand)
  if(!missing(fname)) mcols(ranges)$fname = fname
  ranges
}

get.covered.bp = function(ranges) {
  strand(ranges) = Rle('+')
  coverage = reduce(ranges)
  cov.l = width(coverage)
  cons.bp = sum(cov.l)
  cons.bp
}

findOverlapsDF = function(df1,df2){
  ov = apply(df1,1,function(x) df2[df2$seqname==x['seqname'] & df2$start<=x['end'] & df2$end>=x['start'],])
  ov
}

get.attribute.value = function(df,key) {
  key = paste(key,'=',sep='')
  val = lapply(as.character(df$attribute),function(x)sub(key,'',unlist(strsplit(x,';'))[grep(key,unlist(strsplit(x,';')))]))
  unlist(val)
}

map.on.ensg = function(x, ensg, ense, rfam, trna, rep, bl.ge, bl.ex, bl.rf, bl.tr, bl.re, prc, mir, lnc) {
  overlap = NA
  mapping = NA
  class = NA
  
  ehits = subjectHits(bl.ex[queryHits(bl.ex)==x,])
  ex = ense[ehits,]
  if(nrow(ex)>0) {
    mapping = 'exon'
    overlap = paste(unique(ex$gid),collapse=',')
    class = unique(as.character(ensg[ensg$gid %in% unique(as.character(ex$gid)),'biotype']))
    
    if(sum(class %in% prc)) {
      class = 'coding'
    }
    else if(sum(class %in% mir)) {
      class = 'miRNA'
    }
    else if(sum(class %in% lnc)) {
      class = 'lncRNA'
    }
    else {
      class = 'other'
    }
  }  
  if(is.na(mapping)) {
    ghits = subjectHits(bl.ge[queryHits(bl.ge)==x,])
    ge = ensg[ghits,]
    if(nrow(ge)>0) {
      mapping = 'intron'
      overlap = paste(unique(ge$gid),collapse=',')
      class = unique(as.character(ge$biotype))
      
      if(sum(class %in% prc)) {
        class = 'coding'
      }
      else if(sum(class %in% mir)) {
        class = 'miRNA'
      }
      else if(sum(class %in% lnc)) {
        class = 'lncRNA'
      }
      else {
        class = 'other'
      }      
    }
  }  
  if(is.na(mapping)) {
    if(!is.na(rfam)) {
      rfhits = subjectHits(bl.rf[queryHits(bl.rf)==x,])
      rf = rfam[rfhits,]
      if(nrow(rf)>0) {
        mapping = 'RFAM'
        overlap = 'RFAM'
        class = 'RFAM'
      }
    }
  }
  if(is.na(mapping)) {
    if(!is.na(trna)) {
      trhits = subjectHits(bl.tr[queryHits(bl.tr)==x,])
      tr = trna[trhits,]
      if(nrow(tr)>0) {
        mapping = 'tRNA'
        overlap = 'tRNA'
        class = 'tRNA'
      }
    }
  }
  if(is.na(mapping)) {
    if(!is.na(rep)) {
      rehits = subjectHits(bl.re[queryHits(bl.re)==x,])
      re = rep[rehits,]
      if(nrow(re)>0) {
        mapping = 'repeat'
        overlap = 'repeat'
        class = 'repeat'
      }
    }  
  }
  if(is.na(mapping)) {
    mapping = 'intergenic'
    overlap = 'intergenic'
    class = 'intergenic'
  }  
  c(mapping,overlap,class)
}
