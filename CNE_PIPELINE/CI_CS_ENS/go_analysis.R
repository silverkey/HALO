#------------
# PARAMETERS
#------------

# File containing the annotations from imparanoid
anno.tab = 'cintestinalis_go.txt'

# File containing the CNEs
CNE.tab = 'CI_vs_CS_ensembl_def_bho.txt'

# File containing the mapping of the GO id and definition
gomap = 'go_map'

# Minimum number of counts for a class
min = 10

# To look for enrichments use 'g', for impoverishment use 'l', for both use 't'
prop.alt = 'g'

# Adjusted pvalue to be considered significant
p.filt = 0.1

# The fold the proportion of difference between counts associated to CNEs and universe
mult = 1

get.counts = function(anno) {
  counts = as.data.frame(table(anno$goid))
  colnames(counts) = c('goid','count')
  counts = counts[grep('GO:',counts$goid),]
  counts
}

calculate.enrichments = function(div.sel,div.uni,n.sel,n.uni,go) {
  div = merge(div.sel,div.uni,by.x='goid',by.y='goid',all.y=T)
  div$count.x[is.na(div$count.x)] = 0
  div$pval = apply(div,1,function(x) prop.test(c(as.numeric(x[2]),as.numeric(x[3])),c(n.sel,n.uni),alternative=prop.alt)$p.value)
  div = merge(div,go,by.x='goid',by.y='go_id')
  if(prop.alt == 'g') div = div[(div$count.x/n.sel >= mult*(div$count.y/n.uni)) & div$count.x>=min,]
  if(prop.alt == 'l') div = div[(div$count.x/n.sel <= mult*(div$count.y/n.uni)) & div$count.x>=min,]
  if(prop.alt == 't') div = div[div$count.x>=min,]
  div$padj = p.adjust(div$pval)
  div[order(div$padj,decreasing=T),]
}

t = read.table(file=CNE.tab,sep='\t',head=T,quote='',comment.char='')
# HARDCODED SUBSETTING!!!
sel = subset(t,syn==1 & ((map=='intergenic' & Omap=='intergenic') | (map=='intron' & Omap=='intron')))
selid = unique(unlist(strsplit(as.character(sel$locus),',')))

anno = read.table(file=anno.tab,sep='\t',head=T,comment.char='',quote='')[,c(1,2)]
colnames(anno) = c('gid','goid')

sel.anno = subset(anno,gid %in% selid)

go = read.table(file=gomap,sep='\t',head=T,quote='',comment.char='',stringsAsFactors=F)

nuni = length(unique(anno$gid))
nsel = length(unique(selid))

uni.count = get.counts(anno)
sel.count = get.counts(sel.anno)

res = calculate.enrichments(sel.count,uni.count,nsel,nuni,go)

# sig.bp = subset(res.bp,padj<=p.filt)

