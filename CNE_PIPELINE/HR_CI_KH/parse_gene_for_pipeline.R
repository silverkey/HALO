t = read.table(file='augustusModelGenesWithUTRs_Hr.gff3',sep='\t',head=F,quote='',comment.char='',skip=1)
gff.colnames = c('seqname','source','type','start','end','score','strand','phase','attributes')

get.attribute.value = function(df,key) {
  key = paste(key,'=',sep='')
  val = lapply(as.character(df$attribute),function(x)sub(key,'',unlist(strsplit(x,';'))[grep(key,unlist(strsplit(x,';')))]))
  unlist(val)
}

colnames(t) = gff.colnames

genes = t[t$type=='gene',]
exons = t[t$type=='CDS' | t$type =='five_prime_UTR' | t$type=='three_prime_UTR',]

genes$gid = get.attribute.value(genes,'ID')
exons$gid = get.attribute.value(exons,'Parent')
exons$gid = gsub('.TU.+$','',exons$gid,perl=T)

# ensembl_gene_id	chromosome_name	start_position	end_position	strand	gene_biotype	description
genes$biotype = 'protein_coding'
genes$description = NA
genes = genes[,c('gid','seqname','start','end','strand','biotype','description')]

# ensembl_gene_id	ensembl_exon_id	chromosome_name	exon_chrom_start	exon_chrom_end	strand
exons$eid = NA
exons = exons[,c('gid','eid','seqname','start','end','strand')]

write.table(genes,file='../ANNO/HR_genes.txt',sep="\t",row.names=F,quote=F)
write.table(exons,file='../ANNO/HR_exons.txt',sep="\t",row.names=F,quote=F)


t = read.table(file='KH.KHGene.2012.gff3',sep='\t',head=F,quote='',comment.char='',skip=1)

gff.colnames = c('seqname','source','type','start','end','score','strand','phase','attributes')

get.attribute.value = function(df,key) {
  key = paste(key,'=',sep='')
  val = lapply(as.character(df$attribute),function(x)sub(key,'',unlist(strsplit(x,';'))[grep(key,unlist(strsplit(x,';')))]))
  unlist(val)
}

colnames(t) = gff.colnames

genes = t[t$type=='gene',]
exons = t[t$type=='CDS' | t$type =='five_prime_UTR' | t$type=='three_prime_UTR',]

genes$gid = get.attribute.value(genes,'ID')
exons$gid = get.attribute.value(exons,'Parent')
exons$gid = gsub('.TU.+$','',exons$gid,perl=T)

# ensembl_gene_id	chromosome_name	start_position	end_position	strand	gene_biotype	description
genes$biotype = 'protein_coding'
genes$description = NA
genes = genes[,c('gid','seqname','start','end','strand','biotype','description')]

# ensembl_gene_id	ensembl_exon_id	chromosome_name	exon_chrom_start	exon_chrom_end	strand
exons$eid = NA
exons = exons[,c('gid','eid','seqname','start','end','strand')]

write.table(genes,file='../ANNO/CI_KH_genes.txt',sep="\t",row.names=F,quote=F)
write.table(exons,file='../ANNO/CI_KH_exons.txt',sep="\t",row.names=F,quote=F)

