source('functions3.R')

mc.cores = 24

hom.file = 'HOMOL/webmart.txt'

sp1.blast.file = 'BLAST/CS_vs_CI_def_besthitover_blastn_CI.gff'
sp1.ensg.file = 'ANNO/cintestinalis_gene_coord.txt'
sp1.ense.file = 'ANNO/cintestinalis_exon_coord.txt'
sp1.rfam.file = 'RFAM/CI_Rfam_CI.gff'
sp1.trna.file = 'TRNA/CI_trna.gff'
sp1.rep.file = 'RM/Ciona_intestinalis.KH.71.dna_rm.toplevel.fa.repeatmasker.out.gff'

sp2.blast.file = 'BLAST/CS_vs_CI_def_besthitover_blastn_CS.gff'
sp2.ensg.file = 'ANNO/csavignyi_gene_coord.txt'
sp2.ense.file = 'ANNO/csavignyi_exon_coord.txt'
sp2.rfam.file = NA
sp2.trna.file = NA
sp2.rep.file = NA

check.syn = function(x) {
  fl1 = sp1.res$fl.id[[x]]
  fl2 = sp2.res$fl.id[[x]]
  syn = NA
  if(length(fl1)>0 & length(fl2)>0) {
    n = sum(as.character(unique(hom[hom$sp1 %in% fl1,'sp2'])) %in% fl2)
    syn = ifelse(n>0, 1, 0)
  }
  syn
}

sp1.res = analyze.blast(sp1.blast.file,sp1.ensg.file,sp1.ense.file,sp1.rfam.file,sp1.trna.file,sp1.rep.file)
sp2.res = analyze.blast(sp2.blast.file,sp2.ensg.file,sp2.ense.file,sp2.rfam.file,sp2.trna.file,sp2.rep.file)

gff.colnames = c('seqname','source','type','start','end','score','strand','phase','attributes')

sp1.blast = read.table(file=sp1.blast.file,sep='\t',head=F,quote='',comment.char='')
colnames(sp1.blast) = gff.colnames

sp2.blast = read.table(file=sp2.blast.file,sep='\t',head=F,quote='',comment.char='')
colnames(sp2.blast) = gff.colnames

sp1.ensg = read.table(file=sp1.ensg.file,sep='\t',head=T,quote='',comment.char='')
colnames(sp1.ensg) = c('gid','seqname','start','end','strand','biotype','description')
sp1.ense = read.table(file=sp1.ense.file,sep='\t',head=T,quote='',comment.char='')
colnames(sp1.ense) = c('gid','eid','seqname','start','end','strand')

sp2.ensg = read.table(file=sp2.ensg.file,sep='\t',head=T,quote='',comment.char='')
colnames(sp2.ensg) = c('gid','seqname','start','end','strand','biotype','description')
sp2.ense = read.table(file=sp2.ense.file,sep='\t',head=T,quote='',comment.char='')
colnames(sp2.ense) = c('gid','eid','seqname','start','end','strand')

hom = read.table(file=hom.file,sep='\t',head=T,quote='',comment.char='')
colnames(hom) = c('sp1','sp2')

syn = unlist(mclapply(1:nrow(sp1.blast),function(x) check.syn(x), mc.cores=mc.cores))
