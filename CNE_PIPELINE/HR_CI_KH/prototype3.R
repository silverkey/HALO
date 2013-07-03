source('functions3.R')

mc.cores = 22

hom.file = 'HOMOL/CI_HR_HOMOLOGS_GENE_LEVEL.txt'

sp1.blast.file = 'BLAST/HR_vs_CI_ghost_rm_besthitover_blastn_CI.gff'
sp1.ensg.file = 'ANNO/CI/CI_KH_genes.txt'
sp1.ense.file = 'ANNO/CI/CI_KH_exons.txt'
sp1.rfam.file = 'ANNO/CI/CI_KH_rm_1000bp_Rfam_CI.gff'
sp1.trna.file = 'ANNO/CI/CI_KH_rm_1000bp_trnascan.gff'
sp1.rep.file = 'ANNO/CI/KH-repeatModelerReapeatMasker.gff3'

sp2.blast.file = 'BLAST/HR_vs_CI_ghost_rm_besthitover_blastn_HR.gff'
sp2.ensg.file = 'ANNO/HR/HR_genes.txt'
sp2.ense.file = 'ANNO/HR/HR_exons.txt'
sp2.rfam.file = 'ANNO/HR/HR_rm_1000bp_Rfam_HR.gff'
sp2.trna.file = 'ANNO/HR/HR_rm_1000bp_trnascan.gff'
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
