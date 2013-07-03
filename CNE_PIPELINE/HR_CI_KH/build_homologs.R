# Collect only the homology relationship
tr = read.table(file='../HOMOL/CI_HR_HOMOLOGS_COORDS.txt',sep='\t',quote="",comment.char='',head=T)
sp1.tr = tr[tr$Species=='Ciona_intestinalis',]
sp2.tr = tr[tr$Species=='Halocynthia_roretzi',]
hom = merge(sp1.tr[,c(1,2)],sp2.tr[,c(1,2)],by.x='OrthoID',by.y='OrthoID',sort=F)[,-1]
colnames(hom) = c('sp1','sp2')
hom$simple1 = gsub(".v.+$",'',hom$sp1,perl=T)
hom2 = unique(hom[,c(3,2)])
colnames(hom) = c('sp1','sp2')
write.table(hom2,file='../HOMOL/CI_HR_HOMOLOGS_GENE_LEVEL.txt',sep="\t",row.names=F,quote=F)

