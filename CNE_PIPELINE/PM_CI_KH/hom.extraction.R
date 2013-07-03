h = read.table(file='sqltable.Phallusia_mammillata-Ciona_intestinalis',sep='\t',head=F,comment.char='',quote='')
colnames(h) = c('unkid','ortoid','sp','score','trid')

sp1 = subset(h,sp=='Ciona_intestinalis')
sp2 = subset(h,sp=='Phallusia_mammillata')
m = merge(sp1,sp2,by='ortoid',sort=F)

ht = m[,c('trid.x','trid.y')]
colnames(ht) = c('sp1','sp2')

ht$sp1 = gsub('.v.+$','',ht$sp1)
ht$sp2 = gsub('.TU.+$','',ht$sp2)

write.table(ht,file='CI_KH_PM_HOM.txt',sep='\t',row.names=F,quote=F)
