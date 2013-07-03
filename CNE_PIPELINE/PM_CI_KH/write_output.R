comparison1 = 'CI_vs_PM_def_bho'
comparison2 = 'PM_vs_CI_def_bho'

sp1.blast$attributes = paste(sp1.blast$attributes,sp1.res$map$mapping,sep=';map=')
sp1.blast$attributes = paste(sp1.blast$attributes,sp2.res$map$mapping,sep=';Omap=')
sp1.blast$attributes = paste(sp1.blast$attributes,sp1.res$map$class,sep=';class=')
sp1.blast$attributes = paste(sp1.blast$attributes,sp2.res$map$class,sep=';Oclass=')
sp1.fl.id = unlist(lapply(sp1.res$fl.id,function(x) paste(x,collapse=',')))
sp1.blast$attributes = paste(sp1.blast$attributes,sp1.fl.id,sep=';locus=')
sp1.blast$attributes = paste(sp1.blast$attributes,syn,sep=';syn=')
write.table(sp1.blast,file=paste(comparison1,'gff',sep='.'),sep="\t",row.names=F,quote=F)

sp1.blast$evalue = as.numeric(get.attribute.value(sp1.blast,'E'))
sp1.blast$length = as.numeric(sp1.blast$end-sp1.blast$start+1)
sp1.blast$map = sp1.res$map$mapping
sp1.blast$Omap = sp2.res$map$mapping
sp1.blast$class = sp1.res$map$class
sp1.blast$Oclass = sp2.res$map$class
sp1.blast$locus = sp1.fl.id
sp1.blast$ID = get.attribute.value(sp1.blast,'ID')
sp1.blast$Target = get.attribute.value(sp1.blast,'Target')
sp1.blast$O = get.attribute.value(sp1.blast,'O')
sp1.blast$syn = syn
write.table(sp1.blast[,-9],file=paste(comparison1,'txt',sep='.'),sep="\t",row.names=F,quote=F)


sp2.blast$attributes = paste(sp2.blast$attributes,sp2.res$map$mapping,sep=';map=')
sp2.blast$attributes = paste(sp2.blast$attributes,sp1.res$map$mapping,sep=';Omap=')
sp2.blast$attributes = paste(sp2.blast$attributes,sp2.res$map$class,sep=';class=')
sp2.blast$attributes = paste(sp2.blast$attributes,sp1.res$map$class,sep=';Oclass=')
sp2.fl.id = unlist(lapply(sp2.res$fl.id,function(x) paste(x,collapse=',')))
sp2.blast$attributes = paste(sp2.blast$attributes,sp2.fl.id,sep=';locus=')
sp2.blast$attributes = paste(sp2.blast$attributes,syn,sep=';syn=')
write.table(sp2.blast,file=paste(comparison2,'gff',sep='.'),sep="\t",row.names=F,quote=F)

sp2.blast$evalue = as.numeric(get.attribute.value(sp2.blast,'E'))
sp2.blast$length = as.numeric(sp2.blast$end-sp2.blast$start+1)
sp2.blast$map = sp2.res$map$mapping
sp2.blast$Omap = sp1.res$map$mapping
sp2.blast$class = sp2.res$map$class
sp2.blast$Oclass = sp1.res$map$class
sp2.blast$locus = sp2.fl.id
sp2.blast$ID = get.attribute.value(sp2.blast,'ID')
sp2.blast$Target = get.attribute.value(sp2.blast,'Target')
sp2.blast$O = get.attribute.value(sp2.blast,'O')
sp2.blast$syn = syn
write.table(sp2.blast[,-9],file=paste(comparison2,'txt',sep='.'),sep="\t",row.names=F,quote=F)
