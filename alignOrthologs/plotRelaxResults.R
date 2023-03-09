library(tidyverse)
library(rjson)
library(grid)
library(gridExtra)

g_legend<-function(a.gplot){
	a.gplot=a.gplot+theme(legend.position="right")
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]
	return(legend)
}


args=commandArgs(T)


cutoff_PValue=as.numeric(args[1])
colorFile=args[2] %>% as.character()
colorEq=read.table(colorFile, sep="\t", col.names=c("specie","colour"),comment.char="@") %>% as_tibble() %>% pull(colour,specie)

pantherTab=as.character(args[3])
speciesOrder=as.character(args[4]) %>% strsplit(",") %>% .[[1]]


all_files=args[5:length(args)]

all_files_fields=all_files %>% strsplit(split=".",fixed=T)
pfs.lst=vector(mode="list",length(all_files_fields))
for(i in seq(1,length(all_files_fields))){
	pf=as.character(all_files_fields[[i]][2])
	pfs.lst[[i]]=pf
}
pfs=unlist(pfs.lst)

r.lst=vector(mode="list",length(pfs))
names(r.lst)=pfs
dnds.lst=vector(mode="list",length(pfs))
names(dnds.lst)=pfs

#allowedSpecies=c("Austrolebias_charrua","Cynopoecilus_melanotaenia","Austrofundulus_limnaeus","Nothobranchius_furzeri","Oryzias_latipes")
#allowedSpecies=names(colorEq)
allowedSpecies=speciesOrder

# load all data
for(f in all_files){

	print(paste0("working on ",f,"..."))
	d=fromJSON(file=f)
	panf=strsplit(f,split=".",fixed=T)[[1]][2]
	#d=fromJSON(file="panther.PTHR48177.out-gb.noSpaces.uniqueId.noGapSeqs.noStopCodon.RELAX.json")

	dnds.v=d$fits$`MG94xREV with separate rates for branch sets`$`Rate Distributions`$`non-synonymous/synonymous rate ratio for *test*` %>% unlist
	dnds=dnds.v[1]/dnds.v[2]
	dnds.lst[[panf]]=dnds %>% as.numeric()
	species=names(d$`test results`)
	species=species[species %in% allowedSpecies]
	d.lst=vector(mode="list",length(species))
	names(d.lst)=species
	for(specie in species){
		print(paste0(" -> in the specie ",specie,"..."))
		d.testResults=d$`test results`
		d.sp=d.testResults[[specie]] %>% unlist
		d.sp["specie"]=as.character(specie)
		d.sp["Panther family"]=panf %>% as.character() #d$input$`file name` %>% strsplit(.,split=".",fixed=T) %>% .[[1]] %>% .[2]
		
		d.lst[[specie]]=d.sp
	}

	r=do.call("rbind",d.lst)
	r.lst[[f]]=r
}
species=allowedSpecies

# Generate CONSOLIDATED TABLE
consolidatedTable=do.call("rbind",r.lst)
consolidatedTable=consolidatedTable %>% as_tibble() %>% mutate(	`Corrected p-value`=as.numeric(`Corrected p-value`),
					     	`LRT`=as.numeric(`LRT`),
						`MLE`=as.numeric(`MLE`),
						`p-value`=as.numeric(`p-value`))







# Obtain stats and KaKs
pfDnds.v=unlist(dnds.lst)
pfDnds.tibble=tibble(`Panther family`=names(pfDnds.v),`Ka/Ks`=pfDnds.v)
consolidatedTable=consolidatedTable %>% mutate(zscore = (MLE - mean(MLE))/sd(MLE))

n_sig=consolidatedTable %>% filter(`Corrected p-value`< 0.05) %>% pull(`Panther family`) %>% unique %>% length
n_total=consolidatedTable %>% pull(`Panther family`) %>% unique %>% length



# Plotting time

#pdf("plotRelaxResults.pdf",width=6,height=5)
pdf("plotRelaxResults.pdf",width=6,height=6)


# first exploration:
# Scatterplot
consolidatedTable %>% ggplot(aes(x=log10(`p-value`),y=`MLE`,colour=specie))+geom_point(shape=1,size=2)+geom_vline(xintercept = log10(as.numeric(cutoff_PValue)), linetype="dotted",colour="red")+ggtitle(paste0("k MLE and significance"), subtitle=paste0("Total num. of families = ",n_total,"\nNum. of families with significant records = ",n_sig))+theme_minimal()+scale_colour_manual(values=colorEq)
# histogram k MLE by specie (all data)
consolidatedTable %>% ggplot(aes(x=`MLE`,fill=specie))+geom_histogram(position="dodge")+geom_vline(xintercept=1,linetype="dotted",colour="red")+ggtitle("k MLE distribution")+theme_minimal()+scale_fill_manual(values=colorEq)
# histogram k MLE by specie (significant)
consolidatedTable %>% filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% ggplot(aes(x=`MLE`,fill=specie))+geom_histogram(position="dodge")+geom_vline(xintercept=1,linetype="dotted",colour="red")+ggtitle(paste0("k MLE distribution (Adj. P-value < ",cutoff_PValue,")"))+theme_minimal()+scale_fill_manual(values=colorEq)
# histogram LTR (all data)
consolidatedTable %>% ggplot(aes(x=`LRT`,fill=specie))+geom_histogram(position="dodge")+ggtitle("LRT distribution")+theme_minimal()+scale_fill_manual(values=colorEq)
# histogram LTR (significant)
consolidatedTable %>% filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% ggplot(aes(x=`LRT`,fill=specie))+geom_histogram(position="dodge")+ggtitle(paste0("LRT distribution (Adj. P-value < ",cutoff_PValue,")"))+theme_minimal()+scale_fill_manual(values=colorEq)



# KaKs plot
consolidatedTable.onlySig=consolidatedTable %>% mutate(`k MLE`=`MLE`, `Significant`=ifelse(`Corrected p-value`< cutoff_PValue,T,F)) %>% filter(`Significant`==T)
consolidatedTable=consolidatedTable %>% mutate(`k MLE`=`MLE`, `Significant`=ifelse(`Corrected p-value`< cutoff_PValue,T,F))
pfDnds.tibble %>% ggplot(aes(x=`Ka/Ks`)) + geom_density()+ geom_vline(xintercept=1,colour="red", linetype="dashed")+theme_classic() + ggtitle("Rate of non-synonymous/Rate of synonymous subtitutions (Ka/Ks)")

# agregamos info de DnDs por pantherfam
consolidatedTable=consolidatedTable %>% left_join(pfDnds.tibble,by="Panther family")
# configuramos la var specie como factor, para que queden ordenados
consolidatedTable$specie=factor(consolidatedTable$specie, levels=speciesOrder)


# plot density KaKs per specie. (all together version)
consolidatedTable %>% filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% ggplot(aes(x=`Ka/Ks`,col=specie))+ geom_density()+ geom_vline(xintercept=1,colour="red", linetype="dashed")+theme_classic() + ggtitle("Rate of non-synonymous/Rate of synonymous subtitutions (Ka/Ks)",subtitle="per specie (only significant records)")+scale_colour_manual(values=colorEq)+xlim(0,1.5)


























# Comparamos All-vs-all
# y calculamos los tests estadisticos para las distribuciones de k MLE y KaKs

#d.wilcoxTest.KaKsOver1.lst=vector(mode="list",length(species))
#names(d.wilcoxTest.KaKsOver1.lst)=species
#d.tTest.KaKsOver1.lst=vector(mode="list",length(species))
#names(d.tTest.KaKsOver1.lst)=species

pvaluesRef.t.KaKsOver1=vector(mode="list",length(species))
names(pvaluesRef.t.KaKsOver1)=species
pvaluesRef.wilcox.KaKsOver1=vector(mode="list",length(species))
names(pvaluesRef.wilcox.KaKsOver1)=species


pvaluesRef.k.t=vector(mode="list",length(species))
names(pvaluesRef.k.t)=species
pvaluesRef.k.wilcox=vector(mode="list",length(species))
names(pvaluesRef.k.wilcox)=species

#d.wilcoxTest.KaKsAll.lst=vector(mode="list",length(species))
#names(d.wilcoxTest.KaKsAll.lst)=species
#d.tTest.KaKsAll.lst=vector(mode="list",length(species))
#names(d.tTest.KaKsAll.lst)=species
pvaluesRef.KaKsAll.t=vector(mode="list",length(species))
names(pvaluesRef.KaKsAll.t)=species
pvaluesRef.KaKsAll.wilcox=vector(mode="list",length(species))
names(pvaluesRef.KaKsAll.wilcox)=species

for(refSpecie in species){

	# K (all significant) # # #
	
	# reference data
	kMLE.ref=	consolidatedTable %>% 
		filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% 
		filter(specie==refSpecie) %>% 
		pull(`k MLE`)

	d.wilcoxTestK.lst=vector(mode="list",length(species))
	names(d.wilcoxTestK.lst)=species
	d.tTestK.lst=vector(mode="list",length(species))
	names(d.tTestK.lst)=species

	pvalues.k.t=rep(0,length(species))
	names(pvalues.k.t)=species

	pvalues.k.wilcox=rep(0,length(species))
	names(pvalues.k.wilcox)=species

	for(sp in species){
		print(paste0("==> k MLE : ",refSpecie," vs ",sp," <=="))
		# specie data
		kMLE.sp=	consolidatedTable %>% 
				filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% 
				filter(specie==as.character(sp)) %>% 
				pull(`k MLE`)

		# alternative
		al="two.sided"

		# wilcoxon rank-sum test (K)
		d.wilcoxTestK.lst[[sp]]=wilcox.test(kMLE.ref,
						   kMLE.sp,
						   alternative=al)
		# t-student test (K)
		d.tTestK.lst[[sp]]=t.test(kMLE.ref,
					 kMLE.sp,
					 alternative=al)
		
		# obtenemos y guardamos los pvalues
		pvalues.k.t[sp]=d.tTestK.lst[sp] %>% unlist %>% .[3] %>% as.numeric()    # T-test version
		pvalues.k.wilcox[sp]=d.wilcoxTestK.lst[sp] %>% unlist %>% .[2] %>% as.numeric() # Wilcox version
	}



	# KaKs (K > 1) # # #

	# reference data
	KaKs.ref=	consolidatedTable %>% 
			filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% 
			filter(`k MLE`>1) %>% 
			filter(specie==refSpecie) %>% 
			pull(`Ka/Ks`)

	d.wilcoxTestKaKs.lst=vector(mode="list",length(species))
	names(d.wilcoxTestKaKs.lst)=species
	d.tTestKaKs.lst=vector(mode="list",length(species))
	names(d.tTestKaKs.lst)=species

	pvalues.t.KaKsOver1=rep(0,length(species))
	names(pvalues.t.KaKsOver1)=species
	pvalues.wilcox.KaKsOver1=rep(0,length(species))
	names(pvalues.wilcox.KaKsOver1)=species

	for(sp in species){
		print(paste0("==> Ka/Ks : ",refSpecie," vs ",sp," <=="))

		# specie data
		KaKs.sp=	consolidatedTable %>% 
				filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% 
				filter(`k MLE`>1) %>% 
				filter(specie==as.character(sp)) %>% 
				pull(`Ka/Ks`)

		al="two.sided"
		# wilcoxon rank-sum test (KaKs; K>1)
		d.wilcoxTestKaKs.lst[[sp]]=wilcox.test(KaKs.ref,
						       KaKs.sp,
						       alternative=al)
		# t-student test (KaKs; K>1)
		d.tTestKaKs.lst[[sp]]=t.test(KaKs.ref,
					     KaKs.sp,
					     alternative=al)
		
		# obtenemos y guardamos los pvalues
		pvalues.t.KaKsOver1[sp]=d.tTestKaKs.lst[sp] %>% unlist %>% .[3] %>% as.numeric()    # T-test version
		pvalues.wilcox.KaKsOver1[sp]=d.wilcoxTestKaKs.lst[sp] %>% unlist %>% .[2] %>% as.numeric() # Wilcox version
		
	}
	# KaKsAll # # #

	# reference data
	KaKsAll.ref=	consolidatedTable %>% 
			filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% 
			filter(specie==refSpecie) %>% 
			pull(`Ka/Ks`)

	d.wilcoxTestKaKsAll.lst=vector(mode="list",length(species))
	names(d.wilcoxTestKaKsAll.lst)=species
	d.tTestKaKsAll.lst=vector(mode="list",length(species))
	names(d.tTestKaKsAll.lst)=species
	pvalues.KaKsAll.t=rep(0,length(species))
	names(pvalues.t.KaKsOver1)=species
	pvalues.KaKsAll.wilcox=rep(0,length(species))
	names(pvalues.wilcox.KaKsOver1)=species

	for(sp in species){
		print(paste0("==> Ka/Ks : ",refSpecie," vs ",sp," <=="))

		# specie data
		KaKsAll.sp=	consolidatedTable %>% 
				filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% 
				filter(`k MLE`>1) %>% 
				filter(specie==as.character(sp)) %>% 
				pull(`Ka/Ks`)

		al="two.sided"
		# wilcoxon rank-sum test (KaKsAll; K>1)
		d.wilcoxTestKaKsAll.lst[[sp]]=wilcox.test(KaKsAll.ref,
						       KaKsAll.sp,
						       alternative=al)
		# t-student test (KaKsAll; K>1)
		d.tTestKaKsAll.lst[[sp]]=t.test(KaKsAll.ref,
					     KaKsAll.sp,
					     alternative=al)
		
		# obtenemos y guardamos los pvalues
		pvalues.KaKsAll.t[sp]=d.tTestKaKsAll.lst[sp] %>% unlist %>% .[3] %>% as.numeric()    # T-test version
		pvalues.KaKsAll.wilcox[sp]=d.wilcoxTestKaKsAll.lst[sp] %>% unlist %>% .[2] %>% as.numeric() # Wilcox version
		
	}
	# contruimos listas de vectores
	pvaluesRef.t.KaKsOver1[[refSpecie]]=pvalues.t.KaKsOver1
	pvaluesRef.wilcox.KaKsOver1[[refSpecie]]=pvalues.wilcox.KaKsOver1
	pvaluesRef.k.t[[refSpecie]]=pvalues.k.t
	pvaluesRef.k.wilcox[[refSpecie]]=pvalues.k.wilcox
	pvaluesRef.KaKsAll.t[[refSpecie]]=pvalues.KaKsAll.t
	pvaluesRef.KaKsAll.wilcox[[refSpecie]]=pvalues.KaKsAll.wilcox
}

# consolidate data in matrix format
pvaluesMatrix.KaKsOver1.t=do.call("cbind",pvaluesRef.t.KaKsOver1)
pvaluesMatrix.KaKsOver1.wilcox=do.call("cbind",pvaluesRef.wilcox.KaKsOver1)
pvaluesMatrix.k.t=do.call("cbind",pvaluesRef.k.t)
pvaluesMatrix.k.wilcox=do.call("cbind",pvaluesRef.k.wilcox)
pvaluesMatrix.KaKsAll.t=do.call("cbind",pvaluesRef.KaKsAll.t)
pvaluesMatrix.KaKsAll.wilcox=do.call("cbind",pvaluesRef.KaKsAll.wilcox)

# and then in long table format
library(reshape2)
pvaluesTable.KaKsOver1.t=pvaluesMatrix.KaKsOver1.t %>% melt %>% tibble
colnames(pvaluesTable.KaKsOver1.t)=c("Specie","Ref. specie","p-Value")
pvaluesTable.KaKsOver1.wilcox=pvaluesMatrix.KaKsOver1.wilcox %>% melt %>% tibble
colnames(pvaluesTable.KaKsOver1.wilcox)=c("Specie","Ref. specie","p-Value")
pvaluesTable.k.t=pvaluesMatrix.k.t %>% melt %>% tibble
colnames(pvaluesTable.k.t)=c("Specie","Ref. specie","p-Value")
pvaluesTable.k.wilcox=pvaluesMatrix.k.wilcox %>% melt %>% tibble
colnames(pvaluesTable.k.wilcox)=c("Specie","Ref. specie","p-Value")
pvaluesTable.KaKsAll.t=pvaluesMatrix.KaKsAll.t %>% melt %>% tibble
colnames(pvaluesTable.KaKsAll.t)=c("Specie","Ref. specie","p-Value")
pvaluesTable.KaKsAll.wilcox=pvaluesMatrix.KaKsAll.wilcox %>% melt %>% tibble
colnames(pvaluesTable.KaKsAll.wilcox)=c("Specie","Ref. specie","p-Value")












save.image("test.RData")
# plot histogram KaKs per specie. (facets version)

# generate text data for plotting facets
mode <- function(x) {
   return(as.numeric(names(which.max(table(x)))))
}
data.text.kaks=	consolidatedTable %>% 
		filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% 
		group_by(specie) %>% 
		summarise(Mean=mean(`Ka/Ks`),Median=median(`Ka/Ks`), Mode=mode(`Ka/Ks`), SD=sd(`Ka/Ks`))

numSigPerSp.wilcox=pvaluesTable.KaKsAll.wilcox %>% filter(`p-Value` < cutoff_PValue) %>% group_by(Specie) %>% summarise(`N. sig. W.`=n()) %>% mutate(specie=Specie) %>% dplyr::select(-Specie)
numSigPerSp.t=pvaluesTable.KaKsAll.t %>% filter(`p-Value` < cutoff_PValue) %>% group_by(Specie) %>% summarise(`N. sig. T.`=n()) %>% mutate(specie=Specie) %>% dplyr::select(-Specie)
data.text.kaks=data.text.kaks %>% left_join(numSigPerSp.wilcox,by="specie")
data.text.kaks=data.text.kaks %>% left_join(numSigPerSp.t,by="specie")
nSigPantherfams=consolidatedTable %>% filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% dplyr::select(`specie`,`Panther family`) %>% group_by(specie) %>% summarise(N=n())
data.text.kaks=data.text.kaks %>% left_join(nSigPantherfams,by="specie")
data.text.kaks=data.text.kaks %>% mutate(label=paste0("N = ",`N`,"\nMean = ",format(Mean,digits=3),"\nSD = ",format(SD,digits=3),"\nMode = ",format(Mode,digits=3)), label2=paste0(" W = ",`N. sig. W.`,"\n T = ",`N. sig. T.`))

# plot facets
consolidatedTable %>% filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% ggplot(aes(x=`Ka/Ks`,fill=specie))+ geom_histogram(position="dodge",bins=50)+
	facet_wrap(~specie,  ncol=2)+
	geom_vline(xintercept=1,colour="red", linetype="dashed")+
	geom_vline(aes(xintercept=`Mean`),colour="grey", linetype="dashed",data=data.text.kaks)+
	theme_classic() +
	ggtitle("Ka/Ks distribution per specie",subtitle="(only significant records)")+
	scale_fill_manual(values=colorEq)+xlim(0,1.5)+
	ylab("Number of panther families")+
	geom_text(data=data.text.kaks,aes(x=0.9,y=35,label=label,group=specie),size=2.2,hjust=1) #+
#	geom_text(data=data.text.kaks,aes(x=1.5,y=35,label=label2,group=specie),size=2.2,hjust=1)



# plot facets (only k>1 (intensified selection))

data.text.kaks_Kover1=	consolidatedTable %>% 
			filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% 
			filter(`k MLE`>1) %>% 
			group_by(specie) %>% 
			summarise(
				  #N=n(`Ka/Ks`),
				  Mean=mean(`Ka/Ks`), 
				  Median=median(`Ka/Ks`), 
				  Mode=mode(`Ka/Ks`), 
				  SD=sd(`Ka/Ks`))


AchCompPval.wilcox=pvaluesTable.KaKsOver1.wilcox %>% filter(`Ref. specie`=="Austrolebias_charrua") %>% select(Specie,`p-Value`) %>% mutate(specie=`Specie`, `p-Value(W)`=`p-Value`) %>% dplyr::select(-Specie,-`p-Value`)
data.text.kaks_Kover1=data.text.kaks_Kover1 %>% left_join(AchCompPval.wilcox,by="specie")
AchCompPval.t=pvaluesTable.KaKsOver1.t %>% filter(`Ref. specie`=="Austrolebias_charrua") %>% select(Specie,`p-Value`) %>% mutate(specie=`Specie`, `p-Value(T)`=`p-Value`) %>% dplyr::select(-Specie,-`p-Value`)
data.text.kaks_Kover1=data.text.kaks_Kover1 %>% left_join(AchCompPval.t,by="specie")

numSigPerSp.wilcox=pvaluesTable.KaKsOver1.wilcox %>% filter(`p-Value` < cutoff_PValue) %>% group_by(Specie) %>% summarise(`N. sig. W.`=n()) %>% mutate(specie=Specie) %>% dplyr::select(-Specie)
numSigPerSp.t=pvaluesTable.KaKsOver1.t %>% filter(`p-Value` < cutoff_PValue) %>% group_by(Specie) %>% summarise(`N. sig. T.`=n()) %>% mutate(specie=Specie) %>% dplyr::select(-Specie)

data.text.kaks_Kover1=data.text.kaks_Kover1 %>% left_join(numSigPerSp.wilcox,by="specie")
data.text.kaks_Kover1=data.text.kaks_Kover1 %>% left_join(numSigPerSp.t,by="specie")
nSigPantherfams=consolidatedTable %>% filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% filter(`k MLE`>1) %>% dplyr::select(`specie`,`Panther family`) %>% group_by(specie) %>% summarise(N=n())
data.text.kaks_Kover1=data.text.kaks_Kover1 %>% left_join(nSigPantherfams,by="specie")

data.text.kaks_Kover1=data.text.kaks_Kover1 %>% mutate(label=paste0("N = ",`N`,"\nMean = ",format(Mean,digits=3),"\nSD = ",format(SD,digits=3),"\nMode = ",format(Mode,digits=3)), label2=paste0(" W = ",`N. sig. W.`,"\n  T = ",`N. sig. T.`,"\np-Value(W) = ",signif(`p-Value(W)`,digits=1),"\np-Value(T) = ",signif(`p-Value(T)`,digits=1)))

consolidatedTable %>% filter(`Corrected p-value` < as.numeric(cutoff_PValue)) %>% filter(`k MLE`>1) %>% 
	ggplot(aes(x=`Ka/Ks`,fill=specie))+ geom_histogram(position="dodge",bins=50)+
	facet_wrap(~specie,  ncol=2)+
	geom_vline(xintercept=1,colour="red", linetype="dashed")+
	geom_vline(aes(xintercept=`Mean`),colour="grey", linetype="dashed",data=data.text.kaks_Kover1)+
	theme_classic() +
	ggtitle("Ka/Ks distribution per specie",subtitle="(only significant records having K > 1)")+
	ylab("Number of panther families")+
	scale_fill_manual(values=colorEq)+xlim(0,1.5)+
	geom_text(data=data.text.kaks_Kover1,aes(x=0.9,y=10,label=label,group=specie),size=2.2,hjust=1)+
	geom_text(data=data.text.kaks_Kover1,aes(x=1.5,y=10,label=label2,group=specie),size=2.2,hjust=1)

















# plot pvalue heatmaps

print("ready to print heatmaps")
#  KaKs

pvaluesTable.KaKsOver1.wilcox %>% ggplot(aes(x=`Specie`,y=`Ref. specie`,fill=`p-Value`))+geom_tile()+geom_text(aes(label=format(signif(`p-Value`,digits=3),digits=3)), size=2)+theme_minimal()+scale_fill_gradient(low = "#ffeda0", high = "#f03b20")+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),legend.position="bottom")+ggtitle("All-vs-all comparison of Ka/Ks distributions (k > 1)",subtitle="Two-sided Wilcoxon rank-sum test")
pvaluesTable.KaKsOver1.t %>% ggplot(aes(x=`Specie`,y=`Ref. specie`,fill=`p-Value`))+geom_tile()+geom_text(aes(label=format(signif(`p-Value`,digits=3),digits=3)), size=2)+theme_minimal()+scale_fill_gradient(low = "#ffeda0", high = "#f03b20")+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),legend.position="bottom")+ggtitle("All-vs-all comparison of Ka/Ks distributions (k > 1)",subtitle="Two-sided t-student test")

#  k MLE

pvaluesTable.k.wilcox %>% ggplot(aes(x=`Specie`,y=`Ref. specie`,fill=`p-Value`))+geom_tile()+geom_text(aes(label=format(signif(`p-Value`,digits=3),digits=3)), size=2)+theme_minimal()+scale_fill_gradient(low = "#ffeda0", high = "#f03b20")+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),legend.position="bottom")+ggtitle("All-vs-all comparison of k distributions",subtitle="Two-sided Wilcoxon rank-sum test")
pvaluesTable.k.t %>% ggplot(aes(x=`Specie`,y=`Ref. specie`,fill=`p-Value`))+geom_tile()+geom_text(aes(label=format(signif(`p-Value`,digits=3),digits=3)), size=2)+theme_minimal()+scale_fill_gradient(low = "#ffeda0", high = "#f03b20")+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),legend.position="bottom")+ggtitle("All-vs-all comparison of k distributions",subtitle="Two-sided t-student test")

numSigPerSp.wilcox.kaks=pvaluesTable.KaKsOver1.wilcox %>% filter(`p-Value` < cutoff_PValue) %>% group_by(Specie) %>% summarise(`N. sig.`=n()) %>% mutate(specie=Specie) %>% dplyr::select(-Specie)
numSigPerSp.t.kaks=pvaluesTable.KaKsOver1.t %>% filter(`p-Value` < cutoff_PValue) %>% group_by(Specie) %>% summarise(`N. sig.`=n()) %>% mutate(specie=Specie) %>% dplyr::select(-Specie)















# violin plot
p.v=consolidatedTable %>% ggplot(aes(x=`k MLE`,y=specie, fill=specie))+
	geom_violin(aes(linetype=NA,alpha=0.2))+
	geom_point(data=consolidatedTable.onlySig,
		   aes(colour=`Significant`),
		   position="jitter",
		   size=2,
		   alpha=0.15,
		   colour="black",
		   #linetype=NA,
		   stroke=0,
		   shape=20)+
	geom_vline(xintercept=1, colour="red",linetype="dotted",alpha=1)+
	#theme_minimal()+
	theme_classic()+
	ggtitle("k MLE distribution",)+
	scale_fill_manual(values=colorEq)+
	theme(	legend.position="none",
		axis.text.y = element_text(angle = 45, vjust = -1, hjust=1),
		axis.title.y= element_blank())

p.v
p.legend=p.v+theme(legend.position="right")
p.legend=g_legend(p.legend)
grid.arrange(p.legend)

	# sacamos algunos parametros estadisticos:
	avgMLE_excluding_neutral=consolidatedTable %>% filter(`k MLE` > 1) %>% group_by(specie) %>% summarise(avgMLE=mean(MLE))
	avgMLE=consolidatedTable %>% group_by(specie) %>% summarise(avgMLE=mean(MLE))

dev.off()

























# obtenemos listas

# en formato de archivos de text (lst)
for(sp in species){
	consolidatedTable.onlySig %>% filter(`k MLE`>1) %>% filter(specie==sp) %>% dplyr::select(`Panther family`) %>% write.table(paste0("pantherfams_k_over_1.",sp,".lst"), row.names=F,col.names=F, quote=F)
}


# y cruzamos con panther db (version PANTHER.db, funciona, pero hay un problema al correrlo en el docker)
library(PANTHER.db)
panelPanther=consolidatedTable.onlySig %>% pull(`Panther family`) %>% unique
panelPanther=panelPanther %>% gsub("_",":",.)
panelPanther.tibble=PANTHER.db::select(PANTHER.db, keys=panelPanther, columns=c("CLASS_TERM","FAMILY_TERM","SUBFAMILY_TERM","GOSLIM_ID"),keytype="FAMILY_ID")  %>% tibble
#panelPanther.tibble=panelPanther.tibble %>% mutate(`Panther family`=gsub(":","_",FAMILY_ID)) %>% dplyr::select(-FAMILY_ID)
#consolidatedTable.onlySig.wPantherData=consolidatedTable.onlySig %>% left_join(panelPanther.tibble,by="Panther family")

#consolidatedTable.onlySig.wPantherData.nr=consolidatedTable.onlySig.wPantherData %>% distinct(specie,`Panther family`,GOSLIM_ID,FAMILY_TERM,`k MLE`)
#for(sp in species){
#	consolidatedTable.onlySig.wPantherData.nr %>% filter(`k MLE`>1) %>% filter(specie==sp) %>% dplyr::select(GOSLIM_ID) %>% write.table(paste0("pantherfams_k_over_1.",sp,".go"), row.names=F,col.names=F, quote=F)
#}


# version con anotacion provista como .tab
pantherAnn=read.table(pantherTab,sep="\t") %>% tibble()
colnames(pantherAnn)=c("specie","transcriptID","Panther family","Name","HMM evalue","HMM bitscore","Range")
pantherAnn.simple=pantherAnn %>% distinct(`Panther family`,Name) %>% mutate(`Panther family`=gsub(":","_",`Panther family`))
consolidatedTable=consolidatedTable %>% left_join(.,pantherAnn.simple,by="Panther family")


save.image("plotRelaxResults.RData")

