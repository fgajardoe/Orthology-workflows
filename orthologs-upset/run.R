library(ggupset)
library(tidyverse)

args=commandArgs(T)
outprefix=as.character(args[1])
argsLsts=args[2:length(args)]

d.lst=vector(mode="list", length(argsLsts))

species=argsLsts %>% gsub(".PANTHER","",x=.)

names(d.lst)=species
for(sp in species){
	d.lst[[sp]]=read.table(paste0(sp,".PANTHER"),sep="\t",col.names=c("ID","Panther superfamily","Description","evalue","score","unknown")) %>%
		tibble %>%
		mutate(specie=sp)

}
d=do.call("rbind",d.lst)
d.2=d %>% select(`Panther.superfamily`,specie) %>% distinct(Panther.superfamily,specie) %>% mutate(`Panther family`=Panther.superfamily) %>% select(-Panther.superfamily) %>% group_by(`Panther family`) %>% summarise(speciesArr=list(`specie`))


intersections=as.list(species)
intersections[[length(intersections)+1]]=species
pal.v=read.table("colourFileAustrolebiasPanel.tab",comment.char="@") %>% tibble %>% pull(V2,V1)
pdf(file=paste0(outprefix,".upset.pdf"),width=21,height=7)
d.2 %>% distinct(`Panther family`,speciesArr) %>%
	ggplot(aes(x=speciesArr))+
	geom_bar()+
	geom_text(stat='count', aes(label=after_stat(count),angle=90),size=3, hjust=-0.5) +
	scale_x_upset(order_by="degree",reverse=T,sets=species)+ #,intersections=intersections)+
	theme_minimal()+
	theme_combmatrix(combmatrix.panel.line.size = 0, combmatrix.label.make_space = FALSE)+
	xlab("Intersections")+
	theme(axis.text.y=element_text(size=12))

dev.off()

save.image(paste0(outprefix,".RData"))
