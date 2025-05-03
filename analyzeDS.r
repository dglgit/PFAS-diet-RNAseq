suppressPackageStartupMessages({
library(readr)
library(magrittr)
library(tidyverse)
library(org.Mm.eg.db)
library(fgsea)
library(qusage)
library(dplyr)
library(DT)
})
df=read_csv("deseq2-results-hg-star.csv")
df=filter(df,df[['baseMean']]>0)
df=(df %>% filter(!is.na(df[['pvalue']]) & !is.na(df[['stat']])))
df[is.infinite(df[['stat']]) & df[['stat']]>0,'stat']=max(df[is.finite(df[['stat']]),'stat'])*1000
df[is.infinite(df[['stat']]) & df[['stat']]<0,'stat']=min(df[is.finite(df[['stat']]),'stat'])*1000


#print(nrows(df))
mmgmt="m2.cp.v2024.1.Mm.symbols.gmt"
hggmt="c1.all.v2024.1.Hs.symbols.gmt"
pathways=read.gmt(hggmt)
rs2symbol=AnnotationDbi::select(org.Mm.eg.db, key=df[["row"]], columns="SYMBOL",keytype="REFSEQ")
res=merge(df,rs2symbol,by.x="row",by.y="REFSEQ")
print('head res')
print(head(res))
selected=res[c('SYMBOL','stat')]
selected=selected[!is.na(selected[['SYMBOL']]),]
selected=selected[selected[['SYMBOL']]!='',]
print(nrow(selected))
selected=(selected %>% group_by(SYMBOL) %>% summarise(stat=mean(stat)))
print(nrow(selected))
print(head(selected))
ranks=deframe(selected)
print(c(length(unique(names(ranks))),length(names(ranks))))
#print(ranks)
print('rank crap')
ranks=ranks[!is.na(names(ranks))]
ranks=ranks[names(ranks)!='']
print('rank crap over')
#ranks=ranks[is.finite(ranks)]
if(FALSE){
sym2num=c()
num2sym=c()
count=1
pathwaysNum=c()
newRanks=c()
for(sym in names(ranks)){
	if(!(sym %in% names(sym2num))){
		sym2num[sym]=count
		num2sym[count]=sym
		#pathwaysNum[sym]=count
		
		count=count+1
	}
	print(paste(sym2num[sym],sym,sep="|"))
	newRanks[sym2num[sym]]=ranks[sym]
	print(newRanks)
	if(count>4){
		stopifnot(FALSE)
	}
	#readline(prompt="Press [enter] to continue")

	
}

for (p in pathways){
	pathwaysNum[p]=list()
	for (sym in p){
		if(!(sym %in% names(sym2num))){
				sym2num[sym]=count
				num2sym[count]=sym
				count=count+1
		}
		append(pathwaysNum[[sym]],sym2num[sym])
	}
}
getSyms=function(nums){
	res=list()
	for(n in nums){
		append(res,num2sym[n])
	}
	return (res)
}
getSym=function(num) num2sym[num]
getNum=function(sym) sym2num[sym]
}
print('ball1')
fgseaRes=fgsea(pathways=pathways,stats=ranks,nperm=1000)
print('ball2')
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
print(fgseaResTidy)
print(fgseaResTidy %>% arrange(padj))
fgseaResTidy[['leadingEdge']]=vapply(fgseaResTidy[['leadingEdge']],paste, collapse=', ',character(1L))
#write.csv(as.data.frame(fgseaResTidy),"fgsea-results1.csv")
# Show in a nice table:
write.table(fgseaResTidy %>% arrange(padj),file="fgsea-results1.csv")
#fgseaResTidy %>%
#  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
#  arrange(padj) %>%
#  DT::datatable()
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
