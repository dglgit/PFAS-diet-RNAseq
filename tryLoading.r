library(readr)
suppressMessages(library("DESeq2"))

# Testing out how to read TSVs
#data=read_tsv("Galaxy322-[query results on data 320].tabular")
#df_data <- as.data.frame(data)
#data2=read_tsv("Galaxy188-[featureCounts on data 112 and data 44_ Counts].tabular")
#df_data2 = as.data.frame(data2)
#dsData=as.data.frame(read_tsv("Galaxy320-[Join on data 202 and data 197].tabular"))

metaData=as.data.frame(read_tsv("metadata.tsv"))
#data.frame(id=c("SRR16155199","SRR16155211"),dex=c("control", "treated"))
print(metaData)
bigList=list()
for (i in 1:nrow(metaData)){
	if (metaData[i,'dex']=='treated'){
		a=as.data.frame(read_tsv(paste("pfos-mm/",metaData[i,'id'],".tabular",sep='')))
		bigList=append(bigList,list(a))
	}else{
		a=as.data.frame(read_tsv(paste("control-mm/",metaData[i,'id'],".tabular",sep='')))
		bigList=append(bigList,list(a))
	}
}
print('read')
#print(bigList[1])
print(colnames(bigList[[1]]))
print(colnames(bigList[[2]]))
print(length(bigList))
bigDf=merge(bigList[[1]],bigList[[2]],by='Geneid')
for(i in 3:length(bigList)){
	bigDf=merge(x=bigDf,y=bigList[[i]],by='Geneid')
}
print('merged')
dds=DESeqDataSetFromMatrix(countData=bigDf,colData=metaData,design=~dex, tidy=TRUE)
dds=DESeq(dds)
res=as.data.frame(results(dds,tidy=TRUE))
res=res[!is.na(res[['padj']]),]
head(res)
#res=filter(res, !is.na(res[['padj']]))
write.csv(res,"deseq2-results1.csv")
