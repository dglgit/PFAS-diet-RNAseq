library(readr)
suppressMessages(library("DESeq2"))
fname="pfos-hg-star/SRR16155213.tabular"
#headers=strsplit(readLines(fname,n=1),'\t')
headers=c("geneID","counts_unstrand","counts_firstStrand","counts_secondStrand")
print(headers)
data = read.table(fname,header=F,sep='\t',skip=4,col.names=headers)
data = data[c('geneID',"counts_unstrand")]
print(head(data))
readStar=function(fname){
    print(fname)
    #print(strsplit(fname,'[.]')[[2]])
    entryName=strsplit(strsplit(fname,'[.]')[[1]][1],'[/]')[[1]][2]
    headers=c('gene_id',"counts_unstrand","counts_firstStrand","counts_secondStrand")
    data = read.table(fname,header=F,sep='\t',skip=4,col.names=headers)
    data = data[c('gene_id',"counts_unstrand")]
    names(data)[names(data)=='counts_unstrand']=entryName
    return (data)
}

metaData=as.data.frame(read_tsv("metadata.tsv"))
bigData=list()
for(i in 1:nrow(metaData)){
    
    if(metaData[i,'dex']=='treated'){
        a=readStar(paste("pfos-hg-star/",metaData[i,'id'],'.tabular',sep=''))#read.table(paste("pfos-hg-star/",metaData[i,'id'],'.tabular'),header=F,sep='\t',skip=4,col.names=headers)
        bigData=append(bigData,list(a))
    } else{
        a=readStar(paste("control-hg-star/",metaData[i,'id'],'.tabular',sep=''))
        bigData=append(bigData,list(a))
    }
}
print(length(bigData))
bigDf=merge(bigData[[1]],bigData[[2]],by='gene_id')
for(i in 3:length(bigData)){
	bigDf=merge(x=bigDf,y=bigData[[i]],by='gene_id')
}
print(head(as.data.frame(bigData[[1]])))
print(head(bigDf))
dds=DESeqDataSetFromMatrix(countData=bigDf,colData=metaData,design=~dex, tidy=TRUE)
dds=DESeq(dds)
res=as.data.frame(results(dds,tidy=TRUE))
res=res[!is.na(res[['padj']]),]
head(res)
write.csv(res,"deseq2-results-hg-star.csv")
