#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# ##To execute script at the shell. 
# Rscript --vanilla REGCAP_RC2_merge.R <directory name>

setwd(args[1])

rc2<-list.files()[grep(".rc2", list.files())] #lists all the files in the directory

cc<-read.table(rc2[1],sep="\t", header=TRUE) #for use identifying the col names

summaries<-matrix(NA,ncol=ncol(cc),nrow=length(rc2))
rownames(summaries)<-c(rc2)
colnames(summaries)<-colnames(cc)

#Creates a data frame with the exploded file names as separate columns for analysis. Amended to data matrixes prior to writing to file. 
rows<-as.data.frame(unlist(strsplit(rc2,'_')))
ncol<-length(strsplit(rc2[1], '_')[[1]])
rowsn<-data.frame(matrix(rows[,1], ncol=ncol, nrow=length(rc2), byrow=TRUE), stringsAsFactors=TRUE)
rowsn[,ncol]<-data.frame(matrix(unlist(strsplit(as.character(rowsn[,ncol]), ".", fixed=TRUE)), ncol=2, nrow=length(rc2), byrow=TRUE), stringsAsFactors=TRUE)[,1]
rowsn<-cbind(as.factor(rc2), rowsn)

#fills in the empty matrix "summaries" with the rc2 file contents

for(i in 1:length(rc2)){
summaries[i,]<-as.matrix(read.table(rc2[i],sep="\t",skip=1))
}

name<-strsplit(getwd(), '/')[[1]][length(strsplit(getwd(), '/')[[1]])]
name<-paste(name, '_RC2_merge.csv', sep="")

write.csv(cbind(rowsn, summaries),name, row.names=FALSE)