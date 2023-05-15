#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 3rd 2022 ###
### Merge histone modification annotations ###
library(optparse)

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input list of histone modifications", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
		      help="output file", metavar="character"),
	make_option(c("-w", "--workdir"), type="character", default=NULL, 
		  	  help="directory containing HM_hmname.txt", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

hm_list<-read.table(opt$input,header = T,stringsAsFactors = F)
data.table<-data.frame(apply(data.frame(hm_list[,1]),1,function(x){
	indata<-read.table(paste0(opt$workdir,"/","HM_",x,".txt"),header =T,stringsAsFactors = F)
	indata[,2]
}))
indata<-read.table(paste0(opt$workdir,"/","HM_",hm_list[1,1],".txt"),header =T,stringsAsFactors = F)
data.table=cbind(indata[,1],data.table)
names(data.table)<-c("id",hm_list[,1])
write.table(data.table,opt$output,row.names=F,col.names=T,quote=F,sep="\t")
