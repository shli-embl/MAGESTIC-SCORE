#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 10 2022 ###
### Merge bed files into one ###
library(optparse)

option_list = list(
    make_option(c("-f", "--folder"), type="character", default=NULL, 
              help="input folder", metavar="character"),
    make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input file containing sample names", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
		          help="output file", metavar="character"),
    make_option(c("-w", "--window"), type="integer", default=NULL, 
                help="window size", metavar="character"),
    make_option(c("-s", "--step"), type="integer", default=NULL, 
                help="step size", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sample.list<-read.table(opt$input,header = T, stringsAsFactors = F,sep="\t")

init<-0
for(sample in sample.list[,1])
{
  tmp.table<-read.table(paste0(opt$folder,"/",sample,".sample.w",opt$window,".s",opt$step,".bed"),comment.char = "",sep = "\t", stringsAsFactors = F,header = F)
  tmp.table$index<-paste(tmp.table[,1],tmp.table[,2],tmp.table[,3],sep = "_")
  if(init==0)
  {
    merge.table <- tmp.table
    init<-1
    names(merge.table)[ncol(merge.table)-1]<-sample
  }
  else
  {
    merge.table <- merge(merge.table,tmp.table[,4:5],by = "index")
    names(merge.table)[ncol(merge.table)]<-sample
  }
}
write.table(merge.table[,!names(merge.table) %in% "index"],opt$output,col.names = T, row.names = F, quote = F,sep = "\t")