#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 10 2022 ###
### Merge annotation tables into one ###
library(optparse)

option_list = list(
    make_option(c("-f", "--folder"), type="character", default=NULL, 
              help="input folder", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
		      help="output file", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

init<-0
for(file in list.files(opt$folder))
{
  tmp.table<-read.table(paste0(opt$folder,"/",file),comment.char = "",stringsAsFactors = F,header = T)
  if(init==0)
  {
    merge.table <- tmp.table
    init<-1
  }
  else
  {
    merge.table <- merge(merge.table,tmp.table,by.x = 1, by.y = 1)
  }
}
names(merge.table)[1] <- "id"
write.table(merge.table,opt$output,col.names = T, row.names = F, quote = F,sep = "\t")