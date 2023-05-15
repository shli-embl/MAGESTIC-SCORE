#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 3rd 2022 ###
### Extract histone modification data for a list of target sites ###
library(optparse)

option_list = list(
    make_option(c("-b", "--bed1"), type="character", default=NULL, 
              help="bed file for histone modification signal", metavar="character"),
	make_option(c("-e", "--bed2"), type="character", default=NULL, 
		      help="bed file for background input", metavar="character"),
    make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input list of sites", metavar="character"),
    make_option(c("-o", "--outprefix"), type="character", default=NULL, 
		      help="output files prefix", metavar="character"),
	make_option(c("-l", "--label"), type="character", default=NULL, 
		  	  help="label of the modification", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

hist_mod<-read.table(opt$bed1,header = F,stringsAsFactors = F)
input_background<-read.table(opt$bed2,header = F,stringsAsFactors = F)
hist_mod<-cbind(hist_mod,input_background[,4])
rm(input_background)
names(hist_mod)<-c("chr","start","end","hist_mod","input")
hist_mod <- cbind(hist_mod[,1:3],log2((hist_mod$hist_mod+25)/(hist_mod$input+25)))
names(hist_mod)<-c("chr","start","end",opt$label)
rownames(hist_mod)<-paste(hist_mod$chr,hist_mod$start,sep="_")
#### Output histone modifications of samples ####
inputmeta<-read.table(opt$input,header=T,stringsAsFactors = F,sep="\t")
### get chr_start ###
window1_loc<-paste0(inputmeta[,"chr"],"_",as.integer(inputmeta[,"coord"]-inputmeta[,"coord"]%%50-50))
window2_loc<-paste0(inputmeta[,"chr"],"_",as.integer(inputmeta[,"coord"]-inputmeta[,"coord"]%%50))
### edited to account for edge cases ###
result<-cbind(inputmeta[,"sample"],
data.frame(apply(data.frame(window1_loc,window2_loc),1,function(x){
	if(!(x[1] %in% rownames(hist_mod))){
		out <- hist_mod[x[2],opt$label]
	} else if(!(x[2] %in% rownames(hist_mod))){
		out <- ist_mod[x[1],opt$label]
	} else {
		out <- (hist_mod[x[1],opt$label]+hist_mod[x[2],opt$label])/2
		
	}
	out
})))
names(result)<-c("id",opt$label)
write.table(result,opt$outprefix,row.names=F,col.names=T,quote=F,sep="\t")
