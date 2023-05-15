#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 3rd 2022 ###
### do normalizations to the coverage matrix ###
library(optparse)
library(DESeq2)

option_list = list(
    make_option(c("-b", "--bed"), type="character", default=NULL, 
                help="input bed file", metavar="character"),
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="input list of sites", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
		            help="output file", metavar="character"),
    make_option(c("-w", "--window"), type="integer", default=NULL, 
                help="window", metavar="character"),
    make_option(c("-s", "--step"), type="integer", default=NULL, 
                help="step", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#### DESeq2 ratio of median for normalizing library size difference between samples ####
## plus pseudo count +1 ##
input<-read.table(opt$bed,header = T,stringsAsFactors = F)
input<-input[order(input[,1],input[,2]),]
DESeq_input<-input[4:dim(input)[2]]

## first run of library size calculation ##
#lib_size_factor_pre<-estimateSizeFactorsForMatrix(counts = DESeq_input)

## after inspecting the distribution of lib size, a threshold at 0.125 is decided for filtering out low DNA amount samples ~200 ##
#DESeq_input<-DESeq_input[,lib_size_factor_pre>0.125]
dds <- DESeqDataSetFromMatrix(countData = DESeq_input,
                              colData = data.frame(Name = c(colnames(DESeq_input)),
                                                   Type = c(rep(c("S0","S1"),length.out=dim(DESeq_input)[2]))),
                              design = ~Type )
dds <- estimateSizeFactors(dds)
input_norm<-cbind(input[,1:3],counts(dds,normalized = T))
#lib_size_factor_post<-dds$sizeFactor
#### Normalization END ####
## Round the table to 5 decimal numbers ##
#input_norm[,4:dim(input_norm)[2]]<-round(input_norm[,4:dim(input_norm)[2]],4)
#input_round<-cbind(input_norm[,1:3],round(input_norm[,4:dim(input_norm)[2]]))
## compute rank of rank ##
input_rank<-t(apply(input_norm[,4:dim(input_norm)[2]],1,function(x){as.numeric(sprintf("%.3f",rank(x)/length(x)))}))
input_rank_rank<-cbind(input_norm[,1:3],apply(input_rank,2,function(x){as.numeric(sprintf("%.3f",rank(x)/length(x)))}))
names(input_rank_rank)<-names(input_norm)
samlist<-read.table(opt$input,header = T,stringsAsFactors = F,sep="\t")
result<-t(apply(as.matrix(samlist[,c("sample","chr","PAM_coordinate")]),1,function(x){
  s1<-as.numeric(x[3])-as.numeric(x[3])%%opt$step-opt$step
  s2<-s1+opt$step
 if(x[1] %in% names(input_rank_rank))
	 { 
		if(!(s1 %in% names(table(input_rank_rank[input_rank_rank[,1]==x[2],2]))))
		{
			return(c(NA,input_rank_rank[(input_rank_rank[,1]==x[2])&(input_rank_rank[,2]==s2),x[1]]))
		} else if(!(s2 %in% names(table(input_rank_rank[input_rank_rank[,1]==x[2],2])))) {
			return(c(input_rank_rank[(input_rank_rank[,1]==x[2])&(input_rank_rank[,2]==s1),x[1]],NA))
		} else {
			return(c(input_rank_rank[(input_rank_rank[,1]==x[2])&(input_rank_rank[,2]==s1),x[1]],input_rank_rank[(input_rank_rank[,1]==x[2])&(input_rank_rank[,2]==s2),x[1]]))
		}
	 }
	 else
	 {
		 return(c(NA,NA))
	 }
}))
result<-cbind(samlist[,c("sample","chr","PAM_coordinate")],result)
colnames(result)<-c("id","chr","coord","RR1","RR2")
write.table(result,opt$output,col.names = T, row.names=F, quote = F,sep = "\t")