#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 4rd 2022 ###
### Prepare feature set table for machine learning training ###
library(optparse)

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input matrix of features", metavar="character"),
	make_option(c("-e", "--edoutcome"), type="character", default=NULL, 
		      help="observed editing outcome profiles", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
		      help="output file", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#### Load sub-functions ####
sequence_to_onehot_matrix<-function(input_sequence){
  output<-apply(as.matrix(strsplit(input_sequence,"")[[1]]),1,function(x){
    y<-c(0,0,0,0)
    if(x=="A"){y[1]<-1}
    else if(x == "T"){y[2]<-1}
	else if(x == "C"){y[3]<-1}
    else if(x == "G"){y[4]<-1}
    y
  })
  rownames(output)<-c("A","T","C","G")
  colnames(output)<-1:ncol(output)
  output
}
sequence_to_onehot_matrix_dimer<-function(input_sequence){
  string<-strsplit(input_sequence,"")[[1]]
  nucleotide<-c("A","T","C","G")
  output<-data.frame(t(as.matrix(rep(0,16*(length(string)-1)))))
  j<-1
  for(i in 1:(length(string)-1))
  {
    for(first in nucleotide)
    {
      for(second in nucleotide)
      {
        if((first==string[i])&(second==string[i+1]))
        {
          output[j]<-1
        }
        colnames(output)[j]<-paste0(first,second,i)
        j<-j+1
      }
    }
  }
  as.matrix(output)
}
get_gc_content<-function(seq){
  seq.split<-strsplit(seq,split="")
  seq.split<-seq.split[[1]]
  seq.freqtab<-data.frame(table(seq.split))
  seq.gc<-sum(seq.freqtab[seq.freqtab$seq.split %in% c("G","C"),"Freq"])/sum(seq.freqtab[,"Freq"])
  seq.gc
}
#### END ####

#### Load sample table ####
sample.data<-read.table(opt$input,sep="\t",header=T,stringsAsFactors = F)
### compute derived features form sequences: GC, one-hot matrixes ###
sample.data$guide_GC<-apply(as.matrix(substr(sample.data[,"target_sequence_30mer"],6,25)),1,get_gc_content)
sample.data.guide.onehot<-t(apply(as.matrix(sample.data[,"target_sequence_30mer"]),1,function(x){as.numeric(sequence_to_onehot_matrix(x))}))
colnames(sample.data.guide.onehot)<-paste0(rep(c("A","T","C","G"),30),rep(25:-4,each=4))
## di-nucleotides of 30mer ##
sample.data.guide.onehot.dimer<-data.frame(t(apply(as.matrix(sample.data[,"target_sequence_30mer"]),1,function(x){sequence_to_onehot_matrix_dimer(x)})))
nucleotide<-c("A","T","C","G")
j<-1
for(i in 24:-4)
{
  for(first in nucleotide)
  {
    for(second in nucleotide)
    {
      colnames(sample.data.guide.onehot.dimer)[j]<-paste0(first,second,i)
      j<-j+1
    }
  }
}
editing_outcome<-read.table(opt$edoutcome,header=T,stringsAsFactors = F,sep="\t")

sample.data<-cbind(sample.data,sample.data.guide.onehot,sample.data.guide.onehot.dimer,editing_outcome[,c("Editing_outcome_NE","Editing_outcome_SV","Editing_outcome_LD","Editing_outcome_NRecT")])
sample.data<-sample.data[!is.na(sample.data$Editing_outcome_NE),]
sample.data<-sample.data[!duplicated(sample.data[,c("target_sequence_30mer","Editing_outcome_NE","Editing_outcome_SV","Editing_outcome_LD","Editing_outcome_NRecT"),]),!names(sample.data) %in% c("target_sequence_30mer")]

write.table(sample.data,opt$output,row.names = F,col.names = T,quote = F,sep="\t")
