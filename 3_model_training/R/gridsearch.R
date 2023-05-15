#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 18th 2022 ###
### concatenate results ###
library(optparse)
library(dplyr)

### a string of pseudo seeds to ensure the result can be reproduced ###

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input dir", metavar="character"),
  make_option(c("--features"), type="character", default=NULL, 
              help="feature sets", metavar="character"),
  make_option(c("--outcome"), type="character", default=NULL, 
              help="outcome names", metavar="character"),
  make_option(c("--method"), type="character", default=NULL, 
              help="method for treating data imbalance", metavar="character"),
  make_option(c("--metric"), type="character", default=NULL, 
			  help="metric", metavar="character"),
  make_option(c("--model"), type="character", default=NULL, 
              help="gbm, glmnet or rf", metavar="character"),
  make_option(c("--size"), type="character", default=NULL, 
              help="proportion of major class to be preserved (from 0 to 1)", metavar="character"),
  make_option(c("--shrinkage"), type="character", default=NULL, 
              help="shrinkage", metavar="character"),
  make_option(c("--intdepth"), type="character", default=NULL, 
              help="interaction depth", metavar="character"),
  make_option(c("--ntrees"), type="character", default=NULL, 
              help="ntrees", metavar="character"),
  make_option(c("--alpha"), type="character", default=NULL, 
              help="alpha", metavar="character"),
  make_option(c("--lambda"), type="character", default=NULL, 
              help="lambda", metavar="character"),
  make_option(c("--rfntrees"), type="character", default=NULL, 
              help="rfntrees", metavar="character"),
  make_option(c("--mtry"), type="character", default=NULL, 
			  help="mtry", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

result<-data.frame(method_training="gbm",predict_class="X",feature_number=0,method_dataimbalance="undersample",undersampleSize=0,metric = 0,shrinkage=NA,inter_depth=NA,ntrees=NA,alpha=NA,lambda=NA,rfntrees=NA,mtry=NA,
                   Sensitivity=0,Specificity=0,Pos_pred_value=0,Neg_pred_value=0,
                   Precision=0,Recall=0,F1=0,Prevalance=0,
                   Detection_rate=0,Detection_prevalance=0,Balanced_accuracy=0,AUROC=0,AUPRC=0)
for(outcome in strsplit(opt$outcome,split=",")[[1]])
{
  for(fset in strsplit(opt$features,split=",")[[1]])
  {
    for(methoddb in strsplit(opt$method,split=",")[[1]])
    {
      for(size in strsplit(opt$size,split=",")[[1]])
      {
	  	  for(metric in strsplit(opt$metric,split=",")[[1]])
	      {
	  	    if(opt$model=="gbm")
	  	    {
	  	      for(shrinkage in strsplit(opt$shrinkage,split=",")[[1]])
            {
        	    for(intdepth in strsplit(opt$intdepth,split=",")[[1]])
              {
                for(ntrees in strsplit(opt$ntrees,split=",")[[1]])
                {
                  tmp.table<-read.table(paste0(opt$input,"/","gbm_",metric,".",outcome,".",fset,".",methoddb,".sz",size,".shr",shrinkage,".int",intdepth,".ntr",ntrees,".txt"),header=T,stringsAsFactors = F,sep="\t")
                  tmp.table<-tmp.table[,!names(tmp.table) %in% c("cv","repeats")]
                  tmp.result<-result[1,]
			            tmp.result[,c("Sensitivity","Specificity","Pos_pred_value","Neg_pred_value","Precision","Recall","F1","Prevalance","Detection_rate","Detection_prevalance","Balanced_accuracy","AUROC","AUPRC")]<-apply(tmp.table[,c("Sensitivity","Specificity","Pos_pred_value","Neg_pred_value","Precision","Recall","F1","Prevalance","Detection_rate","Detection_prevalance","Balanced_accuracy","AUROC","AUPRC")],2,function(x){
                    median(x,na.rm = T)
                  })
                  tmp.result[1,c("method_training","predict_class","feature_number","method_dataimbalance","undersampleSize","metric","shrinkage","inter_depth","ntrees","alpha","lambda","rfntrees","mtry")]<-tmp.table[1,c("method_training","predict_class","feature_number","method_dataimbalance","undersampleSize","metric","shrinkage","inter_depth","ntrees","alpha","lambda","rfntrees","mtry")]
                  result<-rbind(result,tmp.result)
                }
		          }
        	  }
	  	    } else if(opt$model=="glmnet"){
	  	      for(alpha in strsplit(opt$alpha,split=",")[[1]])
	  	      {
	  	        for(lambda in strsplit(opt$lambda,split=",")[[1]])
	  	        {
	  	          tmp.table<-read.table(paste0(opt$input,"/","glmnet_",metric,".",outcome,".",fset,".",methoddb,".sz",size,".alp",alpha,".lam",lambda,".txt"),header=T,stringsAsFactors = F,sep="\t")
	  	          tmp.table<-tmp.table[,!names(tmp.table) %in% c("cv","repeats")]
	  	          tmp.result<-result[1,]
	  	          tmp.result[,c("Sensitivity","Specificity","Pos_pred_value","Neg_pred_value","Precision","Recall","F1","Prevalance","Detection_rate","Detection_prevalance","Balanced_accuracy","AUROC","AUPRC")]<-apply(tmp.table[,c("Sensitivity","Specificity","Pos_pred_value","Neg_pred_value","Precision","Recall","F1","Prevalance","Detection_rate","Detection_prevalance","Balanced_accuracy","AUROC","AUPRC")],2,function(x){
	  	                median(x,na.rm = T)
	  	          })
	  	          tmp.result[1,c("method_training","predict_class","feature_number","method_dataimbalance","undersampleSize","metric","shrinkage","inter_depth","ntrees","alpha","lambda","rfntrees","mtry")]<-tmp.table[1,c("method_training","predict_class","feature_number","method_dataimbalance","undersampleSize","metric","shrinkage","inter_depth","ntrees","alpha","lambda","rfntrees","mtry")]
	  	          result<-rbind(result,tmp.result)
	  	        }
	  	      } 
	  	    } else if(opt$model=="rf"){
	  	      for(rfntrees in strsplit(opt$rfntrees,split=",")[[1]])
	  	      {
	  	        for(mtry in strsplit(opt$mtry,split=",")[[1]])
	  	        {
	  	          tmp.table<-read.table(paste0(opt$input,"/","rf_",metric,".",outcome,".",fset,".",methoddb,".sz",size,".ntr",rfntrees,".mtr",mtry,".txt"),header=T,stringsAsFactors = F,sep="\t")
	  	          tmp.table<-tmp.table[,!names(tmp.table) %in% c("cv","repeats")]
	  	          tmp.result<-result[1,]
	  	          tmp.result[,c("Sensitivity","Specificity","Pos_pred_value","Neg_pred_value","Precision","Recall","F1","Prevalance","Detection_rate","Detection_prevalance","Balanced_accuracy","AUROC","AUPRC")]<-apply(tmp.table[,c("Sensitivity","Specificity","Pos_pred_value","Neg_pred_value","Precision","Recall","F1","Prevalance","Detection_rate","Detection_prevalance","Balanced_accuracy","AUROC","AUPRC")],2,function(x){
	  	              median(x,na.rm = T)
	  	          })
	  	          tmp.result[1,c("method_training","predict_class","feature_number","method_dataimbalance","undersampleSize","metric","shrinkage","inter_depth","ntrees","alpha","lambda","rfntrees","mtry")]<-tmp.table[1,c("method_training","predict_class","feature_number","method_dataimbalance","undersampleSize","metric","shrinkage","inter_depth","ntrees","alpha","lambda","rfntrees","mtry")]
	  	          result<-rbind(result,tmp.result)
	  	        }
	  	      } 
	  	    }
	  	  }
      }
    }
  }
}
write.table(result[-1,],opt$output,col.names=T,row.names=F,quote=F,sep="\t")