#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 4rd 2022 ###
### training individual gbm model ###
library(optparse)
library(caret)
library(gbm)
library(smotefamily)
library(PerfMeas)
library(yardstick)
library(dplyr)
library(gridExtra)

### a string of pseudo seeds to ensure the result can be reproduced ###
set.seed(23)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input matrix of features", metavar="character"),
  make_option(c("-f", "--features"), type="character", default=NULL, 
              help="files containing names of features", metavar="character"),
  make_option(c("--outcome"), type="character", default=NULL, 
              help="outcome name", metavar="character"),
  make_option(c("--model"), type="character", default=NULL, 
              help="gbm, glmnet or rf", metavar="character"),
  make_option(c("--method"), type="character", default=NULL, 
              help="method for treating data imbalance", metavar="character"),
  make_option(c("--size"), type="numeric", default=NULL, 
              help="proportion of major class to be preserved (from 0 to 1)", metavar="character"),
  make_option(c("--metric"), type="character", default=NULL, 
              help="metric used for training (Accuracy or Kappa)", metavar="character"),
  make_option(c("--shrinkage"), type="numeric", default=NULL, 
              help="GBM: shrinkage", metavar="character"),
  make_option(c("--intdepth"), type="integer", default=NULL, 
              help="GBM: interaction depth", metavar="character"),
  make_option(c("--ntrees"), type="integer", default=NULL, 
              help="GBM: ntrees", metavar="character"),
  make_option(c("--alpha"), type="numeric", default=NULL, 
              help="GLMNET: alpha", metavar="character"),
  make_option(c("--lambda"), type="numeric", default=NULL, 
              help="GLMNET: lambda", metavar="character"),
  make_option(c("--mtry"), type="numeric", default=NULL, 
              help="RF: mtry", metavar="character"),
  make_option(c("--rfntrees"), type="numeric", default=NULL, 
              help="RF: ntrees", metavar="character"),
  make_option(c("--cv"), type="integer", default=NULL, 
              help="cross validation", metavar="character"),
  make_option(c("--repeats"), type="integer", default=NULL, 
              help="repeats", metavar="character"),
  make_option(c("--ROC"), type="character", default=NULL, 
              help="ROC pdf", metavar="character"),
  make_option(c("--PRC"), type="character", default=NULL, 
              help="PRC pdf", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#### Load sub-functions ####
auroc<-function(score,cls){
  n1<-sum(!cls); sum(cls)->n2;
  U<-sum(rank(score)[!cls])-n1*(n1+1)/2;
  return(1-U/n1/n2);
}
#### END ####
#opt<-data.frame(input="/g/steinmetz/project/yeast_crispr/MAGESTIC_REDI_final/CRISPR_predict_outcome/input_feature_table.txt",
#                features="CRISPR_predict_outcome/inputs/fset3.txt",model="gbm",outcome="Editing_outcome_SV",
#                method="us",size=0.5,shrinkage=0.01,ntrees=500,intdepth=1,metric="Accuracy",
#                cv=3,repeats=5)
#### Load feature table ####
indata<-read.table(opt$input,sep="\t",header=T,stringsAsFactors = F)
fdata<-read.table(opt$features,header = F,stringsAsFactors = F)
#### Assembly feature table and scale columns ####
index.fs<-fdata[,1]
### trim data table ####
indata<-indata[,c(opt$outcome,index.fs)]
names(indata)[1]<-"class"
### remove invariable features ###
var<-apply(indata,2,function(x){
  if(length(names(table(x)))>1)
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
})
indata<-indata[,var]

### start gridsearch model parameters ###
i<-1
result<-data.frame(method_training="gbm",predict_class="X",feature_number=0,method_dataimbalance="undersample",undersampleSize=0,metric="Accuracy",shrinkage=NA,inter_depth=NA,ntrees=NA,alpha=NA,lambda=NA,mtry=NA,rfntrees=NA,repeats=0,cv_fold=0,
                              Sensitivity=0,Specificity=0,Pos_pred_value=0,Neg_pred_value=0,
                              Precision=0,Recall=0,F1=0,Prevalance=0,
                              Detection_rate=0,Detection_prevalance=0,Balanced_accuracy=0,AUROC=0,AUPRC=0)
#### Data augmentation compare #1 downsampling ####
roc_gbm<-list()
pr_gbm<-list()
#AUROC.sum<-0
#AUPRC.sum<-0
#generate input table
for(repeats in 1:opt$repeats)
{ 
  neg_class_table<-indata[indata$class=="0",]
  pos_class_table<-indata[indata$class=="1",]

  ## require caret package ##
  folds_0<-createFolds(y = neg_class_table[,2],k=opt$cv, list = F)
  folds_1<-createFolds(y = pos_class_table[,2],k=opt$cv, list = F)
  for(f in unique(folds_0))
  {
    ind_0 <- which(folds_0 == f) 
    ind_1 <- which(folds_1 == f)
    data_train <- rbind(neg_class_table[-as.numeric(ind_0),],
                        pos_class_table[-as.numeric(ind_1),])
    data_valid <- rbind(neg_class_table[as.numeric(ind_0),],
                        pos_class_table[as.numeric(ind_1),])
    n_neg_class<-dim(data_train[data_train$class==0,])[1]
    n_pos_class<-dim(data_train[data_train$class==1,])[1]
    if(opt$method=="us")
    {
      data_train <- rbind(neg_class_table[-as.numeric(ind_0),][sample(n_neg_class,round(n_neg_class*opt$size)),],
                          pos_class_table[-as.numeric(ind_1),])
    } else if(opt$method=="usos") {
      if(round(n_neg_class*opt$size) > n_pos_class)
      {
        data_train <- rbind(neg_class_table[-as.numeric(ind_0),][sample(n_neg_class,round(n_neg_class*opt$size)),],
                            pos_class_table[-as.numeric(ind_1),],
                            pos_class_table[-as.numeric(ind_1),][sample(n_pos_class,round(n_neg_class*opt$size-n_pos_class),replace=T),])
      } else {
        data_train <- rbind(neg_class_table[-as.numeric(ind_0),][sample(n_neg_class,round(n_neg_class*opt$size)),],
                           pos_class_table[-as.numeric(ind_1),])
      }
    } else if(opt$method=="ussmote") {
      if(round(n_neg_class*opt$size) > n_pos_class)
      {
        data_train_gen_data<-SMOTE(X = data_train[,!names(data_train) %in% c("class")],target = data_train[,"class"], K = 3)
        data_train<-rbind(data_train_gen_data$orig_N[sample(dim(data_train_gen_data$orig_N)[1],round(dim(data_train_gen_data$orig_N)[1]*opt$size)),],
                          data_train_gen_data$orig_P,
                          data_train_gen_data$syn_data[sample(dim(data_train_gen_data$syn_data)[1],min(round(dim(data_train_gen_data$orig_N)[1]*opt$size)-dim(data_train_gen_data$orig_P)[1],dim(data_train_gen_data$syn_data)[1])),])
      } else {
        data_train <- rbind(neg_class_table[-as.numeric(ind_0),][sample(n_neg_class,round(n_neg_class*opt$size)),],
                            pos_class_table[-as.numeric(ind_1),])
        
      }
    } else if(opt$method=="usadasyn") {
      if(round(n_neg_class*opt$size) > n_pos_class)
      {
        data_train_gen_data<-ADAS(X = data_train[,!names(data_train) %in% c("class")],target = data_train[,"class"],K = 3)
        data_train<-rbind(data_train_gen_data$orig_N[sample(dim(data_train_gen_data$orig_N)[1],round(dim(data_train_gen_data$orig_N)[1]*opt$size)),],
                                                          data_train_gen_data$orig_P,
                                                          data_train_gen_data$syn_data[sample(dim(data_train_gen_data$syn_data)[1],min(round(dim(data_train_gen_data$orig_N)[1]*opt$size)-dim(data_train_gen_data$orig_P)[1],dim(data_train_gen_data$syn_data)[1])),])
      } else {
        data_train <- rbind(neg_class_table[-as.numeric(ind_0),][sample(n_neg_class,round(n_neg_class*opt$size)),],
                            pos_class_table[-as.numeric(ind_1),])
      }
    } else {
      stop("unknown argument for --method")
    }
    if(opt$model=="gbm")
    {
      modelGrid <- expand.grid(
        n.trees = opt$ntrees,
        shrinkage= opt$shrinkage,n.minobsinnode = 10,interaction.depth=opt$intdepth)
      ml_model<-train(factor(class) ~ ., data = data_train,method = opt$model, metric = opt$metric, tuneGrid = modelGrid,trControl =  trainControl(method = "none"))
    } else if(opt$model=="glmnet") {
      ### standarization of data prior to glmnet ###
      train_mean<-apply(data_train[,!names(data_train) %in% "class"],2,mean)
      train_sd<-apply(data_train[,!names(data_train) %in% "class"],2,sd)
      data_train[,!names(data_train) %in% "class"]<-apply(as.matrix(names(data_train[,!names(data_train) %in% "class"])),1,function(x){
        (data_train[,x]-train_mean[x])/train_sd[x]
      })
      data_valid[,!names(data_valid) %in% "class"]<-apply(as.matrix(names(data_valid[,!names(data_valid) %in% "class"])),1,function(x){
        (data_valid[,x]-train_mean[x])/train_sd[x]
      })
      ### remove invariable sites
      data_train<-data_train[,!names(data_train) %in% names(train_sd)[train_sd==0]]
      data_valid<-data_valid[,!names(data_valid) %in% names(train_sd)[train_sd==0]]
      modelGrid <- expand.grid(
        alpha = opt$alpha,
        lambda= opt$lambda)
      ml_model<-train(factor(class) ~ ., data = data_train,method = opt$model, metric = opt$metric, tuneGrid = modelGrid,trControl =  trainControl(method = "none"))
    } else if(opt$model=="rf") {
      modelGrid <- expand.grid(
        .mtry = opt$mtry)
      ml_model<-train(factor(class) ~ ., data = data_train,method = opt$model, metric = opt$metric, tuneGrid = modelGrid,trControl =  trainControl(method = "none"),ntrees = opt$rfntrees)
    } else {
      stop("unknown argument for --model")
    }
    
    #ml_model<-gbm(class ~ ., data = data_train,n.minobsinnode = 10, shrinkage = opt$shrinkage, interaction.depth = opt$intdepth, n.trees = opt$ntrees,distribution = "bernoulli")
    valid_prediction<-data.frame(label=data_valid$class,prediction=predict(ml_model,data_valid,type = "raw"))
    ## require caret package ##
    conf_mat<-confusionMatrix(factor(valid_prediction$prediction),factor(valid_prediction$label),positive = "1")
    gridsearch_result_new<-result[1,]
    gridsearch_result_new$method_training<-opt$model
    gridsearch_result_new$predict_class<-opt$outcome
    gridsearch_result_new$feature_number<-length(index.fs)
    gridsearch_result_new$method_dataimbalance<-opt$method
    gridsearch_result_new$metric<-opt$metric
    gridsearch_result_new$undersampleSize<-opt$size
    if(opt$model=="gbm")
    {
      gridsearch_result_new$shrinkage<-opt$shrinkage
      gridsearch_result_new$inter_depth<-opt$intdepth
      gridsearch_result_new$ntrees<-opt$ntrees
    } else if(opt$model=="glmnet") {
      gridsearch_result_new$alpha<-opt$alpha
      gridsearch_result_new$lambda<-opt$lambda
    } else if(opt$model=="rf") {
      gridsearch_result_new$rfntrees<-opt$rfntrees
      gridsearch_result_new$mtry<-opt$mtry
    } else {
      stop("unknown argument for --model")
    }
    gridsearch_result_new$repeats<-repeats
    gridsearch_result_new$cv_fold<-f
    gridsearch_result_new[,16:26]<-conf_mat$byClass

    ### calculate PR curve and AUC ###
    pred_valid_tmp<-predict(ml_model,data_valid,type="prob")[,2]
#    gbm.valid.prediction <- data.frame(prediction = as.matrix(pred_valid_tmp), observed = factor(data_valid$class,levels=c(1,0)))
    gridsearch_result_new$AUROC<-auroc(pred_valid_tmp,ifelse(data_valid$class=="1",TRUE,FALSE))
    gridsearch_result_new$AUPRC<-as.numeric(AUPRC(list(precision.at.all.recall.levels(pred_valid_tmp,data_valid$class))))
#    AUROC.sum<-AUROC.sum + gridsearch_result_new$AUROC
#    AUPRC.sum<-AUPRC.sum + gridsearch_result_new$AUPRC
    
    result<-rbind(result,gridsearch_result_new)
    # computing precision recall values
#    roc_gbm[[paste0(repeats,f)]] <- gbm.valid.prediction %>%
#      roc_curve(observed, prediction)
#    roc_gbm[[paste0(repeats,f)]]$fold<-paste0(repeats,f)
#    pr_gbm[[paste0(repeats,f)]] <- gbm.valid.prediction %>%
#      pr_curve(observed, prediction)
#    pr_gbm[[paste0(repeats,f)]]$fold<-paste0(repeats,f)
                
    i<-i+1
  }
}
#roc_gbm_all<-do.call(rbind,roc_gbm)
#pr_gbm_all<-do.call(rbind,pr_gbm)
#p1<-roc_gbm_all %>%
#  ggplot() +
#  geom_path(aes(1 - specificity, sensitivity,color = factor(fold))) + # connect the points in the order in which they appear in the data to form a curve
#  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # add a reference line by convention
#  coord_fixed(ratio = 1, xlim=c(0,1),ylim=c(0,1))  +
#  annotate("text", x=0.8, y=0.1, label= paste0("AUC = ",sprintf("%.3f",AUROC.sum/(opt$repeats*opt$cv))))+
#  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#p2<-pr_gbm_all %>%
#  ggplot() +
#  geom_path(aes(recall, precision,color = factor(fold))) + # connect the points in the order in which they appear in the data to form a curve
#  coord_fixed(ratio = 1, xlim=c(0,1),ylim=c(0,1)) +
#  annotate("text", x=0.8, y=0.9, label= paste0("AUC = ",sprintf("%.3f",AUPRC.sum/(opt$repeats*opt$cv))))+
#  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#pdf(opt$ROC,height=6,width=6.6,useDingbats = F)
#grid.arrange(p1)
#dev.off()
            
#pdf(opt$PRC,height=6,width=6.6,useDingbats = F)
#grid.arrange(p2)
#dev.off()
write.table(result[-1,],opt$output,row.names = F,col.names = T,quote = F,sep="\t")

#### END ####
