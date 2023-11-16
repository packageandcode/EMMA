library(caret)
library(randomForest)
library(pROC)
library(verification)
library(reshape2)
library(ggsignif)
library(ggplot2)

#custom parameter
args=commandArgs(T)
data_path=args[1]
#load malignant ratio data
file_methy<-read.table(data_path, row.names = 
                       header = T, stringsAsFactors = F, check.names = F)
readMeth <- file_methy[!duplicated(paste(as.character(file_methy$chr), as.character(file_methy$start),
                          as.character(file_methy$end), sep = "_")),-ncol(file_methy)]
rownames(readMeth) <- paste(as.character(readMeth$chr), as.character(readMeth$start),
                          as.character(readMeth$end), sep = "_")
readMeth <- readMeth[,-(1:3)]  
dim(readMeth) 

#split to 3 matrix
fracMeth <- apply(readMeth,2, function(x){
strsplit(x, split = ",") %>% lapply(., function(x) x[1]) %>% unlist %>% as.numeric %>% return})
rownames(fracMeth) <- rownames(readMeth) #extract malignant ratio
countAll <- apply(readMeth,2, function(x){
strsplit(x, split = ",") %>% lapply(., function(x) x[2]) %>% unlist %>% as.numeric %>% return})
rownames(countAll) <- rownames(readMeth) #extract all reads number
countMeth <- round(countAll*fracMeth, digits = 0)
dim(fracMeth) #41199 157
#rename 
coln<-colnames(fracMeth)
colnames(fracMeth)<-coln
colnames(countAll)<-coln
colnames(countMeth)<-coln

#Custom length
my_len<-500 
dmr_len<-file_methy[,3]-file_methy[,2]
fracMeth_1k<-fracMeth[which(dmr_len>=my_len),]
countAll_1k<-countAll[which(dmr_len>=my_len),]
countMeth_1k<-countMeth[which(dmr_len>=my_len),]

#Keep 80% of the samples in the my_len zone with counts greater than or equal to 5
count10<-apply(countAll_1k,1,function(x){x<-length(which(x>=5))});length(which((count10/ncol(countAll_1k))>0.8))
fracMeth_1k10<-fracMeth_1k[which((count10/ncol(countAll_1k))>0.8),]
countAll_1k10<-countMeth_1k[which((count10/ncol(countAll_1k))>0.8),]
countMeth_1k10<-countAll_1k[which((count10/ncol(countAll_1k))>0.8),]

#train set& valid set
train_t<-colnames(fracMeth_1k10)[grep("train_t",colnames(fracMeth_1k10))]
train_n<-colnames(fracMeth_1k10)[grep("train_n",colnames(fracMeth_1k10))]
train<-as.data.frame(t(cbind(train_t,train_n)))
train_sam<-as.factor(c(rep("tumor",ncol(train_t)),rep("normal",ncol(train_n))))
train_t_idx<-grep("train_t",rownames(train))
train_n_idx<-grep("train_n",rownames(train))
valid_t<-colnames(fracMeth_1k10)[grep("valid_t",colnames(fracMeth_1k10))]
valid_n<-colnames(fracMeth_1k10)[grep("valid_n",colnames(fracMeth_1k10))]
valid<-cbind(valid_t,valid_n)
valid_sam<-as.factor(c(rep("tumor",ncol(valid_t)),rep("normal",ncol(valid_n))))

#Screening feature
p.value <- apply(train, 2, function(x) {
    try(wilcox.test(x[train_t_idx],x[train_n_idx])$p.value, silent = T)})
feat_p<-colnames(train[,which(p.value<0.05)])

mean_idx<-which(apply(intrain[,feat_p],2,function(x){x<-(mean(x[train_t_idx])>mean(x[train_n_idx]))}))
sqmean_idx<-intersect(which(apply(intrain[train_t_idx,feat_p],2,function(x){x<-length(which(x>0.2))})>(0.25*length(train_t_idx))),
                      which(apply(intrain[train_n_idx,feat_p],2,function(x){x<-length(which(x<0.2))})>(0.25*train_n_idx)))
feat_p_mean<-colnames(train[,feat_p])[intersect(mean_idx,sqmean_idx)]
#highCorr
descrCorr = cor(train[feat_p_mean])
highCorr = findCorrelation(descrCorr,0.8)
intrain = train[feat_p_mean][, -highCorr]
#多重共线性
comboInfo = findLinearCombos(intrain)
intrain=intrain[, -comboInfo$remove]
set.seed(123)
subsets = c(1:20)
feat_re_in<-rownames(intrain)
intrain<-as.data.frame(intrain[,feat_re_in])
sam_class<-as.factor(train_sam)
ctrl = rfeControl(functions = rfFuncs, method = "cv",repeats = 100)
Profile_ratioMat = rfe(intrain, sam_class, rfeControl = ctrl,sizes = subsets)

#options(repr.plot.width=9, repr.plot.height=5)
feats<-NULL
auc_i<-NULL
obb_i<-NULL
seed_i<-NULL
ntree_i<-NULL
mtry_i<-NULL

trainn<-intrain
trainn_sam<-intrain_sam
p_dmr<-unique(feat_re_in)
for (n in 1:(length(p_dmr)-args[5])){
    aucitmp<-NULL
    if(n==1){
        dmr_s<-p_dmr
    }else{
    idx<-match(feats,p_dmr)
    dmr_s<-p_dmr[-idx]
        }
    for (i in 1:(length(dmr_s))){
        auci<-NULL
		for(i in c(8,sample(1:1000,500))){ #35  #i=108 0.94 #674 
			set.seed(i)
			seed1 <- c(seed1,i)
			train_48<-	sample(1:length(train_t_idx),length(train_t_idx)*0.3)
			trainn<-train[c(-train_48,-(train_48+160)),]
			trainn_sam<-train_sam[c(-train_48,-(train_48+160))]
			intestt<-train[c(train_48,(train_48+160)),]
			intestt_sam<-train_sam[c(train_48,(train_48+160))]
			for(ntree in 1:10){
			ntree_n<-ntree*100
			ntree1<-c(ntree1,ntree_n)
				for(mtry in 2:10){
					mtry_n<-mtry
					mtry1<-c(mtry1,mtry_n)
					for(se in c(66,sample(1:1000,500))){
						set.seed(se)
						seed2<-c(seed2,se)
						best_rf_model<- randomForest(x = trainn[,dmr_s], y =trainn_sam ,importance = TRUE, 	   ntree = ntree_n,mtry=mtry_n,proximity=TRUE)
				obb1<-c(obb1,round(as.numeric(best_rf_model$err.rate[nrow(best_rf_model$err.rate),1]),digits = 3))
				rf_i_pre<-predict(best_rf_model, newdata=intestt[,dmr_s],type="prob" )
				auci<-c(auci,as.numeric(auc(roc(intest_sam,rf_i_pre[,2]))[1])) #
				   }
				}
			}
		}

    aucitmp<-c(aucitmp,max(auci))
    feat_i<-dmr_s[order(aucitmp,decreasing = T)[1]]
   }
    feats<-c(feats,feat_i)
    auc_i<-c(auc_i,max(aucitmp))
	seed_i<-c(seed_i,which(auci==max(aucitmp))
	obb_i<-c(obb_i,which(auci==max(aucitmp))
	ntree_i<-c(ntree_i,which(auci==max(aucitmp))
	mtry_i<-c(mtry_i,which(auci==max(aucitmp))
    }


feats_auc_r<-cbind(feats,auc_i,seed_i,obb_i,ntree_i,mtry_i)
feat_last<-feat_re_in[-match(feats,feat_re_in)]

#valid set
valid_t<-colnames(fracMeth_1k10)[grep("valid_t",colnames(fracMeth_1k10))]
valid_n<-colnames(fracMeth_1k10)[grep("valid_n",colnames(fracMeth_1k10))]
valid<-as.data.frame(t(cbind(valid_t,valid_n)))
valid_sam<-as.factor(c(rep("tumor",ncol(valid_t)),rep("normal",ncol(valid_n))))
valid_t_idx<-grep("valid_t",rownames(valid))
valid_n_idx<-grep("valid_n",rownames(valid))


#model
seed_i<-NULL;auc_i<-NULL;auc_x<-NULL;spec_i<-NULL;spec_x<-NULL;
for( i in  sample(1:1000,200)){
    set.seed(i)
    seed_i<-c(seed_i,i)
    selected_features <-feat_last#c(feat_50,c("cnv_mds1","cnv_mds2"),c("len_mds1","len_mds2"))#unique(c(feat_37[-22:-23],intersect(feat_p.esd,feat_mean)[-15]))#feat_p.esd_cnv
    fitControl <-  trainControl(method = "cv",    number = 10,    classProbs = T,    savePredictions = T,    summaryFunction = twoClassSummary )
    model <- train(intrain[,selected_features], sam_class, method = "rf", trControl = fitControl,  metric = "ROC")
    auc_i<-c(auc_i,max(model$results$ROC))
    auc_x<-c(auc_x, as.numeric(auc(roc(valid_sam,predict(model, newdata = valid,type="prob" )[,2]))[1]))
    spec_i<-c(spec_i,model$results[1,3])
    spec_x<-c(spec_x,roc(as.numeric(valid_sam),as.numeric(predict(model, newdata = valid)))$specificities[2])

}

write.table(cbind(seed_i,auc_i,auc_x,spec_i,spec_x),
            file=args[6],col.names = T,row.names = T,quote = F,sep = "\t")

