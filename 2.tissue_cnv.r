correctGC <- function(data) {
  filtered_data <- data[data$CG > 0.3, ]
  model <- loess(counts ~ CG, data = filtered_data)
  filtered_data$Fit <- predict(model)
  Median1 <- median(filtered_data$counts)
  filtered_data$WeightValue <- Median1 / filtered_data$Fit

  data$loessCounts <- data$counts
  data[data$CG > 0.3, ]$loessCounts <- filtered_data$counts * filtered_data$WeightValue

  return(data)
}


args=commandArgs(T)
normal<-read.table(args[1],header=F) #file_counts[which(file_counts$sample=="015N"),]
tumor<-read.table(args[2],header=F) #file_counts[which(file_counts$sample=="015T"),]

colnames(normal)<-c("chr","start","end","AT","CG","counts")
colnames(tumor)<-c("chr","start","end","AT","CG","counts")

normal<-correctGC(normal)
tumor<-correctGC(tumor)
#plot(normal$CG,normal$loessCounts,cex=0.1,xlab="GC Content",ylab=expression(paste("counts",sep="")),xlim = c(0.3,0.6))
#plot(normal$CG,normal$counts,cex=0.1,xlab="GC Content",ylab=expression(paste("counts",sep="")),xlim = c(0.3,0.6))

#cnv
log2ratio<-as.data.frame(log2((tumor$loessCounts+0.01)/(normal$loessCounts+0.01)))
colnames(log2ratio)<-args[3]
write.table(log2ratio,file=args[4],sep = "\t",row.names = F,col.names = T,quote = F)

#############GMM+HMM计算cnv
library(pheatmap)
library(reshape2)
library(HMM)
library(ggplot2)
library(mclust)
df<-read.table("/share/pub/wangxy/dlj/hg38.chrom.sizes",header = F)[1:24,]
colnames(df) <- c("chromName","chromlength")
chr<-paste0("chr",c(1:22))
chr_idx<-match(chr,df$chromName)
df$chromName
df<-df[chr_idx,]

df$chromNum <- 1:length(df$chromName)

df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度
# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])
# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle
df

scores<-read.table("/share2/pub/dailj/dailj/escc_cfDNA/CNV/tissue/atcg.depth/9.zscore.chrX.CNV/all.mean.sd.zscore.nochrY.bed",sep="\t",header=T,stringsAsFactors = F)

scores<-read.table("/share2/pub/dailj/dailj/escc_cfDNA/CNV/tissue/atcg.depth/5.CNV/new_try/log2R.all.bed",sep="\t",header=T,stringsAsFactors = F)

scores<-read.table("/share2/pub/dailj/dailj/fdd/cnv/009.nodone.log2R.txt",sep="\t",header=T,stringsAsFactors = F)
pre_ncol<-3    
ti_cnv<-c()
for(i in 1:(ncol(scores)-pre_ncol)){
    predicted_states<-c()
    obs_data <-scores[,c(pre_ncol+i)][which(scores$chr!=("chrY")&scores$chr!=("chrX"))]obs_data[is.na(obs_data)]=-0.1
    set.seed(66);prepare<-kmeans(matrix(obs_data), centers = 3)
    summary(prepare$cluster)
    pre_cluster <- prepare$cluster
    pre_cluster[which(prepare$cluster==order(prepare$centers)[1])]<-(1)
    pre_cluster[which(prepare$cluster==order(prepare$centers)[2])]<-(2)
    obs_data[which(prepare$cluster==order(prepare$centers)[3])]<-obs_data[which(prepare$cluster==order(prepare$centers)[3])]*1.5
    obs_data[which(prepare$cluster==order(prepare$centers)[1])]<-obs_data[which(prepare$cluster==order(prepare$centers)[1])]*0.5
    obs_data[which(prepare$cluster==order(prepare$centers)[2])]<-obs_data[which(prepare$cluster==order(prepare$centers)[2])]*1
    summary(pre_cluster)
    set.seed(66);
    gmm_model<-Mclust(c(obs_data), G=3,modelNames = "V")
    center_d<-order(gmm_model$parameters$mean)[1]
    center_a<-order(gmm_model$parameters$mean)[3]
    labels<-ifelse( (apply(gmm_model$z, 1, which.max) == center_d), "DELE",
                             ifelse((apply(gmm_model$z, 1, which.max) == center_a),"AMPL", "noCNV"))#[-1:-100]
    
    table(labels)
    #HMM
    for (c in 1:22) {
    chr_chr<-labels[which(scores$chr[which(scores$chr!=("chrY")&scores$chr!=("chrX"))]==paste0("chr",c))]#[which(obs_data$label==paste0("chr",c))]#
    # 定义状态和初始概率
    states <- c( "noCNV","DELE", "AMPL")
    startProbs <- c(0.8, 0.14, 0.06)  #table(ti_cnv$X002)
    
    # 定义状态转移概率矩阵
    transProbs <- matrix(c(0.9, 0.09, 0.01,  # 从DELE到DELE的概率为0.7，到noCNV的概率为0.2，到AMPL的概率为0.1
                           0.005, 0.99, 0.005,  # 从noCNV到DELE的概率为0.3，到noCNV的概率为0.6，到AMPL的概率为0.1
                           0.001, 0.009, 0.99), # 从AMPL到DELE的概率为0.2，到noCNV的概率为0.2，到AMPL的概率为0.6
                         nrow = length(states), ncol = length(states), byrow = TRUE,
                         dimnames = list(states, states))
    # 定义状态发射概率矩阵
    emissionProbs <- matrix(c(0.8, 0.1, 0.1,  # DELE状态下观测到DELE的概率为0.2，观测到noCNV的概率为0.6，观测到AMPL的概率为0.2
                              0.19, 0.8, 0.01,  # noCNV状态下观测到DELE的概率为0.5，观测到noCNV的概率为0.2，观测到AMPL的概率为0.3
                              0.19, 0.01, 0.8), # AMPL状态下观测到DELE的概率为0.3，观测到noCNV的概率为0.3，观测到AMPL的概率为0.4
                            nrow = length(states), ncol = 3, byrow = TRUE,
                            dimnames = list(states, c( "noCNV","DELE", "AMPL")))
    # 构建HMM模型
    hmm_model <- initHMM(States = states, Symbols = c("noCNV","DELE",  "AMPL"),    
                         startProbs = startProbs, transProbs= transProbs, emissionProbs = emissionProbs)
    # 运用Viterbi算法进行预测
    predicted_states <- c(predicted_states,viterbi(hmm_model, chr_chr))
    
  }
  ti_cnv<-cbind(ti_cnv,predicted_states)
}
head(ti_cnv)
colnames(ti_cnv)<-colnames(scores)[-1:-pre_ncol]
ti_cnv<-cbind(scores[which(scores$chr!=("chrY")&scores$chr!=("chrX")),c(1:3)],ti_cnv)

write.table(ti_cnv,
            file=args[3],col.names = T,row.names = F,quote = F,sep = "\t")
