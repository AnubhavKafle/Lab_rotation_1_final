########################## Differential expression analysis######################################
Tumor_Normal_Match=c("BRCA","COAD","HNSC","KIRC","LIHC","LUAD","LUSC","PRAD","THCA") ## Only for samples that have tumor normal match with 40 or more samples.

Cancer_Diff=list()
Normal_Diff=list()

for ( i in Tumor_Normal_Match)
{
Cancer_Diff[[i]]=read.table(list.files("~/Desktop/stddata__2016_01_28/Tumor_Normal_Match/",pattern=paste(i,"_Cancer",sep=""),full.names = T), header=T,sep="\t")

Normal_Diff[[i]]=read.table(list.files("~/Desktop/stddata__2016_01_28/Tumor_Normal_Match/",pattern=paste(i,"_Normal",sep=""),full.names = T),sep="\t",header=T)

}

######################################################################################
#    Function to calculate differential expression for all the cancer types          #
######################################################################################
#######################################################################
library(DESeq2)  ########Call Deseq2 package for differential analysis#
#######################################################################
Calculation_DiffExpr = function(cancer_case,normal_case) 
{
  sample=c(rep("Cancer",length(colnames(cancer_case[,-1]))),rep("Normal",length(colnames(normal_case[,-1]))))
  p=c(colnames(cancer_case[,-1]),colnames(normal_case[,-1]))
  sample_info=data.frame(sample,row.names=p)
  ##################### 
  #Preparing data for DESEQ2 analysis
  ####################
  cancer_matrix=cancer_case[,-1]
  rownames(cancer_matrix)=cancer_case[,1]
  normal_matrix=normal_case[,-1]
  rownames(normal_matrix)=normal_case[,1]
  Cancer_RNASeq_count=cbind(cancer_matrix,normal_matrix)
  RNASeq_count_Rounding=apply(Cancer_RNASeq_count,MARGIN = 2,FUN = floor)
  ##########################################################
  deseq.matrix=DESeqDataSetFromMatrix(RNASeq_count_Rounding,colData = sample_info,design = formula(~ sample)) #Check with Frank about the accuracy of the linear model
  dds = DESeq(deseq.matrix)
  size.factors=sizeFactors(dds)
  #size.factors
  normalized.counts = t(t(RNASeq_count_Rounding)/size.factors)
  head(RNASeq_count_Rounding)
  head(RNASeq_count_Rounding)
  results = results(dds, contrast=c("sample", "Cancer", "Normal"))
  #pval=results$padj
  #log2fc=results$log2FoldChange
  #fc=2^log2fc
  return(results)
}
############### Differential analysis ends here####################
DeSEq2=list()
for ( number in 1:length(names(Cancer_Diff)))
{
  
  DeSEq2[[names(Cancer_Diff)[number]]]=Calculation_DiffExpr(Cancer_Diff[[number]],Normal_Diff[[number]])
  write.table(DeSEq2[[number]],file=paste("DeSEQ2_analysis/",names(Cancer_Diff[number]),sep=""),sep="\t",row.names= T)
}
  




