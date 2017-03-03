UVM_RNASeq=read.table("UVM.uncv2.mRNAseq_raw_counts.txt",header = T,sep="\t")

UVM_Cancer=UVM_RNASeq[,grep("*\\.01$", colnames(UVM_RNASeq))]

UVM_Cancer=cbind(UVM_RNASeq[,1],UVM_Cancer)
colnames(UVM_Cancer)[1]= "Gene_symbol"
dim(UVM_Cancer)

UVM_Normal=UVM_RNASeq[,grep("*\\.11$", colnames(UVM_RNASeq))]

############# Testing for other normal types
UVM_Normal_test=UVM_RNASeq[,grep("*\\.10$", colnames(UVM_RNASeq))]
dim(UVM_Normal_test)
###########################################

UVM_Normal=cbind(UVM_RNASeq[,1],UVM_Normal)
colnames(UVM_Normal)[1]= "Gene_symbol"
dim(UVM_Normal)
match=c()

ss=c()
for ( i in colnames(UVM_Cancer))
{
  a=paste(unlist(strsplit(i,"[.]"))[c(1,2,3)],sep=".",collapse = ".") 
  for ( j in colnames(UVM_Normal))
  {
    ss=c(ss,paste(unlist(strsplit(j,"[.]"))[c(1,2,3)],sep=".",collapse = "."))
    if (length(intersect(a,ss)== 1))
    {
      match = c(match,i)
      break
    }
    ss=c()
  }
}
UVM_Cancer_Match=UVM_Cancer[,match]
dim(UVM_Cancer_Match)
write.table(UVM_Cancer_Match,file="./Tumor_Normal_Match/UVM_Cancer_Matched.txt",sep="\t",row.names=F)
write.table(UVM_Normal,file="./Tumor_Normal_Match/UVM_Normal_Matched.txt",sep="\t",row.names=F)

################# clear workspace in R 
rm(list=ls())

