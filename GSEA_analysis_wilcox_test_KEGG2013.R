#############################################
#GSEA Analysis on both gene expression and mutation survival output files . Ranking the genes by summing up each analysis rank and then ranking again.
##############################################
expression_survival = read.csv("/media/pathway/48F918113AE39451/New_results/Genes_with_only_Pvalues_new_coxph_expresion.txt",header=T, sep="\t",row.names = 1)#read in dataframe from expression survival analysis
mutation_survival = read.csv("/media/pathway/48F918113AE39451/New_results/Mutation_survival_onlyPval.csv",row.names = 1) #Read in dataframe from mutation survival analysis results
expression_survival = expression_survival[-grep("^\\?",rownames(expression_survival)),] #Removing predicted LOCUS genes starting with '?' for HNSC gene name 
###################################### Loading pathway knowledge
load("/home/pathway/Desktop/labRotation1_AnubhavK/pathwaydata/fwdpwanalyses/pwgenes.RData") #Frank-given data for pathway
pathway_details=read.csv(file="/home/pathway/Desktop/labRotation1_AnubhavK/pathwaydata/fwdpwanalyses/pathways.csv") #Frank-given data about pathway IDs and gene names associated
pw_names = as.character(pathway_details$pw)
pw_description = as.character(pathway_details$desc)
rm(pathway_details)
###################################
featured_genes = rownames(mutation_survival)
for (cancer in colnames(mutation_survival))
{
  ranks_expression = rank(expression_survival[,cancer])
  ranks_mutation = rank(mutation_survival[,cancer])
  sum_ranks_both = apply(cbind(ranks_expression, ranks_mutation),1,sum)
  ranks_acc = rank(sum_ranks_both)
  pw_genesNrs = list()
  pw_genes = list()
  for ( pw in pw_names ) {
  pathway_genes = pwgenes[[pw]] #genes in a particular pathway
  matched_pathway_genes = pathway_genes %in% featured_genes # Matching genes in a particular pathway  
  pw_genes[[pw]] = unique(pathway_genes[matched_pathway_genes])#Filtering genes that dont match  
  pw_genesNrs[[pw]] = length(pathway_genes[matched_pathway_genes]) # maintaining a record of no of genes left
  }

  W_pValue = list()
  for ( pw in pw_names ) {
  reduced_genes = pw_genes[[pw]]
  
  pw_ranks = ranks_acc[match(reduced_genes, featured_genes)] 
  rest_of_genes_ranks = ranks_acc[!ranks_acc %in% pw_ranks]
  if(!(length(pw_ranks) == 0)) {
    wilcox_test = wilcox.test(pw_ranks,rest_of_genes_ranks,
                              alternative="less")
    W_pValue[pw] = wilcox_test$p.value
  } else { W_pValue[1] = 1 }
  
  q_val = p.adjust(W_pValue, method = "fdr")
  }
  write.csv(q_val,file = paste(cancer,"_GSEA_Analysis_Mutation_Expression_survival.csv",sep=""), row.names = T)
}
