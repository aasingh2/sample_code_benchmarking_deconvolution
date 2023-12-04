## Author: Aakanksha Singh
## Date: 2 August 2022
## Version: 1.0
## Aim: To generate pseudobulk samples with known numbers of cell types from Single
##      SingleCellExperiment Object

# For the function generating_matrix_for_deconv
#   Input: (a) A processed SingleCellExperiment dataset (output from processing_sce())
#          (b) The title of the datset (string)
#          (c) The name and number of each celltype (atm 4 types of cells) 
#   Output: (a) A matrix summarising the SingleCellExperiment that can be used as a bulk sample for further analysis
#           (b) Plot of mean log expression vs variance of log expression with a treand line
#           (c) A tSNE plot colored by predicted celltype
#           (d) A table with the number of the predicted cell types
#           (e) The output variable stores the analysed SingleCellExperiment Object with cell type annotations. 



generating_matrix_for_deconv<-function(sce,title,cell1,num1,cell2,num2,cell3,num3,cell4,num4)
{ 
  # Adding an additional column with consolidated cell types, 
  # removing rows with no Hugo Symbol,
  # Setting the Hugo symbol as rownames instead of ensmbl ids
  #sce$cell<-sce$predicted_celltype
  #colData(sce)$cell <- ifelse(colData(sce)$cell %in% c("CD8+ T cells", "CD4+ T cells"), "T cells", colData(sce)$cell)
  #rownames(sce) <- rowData(sce)$Symbol
  #sce<-sce[!(is.na(row.names(sce))),]
  #table(sce$cell)

  # Getting cell types of each type
  cell1<-sce[,sce$cell==cell1]
  cell2<-sce[,sce$cell==cell2]
  cell3<-sce[,sce$cell==cell3]
  cell4<-sce[,sce$cell==cell4]
  
  # Generating a random sample with the number of known cell types
  known_sample<-cbind(cell1[,sample(ncol(cell1),num1)] ,cell2[,sample(ncol(cell2),num2)],
                      cell3[,sample(ncol(cell3),num3)],cell4[,sample(ncol(cell4),num4)])
  
  # Summing the counts of all the cell rowwise to get a pseudobulk sample
  samp<-as.data.frame(Matrix::rowSums(assays(known_sample)$counts))
  samp<-cbind(rownames(known_sample),samp)
  
  # Aggregating rows with common gene symbol
  samp<-aggregate(samp[,-1],by=samp[1],sum)
  #colnames(samp)=c("HUGO",title)
  
  # Getting TPM ( where TPM= (count/total count)* 1000000)
  count_sum<-sum(samp$x)
  samp[,2]<-sapply(samp[,2],function(i) (i/count_sum)*1000000)
  colnames(samp)=c("HUGO",(paste0(title,"_TPM")))
  return(samp)
  
}


generating_raw_pseudobulk_matrix<-function(sce,title,cell1,num1,cell2,num2,cell3,num3,cell4,num4)
{ 
  # Adding an additional column with consolidated cell types, 
  # removing rows with no Hugo Symbol,
  # Setting the Hugo symbol as rownames instead of ensmbl ids
  #sce$cell<-sce$predicted_celltype
  #colData(sce)$cell <- ifelse(colData(sce)$cell %in% c("CD8+ T cells", "CD4+ T cells"), "T cells", colData(sce)$cell)
  #rownames(sce) <- rowData(sce)$Symbol
  #sce<-sce[!(is.na(row.names(sce))),]
  #table(sce$cell)
  
  # Getting cell types of each type
  cell1<-sce[,sce$cell==cell1]
  cell2<-sce[,sce$cell==cell2]
  cell3<-sce[,sce$cell==cell3]
  cell4<-sce[,sce$cell==cell4]
  
  # Generating a random sample with the number of known cell types
  known_sample<-cbind(cell1[,sample(ncol(cell1),num1)] ,cell2[,sample(ncol(cell2),num2)],
                      cell3[,sample(ncol(cell3),num3)],cell4[,sample(ncol(cell4),num4)])
  
  # Summing the counts of all the cell rowwise to get a pseudobulk sample
  samp<-as.data.frame(Matrix::rowSums(assays(known_sample)$counts))
  samp<-cbind(rownames(known_sample),samp)
  
  # Aggregating rows with common gene symbol
  samp<-aggregate(samp[,-1],by=samp[1],sum)
  colnames(samp)=c("HUGO",title)
  
  # Getting TPM ( where TPM= (count/total count)* 1000000)
  # count_sum<-sum(samp$x)
  # samp[,2]<-sapply(samp[,2],function(i) (i/count_sum)*1000000)
  # colnames(samp)=c("HUGO",(paste0(title,"_TPM")))
  return(samp)
  
}
