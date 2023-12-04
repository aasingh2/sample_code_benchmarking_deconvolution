## Author: Aakanksha Singh
## Date: 2 August 2022
## Version: 1.0
## Aim: To dconvolute with SpatialDeconv



# For the function processing_sce
#   Input: (a) A dataframe with the peseudobulk samples (GeneSymbols as rows and colnames as sample)
#          (b) A Processed SingleCellExperiment dataset
#   Output: (a) A dataframe with proportions of all the cell types


spatial_deconv <- function(spatial_data, single_cell_data, matrix_type = c("safeTME", "custom"))

{
  # Data wrangling for Spatial Deconv
  
     ## Coverting col to rownames and converting dataframe to numerical matrix
        rownames(spatial_data)<-spatial_data$HUGO
        spatial_data<-data.matrix(spatial_data[-1])
    
    
      ## Convert the deconv_table to a sce which is then converted to seurat with a spatial component
         sc_deconv_table <- SpatialExperiment(
         assays = list(counts = spatial_data, Spatial=spatial_data),
         rowData = rownames(spatial_data),
         colData = colnames(spatial_data))
    
         logcounts(sc_deconv_table)<-log2(counts(sc_deconv_table))
         sc_deconv_table_seurat<- as.Seurat(sc_deconv_table)
    
      ## Adding the Spatial Assay to the seurat object to be able to run spatial decon
         sc_deconv_table_seurat@assays[["Spatial"]]<-CreateAssayObject(data = spatial_data)
    
  # Running spatialdeconv according to the type of signature matrix
      
         if (matrix_type=="custom")
          {
              # Custom profile matrix 
              single_cell_data$cell<-single_cell_data$predicted_celltype
              colData(single_cell_data)$cell <- ifelse(colData(single_cell_data)$cell %in% c("CD8+ T cells", "CD4+ T cells"), "T cells", colData(single_cell_data)$cell)
              rownames(single_cell_data) <- rowData(single_cell_data)$Symbol
              single_cell_data<-single_cell_data[!(is.na(row.names(single_cell_data))),]
      
              annot<-as.data.frame(cbind(single_cell_data$Barcode,single_cell_data$cell))

              x<-as.data.frame(assay(single_cell_data))
              colnames(x)<-c(single_cell_data$Barcode)

              # The function does not accept a S4 object
              custom_matrix<-create_profile_matrix(mtx = x, 
                                            cellAnnots = annot, 
                                            cellTypeCol = "V2", 
                                            cellNameCol = "V1", 
                                            matrixName = "custom_mini_colon",
                                            outDir = NULL, 
                                            normalize = FALSE, 
                                            minCellNum = 5, 
                                            minGenes = 10,
                                            discardCellTypes = FALSE)
      
              res = runspatialdecon(sc_deconv_table_seurat,
                                    bg = 0.01,
                                    X = custom_matrix,
                                     align_genes = TRUE)
      
              result_df<-as.data.frame(res$prop_of_nontumor)
              result_df$cell_type<-row.names(result_df)
              row.names(result_df)=NULL
      
              result_df$cell_type<-mgsub::mgsub((result_df$cell_type),c("B cells","Monocytes","T cells","NK cells",
                                          "Dendritic cells","Progenitors"), c("B cell","Monocyte","T cell","NK cell", 
                                          "Dendritic cell","Progenitor"))
          }
  
          if ((matrix_type=="safeTME") | (missing(single_cell_data)))
          {
        
              res = runspatialdecon(sc_deconv_table_seurat,
                         bg = 0.01,
                         X = safeTME,
                         align_genes = TRUE)
  
              result_df<-as.data.frame(res$prop_of_nontumor)
              result_df$cell_type<-row.names(result_df)
              row.names(result_df)=NULL


              result_df$cell_type<- mgsub::mgsub((result_df$cell_type), c( "macrophages","mast","B.naive", "B.memory" ,
                                "plasma" ,"T.CD4.naive","T.CD4.memory","T.CD8.naive" ,  "T.CD8.memory" ,"NK", "pDCs",
                                "mDCs","monocytes.C" , "monocytes.NC.I" ,   "neutrophils","Treg", "endothelial.cells", 
                                "fibroblasts"), c("Macropahge","Mast cell","B cell","B cell","Plasma cell","T cell",
                                "T cell","T cell","T cell","NK cell","Dendritic cell","Dendritic cell","Monocyte",
                                "Monocyte","Neutrophil","T cell","Endothelial cell","Fibroblast cell"))

      }

      return(result_df)

}



