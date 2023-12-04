## Author: Aakanksha Singh
## Date: 23 September 2022
## Version: 1.0
## Aim: Processing a datset using spatial deconv which is designed specifically for data from Nanostring GeoMx spatial profiler.
##      In this updated method, I process the data like I would process a normal nanostring dataset
##      How is it different from the spatial_deconv()
##      1. Use the matrix (cut out the commands to convert it to a seurat object)
##      2. Use the cell_counts featute of the spatialdeconv to get cellsper100

# There are two functions in this file, spatial_deconv_updated() when you have the number of cells in the bulk sample and a
# second function called spatial_deconv_updated_without_cell_number() for when we do not have the number of cells in the 
# the pseudobulk sample
#   Input: (a) The deconvolution table in the form of a data frame
#   Output: (a) A dataframe with cell deconvolutions from spatial_deco with tme and spatial decon with your signature matrix

spatial_deconv_updated <- function(spatial_data, single_cell_data, matrix_type = c("safeTME", "custom","LM22","TIL10","Immunostates","Glioma_custom_matrix"))
  
{
    lm22_mat= lm22 %>% column_to_rownames(var = "Gene symbol") %>% as.matrix()
    til10_mat = til10 %>% column_to_rownames(var="GeneSymbol") %>% as.matrix()
    immunostates_mat= immunostates %>% column_to_rownames(var= "Gene") %>% as.matrix()
    glioma_sce_mat= glioma_sce %>% column_to_rownames("GeneSymbol") %>% as.matrix()
  
  # Data wrangling for Spatial Deconv
  
  ## Coverting col to rownames and converting dataframe to numerical matrix
  rownames(spatial_data)<-spatial_data$HUGO
  spatial_data<-data.matrix(spatial_data[-1])
  
  # Running spatialdeconv according to the type of signature matrix
  
  if ((matrix_type=="custom") & (!missing(single_cell_data)))
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
    
    
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = custom_matrix,
                        cell_counts = aoi_table$aoi
                        )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100$cell_type<- mgsub::mgsub((cellsper100$cell_type), c("B cells","Monocytes","T cells","NK cells","Dendritic cells","Progenitors" ), 
                                                                  c("B cell","Monocyte","T cell","NK cell","Dendritic cell","Progenitor cell"))
  # Change the naming if the cell annotation differs in the single cell data
    
  }
  
  else if (matrix_type=="safeTME")
  { 
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = safeTME,
                        cell_counts = aoi_table$aoi
                        )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100$cell_type<- mgsub::mgsub((cellsper100$cell_type), c( "macrophages","mast","B.naive", "B.memory" ,
                                                                     "plasma" ,"T.CD4.naive","T.CD4.memory","T.CD8.naive" ,  "T.CD8.memory" ,"NK", "pDCs",
                                                                     "mDCs","monocytes.C" , "monocytes.NC.I" ,   "neutrophils","Treg", "endothelial.cells", 
                                                                     "fibroblasts"), c("Macrophage","Mast cell","B cell","B cell","Plasma cell","T cell",
                                                                                       "T cell","T cell","T cell","NK cell","Dendritic cell","Dendritic cell","Monocyte",
                                                                                       "Monocyte","Neutrophil","T cell","Endothelial cell","Fibroblast cell"))
      
  }
  else if (matrix_type=="TIL10")
  { 
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = til10_mat,
                        cell_counts = aoi_table$aoi
                        )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100$cell_type<- mgsub::mgsub((cellsper100$cell_type), c("B.cells" ,"Macrophages.M1" , "Macrophages.M2" , "Monocytes", "Neutrophils" ,"NK.cells" ,"T.cells.CD4","T.cells.CD8","Tregs","Dendritic.cells"), 
                                         c("B cell","Macrophage","Macrophage","Monocyte","Neutrophil","NK cell","T cell","T cell","T cell","Dendritic cell"))
    
  }
  
  else if (matrix_type=="LM22")
  { 
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = lm22_mat,
                        cell_counts = aoi_table$aoi
                        )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100$cell_type<- mgsub::mgsub((cellsper100$cell_type), c("B cells naive","B cells memory" ,"Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting" ,"T cells CD4 memory activated",
                                                                    "T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0",              
                                                                    "Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils",                 
                                                                    "Neutrophils","T cells regulatory (Tregs)"), 
                                                                   c("B cell","B cell" ,"Plasma cell","T cell","T cell","T cell" ,"T cell",
                                                                     "T cell","T cell","T cell","NK cell","NK cell","Monocyte","Macrophage",              
                                                                     "Macrophage","Macrophage","Dendritic cell","Dendritic cell","Mast cell","Mast cell","Eosinophil",                 
                                                                      "Neutrophil","T cell"))
  }
  
  else if (matrix_type=="Immunostates")
  { 
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = immunostates_mat,
                        cell_counts = aoi_table$aoi
                        )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100$cell_type<- mgsub::mgsub((cellsper100$cell_type), c("CD14_positive_monocyte","CD16_positive_monocyte","CD4_positive_alpha_beta_T_cell","CD56bright_natural_killer_cell","CD56dim_natural_killer_cell",   
                                                                    "CD8_positive_alpha_beta_T_cell","MAST_cell","basophil","eosinophil","gamma_delta_T_cell","hematopoietic_progenitor",      
                                                                    "macrophage_m0","macrophage_m1","macrophage_m2","memory_B_cell","myeloid_dendritic_cell","naive_B_cell","neutrophil","plasma_cell","plasmacytoid_dendritic_cell"), 
                                         c("Monocyte","Monocyte","T cell","NK cell","NK cell",   
                                           "T cell","Mast cell","Basophil","Eosinophil","T cell","Hematopoietic progenitor",      
                                           "Macrophage","Macrophage","Macrophage","B cell","Dendritic cell","B cell","Neutrophil","Plasma cell","Plasmacytoid Dendritic cell"))
  }
  
  else if (matrix_type=="Glioma_custom_matrix")
  { 
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = glioma_sce_mat,
                        cell_counts = aoi_table$aoi
                        )
    
    cellsper100<-as.data.frame(resu$cell.counts$cells.per.100)
    sum_of_cols <- colSums(cellsper100 [,])
    cellsper100<-rbind(cellsper100, c(100-sum_of_cols))
    rownames(cellsper100)[nrow(cellsper100)] <-"uncharacterized cell"
    cellsper100$cell_type<-rownames(cellsper100)
    cellsper100$cell_type<- mgsub::mgsub((cellsper100$cell_type), c("malignant cell","macrophage","mural cell","dendritic cell","microglial cell","monocyte",                      
                                                                    "oligodendrocyte","endothelial cell","mature T cell","oligodendrocyte precursor cell","mast cell","B cell",                        
                                                                    "plasma cell","natural killer cell","astrocyte","radial glial cell","neuron"), 
                                         c("Malignant cell","Macrophage","Mural cell","Dendritic cell","Microglial cell","Monocyte",                      
                                           "Oligodendrocyte","Endothelial cell","T cell","Oligodendrocyte precursor cell","Mast cell","B cell",                        
                                           "Plasma cell","NK cell","Astrocyte","Radial glial cell","Neuron"))
    
  }
  
  else
  { print("Check your input")
  }
  
  result_of_function= cellsper100 %>% 
    remove_rownames() %>% 
    group_by(cell_type) %>% 
    summarise_all(sum, na.rm=TRUE)
  
  return(result_of_function)
  
}


###########################################################################################################################################################################################
spatial_deconv_updated_without_cell_number <- function(spatial_data, single_cell_data, matrix_type = c("safeTME", "custom","LM22","TIL10","Immunostates","Glioma_custom_matrix"))
  
{
  lm22_mat= lm22 %>% column_to_rownames(var = "Gene symbol") %>% as.matrix()
  til10_mat = til10 %>% column_to_rownames(var="GeneSymbol") %>% as.matrix()
  immunostates_mat= immunostates %>% column_to_rownames(var= "Gene") %>% as.matrix()
  glioma_sce_mat= glioma_sce %>% column_to_rownames("GeneSymbol") %>% as.matrix()
  
  # Data wrangling for Spatial Deconv
  
  ## Coverting col to rownames and converting dataframe to numerical matrix
  rownames(spatial_data)<-spatial_data$HUGO
  spatial_data<-data.matrix(spatial_data[-1])
  
  # Running spatialdeconv according to the type of signature matrix
  
  if ((matrix_type=="custom") & (!missing(single_cell_data)))
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
    
    
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = custom_matrix
                        )
    
    res_df <-as.data.frame(resu$prop_of_all)
       res_df$cell_type<-row.names(res_df)
       row.names(res_df)=NULL
     
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c("B cells","Monocytes","T cells","NK cells","Dendritic cells","Progenitors" ), 
                                         c("B cell","Monocyte","T cell","NK cell","Dendritic cell","Progenitor cell"))
    # Change the naming if the cell annotation differs in the single cell data
    
  }
  
  else if (matrix_type=="safeTME")
  { 
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = safeTME,
                        cell_counts = aoi_table$aoi
    )
    
    res_df <-as.data.frame(resu$prop_of_all)
    res_df$cell_type<-row.names(res_df)
    row.names(res_df)=NULL
    
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c( "macrophages","mast","B.naive", "B.memory" ,
                                                                     "plasma" ,"T.CD4.naive","T.CD4.memory","T.CD8.naive" ,  "T.CD8.memory" ,"NK", "pDCs",
                                                                     "mDCs","monocytes.C" , "monocytes.NC.I" ,   "neutrophils","Treg", "endothelial.cells", 
                                                                     "fibroblasts"), c("Macrophage","Mast cell","B cell","B cell","Plasma cell","T cell",
                                                                                       "T cell","T cell","T cell","NK cell","Dendritic cell","Dendritic cell","Monocyte",
                                                                                       "Monocyte","Neutrophil","T cell","Endothelial cell","Fibroblast cell"))
    
  }
  else if (matrix_type=="TIL10")
  { 
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = til10_mat,
                        cell_counts = aoi_table$aoi
    )
    
    res_df <-as.data.frame(resu$prop_of_all)
    res_df$cell_type<-row.names(res_df)
    row.names(res_df)=NULL
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c("B.cells" ,"Macrophages.M1" , "Macrophages.M2" , "Monocytes", "Neutrophils" ,"NK.cells" ,"T.cells.CD4","T.cells.CD8","Tregs","Dendritic.cells"), 
                                         c("B cell","Macrophage","Macrophage","Monocyte","Neutrophil","NK cell","T cell","T cell","T cell","Dendritic cell"))
    
  }
  
  else if (matrix_type=="LM22")
  { 
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = lm22_mat,
                        cell_counts = aoi_table$aoi
    )
    
    res_df <-as.data.frame(resu$prop_of_all)
    res_df$cell_type<-row.names(res_df)
    row.names(res_df)=NULL
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c("B cells naive","B cells memory" ,"Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting" ,"T cells CD4 memory activated",
                                                                    "T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0",              
                                                                    "Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils",                 
                                                                    "Neutrophils","T cells regulatory (Tregs)"), 
                                         c("B cell","B cell" ,"Plasma cell","T cell","T cell","T cell" ,"T cell",
                                           "T cell","T cell","T cell","NK cell","NK cell","Monocyte","Macrophage",              
                                           "Macrophage","Macrophage","Dendritic cell","Dendritic cell","Mast cell","Mast cell","Eosinophil",                 
                                           "Neutrophil","T cell"))
  }
  
  else if (matrix_type=="Immunostates")
  { 
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = immunostates_mat,
                        cell_counts = aoi_table$aoi
    )
    
    res_df <-as.data.frame(resu$prop_of_all)
    res_df$cell_type<-row.names(res_df)
    row.names(res_df)=NULL
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c("CD14_positive_monocyte","CD16_positive_monocyte","CD4_positive_alpha_beta_T_cell","CD56bright_natural_killer_cell","CD56dim_natural_killer_cell",   
                                                                    "CD8_positive_alpha_beta_T_cell","MAST_cell","basophil","eosinophil","gamma_delta_T_cell","hematopoietic_progenitor",      
                                                                    "macrophage_m0","macrophage_m1","macrophage_m2","memory_B_cell","myeloid_dendritic_cell","naive_B_cell","neutrophil","plasma_cell","plasmacytoid_dendritic_cell"), 
                                         c("Monocyte","Monocyte","T cell","NK cell","NK cell",   
                                           "T cell","Mast cell","Basophil","Eosinophil","T cell","Hematopoietic progenitor",      
                                           "Macrophage","Macrophage","Macrophage","B cell","Dendritic cell","B cell","Neutrophil","Plasma cell","Plasmacytoid Dendritic cell"))
  }
  
  else if (matrix_type=="Glioma_custom_matrix")
  { 
    a <- as.matrix(spatial_data)
    aoi_table <- data.frame("Sample"=c("liver","pbmc8k","pbmc6k"),
                            "aoi"=c(sum(liver_cell_table$inDeconv),sum(pbmc8k_cell_table$inDeconv),sum(pbmc6k_cell_table$inDeconv)))
    resu = spatialdecon(a,
                        bg = 0.1,
                        X = glioma_sce_mat,
                        cell_counts = aoi_table$aoi
    )
    
    res_df <-as.data.frame(resu$prop_of_all)
    res_df$cell_type<-row.names(res_df)
    row.names(res_df)=NULL
    res_df$cell_type<- mgsub::mgsub((res_df$cell_type), c("malignant cell","macrophage","mural cell","dendritic cell","microglial cell","monocyte",                      
                                                                    "oligodendrocyte","endothelial cell","mature T cell","oligodendrocyte precursor cell","mast cell","B cell",                        
                                                                    "plasma cell","natural killer cell","astrocyte","radial glial cell","neuron"), 
                                         c("Malignant cell","Macrophage","Mural cell","Dendritic cell","Microglial cell","Monocyte",                      
                                           "Oligodendrocyte","Endothelial cell","T cell","Oligodendrocyte precursor cell","Mast cell","B cell",                        
                                           "Plasma cell","NK cell","Astrocyte","Radial glial cell","Neuron"))
    
  }
  
  else
  { print("Check your input")
  }
  
  result_of_function= res_df %>% 
    remove_rownames() %>% 
    group_by(cell_type) %>% 
    summarise_all(sum, na.rm=TRUE)
  
  return(result_of_function)
  
}