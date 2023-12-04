## Author: Aakanksha Singh
## Date: 17 August 2023
## Version: 1.0
## Aim: Bootstrapping the deconvolution without the sample sizes for the  different signature matrices

# For the function processing_sce
#   Input: (a) The deconvolution table in the form of a data frame
#           (b) Signature refernce to be used when deconvoluting while using sce
#          
#   Output: (a) A dataframe with cell deconvolutions from SpatialDecon with differnt signature matrices
#                   1. safeTME
#                   2. Custom pbmc_8k matrix
#                   3. TIL10
#                   4. LM22
#                   5. Immunostates
#                   6. Glioma_sce_matrix

 signature_matrix_bootstrapping_without_subsampling <- function(deconv_table, signature_reference)

 { 
 
  # Converting mcp to fraction snad adding the deconv type column
  sd_tme <- spatial_deconv_updated(deconv_table,signature_reference,matrix_type = "safeTME")
  sd_tme$deconv_type <- "Spatial_decon_tme"
  
  # Spatial deconv (with custom matrix)
  sd_cus <- spatial_deconv_updated(deconv_table,signature_reference,matrix_type = "custom")
  sd_cus$deconv_type <- "Custom_pbmc8k"
  
  #Spatial decon with TIL10 signature
  sd_til10 <- spatial_deconv_updated(deconv_table,signature_reference,matrix_type = "TIL10")
  sd_til10$deconv_type <- "TIL10"
  
  # Spatial decon with LM22 signature
  sd_lm22 <- spatial_deconv_updated(deconv_table,signature_reference,matrix_type = "LM22")
  sd_lm22$deconv_type <- "LM22"
  
  # Spatial decon with Immunostates matrix
  sd_immunostates <- spatial_deconv_updated(deconv_table,signature_reference,matrix_type = "Immunostates")
  sd_immunostates$deconv_type <- "Immunostates"
  
  # Spatil decon with glioma sce matrix
  sd_glioma_sce <- spatial_deconv_updated(deconv_table,signature_reference,matrix_type = "Glioma_custom_matrix")
  sd_glioma_sce$deconv_type <- "Glioma_sce"
  
  print(" All the deconvolutions with SpatialDecon with different signature matrices completed successfully")
  
  # Making a singular qunatitative table
  
  # Joining all the tables
  quantitative_tables <- rbind(sd_tme, sd_cus, sd_til10, sd_lm22, sd_immunostates, sd_glioma_sce)
  
  print("All quantitative tables joined succesfully")
  
  # Grouping one of cells
   quantitative_tables$cell_type<-str_replace_all(quantitative_tables$cell_type, fixed(c("T cell CD4+ (non-regulatory)" = "T cell" , "T cell CD8+" = "T cell","T cell regulatory (Tregs)" = "T cell", "T cell CD4+" = "T cell", "Macrophage M1" = "Macrophage",    "Macrophage M2" = "Macrophage")))
  
  # Summarising all the T cells together
  summarised_quantitative_table<-aggregate(cbind(liver_TPM,pbmc8k_TPM,pbmc6k_TPM) ~ cell_type + deconv_type , data = quantitative_tables, FUN= sum , na.rm=TRUE)
  
  # Initializing the empty columns
  summarised_quantitative_table[,c("pbmc6k_actual", "pbmc8k_actual","liver_actual")] <- NA
  
  
  print("Summarised quantitative table is being processed")
  
  # Qualitative Plots
  
  # Adding actual proportions to the tables
  for (i in 1:nrow(summarised_quantitative_table)) 
  {
    if (summarised_quantitative_table$cell_type[i]=="B cell") 
    {
      summarised_quantitative_table[i,"pbmc6k_actual"]<- pbmc6k_cell_table[pbmc6k_cell_table$cell_type=="B cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"pbmc8k_actual"] <- pbmc8k_cell_table[pbmc8k_cell_table$cell_type=="B cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"liver_actual"] <- liver_cell_table[liver_cell_table$cell_type=="B cells","inDeconv_Frac"] 
    } 
    
    else if (summarised_quantitative_table$cell_type[i]=="T cell")
    {
      summarised_quantitative_table[i,"pbmc6k_actual"]<- pbmc6k_cell_table[pbmc6k_cell_table$cell_type=="T cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"pbmc8k_actual"] <- pbmc8k_cell_table[pbmc8k_cell_table$cell_type=="T cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"liver_actual"] <- liver_cell_table[liver_cell_table$cell_type=="T cells","inDeconv_Frac"] 
    }
    
    else if (summarised_quantitative_table$cell_type[i]=="Monocyte")
    {
      summarised_quantitative_table[i,"pbmc6k_actual"]<- pbmc6k_cell_table[pbmc6k_cell_table$cell_type=="Monocytes","inDeconv_Frac"]
      summarised_quantitative_table[i,"pbmc8k_actual"] <- pbmc8k_cell_table[pbmc8k_cell_table$cell_type=="Monocytes","inDeconv_Frac"]
      summarised_quantitative_table[i,"liver_actual"] <- liver_cell_table[liver_cell_table$cell_type=="Monocytes","inDeconv_Frac"] 
    }
    
    else if (summarised_quantitative_table$cell_type[i]=="NK cell")
    {
      summarised_quantitative_table[i,"pbmc6k_actual"]<- pbmc6k_cell_table[pbmc6k_cell_table$cell_type=="NK cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"pbmc8k_actual"] <- pbmc8k_cell_table[pbmc8k_cell_table$cell_type=="NK cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"liver_actual"] <- liver_cell_table[liver_cell_table$cell_type=="NK cells","inDeconv_Frac"] 
    }
    
    else
    {
      summarised_quantitative_table[i,"pbmc6k_actual"]<- 0
      summarised_quantitative_table[i,"pbmc8k_actual"] <- 0
      summarised_quantitative_table[i,"liver_actual"] <- 0 
    }
  }
  
  return(summarised_quantitative_table)
}

 
##########################################################################################################################################################
 signature_matrix_bootstrapping_without_subsampling_without_cell_numbers <- function(deconv_table, signature_reference)
   
 { 
   
   # Converting mcp to fraction snad adding the deconv type column
   sd_tme <- spatial_deconv_updated_without_cell_number(deconv_table,signature_reference,matrix_type = "safeTME")
   sd_tme$deconv_type <- "Spatial_decon_tme"
   
   # Spatial deconv (with custom matrix)
   sd_cus <- spatial_deconv_updated_without_cell_number(deconv_table,signature_reference,matrix_type = "custom")
   sd_cus$deconv_type <- "Custom_pbmc8k"
   
   #Spatial decon with TIL10 signature
   sd_til10 <- spatial_deconv_updated_without_cell_number(deconv_table,signature_reference,matrix_type = "TIL10")
   sd_til10$deconv_type <- "TIL10"
   
   # Spatial decon with LM22 signature
   sd_lm22 <- spatial_deconv_updated_without_cell_number(deconv_table,signature_reference,matrix_type = "LM22")
   sd_lm22$deconv_type <- "LM22"
   
   # Spatial decon with Immunostates matrix
   sd_immunostates <- spatial_deconv_updated_without_cell_number(deconv_table,signature_reference,matrix_type = "Immunostates")
   sd_immunostates$deconv_type <- "Immunostates"
   
   # Spatil decon with glioma sce matrix
   sd_glioma_sce <- spatial_deconv_updated_without_cell_number(deconv_table,signature_reference,matrix_type = "Glioma_custom_matrix")
   sd_glioma_sce$deconv_type <- "Glioma_sce"
   
   print(" All the deconvolutions with SpatialDecon with different signature matrices completed successfully")
   
   # Making a singular qunatitative table
   
   # Joining all the tables
   quantitative_tables <- rbind(sd_tme, sd_cus, sd_til10, sd_lm22, sd_immunostates, sd_glioma_sce)
   
   print("All quantitative tables joined succesfully")
   
   # Grouping one of cells
   quantitative_tables$cell_type<-str_replace_all(quantitative_tables$cell_type, fixed(c("T cell CD4+ (non-regulatory)" = "T cell" , "T cell CD8+" = "T cell","T cell regulatory (Tregs)" = "T cell", "T cell CD4+" = "T cell", "Macrophage M1" = "Macrophage",    "Macrophage M2" = "Macrophage")))
   
   # Summarising all the T cells together
   summarised_quantitative_table<-aggregate(cbind(liver_TPM,pbmc8k_TPM,pbmc6k_TPM) ~ cell_type + deconv_type , data = quantitative_tables, FUN= sum , na.rm=TRUE)
   
   # Initializing the empty columns
   summarised_quantitative_table[,c("pbmc6k_actual", "pbmc8k_actual","liver_actual")] <- NA
   
   
   print("Summarised quantitative table is being processed")
   
   # Qualitative Plots
   
   # Adding actual proportions to the tables
   for (i in 1:nrow(summarised_quantitative_table)) 
   {
     if (summarised_quantitative_table$cell_type[i]=="B cell") 
     {
       summarised_quantitative_table[i,"pbmc6k_actual"]<- pbmc6k_cell_table[pbmc6k_cell_table$cell_type=="B cells","inDeconv_Frac"]
       summarised_quantitative_table[i,"pbmc8k_actual"] <- pbmc8k_cell_table[pbmc8k_cell_table$cell_type=="B cells","inDeconv_Frac"]
       summarised_quantitative_table[i,"liver_actual"] <- liver_cell_table[liver_cell_table$cell_type=="B cells","inDeconv_Frac"] 
     } 
     
     else if (summarised_quantitative_table$cell_type[i]=="T cell")
     {
       summarised_quantitative_table[i,"pbmc6k_actual"]<- pbmc6k_cell_table[pbmc6k_cell_table$cell_type=="T cells","inDeconv_Frac"]
       summarised_quantitative_table[i,"pbmc8k_actual"] <- pbmc8k_cell_table[pbmc8k_cell_table$cell_type=="T cells","inDeconv_Frac"]
       summarised_quantitative_table[i,"liver_actual"] <- liver_cell_table[liver_cell_table$cell_type=="T cells","inDeconv_Frac"] 
     }
     
     else if (summarised_quantitative_table$cell_type[i]=="Monocyte")
     {
       summarised_quantitative_table[i,"pbmc6k_actual"]<- pbmc6k_cell_table[pbmc6k_cell_table$cell_type=="Monocytes","inDeconv_Frac"]
       summarised_quantitative_table[i,"pbmc8k_actual"] <- pbmc8k_cell_table[pbmc8k_cell_table$cell_type=="Monocytes","inDeconv_Frac"]
       summarised_quantitative_table[i,"liver_actual"] <- liver_cell_table[liver_cell_table$cell_type=="Monocytes","inDeconv_Frac"] 
     }
     
     else if (summarised_quantitative_table$cell_type[i]=="NK cell")
     {
       summarised_quantitative_table[i,"pbmc6k_actual"]<- pbmc6k_cell_table[pbmc6k_cell_table$cell_type=="NK cells","inDeconv_Frac"]
       summarised_quantitative_table[i,"pbmc8k_actual"] <- pbmc8k_cell_table[pbmc8k_cell_table$cell_type=="NK cells","inDeconv_Frac"]
       summarised_quantitative_table[i,"liver_actual"] <- liver_cell_table[liver_cell_table$cell_type=="NK cells","inDeconv_Frac"] 
     }
     
     else
     {
       summarised_quantitative_table[i,"pbmc6k_actual"]<- 0
       summarised_quantitative_table[i,"pbmc8k_actual"] <- 0
       summarised_quantitative_table[i,"liver_actual"] <- 0 
     }
   }
   
   return(summarised_quantitative_table)
 }
 