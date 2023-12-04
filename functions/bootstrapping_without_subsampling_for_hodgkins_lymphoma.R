## Author: Aakanksha Singh
## Date: 4 October 2022
## Version: 1.0
## Aim: Bootstrapping the deconvolution without the sample sizes ( Updated for lymphoma)

# For the function processing_sce
#   Input: (a) The deconvolution table in the form of a data frame
#          
#   Output: (a) A dataframe with cell deconvolutions from Quantiseq,epic,mcp,spotlight, spatial_deco with tme and spatial decon with your signature matrix


bootstrapping_without_subsampling_for_hodgkins_lymphoma <- function(deconv_table)
  
{
  sample1<-deconv_table
  
  ## Calling the deconvolution methods
  
  # Spotlight
  #spl <-spotlight_deconv(sample1,signature_reference)
  #spl$deconv_type <- "Spotlight"
  
  # Spatial deconv (with TME)
  sd_tme <- spatial_deconv(sample1,signature_reference,matrix_type = "safeTME")
  sd_tme$deconv_type <- "Spatial_deconv_tme"
  
  # Spatial deconv (with custom matrix)
  sd_cus <- spatial_deconv(sample1,signature_reference,matrix_type = "custom")
  sd_cus$deconv_type <- "Spatial_deconv_cus"
  
  # To be put in the generating matrix function, used in all of them
  rownames(sample1) <- sample1[,"HUGO"] 
  sample1 <- sample1[,-1]
  
  # Quantiseq
  sa<-na.omit(sample1)
  quant <- deconvolute(sa, "quantiseq",tumor = FALSE)
  quant$deconv_type<-"Quantiseq"
  
  # Epic
  epic <- deconvolute(sa,"epic",tumor=FALSE)
  epic$deconv_type<-"EPIC"
  
  # MCP
  mcp <- deconvolute(sa, "mcp_counter")
  
  # Converting mcp to fraction snad adding the deconv type column
  mcp <- mcp %>%
    mutate_if(is.numeric, funs(./sum(.))) 
  
  mcp$deconv_type<- "MCP"
  
  
  print(" All the deconvolutions completed successfully")
  
  # Making a singular qunatitative table
  
  # Joining all the tables
  quantitative_tables <- rbind(quant,epic,mcp,sd_tme, sd_cus)
  
  print("All quantitative tables joined succesfully")
  
  # Grouping one of cells
  quantitative_tables$cell_type<-str_replace_all(quantitative_tables$cell_type, fixed(c("T cell CD4+ (non-regulatory)" = "T cell" , "T cell CD8+" = "T cell","T cell regulatory (Tregs)" = "T cell", "T cell CD4+" = "T cell", "Macrophage M1" = "Macrophage",    "Macrophage M2" = "Macrophage")))
  
  # Summarising all the T cells together
  summarised_quantitative_table<-aggregate(cbind(targeted_pan_cancer_TPM,targeted_immunology_TPM,wta_TPM) ~ cell_type + deconv_type , data = quantitative_tables, FUN= sum , na.rm=TRUE)
  
  # Initializing the empty columns
  summarised_quantitative_table[,c("wta_actual", "targeted_immunology_actual","targeted_pan_cancer_actual")] <- NA
  
  
  print("Summarised quantitative table is being processed")
  
  # Qualitative Plots
  
  # Adding actual proportions to the tables
  for (i in 1:nrow(summarised_quantitative_table)) 
  {
    if (summarised_quantitative_table$cell_type[i]=="B cell") 
    {
      summarised_quantitative_table[i,"wta_actual"]<- wta_cell_table[wta_cell_table$cell_type=="B cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"targeted_immunology_actual"] <- targeted_immunology_cell_table[targeted_immunology_cell_table$cell_type=="B cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"targeted_pan_cancer_actual"] <- targeted_pan_cancer_cell_table[targeted_pan_cancer_cell_table$cell_type=="B cells","inDeconv_Frac"] 
    } 
    
    else if (summarised_quantitative_table$cell_type[i]=="T cell")
    {
      summarised_quantitative_table[i,"wta_actual"]<- wta_cell_table[wta_cell_table$cell_type=="T cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"targeted_immunology_actual"] <- targeted_immunology_cell_table[targeted_immunology_cell_table$cell_type=="T cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"targeted_pan_cancer_actual"] <- targeted_pan_cancer_cell_table[targeted_pan_cancer_cell_table$cell_type=="T cells","inDeconv_Frac"] 
    }
    
    else if (summarised_quantitative_table$cell_type[i]=="Monocyte")
    {
      summarised_quantitative_table[i,"wta_actual"]<- wta_cell_table[wta_cell_table$cell_type=="Monocytes","inDeconv_Frac"]
      summarised_quantitative_table[i,"targeted_immunology_actual"] <- targeted_immunology_cell_table[targeted_immunology_cell_table$cell_type=="Monocytes","inDeconv_Frac"]
      summarised_quantitative_table[i,"targeted_pan_cancer_actual"] <- targeted_pan_cancer_cell_table[targeted_pan_cancer_cell_table$cell_type=="Monocytes","inDeconv_Frac"] 
    }
    
    else if (summarised_quantitative_table$cell_type[i]=="NK cell")
    {
      summarised_quantitative_table[i,"wta_actual"]<- wta_cell_table[wta_cell_table$cell_type=="NK cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"targeted_immunology_actual"] <- targeted_immunology_cell_table[targeted_immunology_cell_table$cell_type=="NK cells","inDeconv_Frac"]
      summarised_quantitative_table[i,"targeted_pan_cancer_actual"] <- targeted_pan_cancer_cell_table[targeted_pan_cancer_cell_table$cell_type=="NK cells","inDeconv_Frac"] 
    }
    
    else
    {
      summarised_quantitative_table[i,"wta_actual"]<- 0
      summarised_quantitative_table[i,"targeted_immunology_actual"] <- 0
      summarised_quantitative_table[i,"targeted_pan_cancer_actual"] <- 0 
    }
  }
  
  return(summarised_quantitative_table)
}