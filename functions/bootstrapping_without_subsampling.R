## Author: Aakanksha Singh
## Date: 23 September 2022
## Version: 1.0
## Aim: Bootstrapping the deconvolution without the sample sizes

# For the function processing_sce
#   Input: (a) The deconvolution table in the form of a data frame
#          
#   Output: (a) A dataframe with cell deconvolutions from Quantiseq,epic,mcp,spotlight, spatial_deco with tme and spatial decon with your signature matrix


bootstrapping_without_subsampling <- function(deconv_table)
  
{
  
  ## Calling the deconvolution methods
  
  # Spotlight
  spl <-spotlight_deconv(deconv_table,signature_reference)
  spl$deconv_type <- "Spotlight"
  
  # Spatial deconv (with TME)
  sd_tme <- spatial_deconv(deconv_table,signature_reference,matrix_type = "safeTME")
  sd_tme$deconv_type <- "Spatial_deconv_tme"

  # Spatial deconv (with custom matrix)
  sd_cus <- spatial_deconv(deconv_table,signature_reference,matrix_type = "custom")
  sd_cus$deconv_type <- "Spatial_deconv_cus"
  
  # To be put in the generating matrix function, used in all of them
  rownames(deconv_table) <- deconv_table[,"HUGO"] 
  deconv_table <- deconv_table[,-1]
  
  # Quantiseq
  sa<-na.omit(deconv_table)
  quant <- deconvolute(sa, "quantiseq",tumor = FALSE)
  quant$deconv_type<-"Quantiseq"
  
  # Epic
  epic <- deconvolute(sa,"epic",tumor=FALSE)
  epic$deconv_type<-"EPIC"
  
  # MCP
  mcp <- deconvolute(sa, "mcp_counter")
  
  # Converting mcp to fraction snad adding the deconv type column
  mcp$pbmc8k_TPM=mcp$pbmc8k_TPM/sum(mcp$pbmc8k_TPM)
  mcp$pbmc6k_TPM=mcp$pbmc6k_TPM/sum(mcp$pbmc6k_TPM)
  mcp$liver_TPM=mcp$liver_TPM/sum(mcp$liver_TPM)
  mcp$deconv_type<- "MCP"
 
  
  print(" All the deconvolutions completed successfully")
  
  # Making a singular qunatitative table
  
  # Joining all the tables
  quantitative_tables <- rbind(quant,epic,mcp,spl,sd_tme, sd_cus)
  
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


# Plot function
#quantitative_plots <- function(column_name, plot_title)
#{
#  nb.cols <- 18
#  mycolors <- colorRampPalette(brewer.pal(11, "Set3"))(nb.cols)
  
#  ggplot(final_table,aes(x=actual_prop,y=column_name, color,group=cell_type)) + geom_point(aes(shape=deconv_type, color=cell_type), size=3) +
#    scale_color_manual (values = c ( "B cell"  =  "#8DD3C7" , "Macrophage"  =  "#FFFFB3","Monocyte"  =  "#BEBADA", "Neutrophil"  =  "#FB8072", "NK cell"  =  "#80B1D3" , "T cell "  =  "#FDB462" ,"Myeloid dendritic cell"  =  "#B3DE69" ,      "uncharacterized cell"  =  "#FCCDE5",  "cytotoxicity score"  =  "#C9B4A7" , "Macrophage/Monocyte"  =  "#D9D9D9",       "Endothelial cell"  =  "#BC80BD" ,  "Cancer associated fibroblast"  =  "#CCEBC5" ,"Progenitor" = "#573F94","Dendritic cell" = "#B3DE69", "Mast cell"="#E6B9A1","Plasma cell"="#E6E6CA","Fibroblast cell"="#CCEBC5")) +
#   geom_smooth(method = "lm", formula = y~x, col = "red") +
#   scale_shape_manual(values=c(0,1,2,3, 16, 17,11,4,8))  + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5)) +   
#   theme(axis.text.y   = element_text(size=14),
#         axis.text.x   = element_text(size=14),
#         axis.title.y  = element_text(size=14),
#         axis.title.x  = element_text(size=14),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)) +                                                                                                                                      
#   xlab("Actual proportions") +
#   ylab("Predicted proportions")
#}

#quantitative_graphs <- ggarrange(
#  quantitative_plots(final_table$liver_TPM, "A. Liver"),
#  quantitative_plots(final_table$pbmc8k_TPM,"B. PBMC_8K"),
#  quantitative_plots(final_table$pbmc6k_TPM,"C. PBMC_6K"),
#  widths = c(1,1,1),  align = "h", nrow=1,ncol = 3,common.legend = TRUE,legend = "bottom"
#)




#pdf("Results/Plots/bootstrapping_20.pdf", height = 8, width = 12)
#quantitative_graphs
#dev.off()

bootstrapping_without_subsampling_updated_for_testing_normalisation <- function(deconv_table,pbmc6k_cell_table,pbmc8k_cell_table,liver_cell_table)
  
{
  
  ## Calling the deconvolution methods
  
  # Spotlight
  spl <-spotlight_deconv(deconv_table,signature_reference)
  spl$deconv_type <- "Spotlight"
  
  # Spatial deconv (with TME)
  sd_tme <- spatial_deconv(deconv_table,signature_reference,matrix_type = "safeTME")
  sd_tme$deconv_type <- "Spatial_deconv_tme"
  
  # Spatial deconv (with custom matrix)
  sd_cus <- spatial_deconv(deconv_table,signature_reference,matrix_type = "custom")
  sd_cus$deconv_type <- "Spatial_deconv_cus"
  
  # To be put in the generating matrix function, used in all of them
  rownames(deconv_table) <- deconv_table[,"HUGO"] 
  deconv_table <- deconv_table[,-1]
  
  # Quantiseq
  sa<-na.omit(deconv_table)
  quant <- deconvolute(sa, "quantiseq",tumor = FALSE)
  quant$deconv_type<-"Quantiseq"
  
  # Epic
  epic <- deconvolute(sa,"epic",tumor=FALSE)
  epic$deconv_type<-"EPIC"
  
  # MCP
  mcp <- deconvolute(sa, "mcp_counter")
  
  # Converting mcp to fraction snad adding the deconv type column
  mcp$pbmc8k=mcp$pbmc8k/sum(mcp$pbmc8k)
  mcp$pbmc6k=mcp$pbmc6k/sum(mcp$pbmc6k)
  mcp$liver=mcp$liver/sum(mcp$liver)
  mcp$deconv_type<- "MCP"
  
  
  print(" All the deconvolutions completed successfully")
  
  # Making a singular qunatitative table
  
  # Joining all the tables
  quantitative_tables <- rbind(quant,epic,mcp,spl,sd_tme, sd_cus)
  
  print("All quantitative tables joined succesfully")
  
  # Grouping one of cells
  quantitative_tables$cell_type<-str_replace_all(quantitative_tables$cell_type, fixed(c("T cell CD4+ (non-regulatory)" = "T cell" , "T cell CD8+" = "T cell","T cell regulatory (Tregs)" = "T cell", "T cell CD4+" = "T cell", "Macrophage M1" = "Macrophage",    "Macrophage M2" = "Macrophage")))
  
  # Summarising all the T cells together
  summarised_quantitative_table<-aggregate(cbind(liver,pbmc8k,pbmc6k) ~ cell_type + deconv_type , data = quantitative_tables, FUN= sum , na.rm=TRUE)
  
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