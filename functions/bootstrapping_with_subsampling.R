## Author: Aakanksha Singh
## Date: 31 August 2022
## Version: 1.0
## Aim:  Bootstrapping the deconvolution with sample sizes

## Input: (a) A deconvolution dataframe
## Output: (b) A dataframe with cell deconvolutions from Quantiseq,epic,mcp,spotlight, spatial_deco with tme and spatial decon with your signature matrix

bootstrapping_with_subsampling <- function(deconv_table)

{
  sample1<-deconv_table

## Creating sample replicates
set.seed(245)
sample2<-sample1[sample(c(1:nrow(sample1)),size = 0.9*nrow(sample1)),]

set.seed(282)
sample3<-sample2[sample(c(1:nrow(sample2)),size = 0.9*nrow(sample2)),]

## Calling the deconvolution methods

# Spotlight

spl1 <-spotlight_deconv(sample1,signature_reference)
spl1$sample_number<- c(rep("sample1",nrow(spl1)))

spl2 <-spotlight_deconv(sample2,signature_reference)
spl2$sample_number<- c(rep("sample2",nrow(spl2)))

spl3 <-spotlight_deconv(sample3,signature_reference)
spl3$sample_number<- c(rep("sample3",nrow(spl3)))

# Spatial deconv (with TME)
sd_tme1 <- spatial_deconv(sample1,signature_reference,matrix_type = "safeTME")
sd_tme1$sample_number<- c(rep("sample1",nrow(sd_tme1)))

sd_tme2 <- spatial_deconv(sample2,signature_reference,matrix_type = "safeTME")
sd_tme2$sample_number<- c(rep("sample2",nrow(sd_tme2)))

sd_tme3 <- spatial_deconv(sample3,signature_reference,matrix_type = "safeTME")
sd_tme3$sample_number<- c(rep("sample3",nrow(sd_tme3)))


# Spatial deconv (with custom matrix)
sd_cus1 <- spatial_deconv(sample1,signature_reference,matrix_type = "custom")
sd_cus1$sample_number<- c(rep("sample1",nrow(sd_cus1)))

sd_cus2 <- spatial_deconv(sample2,signature_reference,matrix_type = "custom")
sd_cus2$sample_number<- c(rep("sample2",nrow(sd_cus2)))

sd_cus3 <- spatial_deconv(sample3,signature_reference,matrix_type = "custom")
sd_cus3$sample_number<- c(rep("sample3",nrow(sd_cus3)))


# To be put in the generating matrix function, used in all of them

restest<-lapply(list(sample1, sample2, sample3), function(w) { rownames(w)<-w[,"HUGO"] ; w<-w[,-1]}) 
sample1<-restest[[1]]
sample2<-restest[[2]]
sample3<-restest[[3]]


# Quantiseq
sa1<-na.omit(sample1)
sa2 <-na.omit(sample2)
sa3 <- na.omit(sample3)

quant1 <- deconvolute(sa1, "quantiseq",tumor = FALSE)
quant1$sample_number<- c(rep("sample1",nrow(quant1)))

quant2 <- deconvolute(sa2, "quantiseq",tumor = FALSE)
quant2$sample_number<- c(rep("sample2",nrow(quant2)))

quant3 <- deconvolute(sa3, "quantiseq",tumor = FALSE)
quant3$sample_number<- c(rep("sample3",nrow(quant3)))


# Epic
epic1 <- deconvolute(sa1,"epic",tumor=FALSE)
epic1$sample_number<- c(rep("sample1",nrow(epic1)))

epic2 <- deconvolute(sa2,"epic",tumor=FALSE)
epic2$sample_number<- c(rep("sample2",nrow(epic2)))

epic3 <- deconvolute(sa3,"epic",tumor=FALSE)
epic3$sample_number<- c(rep("sample3",nrow(epic3)))


# MCP
mcp1 <- deconvolute(sa1, "mcp_counter")
mcp1$sample_number<- c(rep("sample1",nrow(mcp1)))

mcp2 <- deconvolute(sa2, "mcp_counter")
mcp2$sample_number<- c(rep("sample2",nrow(mcp2)))

mcp3 <- deconvolute(sa3, "mcp_counter")
mcp3$sample_number<- c(rep("sample3",nrow(mcp3)))


# Converting mcp to fraction snad adding the deconv type column
mcp_tables <-lapply(list(mcp1,mcp2,mcp3), function(i) 
{i$pbmc8k_TPM=i$pbmc8k_TPM/sum(i$pbmc8k_TPM)
i$pbmc6k_TPM=i$pbmc6k_TPM/sum(i$pbmc6k_TPM)
i$liver_TPM=i$liver_TPM/sum(i$liver_TPM)
i$deconv_type<-c(rep("MCP",nrow(i)))
return (i)})

print(" All the deconvolutions completed successfully")

# Making a singular qunatitative table

# Quantiseq tables
quantiseq_tables <- lapply(list(quant1,quant2,quant3), function(i) {i$deconv_type<-c(rep("Quantiseq",nrow(i)));return (i)})

# Epic tables
epic_tables <- lapply(list(epic1,epic2,epic3), function(i) {i$deconv_type<-c(rep("Epic",nrow(i)));return (i)})

# Mcp Tables
#mcp_tables <-lapply(list(mcp1,mcp2,mcp3), function(i) {i$deconv_type<-c(rep("MCP",nrow(i)));return (i)})

# Spotlight tables
spotlight_tables <- lapply(list(spl1,spl2,spl3), function(i) {i$deconv_type<-c(rep("Splotlight",nrow(i)));return (i)})

# Spatial Decov (TME) tables
spatial_deconv_tme_tables <- lapply(list(sd_tme1,sd_tme2,sd_tme3), function(i) {i$deconv_type<-c(rep("Spatial_Deconv_tme",nrow(i)));return (i)})

# Spatial Deconv Customised tables
spatial_deconv_cus_tables <- lapply(list(sd_cus1,sd_cus2,sd_cus3), function(i) {i$deconv_type<-c(rep("Spatial_Deconv_cus",nrow(i)));return (i)})

# Joining all the tables
quantitative_tables <- rbind(rbindlist(quantiseq_tables),rbindlist(epic_tables), rbindlist(mcp_tables),rbindlist(spotlight_tables), rbindlist(spatial_deconv_tme_tables), rbindlist(spatial_deconv_cus_tables))

print("All quantitative tables joined succesfully")

# Grouping one of cells
quantitative_tables$cell_type<-str_replace_all(quantitative_tables$cell_type, fixed(c("T cell CD4+ (non-regulatory)" = "T cell" , "T cell CD8+" = "T cell","T cell regulatory (Tregs)" = "T cell", "T cell CD4+" = "T cell", "Macrophage M1" = "Macrophage",    "Macrophage M2" = "Macrophage")))

# Summarising all the T cells together
summarised_quantitative_table<-aggregate(cbind(liver_TPM,pbmc8k_TPM,pbmc6k_TPM) ~ cell_type + deconv_type + sample_number , data = quantitative_tables, FUN= sum , na.rm=TRUE)

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
#    geom_smooth(method = "lm", formula = y~x, col = "red") +
#    scale_shape_manual(values=c(0,1,2,3, 16, 17,11,4,8))  + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5)) +   
#    theme(axis.text.y   = element_text(size=14),
#          axis.text.x   = element_text(size=14),
#          axis.title.y  = element_text(size=14),
#          axis.title.x  = element_text(size=14),
#          panel.background = element_blank(),
#          panel.grid.major = element_blank(), 
#          panel.grid.minor = element_blank(),
#          axis.line = element_line(colour = "black"),
#          panel.border = element_rect(colour = "black", fill=NA, size=1)) +                                                                                                                                      
#    xlab("Actual proportions") +
#    ylab("Predicted proportions")
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