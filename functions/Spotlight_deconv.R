## Author: Aakanksha Singh
## Date: 2 August 2022
## Version: 1.0
## Aim: To dconvolute with spotlight 




# For the function processing_sce
#   Input: (a) A dataframe with the peseudobulk samples (GeneSymbols as rows and colnames as sample)
#          (b) A Processed SingleCellExperiment dataset
#   Output: (a) A dataframe with proportions of all the cell types



#spatial_data_NS<- deconv_table[,c(1,2,3)]
#sce_NS <-wta

spotlight_deconv<-function(spatial_data_NS, sce_NS)
  
{ rownames(spatial_data_NS) <- spatial_data_NS$HUGO
  spatial_data_NS <- spatial_data_NS[,-1]

  ## Converting the spatial count matrix to a SingleCellExperiment Class
  spe_NS <- SpatialExperiment(
  assays = list(counts = spatial_data_NS),
  rowData = rownames(spatial_data_NS),
  colData = colnames(spatial_data_NS))

  ## Variance modelling
  dec <- modelGeneVar(sce_NS)
  plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
  curve(metadata(dec)$trend(x), col = "blue", add = TRUE)

  ## Getting the top 3000 genes
  hvg <- getTopHVGs(dec, n = 3000)
  colLabels(sce_NS) <- colData(sce_NS)$cell
  
  # Get vector indicating which genes are neither ribosomal or mitochondrial
  genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce_NS))
  
  # Compute marker genes
  mgs <- scoreMarkers(sce_NS,subset.row = genes)
  mgs_fil <- lapply(names(mgs),  function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ] ## Changed for testing the hodgkins lymphoma case
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_fil)

  idx <- split(seq(ncol(sce_NS)), sce_NS$cell)
  # downsample to at most 20 per identity & subset
  n_cells <- 20
  cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
  })
  sce_NS <- sce_NS[, unlist(cs_keep)]

  ## Running the SPOTLight function (Training the NMF model)
  library(SPOTlight)
  colnames(sce_NS)<-sce_NS$cell
  #mgs_df
  res <- SPOTlight(
  x = sce_NS,
  y = spe_NS,
  groups = sce_NS$cell,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

  mat<-res$mat
  mod<-res$NMF

  plotTopicProfiles(
    x = mod,
    y = sce_NS$cell,
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)

  
  plotTopicProfiles(
    x = mod,
    y = sce_NS$cell,
    facet = TRUE,
    min_prop = 0.01,
    ncol = 6)


  sign <- basis(mod)
  colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))

  mat_plot<-t(mat)
  mat_plot<-as.data.frame(mat_plot)
  mat_plot <- tibble::rownames_to_column(mat_plot, "cell_type")
  
  # Making the cell type singular
  mat_plot$cell_type=gsub("*s","",mat_plot$cell_type)
  return(mat_plot)

  }


