## Author: Aakanksha Singh
## Date: 29 July 2022
## Version: 1.0
## Aim: To process single cell data and return a SingleCellExperiment Object 
##     with cell type annotation



# For the function processing_sce
#   Input: (a) A SingleCellExperiment dataset
#          (b) The title of the datset (string)
#   Output: (a) Plot highlight the removed mitochondrial gens
#           (b) Plot of mean log expression vs variance of log expression with a treand line
#           (c) A tSNE plot colored by predicted celltype
#           (d) A table with the number of the predicted cell types
#           (e) The output variable stores the analysed SingleCellExperiment Object with cell type annotations. 




processing_sce <- function(dataset,title) 
{   
  # Step 1: Quality Control
  stats <- high.mito <- list()
  unfiltered<-dataset
  
  current <- dataset
  is.mito <- grep("MT", rowData(dataset)$Symbol_TENx)
  stats <- perCellQCMetrics(dataset, subsets=list(Mito=is.mito))
  high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
  dataset <- dataset[,!high.mito]
  
  current <- unfiltered 
  colData(current) <- cbind(colData(current), stats)
  current$discard <- high.mito
  
  # Note: The code line below generated a plot highlighting the mito genes as outliers
  qcplots<- plotColData(current, x="sum", y="subsets_Mito_percent",
                        colour_by="discard") + scale_x_log10()
  print(qcplots)
  print("Processing sce : Step 1/7 Done")
  
  # Step 2: Log Normalisation
  dataset <- logNormCounts(dataset)
  print("Processing sce : Step 2/7 Done")
  
  # Step 3: Variance Modelling
  dec <- modelGeneVar(dataset)
  hvgs <-getTopHVGs(dataset, prop=0.1)
  
  curdec <- dec
  print(plot(curdec$mean, curdec$total, pch=16, cex=0.5, main= title,
             xlab="Mean of log-expression", ylab="Variance of log-expression"))
  curfit <- metadata(curdec)
  curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
  print("Processing sce : Step 3/7 Done")
  
  # Step 4: Dimensionality Reduction (PCA, tSNE, UMAP)
  
  set.seed(10000)
  dataset <- runPCA(dataset)
  
  set.seed(100000)
  dataset <- runTSNE(dataset)
  
  set.seed(1000000)
  dataset <- runUMAP(dataset)
  print("Processing sce : Step 4/7 Done")
  
  # Step 5: Clustering
  g <- buildSNNGraph(dataset, k=10)
  clust <- igraph::cluster_walktrap(g)$membership
  colLabels(dataset)  <- factor(clust)
  print("Processing sce : Step 5/7 Done")
  # Plot tsne with the cell type labels main
  #tsne <- plotTSNE(dataset, colour_by="label") + ggtitle(dataset)
  
  
  # Step 6: Cell type annotation
  ref<-celldex::MonacoImmuneData(ensembl=TRUE, cell.ont="all")
  predicted_celltype <- SingleR(dataset, ref=ref, labels=ref$label.main)
  
  colData(dataset)$predicted_celltype <- predicted_celltype$labels
  tsne <- plotTSNE(dataset, colour_by="predicted_celltype") + ggtitle(title)
  print(tsne)
  print("Processing sce : Step 6/7 Done")
  
  # Step 7: Adding an addition column that annotes all T cells the same
            # Adding an additional column with consolidated cell types, 
            # removing rows with no Hugo Symbol,
            # Setting the Hugo symbol as rownames instead of ensmbl ids
  
  dataset$cell<-dataset$predicted_celltype
  colData(dataset)$cell <- ifelse(colData(dataset)$cell %in% c("CD8+ T cells", "CD4+ T cells"), "T cells", colData(dataset)$cell)
  rownames(dataset) <- rowData(dataset)$Symbol
  dataset<-dataset[!(is.na(row.names(dataset))),]
  print("Processing sce : Step 7/7 Done")
  print(table(dataset$cell))
  
  return(dataset)
  
}




