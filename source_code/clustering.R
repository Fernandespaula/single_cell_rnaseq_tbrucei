
library(optparse)
library(Seurat)
library(ggplot2)

ClusteringAnalysis = function(inputdir = inputdir,
                              proj_name = proj_name,
                              npcs_dim = npcs_dim,
                              outdir = outdir) {
  
  integrated_data = readRDS(inputdir)
  clust_resolutions = c(0.2, 0.4, 0.5, 0.6, 0.8, 1)
  integrated_data@project.name = proj_name
  message(paste0("Run PCA 1:" ,length(npcs_dim), "..."))
  integrated_data = RunPCA(integrated_data, npcs = 50, verbose = FALSE)
  
  message("RunUMAP, RuntSNE...")
  integrated_data = RunUMAP(integrated_data, reduction = "pca", dims = npcs_dim)
  integrated_data = RunTSNE(integrated_data, reduction = "pca", dims = npcs_dim)
  
  message("Find neighbors...")
  integrated_data = FindNeighbors(integrated_data, reduction = "pca", dims = npcs_dim)
  
  message("Find clusters...")
  integrated_data = FindClusters(integrated_data, reduction = "pca", resolution = clust_resolutions)
  
  saveRDS(integrated_data, paste0(outdir, "/", integrated_data@project.name,"_", length(npcs_dim),"pcas_clustering.RDS"))
  
  return(integrated_data)
}

 integrated_clustering_output = ClusteringAnalysis(inputdir  = "/Users/paulafernandes/sc_projects_output/fly_integration/sg_rbp6_d2/sg_rbp6_d2_integrated_sct_2000.RDS",
                                                   npcs_dim = 1:10,
                                                   proj_name = "sg_rbp6d2_extra",
                                                   outdir = "/Users/paulafernandes/sc_projects_output/fly_integration/sg_rbp6_d2")




