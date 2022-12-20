
print("scRNA-seq analysis: 1- Quality control and normalization")
format(Sys.time(), '%d.%m.%Y')

#libraries 
library(optparse)
library(Seurat)
library(ggplot2)
library(glmGamPoi)
library(gridExtra)

# add ribosomal and mito genes
ribosomal_RNA = c("tmp.1.10", "tmp.1.100", "tmp.1.40", "Tb927.2.1931", "Tb927.2.1942", "Tb927.2.1953" , "Tb927.2.1964", "Tb927.2.1975", "Tb927.2.1986", "Tb927.2.1997",
                   "Tb927.2.2008", "Tb927.2.1407", "Tb927.2.1416", "Tb927.2.1434" , "Tb927.2.1443", "Tb927.2.1452", "Tb927.2.1500", "Tb927.2.1510", "Tb927.2.1530", "Tb927.2.1540", "Tb927.2.1550",
                   "Tb927.3.3421", "Tb927.3.3422", "Tb927.3.3423" , "Tb927.3.3424", "Tb927.3.3425", "Tb927.3.3427", "Tb927.3.3429" , "Tb927.3.3431", "Tb927.3.3432" , "Tb927.3.3434" , "Tb927.3.3436", 
                   "Tb927.3.3438", "Tb927.3.3439", "Tb927.3.3441", "Tb927.3.3444", "Tb927.3.3445", "Tb927.3.3447" , "Tb927.3.3448",  "Tb927.3.3449", "Tb927.3.3451","Tb927.3.3452" , "Tb927.3.3454", 
                   "Tb927.6.184", "Tb927.6.185", "Tb927.7.6866" , "Tb927.7.6881" ,"Tb927.7.6882" , "Tb927.7.6883", "Tb927.7.6884" , "Tb927.7.6885" , "Tb09_rRNA_4", "Tb09_rRNA_1", "Tb09_rRNA_2", 
                   "Tb09_rRNA_3"  , "Tb10_rRNA_1:rRNA", "Tb10_rRNA_2:rRNA", "Tb927.11.rRNA_1", "Tb927.11.rRNA_2")


mito_genes = c("MT-12SrRNA","MT-9SrRNA","MT-ND8","MT-ND9", "MT-ND7", "MT-CO3","MT-CYTB", "MT-ATP6","MT-MURF1","MT-CR3","MT-ND1", 
               "MT-CO2","MT-MURF2","MT-CO1","MT-CR4", "MT-ND4", "MT-ND3","MT-RPS12")

option_list = list(make_option(c("-d", "--inpdir"), default=NULL, help="path to input directory"),
                   
                   make_option(c("-p", "--projname"), default=NULL, help="project name"),
                   
                   make_option(c("-e", "--mincells"), default=3, help="minimum cells"),
                   
                   make_option(c("-g", "--mingenes"), default=10, help="minimum genes"),
                   
                   make_option(c("-f", "--featuremax"), type = "integer",  default=NULL, help="feature max cutoff"),
                   
                   make_option(c("-n", "--featuremin"), type = "integer",  default=NULL, help="feature min cutoff"),
                   
                   make_option(c("-i", "--countmin"), type = "integer", default=NULL, help="count min cutoff"),
                   
                   make_option(c("-c", "--countmax"), type = "integer", default=NULL, help="count max cutoff"),
                   
                   make_option(c("-m", "--mito"), type = "integer", default=NULL, help="mito percent cutoff"),
                   
                   make_option(c("-r", "--rrna"), type = "integer", default=NULL, help="rrna percent cutoff"),
                   
                   make_option(c("-o", "--outdir"), default=NULL, help="output directory"))

opt_parser = OptionParser(option_list=option_list, add_help_option = FALSE)
opt = parse_args(opt_parser)
print(opt)
if(
  is.null(opt$inpdir) || 
  is.null(opt$projname) || 
  is.null(opt$featuremax) ||
  is.null(opt$featuremin) ||
  is.null(opt$countmin) ||
  is.null(opt$countmax) ||
  is.null(opt$mito) ||
  is.null(opt$rrna) ||
  is.null(opt$outdir)
){
  print_help(opt_parser)
  stop("A parameter is missing", call.=FALSE)
}

runQualityControlSeurat = function(input_matrix, 
                                    project_name = "project_name",
                                    min.c = 5,
                                    min.g = 10,
                                    feature_cutoffmax = 1000,
                                    feature_cutoffmin = 20,
                                    count_cutoffmin = 30,
                                    count_cutoffmax = 1000, 
                                    mito_cutoff = 5,
                                    rrna_cutoff = 10,
                                    outdir = "output_analysis")
{
  input_matrix = Read10X(input_matrix)  
  seurat_obj = CreateSeuratObject(input_matrix,  project = project_name, min.cells = min.c, min.features = min.g)
  
  genes = seurat_obj@assays$RNA@counts@Dimnames[[1]] 
  rRNA_genes = subset(ribosomal_RNA, ribosomal_RNA %in% genes)
  
  seurat_obj[["percent.rRNA"]] = PercentageFeatureSet(seurat_obj, features = rRNA_genes)
  seurat_obj[["percent.mito"]] = PercentageFeatureSet(seurat_obj, features = mito_genes)
  
  vlnplot_before_qc = VlnPlot(seurat_obj, 
                               features = c("nFeature_RNA",  "nCount_RNA", "percent.rRNA", "percent.mito"), 
                               ncol = 4, 
                               cols = c("lightblue", "lightyellow", "pink", "lightgreen")) + patchwork::plot_annotation(title = "Number of genes (features), UMIs (counts), percentage of ribosomal RNA and percentage of kDNA per cell:",
                                                                                                                        theme = theme(plot.title = element_text(size = 8))) 
  
  ggsave(file.path(outdir, paste0(seurat_obj@project.name, ".vlnplot_before_qc.png")), vlnplot_before_qc,  width=10, height=6)
  
  hplot_umis_before = ggplot(seurat_obj@meta.data, aes(x = nCount_RNA)) +
    geom_histogram(bins = 100, colour="darkorchid2", fill="darkorchid2") + 
    geom_vline(aes(xintercept=median(nCount_RNA)), color="orangered4", linetype="dashed", size=0.5) +
    labs(title = "UMI Number Before Filtering", x = "number of UMIs") 
  
  hplot_genes_before = ggplot(seurat_obj@meta.data, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 100, colour="goldenrod2", fill= "goldenrod2") + 
    geom_vline(aes(xintercept=median(nFeature_RNA)), color="orangered4", linetype="dashed", size=0.5) +
    labs(title = "Gene Number Before Filtering", x = "number of genes") 
  
  rrna_perct_plot =  ggplot(seurat_obj@meta.data, aes(x = nCount_RNA, y = percent.rRNA)) +
    geom_point(aes(colour = orig.ident), size = 0.1) + labs(title = "") +
    labs(title = "Relationship of ribosomal RNA percentage and UMI count per cell:", x =  "UMI Count", y = "% rRNA genes") 
  ggsave(file.path(outdir, paste0(seurat_obj@project.name, ".rrna_perct_plot_before_qc.png")),  rrna_perct_plot,  width=10, height=6)
  
  message("subsetting data ...")
  message(paste0("feature_cutoffmax", " ", feature_cutoffmax))
  message(paste0("feature_cutoffmin", " ", feature_cutoffmin))
  message(paste0("count_cutoffmax", " ",count_cutoffmax))
  message(paste0("count_cutoffmin", " ",count_cutoffmin))
  message(paste0("mito_cutoff", " ", mito_cutoff))
  message(paste0("rrna_cutoff", " ",rrna_cutoff))
  
  seurat_obj = subset(seurat_obj, 
                       subset <- nFeature_RNA < feature_cutoffmax & 
                         nFeature_RNA > feature_cutoffmin &
                         nCount_RNA < count_cutoffmax &
                         nCount_RNA > count_cutoffmin &
                         percent.mito < mito_cutoff &
                         percent.rRNA < rrna_cutoff)
  
  vlnplot_after_qc = VlnPlot(seurat_obj, 
                              features = c("nFeature_RNA",  "nCount_RNA", "percent.rRNA", "percent.mito"), 
                              ncol = 4, 
                              cols = c("pink", "lightcoral", "orange", "lightblue")) + patchwork::plot_annotation(title = "Number of genes (features), UMIs (counts), percentage of ribosomal RNA and percentage of kDNA per cell:",
                                                                                                                  theme = theme(plot.title = element_text(size = 8))) 
  ggsave(file.path(outdir, paste0(seurat_obj@project.name, ".vlnplot_after_qc.png")), vlnplot_after_qc,  width=10, height=6)
  
  hplot_umis_after = ggplot(seurat_obj@meta.data, aes(x = nCount_RNA)) +
    geom_histogram(bins = 100, colour="darkorchid2", fill="darkorchid2") + 
    geom_vline(aes(xintercept=median(nCount_RNA)), color="orangered4", linetype="dashed", size=0.5) +
    labs(title = "UMI Number After Filtering", x = "number of UMIs") 
  
  hplot_genes_after = ggplot(seurat_obj@meta.data, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 100, colour="goldenrod2", fill= "goldenrod2") + 
    geom_vline(aes(xintercept=median(nFeature_RNA)), color="orangered4", linetype="dashed", size=0.5) +
    labs(title = "Gene Number After Filtering", x = "number of genes") 
  
  hist_plots = grid.arrange(hplot_umis_before, hplot_genes_before, hplot_umis_after, hplot_genes_after, ncol = 1, nrow = 4)
  ggsave(file.path(outdir, paste0(seurat_obj@project.name, ".histplot_qc.png")), hist_plots,  width=10, height=6)
  
  message("Normalization through SCtransform v2...")
  seurat_obj = SCTransform(seurat_obj, verbose = FALSE, vst.flavor="v2")
  
  saveRDS(seurat_obj, paste0(outdir,"/" ,seurat_obj@project.name, "_qc_sctransform.RDS"))
  
  return(seurat_obj)
}

seuobj_filtered = runQualityControlSeurat(input_matrix = opt$inpdir, 
                                          project_name = opt$projname,
                                          feature_cutoffmax = opt$featuremax,
                                          feature_cutoffmin = opt$featuremin,
                                          count_cutoffmin = opt$countmin,
                                          count_cutoffmax = opt$countmax,
                                          mito_cutoff = opt$mito,
                                          rrna_cutoff = opt$rrna,
                                          outdir = opt$outdir)

#Rscript quality_control_norm.R -d /Users/paulafernandes/projects_2022/test/raw_data/filtered_em_multi_d0 -p rbp6d0 -f 1000 -i 30 -c 1000 -m 5 -r 10 -o /Users/paulafernandes/projects_2022/test/output_analysis
# 
# 
# make_option(c("-d", "--inpdir"), default=NULL, help="path to input directory"),
# 
# make_option(c("-p", "--projname"), default=NULL, help="project name"),
# 
# make_option(c("-e", "--mincells"), default=3, help="minimum cells"),
# 
# make_option(c("-g", "--mingenes"), default=10, help="minimum genes"),
# 
# make_option(c("-f", "--featuremax"),  default=NULL, help="feature max cutoff"),
# 
# make_option(c("-i", "--featuremin"),  default=NULL, help="feature min cutoff"),
# 
# make_option(c("-c", "--countmax"),  default=NULL, help="count max cutoff"),
# 
# make_option(c("-m", "--mito"),  default=NULL, help="mito percent cutoff"),
# 
# make_option(c("-r", "--rrna"),  default=NULL, help="rrna percent cutoff"),
# 
# make_option(c("-o", "--outdir"), default=NULL, help="output directory"))

