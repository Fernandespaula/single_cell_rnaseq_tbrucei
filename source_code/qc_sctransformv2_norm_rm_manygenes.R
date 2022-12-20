
print("scRNA-seq analysis: 1- Quality control and normalization")
format(Sys.time(), '%d.%m.%Y')

#libraries 
library(optparse)
library(Seurat)
library(ggplot2)
library(glmGamPoi)
library(gridExtra)

# add genes

histones = c('Tb927.11.1790',
             'Tb927.11.1800',
             'Tb927.11.1860',
             'Tb927.11.1870',
             'Tb927.11.1880',
             'Tb927.10.15350',
             'Tb927.11.7350',
             'Tb927.7.6360',
             'Tb927.2.2670',
             'Tb927.1.2430',
             'Tb927.1.2450',
             'Tb927.1.2470',
             'Tb927.1.2490',
             'Tb927.1.2510',
             'Tb927.1.2530', 
             'Tb927.1.2550',
             'Tb927.10.10460',
             'Tb927.10.10470',
             'Tb927.10.10480',
             'Tb927.10.10490',
             'Tb927.10.10500',
             'Tb927.10.10510',
             'Tb927.10.10520',
             'Tb927.10.10530',
             'Tb927.10.10540',
             'Tb927.10.10550',
             'Tb927.10.10560',
             'Tb927.10.10570',
             'Tb927.10.10580',
             'Tb927.10.10590',
             'Tb927.5.4170',
             'Tb927.5.4180',
             'Tb927.5.4190',
             'Tb927.5.4200',
             'Tb927.5.4210',
             'Tb927.5.4220',
             'Tb927.5.4230',
             'Tb927.5.4240',
             'Tb927.5.4250',
             'Tb927.5.4260',
             'Tb927.7.2820',
             'Tb927.7.2830',
             'Tb927.7.2840',
             'Tb927.7.2850',
             'Tb927.7.2860',
             'Tb927.7.2870',
             'Tb927.7.2880',
             'Tb927.7.2890',
             'Tb927.7.2900',
             'Tb927.7.2910',
             'Tb927.7.2920',
             'Tb927.7.2930',
             'Tb927.7.2940')

elong_factor = c("Tb927.10.2090", "Tb927.10.2100", "Tb927.10.2110")


ribosomal_subunits = c("Tb11.v5.0243",   "Tb927.1.3180"  , "Tb927.10.1080" , "Tb927.10.1090"  ,"Tb927.10.1100" , 
                       "Tb927.10.11390" ,"Tb927.10.11540", "Tb927.10.12680","Tb927.10.13500" ,"Tb927.10.13517", 
                       "Tb927.10.14580","Tb927.10.14600", "Tb927.10.14710", "Tb927.10.15120", "Tb927.10.190" ,   "Tb927.10.3280" 
                       ,"Tb927.10.3930" , "Tb927.10.3940" , "Tb927.10.4110" , "Tb927.10.4120" , "Tb927.10.4120" , "Tb927.10.5030" , "Tb927.10.5330" , "Tb927.10.5480" 
                       ,"Tb927.10.560"  , "Tb927.10.5610" , "Tb927.10.7340" , "Tb927.10.8430" , "Tb927.10.9800" , "Tb927.11.10160", "Tb927.11.11230" ,"Tb927.11.11725"
                       , "Tb927.11.11820", "Tb927.11.11830", "Tb927.11.13007", "Tb927.11.14130" ,"Tb927.11.15880", "Tb927.11.15900", "Tb927.11.16280" ,"Tb927.11.2060" 
                       , "Tb927.11.3000" , "Tb927.11.3230" , "Tb927.11.3590" , "Tb927.11.4300"  ,"Tb927.11.4820" , "Tb927.11.6140" , "Tb927.11.6180" , "Tb927.11.6500" 
                       , "Tb927.11.6510" , "Tb927.11.8200",  "Tb927.11.9710" , "Tb927.11.9720",  "Tb927.11.9730",  "Tb927.2.5910" ,  "Tb927.3.1370"  , "Tb927.3.3310"  
                       ,"Tb927.3.3320"  , "Tb927.3.5050" ,  "Tb927.4.1800"  , "Tb927.4.1860" ,  "Tb927.4.2180"  , "Tb927.4.3550" ,  "Tb927.5.1110"  , "Tb927.5.1610"  
                       , "Tb927.5.1820"  , "Tb927.6.4690" ,  "Tb927.6.4980" ,  "Tb927.6.5040"  , "Tb927.6.5120" ,  "Tb927.6.5130" ,  "Tb927.6.720"  ,  "Tb927.7.1040"  
                       , "Tb927.7.1050"  , "Tb927.7.1730" ,  "Tb927.7.1750" ,  "Tb927.7.2340" ,  "Tb927.7.3680" ,  "Tb927.7.5180"  , "Tb927.8.1110" ,  "Tb927.8.1330"  
                       , "Tb927.8.1340" ,  "Tb927.8.5260"  , "Tb927.8.6030" ,  "Tb927.8.6038" ,  "Tb927.8.6150"  , "Tb927.8.6160"  , "Tb927.8.6180"  , "Tb927.9.11380" 
                       ,  "Tb927.9.11410" , "Tb927.9.11470",  "Tb927.9.11490"  ,"Tb927.9.12200" , "Tb927.9.12240" , "Tb927.9.15110" , "Tb927.9.15170"  ,"Tb927.9.15190" 
                       ,"Tb927.9.15210"  ,"Tb927.9.15360" , "Tb927.9.15420" , "Tb927.9.1810"  , "Tb927.9.1850" ,  "Tb927.9.2020"  , "Tb927.9.3920" ,  "Tb927.9.3990"  
                       ,"Tb927.9.5690"  , "Tb927.9.6070")

ribosomal_RNA = c("tmp.1.10", "tmp.1.100", "tmp.1.40", "tmp.1.100", "tmp.1.100" , "Tb927.2.1931", "Tb927.2.1942", "Tb927.2.1953" , "Tb927.2.1964", "Tb927.2.1975", "Tb927.2.1986", "Tb927.2.1997",
                   "Tb927.2.2008", "Tb927.2.1407", "Tb927.2.1416", "Tb927.2.1434" , "Tb927.2.1443", "Tb927.2.1452", "Tb927.2.1500", "Tb927.2.1510", "Tb927.2.1530", "Tb927.2.1540", "Tb927.2.1550",
                   "Tb927.3.3421", "Tb927.3.3422", "Tb927.3.3423" , "Tb927.3.3424", "Tb927.3.3425", "Tb927.3.3427", "Tb927.3.3429" , "Tb927.3.3431", "Tb927.3.3432" , "Tb927.3.3434" , "Tb927.3.3436", 
                   "Tb927.3.3438", "Tb927.3.3439", "Tb927.3.3441", "Tb927.3.3444", "Tb927.3.3445", "Tb927.3.3447" , "Tb927.3.3448",  "Tb927.3.3449", "Tb927.3.3451","Tb927.3.3452" , "Tb927.3.3454", 
                   "Tb927.6.184", "Tb927.6.185", "Tb927.7.6866" , "Tb927.7.6881" ,"Tb927.7.6882" , "Tb927.7.6883", "Tb927.7.6884" , "Tb927.7.6885" , "Tb09_rRNA_4", "Tb09_rRNA_1", "Tb09_rRNA_2", 
                   "Tb09_rRNA_3"  , "Tb10_rRNA_1:rRNA", "Tb10_rRNA_2:rRNA", "Tb10_rRNA_2:rRNA" , "Tb927.11.rRNA_1", "Tb927.11.rRNA_2", "Tb927.11.rRNA_2",  "Tb10-rRNA-1:rRNA")

tubulin = c("Tb927.1.2330", "Tb927.1.2340" ,"Tb927.1.2350", "Tb927.1.2360", "Tb927.1.2370", "Tb927.1.2380", "Tb927.1.2390", "Tb11.v5.0469")

hyp_protein = c("Tb927.11.18710")

mito_genes = c("MT-12SrRNA","MT-9SrRNA","MT-ND8","MT-ND9","MT-ND7", "MT-CO3","MT-CYTB", "MT-ATP6","MT-MURF1","MT-CR3","MT-ND1", 
               "MT-CO2","MT-MURF2","MT-CO1","MT-CR4", "MT-ND4", "MT-ND3","MT-RPS12","MT-ND5")

option_list = list(make_option(c("-d", "--inpdir"), default=NULL, help="path to input directory"),
                   
                   make_option(c("-p", "--projname"), default=NULL, help="project name"),
                   
                   make_option(c("-e", "--mincells"), default=3, help="minimum cells"),
                   
                   make_option(c("-g", "--mingenes"), default=10, help="minimum genes"),
                   
                   make_option(c("-f", "--featuremax"), type = "integer",  default=NULL, help="feature max cutoff"),
                   
                   make_option(c("-n", "--featuremin"), type = "integer",  default=NULL, help="feature min cutoff"),
                   
                   make_option(c("-i", "--countmin"), type = "integer", default=NULL, help="count min cutoff"),
                   
                   make_option(c("-c", "--countmax"), type = "integer", default=NULL, help="count max cutoff"),
                   
                   make_option(c("-m", "--mito"), type = "integer", default=NULL, help="mito percent cutoff"),
                   
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
                                    outdir = "output_analysis")
{
  input_matrix = Read10X(input_matrix)  
  seurat_obj = CreateSeuratObject(input_matrix,  project = project_name, min.cells = min.c, min.features = min.g)
  
  genes = seurat_obj@assays$RNA@counts@Dimnames[[1]] 
  
  seurat_obj[["percent.mito"]] = PercentageFeatureSet(seurat_obj, features = mito_genes)
  
  vlnplot_before_qc = VlnPlot(seurat_obj, 
                               features = c("nFeature_RNA",  "nCount_RNA", "percent.mito"), 
                               ncol = 3, 
                               cols = c("lightblue", "lightyellow", "lightgreen")) + patchwork::plot_annotation(title = "Number of genes (features), UMIs (counts), percentage of kDNA per cell:",
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
  
  message("subsetting data ...")
  message(paste0("feature_cutoffmax", " ", feature_cutoffmax))
  message(paste0("feature_cutoffmin", " ", feature_cutoffmin))
  message(paste0("count_cutoffmax", " ",count_cutoffmax))
  message(paste0("count_cutoffmin", " ",count_cutoffmin))
  message(paste0("mito_cutoff", " ", mito_cutoff))
  
  seurat_obj = subset(seurat_obj, 
                       subset <- nFeature_RNA < feature_cutoffmax & 
                         nFeature_RNA > feature_cutoffmin &
                         nCount_RNA < count_cutoffmax &
                         nCount_RNA > count_cutoffmin &
                         percent.mito < mito_cutoff)
  
  counts = GetAssayData(seurat_obj, assay = "RNA")
  counts = counts[-(which(rownames(counts) %in% ribosomal_RNA)),]
  counts = counts[-(which(rownames(counts) %in% ribosomal_subunits)),]
  counts = counts[-(which(rownames(counts) %in% tubulin)),]
  counts = counts[-(which(rownames(counts) %in% hyp_protein)),]
  counts = counts[-(which(rownames(counts) %in% elong_factor)),]
  counts = counts[-(which(rownames(counts) %in% histones)),]
  
  seurat_obj = subset(seurat_obj, features = rownames(counts))

  vlnplot_after_qc = VlnPlot(seurat_obj, 
                              features = c("nFeature_RNA",  "nCount_RNA", "percent.mito"), 
                              ncol = 3, 
                              cols = c("pink", "lightcoral", "lightblue")) + patchwork::plot_annotation(title = "Number of genes (features), UMIs (counts),RNA and percentage of kDNA per cell:",
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

