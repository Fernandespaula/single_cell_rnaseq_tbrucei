
library(optparse)
library(Seurat)
library(ggplot2)


option_list = list(make_option(c("-d", "--inpdir"), default=NULL, help="path to input directory"),
                   
                   make_option(c("-p", "--projname"), type = "character", default=NULL, help="give a project name"),
                   
                   make_option(c("-f", "--numfeatures"), type = "integer", default=NULL, help="number of features for integration"),
                   
                   make_option(c("-o", "--outdir"), default=NULL, help="output directory"))

opt_parser = OptionParser(option_list=option_list, add_help_option = FALSE)
opt = parse_args(opt_parser)
print(opt)
if(
  is.null(opt$inpdir) || 
  is.null(opt$projname) || 
  is.null(opt$numfeatures) ||
  is.null(opt$outdir)
){
  print_help(opt_parser)
  stop("A parameter is missing", call.=FALSE)
}

IntegrationSTC = function(inputdir, projectname, numfeatures, outdir) {
  
  message("Loading data...") 
  files_path = list.files(inputdir, pattern = "\\.RDS", full.names = TRUE)
  
  list_datasets = lapply(files_path, readRDS)
  
  features = SelectIntegrationFeatures(object.list = list_datasets, nfeatures = numfeatures)
  
  message("Perform integration...") 
  data_anchors = FindIntegrationAnchors(object.list = list_datasets, normalization.method = "LogNormalize",
                                        anchor.features = features)
  data_integrated = IntegrateData(anchorset = data_anchors, normalization.method = "LogNormalize")
  
  data_integrated = ScaleData(data_integrated, verbose = FALSE)
  
  saveRDS(data_integrated, paste0(outdir, "/", projectname, "_integrated_lognorm_", numfeatures,".RDS"))
  
  return(data_integrated)
}

integrated_data = IntegrationSTC(inputdir = opt$inpdir, 
                                 projectname = opt$projname,
                                 numfeatures = opt$numfeatures,
                                 outdir = opt$outdir)

#Rscript Integration_sct.R --inpdir /Users/paulafernandes/scrnaseq_projects/output_analysis_v811_2_rbp6 --projname rbp6 --numfeatures 2000 --outdir /Users/paulafernandes/scrnaseq_projects/output_analysis_v811_2_rbp6/test 
