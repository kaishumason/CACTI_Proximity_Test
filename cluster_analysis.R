
#c1 = macro
#c2 = rest
#c3 = fib
source("CACTI_PROXIMITY_Functions.R")
#insert directory that contains data
wd = "C:\\Users\\Kaishu\\Dropbox (Penn)\\Visium\\Codex\\aggr\\data"
setwd(wd)
#read in cell types
cell_type = read.csv('corrected_celltypes.csv')
#read in similarity matrix
snn = readRDS('snn.rda')
#read in data
seuratobj = readRDS("aggr.cacti_revised_celltypes.rds")
#remove artifacts
cells = which(seuratobj@meta.data$revised_celltype != "artifact")
barcodes = rownames(seuratobj@meta.data)
art_barcodes = barcodes[cells]
seuratobj =subset(seuratobj,cells = art_barcodes )
#filter similarity matrix and cell type matrix as well (remove artifacts)
snn = snn[cells,cells]
cell_type = cell_type[cells,]

#get patient list
pats = unique(substr(rownames(snn),1,5))
#initialize pvalue vector
p_val = rep(NA,length(pats))
#perform fibroblast proximity test for macrophages vs rest in each patient
counter = 1
for(pat in pats){
  #get cells belonging to index patient
  cells = substr(rownames(snn),1,5) == pat
  temp = subset(seuratobj,cells = rownames(snn)[cells])
  snn_temp = snn[cells,cells]
  CT_temp = cell_type[cells,2]
  #get cells belonging to cell types of macrophages
  c1 = which(substr(CT_temp,1,4)=='macr')
  #fibroblasts
  c3 = which(substr(CT_temp,1,4)=='fibr')
  #everything else
  c2 = c(1:length(CT_temp))[-c(c1,c3)]
  #perform proximity test
  p_val[counter] = proximity_test(c1,c2,c3,snn_temp)
  counter = counter + 1
}









