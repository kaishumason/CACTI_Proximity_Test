#read in data
#Get relevant libraries
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(ggvis)
library(FNN)
library(geometry)
library(Rtsne)

print("Doing Joint Clustering")
setwd("C:\\Users\\Kaishu\\Dropbox (Penn)\\Visium\\Codex")
source("CACTI_PROXIMITY_Functions.R")

#insert directory that contains data
wd = "C:\\Users\\Kaishu\\Dropbox (Penn)\\Visium\\Codex\\aggr\\data"
setwd(wd)
#read in data
seuratobj = readRDS(paste0('aggr',".seurat_obj.rds"))
#remove artifacts
barcodes = rownames(seuratobj@meta.data)
art_barcodes = barcodes[seuratobj@meta.data$high_level_celltype != "artifact"]
seuratobj =subset(seuratobj,cells = art_barcodes )
counts = seuratobj@assays$Spatial@scale.data
#transpose data
data = t(counts)
#get list of coordinates for each patient
coord = list()
n_pat = length(seuratobj@images)
for(j in c(1:n_pat)){
  coord[[j]] = seuratobj@images[[j]]@coordinates[,1:2]
}

#Low resolution Clustering
connect_low = seuratobj@graphs$Spatial_snn
X_low = Idents(seuratobj)

#get umap
seuratobj <- RunUMAP(seuratobj,reduction = "harmony",dims = 1:10)
#umap_emb = Embeddings(seuratobj[["umap"]])


print("Constructing Niche Features")
#get patient numbers
patients = unique(substr(rownames(data),1,5))
#get length of patient datasets
pat_num = rep(NA,n_pat)
X_harm = seuratobj@reductions$harmony@cell.embeddings
for(j in c(1:n_pat)){
  pat_num[j] = nrow(seuratobj@images[[j]]@coordinates)
}
pat_num = cumsum(pat_num)
pat_num = c(0,pat_num)

ptm <- proc.time()
Y_niche=matrix(NA,nrow = nrow(X_harm),ncol = ncol(X_harm))
for(j in c(1:n_pat)){
  print(patients[j])
  #get patient spatial data
  pat_coord = coord[[j]]
  #get where this patient is in total dataset
  pat_index = c((pat_num[j]+1):pat_num[j+1])
  #get the data
  pat_data = X_harm[pat_index,]
  #perform triangulation
  spatial_neigh = delaunayn(pat_coord)
  spatial_neigh_tc = t(spatial_neigh)
  for(i in c(1:nrow(pat_coord))){
    if(i%%10000 == 0){
      print(i)
    }
    tmp=unique(as.vector(spatial_neigh_tc[,colSums(spatial_neigh_tc==i)>0]))
    tmp=tmp[tmp!=i]
    tmp1=pat_data[tmp,]
    tr1=colMeans(tmp1)
    Y_niche[pat_index[i],] = tr1
  }
}
colnames(Y_niche)=paste0('n1avg_',c(1:ncol(Y_niche)))

print("Get joint matrix")
X_niche = as.matrix(X_harm)
Y_niche = Y_niche


#Do joint clustering
#get lambda
lambda_star = lambda_calc(X_low,X_niche,Y_niche,alpha = 0.15,lambda_max = 4,k_param = 20,resolution = 0.6)
#get the optimal resolution and resulting clustering
res_star = res_calc(res_max = 5,X_niche,Y_niche,lambda = lambda_star, clust_max = 60, clust_min = 45, k_param = 20)
#prune clusters
cluster_xy = prune_clust(connect_low, Idents(res_star))
cluster_xy = Idents(res_star)
print("Almost there")
print("Making UMAP Plots")
res_star <- RunUMAP(res_star,dims = 1:48)
tsne_emb_joint = Embeddings(res_star[["umap"]])

print("Exporting Everything")
#write.table(cluster_xy,file = "cluster_xy")
#write.table(Idents(seuratobj),file = "cluster_x")
#write.table(X_low,file = "cluster_x_low")
#write.table(tsne_emb,file = "X_tsne")
#write.table(tsne_emb_joint,file = "XY_tsne")

#write.table(Y_niche,file = "Y_niche")

saveRDS(res_star,"aggr.seurat_obj.rds")





















