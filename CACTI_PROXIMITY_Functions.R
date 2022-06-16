
#This function calculates the optimal lambda value in CACTI via a binary search 
lambda_calc = function(C_low,X,Y, alpha = 0.15, lambda_max = 8,lambda_min = 0, lambda_wid = 0.05,seurat_obj = NULL,k_param = 15,resolution = 1.5,n_pcs = NCOL(X)+NCOL(Y)){
  ptm <- proc.time()
  #if seuratobject is not supplied make one and perform variable feature selection, PCA, etc..
  if (is.null(seurat_obj) == TRUE){
    X = as.data.frame(X)
    Y = as.data.frame(Y)
    new_data = as.data.frame(t(cbind(X,Y)))
    colnames(new_data) = c(1:NCOL(new_data)) 
    rownames(new_data) = c(paste0("X:",c(1:ncol(X))),paste0("Y:",c(1:ncol(Y))))
    #create seurat object
    seurat_obj = CreateSeuratObject(counts = new_data,min.cells = -Inf, min.genes = -Inf)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",nfeatures = 20)
    VariableFeatures(object = seurat_obj) = rownames(new_data)
    seurat_obj <- ScaleData(seurat_obj,do.scale = F,do.center = F,scale.max = Inf)
    data = GetAssayData(object = seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, pcs.print = 0,pcs.compute = NROW(new_data))
  }
  data = GetAssayData(object = seurat_obj)
  #perform binary search
  while(lambda_max  -lambda_min > lambda_wid){
    class_err = 0
    #Get test lambda value
    lambda_test = mean(c(lambda_max,lambda_min))
    print(lambda_test)
    #Change pca embeddings based on lambda.
    seurat_obj@reductions$pca@cell.embeddings <- as.matrix(cbind(X,lambda_test*Y))
    rownames(seurat_obj@reductions$pca@cell.embeddings) = c(1:NROW(seurat_obj@reductions$pca@cell.embeddings))
    #Cluster cells. This requires us to recalculate the SNN graph on each iteration
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs,k.param = k_param,force.recalc = T)
    seurat_obj <- FindClusters(seurat_obj,resolution=resolution,n.start = 20)
    #Calculate the Classification Error
    cluster = Idents(seurat_obj)
    for (val in unique(cluster)){
      cells = which(cluster == val)
      class_err = class_err + length(cells) - max(table(C_low[cells]))
    }
    class_err = class_err/NCOL(data)
    print(class_err)
    print(length(unique(Idents(seurat_obj))))
    # If the classification error is below alpha, we know that lambda_star > = lambda_t
    if (class_err < alpha){
      lambda_min = lambda_test
    } else{
      lambda_max = lambda_test
    }
  }
  lambda_star = lambda_min
  print(paste("Total function time:", proc.time() - ptm))
  print(paste("Optimal Lambda Calculated to be:", lambda_star))
  return(lambda_star)
}

###### Resolution calculation also binary search to make sure that number of clusters found is in acceptable range
res_calc = function(res_max,X,Y,lambda, clust_max, clust_min,seurat_obj = NULL,k_param = 15,res_wid = 0.05,res_min = 0,n_pcs = NCOL(X)+NCOL(Y)){
  ptm <- proc.time()
  res_star = 0
  #if seurat object is not supplied, make one
  if (is.null(seurat_obj) == TRUE){
    X = as.data.frame(X)
    Y = as.data.frame(Y)
    new_data = as.data.frame(t(cbind(X,lambda*Y)))
    colnames(new_data) = c(1:NCOL(new_data)) 
    rownames(new_data) = c(paste0("X:",c(1:ncol(X))),paste0("Y:",c(1:ncol(Y))))
    #create seurat object
    seurat_obj = CreateSeuratObject(counts = new_data,min.cells = -Inf, min.genes = -Inf)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 100)
    VariableFeatures(object = seurat_obj) = rownames(new_data)
    seurat_obj <- ScaleData(seurat_obj,do.scale = F,do.center = F,scale.max = Inf)
    data = GetAssayData(object = seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, pcs.print = 0,pcs.compute = NROW(new_data))
  }
  seurat_obj@reductions$pca@cell.embeddings <- as.matrix(cbind(X,lambda*Y))
  rownames(seurat_obj@reductions$pca@cell.embeddings) = c(1:NROW(seurat_obj@reductions$pca@cell.embeddings))
  n = 0
  while(res_max  - res_min > res_wid){
    res_test = mean(c(res_max,res_min))
    #cluster with resolution = res_t
    if(n == 0){
      seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs,k.param = k_param)
    }
    seurat_obj <- FindClusters(seurat_obj,resolution=res_test,n.start = 100)
    #seurat_obj <- FindClusters(seurat_obj, dims.use = 1:(n_pcs), print.output = FALSE,k.param = k_param,resolution = res_test,force.recalc = (n==0)) 
    n = n+1
    cluster = Idents(seurat_obj)
    #Check if clustering is sufficient
    if (length(unique(cluster)) <= clust_max & length(unique(cluster)) >= clust_min){
      res_star  = res_test
      break
    } else {
      res_max = res_test*(length(unique(cluster)) > clust_max) + res_max*(1-(length(unique(cluster)) > clust_max))
      res_min = res_test*(1-(length(unique(cluster)) > clust_max)) + res_min*(length(unique(cluster)) > clust_max)
    } 
  }
  print(paste("Total function time:", proc.time() - ptm))
  print(paste0("Optimal Resolution Calculated to be:", res_star))
  return(seurat_obj)
  
}

##### Prune clusters: Perform hypothesis test to see if clusters are too similar to one another
prune_clust = function(connect, XY_labels, p_val = 0.05,n_sim = 1000){
  ptm <- proc.time()
  N = length(unique(XY_labels))
  clusters = unique(XY_labels)
  for (i in 0:(N-1)){
    P_val = matrix(-1,N-i,N-i)
    for (j in 1:(N-i)){
      #Get indices of cells that belong to cluster j
      j_val = clusters[j]
      j_index = which(XY_labels == j_val)
      for (k in 1:j){
        #Get indices of cells that belong to cluster k
        k_val = clusters[k]
        k_index = which(XY_labels == k_val)
        #Calculate the connectivity of the two clusters
        T_val = sum(connect[j_index,k_index])
        T_wm = c()
        #Simulate the connectivity of the two clusters under the well-mixed null hypothesis
        for (sim in 1:n_sim){
          samp_index = sample(c(j_index,k_index),length(c(j_index,k_index)),replace = F)
          j_samp = samp_index[1:length(j_index)]
          k_samp = samp_index[(length(j_index)+1):length(samp_index)]
          T_wm = append(T_wm,sum(connect[j_samp,k_samp]))
        }
        #Calculate our test statistic
        sigma_wm = sd(T_wm)
        mu_wm = mean(T_wm)
        T_val = (T_val - mu_wm)/sigma_wm
        #Determine the P-value
        P_val[j,k] = 1-pnorm(-T_val)
      }
    }
    diag(P_val) = -1
    #Find the maximum P-value
    max_p = max(P_val)
    P_index = which(P_val == max_p, arr.ind = T)[1,]
    print(max_p)
    #Calculate the bonferroni factor
    n_bon = (N-i)*(N-i-1)/2
    if(max_p > 0.05/n_bon){
      #Update labels if maximum p-value is insignificant
      XY_labels[which(XY_labels == clusters[P_index[2]])] = clusters[P_index[1]]
      clusters = unique(XY_labels)
    }
    if (max_p < 0.05/n_bon){
      break
    }
    
  }
  print(paste("Total function time:", proc.time() - ptm))
  return(XY_labels)
}




proximity_test = function(c1,c2,c3,snn,nsim = 10000){
  #calculate mean similarity
  T_stat = mean(snn[c1,c3])
  n1 = length(c1)
  T_wm = c()
  for(j in c(1:nsim)){
    #shuffle cell type labels
    z1 = sample(c(c1,c2),n1)
    #recalculate avergae similarity
    T_wm = append(T_wm,mean(snn[z1,c3]))
  }
  sigma_wm = sd(T_wm)
  mu_wm = mean(T_wm)
  #get z score of observed test statistic 
  T_val = (T_stat - mu_wm)/sigma_wm
  #get pvalue, low pvalues correspond to high similarity
  p_val = 1-pnorm(T_val)
  return(p_val)
}
