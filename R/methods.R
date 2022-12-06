#' @title 找出基因表达值高于某阈值的细胞互作
#'
#' @description sender_celltypes和receiver_celltypes的基因表达值均高于设定阈值的细胞互作
#'
#' @param mtx_lr 各细胞类型的基因表达值
#' @param cell.named.types 由细胞类型组成的向量，向量名为细胞名
#' @param LR 先验的配体受体对
#' @param threshold a number,基因表达值应高于的阈值
#' @param sender_celltypes character,发送ligand的细胞类型
#' @param receiver_celltypes character,发送receptor的细胞类型
#' @return interact_threshold 基因表达值高于某阈值的细胞互作
#' @examples
#' mean_lr = matrix(abs(rnorm(60)), 20, 3)
#' colnames(mean_lr) = c("Bcell","CAF","Tcell")
#' rownames(mean_lr) = paste("Gene", 1:20, sep = "")
#' LR = data.frame(ligand=paste("Gene", 1:10, sep = ""),receptor=paste("Gene", 11:20, sep = ""))
#' sender_celltypes = c("Bcell","Tcell")
#' receiver_celltypes = c("CAF","Tcell")
#' threshold = 0.2
#' interact_threshold(mean_lr,LR,threshold,sender_celltypes,receiver_celltypes)
#'
#' @export
#'
interact_threshold <- function(mtx_lr,cell.named.types,LR,sender_celltypes,receiver_celltypes,threshold){

  #取出sender和receiver对应的配体受体对的表达值
  mean_lr <- mean_celltype(mtx_lr,cell.named.types)
  list_mean <- list_mean_lr(mean_lr,LR,sender_celltypes,receiver_celltypes)
  mean.ligand <- list_mean$ligand
  mean.receptor <-  list_mean$receptor
  if(!threshold==0){
    binary_ligand = mean.ligand
    binary_ligand[binary_ligand < threshold] <- 0
    binary_receptor=mean.receptor
    binary_receptor[binary_receptor < threshold] <- 0
  #筛选某基因是否在至少一个细胞类型中表达高于阈值
    judge=rowSums(binary_ligand)>0 & rowSums(binary_receptor) > 0
    binary_ligand=binary_ligand[judge,,drop=FALSE]
    binary_receptor=binary_receptor[judge,,drop=FALSE]
    if(nrow(binary_ligand)==0 | nrow(binary_receptor)==0){
      stop("Cannot found LR Pairs,a smaller threshold is recommended.")
    }
  } else{
    binary_ligand=mean.ligand
    binary_receptor=mean.receptor
  }
  res=data.frame()
  for(cell1 in colnames(binary_ligand)){
    for(cell2 in colnames(binary_receptor)){
      col_res=as.data.frame(cbind(binary_ligand[,cell1],binary_receptor[,cell2]))
      col_res$ligand=rownames(binary_ligand)
      col_res$receptor=rownames(binary_receptor)
      col_res$from=cell1
      col_res$to=cell2
      res=rbind(res,col_res)}
  }
  interaction_threshold=res[res$V1 >= threshold & res$V2 >= threshold,]
  names(interaction_threshold)[1:2]=c("mean_ligand","mean_receptor")
  interaction_threshold$mean_lr <- rowMeans(interaction_threshold[,c("mean_ligand","mean_receptor")])
  interaction_threshold <- interaction_threshold[order(interaction_threshold$mean_lr,decreasing = T),]
  interaction_threshold <- interaction_threshold[,c("ligand","receptor","from","to","mean_ligand","mean_receptor","mean_lr")]
  rownames(interaction_threshold) <- NULL
  return(interaction_threshold)
}
#######
#' @title 找出配体受体基因表达乘积top(product_top)的细胞互作
#'
#' @description sender_celltypes和receiver_celltypes的基因表达值均高于设定阈值的细胞互作
#'
#' @param mtx_lr 各细胞类型的基因表达值
#' @param LR 先验的配体受体对
#' @param cell.named.types 由细胞类型组成的向量，向量名为细胞名
#' @param product_top a number,percentage
#' @param sender_celltypes character,发送ligand的细胞类型
#' @param receiver_celltypes character,产生receptor的细胞类型
#' @return interaction_product 配体受体基因表达乘积高于某阈值的细胞互作
#' @examples
#' mean_lr = matrix(abs(rnorm(60)), 20, 3)
#' colnames(mean_lr) = c("Bcell","CAF","Tcell")
#' rownames(mean_lr) = paste("Gene", 1:20, sep = "")
#' LR = data.frame(ligand=paste("Gene", 1:10, sep = ""),receptor=paste("Gene", 11:20, sep = ""))
#' sender_celltypes = c("Bcell","Tcell")
#' receiver_celltypes = c("CAF","Tcell")
#' product_top = 0.8
#' interact_product(mean_lr,LR,sender_celltypes,receiver_celltypes,product_top)
#'
#' @export
#'
interact_product <- function(mtx_lr,cell.named.types,LR,sender_celltypes,receiver_celltypes,product_top){
  mean_lr <- mean_celltype(mtx_lr,cell.named.types)
  list_mean <- list_mean_lr(mean_lr,LR,sender_celltypes,receiver_celltypes)
  mean.ligand <- list_mean$ligand
  mean.receptor <-  list_mean$receptor

  res=data.frame()
  for(cell1 in colnames(mean.ligand)){
    for(cell2 in colnames(mean.receptor)){
      col_res=as.data.frame(cbind(mean.ligand[,cell1],mean.receptor[,cell2]))
      #  print(head(col_res))
      col_res$ligand=rownames(mean.ligand)
      col_res$receptor=rownames(mean.receptor)
      print(dim(col_res))
      col_res$from=cell1
      col_res$to=cell2
      res=rbind(res,col_res)}}
  lr_product <- res$V1*res$V2
  res <- cbind(res,lr_product)

  #step2.2, 计算LR乘积，并筛选乘积值top的LR Pairs
  res <- res[order(res$lr_product,decreasing = T),]
  res <- res[res$lr_product > quantile(res$lr_product,1-product_top),]
  interaction_product =res
  rownames(interaction_product ) <- NULL
  colnames(interaction_product)[1:2] <- c("mean_ligand","mean_receptor")
  interaction_product <- interaction_product[,c("ligand","receptor","from","to","mean_ligand","mean_receptor","lr_product")]
  return(interaction_product)
}
######
#' @title 找出配体受体marker基因之间的细胞互作
#'
#' @description sender_celltypes和receiver_celltypes的marker基因的细胞互作
#'
#' @param sub_object seurat object
#' @param LR ligand-receptor pairs
#' @param sender_celltypes character,发送ligand的细胞类型
#' @param receiver_celltypes character,产生receptor的细胞类型
#' @param marker_filter logi,True or False
#' @param avg_log2FC Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param padj Limit testing to genes which show,Default is 0.05
#' @return interaction_combination 配体受体marker基因的细胞互作
#' @examples
#' sub_object = readRDS(system.file("extdata","pbmc3k_final.rds", package = "iPhonedev2"))
#' LR = Ramilowski_human
#' sender_celltypes = c("CD8 T","NK")
#' receiver_celltypes = c("DC","Platelet")
#' interact_combination(sub_object,LR,marker_filter=FALSE,
#' avg_log2FC=0.5,padj=0.05,sender_celltypes,receiver_celltypes)
#'
#' @export
#'
interact_combination <- function(sub_object,LR,marker_filter=FALSE,avg_log2FC=1,padj=0.05,sender_celltypes,receiver_celltypes){

  expr.markers <- FindAllMarkers(sub_object,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  expr.markers <- expr.markers[expr.markers$p_val_adj < padj   #参数 0.05
                                & expr.markers$avg_log2FC > avg_log2FC,]

  if(marker_filter == TRUE){
    cluster_mean_expression=data.frame()
    meta=sub_object@meta.data
    for(clust in unique(meta$celltype)){
      cells=rownames(meta)[meta$celltype==clust]
      othercells=rownames(meta)[meta$celltype !=clust]
      exprmat=as.matrix(sub_object@assays$RNA@data)[,cells]
      otherexprmat=as.matrix(sub_object@assays$RNA@data)[,othercells]
      rowmean=rowMeans(exprmat)
      otherrowmean=rowMeans(otherexprmat)
      rowmeandf=data.frame(value1=rowmean,value2=otherrowmean)
      colnames(rowmeandf)=c(clust,paste0("not_",clust))
      if(ncol(cluster_mean_expression)==0){
        cluster_mean_expression=rowmeandf
      }else{
        cluster_mean_expression=cbind(cluster_mean_expression,rowmeandf)
      }
    }
    maxcell=apply(cluster_mean_expression,1,function(x){
      colnames(cluster_mean_expression)[order(x,decreasing = T)[1]]
    })
    expr.markers$maxcell=maxcell[expr.markers$gene] #maxcluster 基因均值哪一类最高
    expr.markers = expr.markers[as.character(expr.markers$cluster)==expr.markers$maxcell,]
  }
  ##############
  expr.markers = expr.markers[expr.markers$cluster %in% c(sender_celltypes,receiver_celltypes),]

  #3.2, 筛选属于marker的LR Pairs
  lr=LR[c(LR$ligand %in% expr.markers$gene & LR$receptor %in% expr.markers$gene),]

  #3.3, 得到LR配对，avg_log2FC和细胞类型注释
  ligand_celltype=merge(lr,expr.markers,by.x="ligand",by.y="gene",sort = F)[,c("ligand","receptor","avg_log2FC","cluster")]
  interaction_combination=merge(ligand_celltype,expr.markers,by.x="receptor",by.y="gene",sort = F)[,c("ligand","receptor","avg_log2FC.x","avg_log2FC.y","cluster.x","cluster.y")]

  colnames(interaction_combination) <- c("ligand","receptor","avg_log2FC.l","avg_log2FC.r","from","to")
  interaction_combination$from <- as.character(interaction_combination$from)
  interaction_combination$to <- as.character(interaction_combination$to)
  interaction_combination=interaction_combination[c(interaction_combination$from %in% sender_celltypes
                                                      & interaction_combination$to %in% receiver_celltypes),]

  return(interaction_combination)
}
#########
#' @title Find ligand-receptor pairs from expression matrix
#'
#' @param object Seurat object
#' @param DB database of ligand-receptor pairs
#' @param celltype Seurat object的细胞类型
#' @param method "expression_threshold", "expression_product", "differential_combinations"
#' @param threshold a number,基因表达均值应高于的阈值
#' @param product_top a number,配体受体对的基因表达均值乘积的percentage
#' @param receiver_celltypes character,产生receptor的细胞类型
#' @param sender_celltypes character,发送ligand的细胞类型
#' @param marker_filter logi,True or False
#' @param avg_log2FC Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param padj Limit testing to genes which show,Default is 0.05
#'
#' @return ligand-receptor pairs from expression matrix
#' @examples
#' load(url("https://figshare.com/ndownloader/files/37663929"))
#' lrDB <- CellChatDB.human
#' receiver_celltypes=c("T cell","CAF","Endothelial","Macrophage","myofibroblast")
#' sender_celltypes="primary"
#' interaction1 <- find_interactions(object=hnscc_primary,DB=lrDB,celltype="celltype",
#'                                   method="expression_threshold",
#'                                   receiver_celltypes=receiver_celltypes,
#'                                   sender_celltypes=sender_celltypes,threshold=3)
#' interaction2 <- find_interactions(object=hnscc_primary,DB=lrDB,celltype="celltype",
#'                                   method="expression_product",
#'                                   receiver_celltypes=receiver_celltypes,
#'                                   sender_celltypes=sender_celltypes,product_top=0.5)
#' interaction3 <- find_interactions(object=hnscc_primary,DB=lrDB,celltype="celltype",
#'                                   method="differential_combinations",
#'                                   receiver_celltypes=receiver_celltypes,
#'                                   sender_celltypes=sender_celltypes,
#'                                   marker_filter=TRUE,avg_log2FC=1,padj=0.05)
#'
#' @export
#'
find_interactions <- function(object,DB,
                              celltype="celltype",
                              method=c("expression_threshold", "expression_product", "differential_combinations"),
                              threshold=0,product_top=0.5,receiver_celltypes="ALL",sender_celltypes="ALL",
                              marker_filter=FALSE,avg_log2FC=1,padj=0.05){

  method <- match.arg(method)
  mtx <- as.matrix(object@assays$RNA@data) #object profile
  CT <- object@meta.data[,celltype] #object$celltype
  names(CT) <- colnames(object)
  AI <- object@active.ident

  # step1,如果存在多亚基，计算complex的geometric.mean,结果保存在mtx profile
  LRDB <- DB[,c("ligand","receptor")]
  LRDB <- LRDB[!duplicated(paste0(LRDB$ligand,LRDB$receptor)),]
  gene <- unique(c(LRDB$ligand,LRDB$receptor))
  gene_name_list <- strsplit(gene,"_") #考虑多亚基
  gene2 <- unlist(gene_name_list)

  if(identical(gene,gene2)==F){
    mtx_complex <- complex(mtx,gene_name_list)
    object=CreateSeuratObject(counts = mtx_complex, project = "complex")
    object$celltype <- CT
    object@active.ident <- AI
    gene <- gsub("_","--",gene)
    mtx <- mtx_complex
  }
  #####
  if (identical(receiver_celltypes,"ALL")) {
    receiver_celltypes=unique(object$celltype)
  } else if(!all(receiver_celltypes %in% object$celltype)){
    stop("receiver_celltypes must in object's celltype")
  }
  if (identical(sender_celltypes,"ALL")) {
    sender_celltypes=unique(object$celltype)
  } else if(!all(sender_celltypes%in% object$celltype)){
    stop("sender_celltypes must in object's celltype")
  }
  cell.named.types=CT[CT %in% (c(receiver_celltypes,sender_celltypes))]
  # step2, 过滤配受体基因,要求至少在一个细胞类型的10%细胞中表达,得到profile
  sub_object=object[gene,]
  sub_object@meta.data=object@meta.data 
  Idents(sub_object)=Idents(object)
  all_expressed_genes = get_expressed_genes(sub_object, ident = c(receiver_celltypes,sender_celltypes),pct = 0.1)
  expressed_genes = all_expressed_genes %>% unlist() %>% unique()
  mtx_lr <- mtx[rownames(mtx) %in% expressed_genes,]
  LRDB$ligand <- gsub("_","--",LRDB$ligand)
  LRDB$receptor <- gsub("_","--",LRDB$receptor)
  ## 对配体受体库进行过滤
  LR=LRDB[LRDB$ligand %in% expressed_genes & LRDB$receptor %in% expressed_genes,]

  #step3,三种方法find interactions
  if (method == "expression_threshold") {   # 1.表达阈值法筛选LR Pairs
    interaction <- interact_threshold(mtx_lr,cell.named.types,LR,sender_celltypes,receiver_celltypes,threshold)
  }
  else if (method == "expression_product") {   # 2.表达乘积法筛选LR Pairs
    interaction <- interact_product(mtx_lr,cell.named.types,LR,sender_celltypes,receiver_celltypes,product_top)
  }
  else if (method == "differential_combinations") {  # 3, 差异组合法匹配LR Pairs
    interaction <- interact_combination(sub_object,LR,receiver_celltypes=receiver_celltypes,sender_celltypes=sender_celltypes,
                                        marker_filter,avg_log2FC,padj)
  }
  return(list(interaction=interaction,object=object))
}
#####
