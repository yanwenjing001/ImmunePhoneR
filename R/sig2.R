#' 利用permutation计算配体受体对的pval,找出有显著意义的配体受体对
#'
#' @param object1 Seurat object1
#' @param object2 Seurat object2
#' @param interaction1 data.frame,必须包含ligand,receptor,from,to
#' @param interaction2 data.frame,必须包含ligand,receptor,from,to
#' @param pval number,默认值为0.05
#' @param sample1 name of sample1
#' @param sample2 name of sample2
#' @return list,pvalue,有显著意义的配体受体对;
#'              significant,data.frame,行是ligand基因名,列是receptor基因名;
#'              receptor_logfc,receptor基因在特定细胞类型中的logfc值;
#'              ligand_logfc,ligand基因在特定细胞类型中的logfc值
#'
#' @export
#'
sig2 <- function(object1,object2,interaction1,interaction2,pval=0.05,sample1="sample1",sample2="sample2"){

  ##ligand,receptor data.frame
  interaction=rbind(interaction1,interaction2)
  lt <- as.data.frame(matrix(0,ncol=length(unique(interaction$receptor)),nrow=length(unique(interaction$ligand))))
  colnames(lt) <- unique(interaction$receptor)
  rownames(lt) <- unique(interaction$ligand)
  for(i in 1:nrow(interaction)){
    lt[interaction$ligand[i],interaction$receptor[i]]=lt[interaction$ligand[i],interaction$receptor[i]]+1
  }

  #step1：求ligand和receptor基因在每个细胞类型中的均值
  mtx1 <- as.matrix(object1@assays$RNA@data)
  mtx_interaction1 <- mtx1[unique(c(interaction$ligand,interaction$receptor)),]
  cell.named.types=object1$celltype
  mean_ct <- mean_celltype(mtx_interaction1,cell.named.types)

  receiver=colnames(mean_ct) %in% unique(interaction$to)
  sender=colnames(mean_ct) %in% unique(interaction$from)
  #step2:求每一对from-to pairs在sender和receiver细胞类型之间的矩阵乘积
  exp_r <- mean_ct[rownames(mean_ct) %in% interaction$receptor,][,receiver,drop=FALSE]
  exp_l <- mean_ct[rownames(mean_ct) %in% interaction$ligand,][,sender,drop=FALSE]

  lt.significant <- lt
  pvalue1 <- data.frame()
  for(i in 1:nrow(interaction)){
    print(i)
    ligand <- interaction$ligand[i]
    receptor <- interaction$receptor[i]
    from <- interaction$from[i]
    to <- interaction$to[i]
    if(from %in% colnames(exp_l) & to %in% colnames(exp_r)){
      l <- exp_l[ligand,,drop=FALSE][,from,drop=FALSE]
      r <- exp_r[receptor,,drop=FALSE][,to,drop=FALSE]
      prod <- production(l,r)
    }else{
      prod <- data.frame(0)
      colnames(prod) <- to
      rownames(prod) <- from
    }
    prod.pvalue <- prod
    for(send in from){
      for(rec in to){
        if(from %in% colnames(exp_l) & to %in% colnames(exp_r)){
          #取sender和receiver的细胞数目
          seu_sender <- subset(object1,idents = send)
          seu_receiver <- subset(object1,idents = rec)
          n_sender <- ncol(seu_sender@assays$RNA@data)
          n_receiver <- ncol(seu_receiver@assays$RNA@data)
          value_to_test <- as.numeric(prod)
          vector_for_test <- permutation(object1,ligand,receptor,n_sender,n_receiver)
          p.value <- sum(vector_for_test > value_to_test)/length(vector_for_test)
          prod.pvalue[send,rec] <- p.value
        } else {
          prod.pvalue[send,rec] <- 1
        }
        prod.pvalue <- as.matrix(prod.pvalue)
      }
    }
    pvalue <- melt(prod.pvalue)
    colnames(pvalue) <- c("sender","receiver","pval")
    pvalue$ligand <- rownames(l)
    pvalue$receptor <- rownames(r)

    pvalue1 <- rbind(pvalue1,pvalue)
  }

  pvalue_order <- pvalue1[order(pvalue1$pval),]
  pvalue_min <- pvalue_order[!duplicated(pvalue_order[,c("ligand","receptor")]),]

  for(i in 1:nrow(pvalue_min)){
    #求一对from-to pairs在sender和receiver细胞类型之间的矩阵乘积
    ligand2 <- pvalue_min$ligand[i]
    receptor2 <- pvalue_min$receptor[i]
    lt.significant[ligand2,receptor2] <- sum(pvalue_min[c(pvalue_min$ligand %in% ligand2 & pvalue_min$receptor %in% receptor2),]$pval < pval)
  }

  pvalue2 <- cbind(interaction,pval=pvalue1[,"pval"])
  pvalue3 <- pvalue2[pvalue2$pval <= pval,]

  sig1 <- list(pvalue1=pvalue2,significant1=lt.significant)

  #####2
  #step2：求ligand和receptor基因在每个细胞类型中的均值
  mtx2 <- as.matrix(object2@assays$RNA@data)
  mtx_interaction2 <- mtx2[unique(c(interaction$ligand,interaction$receptor)),]
  cell.named.types=object2$celltype
  mean_ct <- mean_celltype(mtx_interaction2,cell.named.types)

  receiver=colnames(mean_ct) %in% unique(interaction$to)
  sender=colnames(mean_ct) %in% unique(interaction$from)
  #step2:求每一对from-to pairs在sender和receiver细胞类型之间的矩阵乘积
  exp_r <- mean_ct[rownames(mean_ct) %in% interaction$receptor,][,receiver,drop=FALSE]
  exp_l <- mean_ct[rownames(mean_ct) %in% interaction$ligand,][,sender,drop=FALSE]

  lt.significant <- lt
  pvalue1 <- data.frame()
  for(i in 1:nrow(interaction)){
    print(i)
    ligand <- interaction$ligand[i]
    receptor <- interaction$receptor[i]
    from <- interaction$from[i]
    to <- interaction$to[i]
    if(from %in% colnames(exp_l) & to %in% colnames(exp_r)){
      l <- exp_l[ligand,,drop=FALSE][,from,drop=FALSE]
      r <- exp_r[receptor,,drop=FALSE][,to,drop=FALSE]
      prod <- production(l,r)
    }else{
      prod <- data.frame(0)
      colnames(prod) <- to
      rownames(prod) <- from
    }
    prod.pvalue <- prod
    for(send in from){
      for(rec in to){
        if(from %in% colnames(exp_l) & to %in% colnames(exp_r)){
          #取sender和receiver的细胞数目
          seu_sender <- subset(object2,idents = send)
          seu_receiver <- subset(object2,idents = rec)
          n_sender <- ncol(seu_sender@assays$RNA@data)
          n_receiver <- ncol(seu_receiver@assays$RNA@data)
          value_to_test <- as.numeric(prod)
          vector_for_test <- permutation(object2,ligand,receptor,n_sender,n_receiver)
          p.value <- sum(vector_for_test > value_to_test)/length(vector_for_test)
          prod.pvalue[send,rec] <- p.value
        } else {
          prod.pvalue[send,rec] <- 1
        }
        prod.pvalue <- as.matrix(prod.pvalue)
      }
    }
    pvalue <- melt(prod.pvalue)
    colnames(pvalue) <- c("sender","receiver","pval")
    pvalue$ligand <- rownames(l)
    pvalue$receptor <- rownames(r)

    pvalue1 <- rbind(pvalue1,pvalue)
  }

  pvalue_order <- pvalue1[order(pvalue1$pval),]
  pvalue_min <- pvalue_order[!duplicated(pvalue_order[,c("ligand","receptor")]),]

  for(i in 1:nrow(pvalue_min)){
    #求一对from-to pairs在sender和receiver细胞类型之间的矩阵乘积
    ligand2 <- pvalue_min$ligand[i]
    receptor2 <- pvalue_min$receptor[i]
    lt.significant[ligand2,receptor2] <- sum(pvalue_min[c(pvalue_min$ligand %in% ligand2 & pvalue_min$receptor %in% receptor2),]$pval < pval)
  }

  pvalue2 <- cbind(interaction,pval=pvalue1[,"pval"])
  pvalue3 <- pvalue2[pvalue2$pval <= pval,]

  sig2 <- list(pvalue2=pvalue2,significant2=lt.significant)

  for_map <- logfc2(object1,object2,sig1,sig2,sample1,sample2)
  result2 <- c(list(pvaule1=sig1$pvalue1,pvalue2=sig2$pvalue2),for_map)
}
###############
#' 计算ligand和receptor在指定细胞类型中的logfc
#'
#' @param object1 Seurat object1
#' @param object2 Seurat object2
#' @param sig1 significant result of object1
#' @param sig2 significant result of object2
#' @param sample1 name of sample1
#' @param sample2 name of sample2
#' @return list,receptor_logfc,receptor基因在特定细胞类型中的logfc值;
#'              ligand_logfc,ligand基因在特定细胞类型中的logfc值
#' @export
#'
logfc2 <- function(object1,object2,sig1,sig2,sample1="sample1",sample2="sample2"){

  interaction1 <- sig1$pvalue1
  middle1 <- sig1$significant1
  lr_net1 <- interaction1[,c("ligand","receptor","from","to")]

  interaction2 <- sig2$pvalue2
  middle2 <- sig2$significant2
  lr_net2 <- interaction2[,c("ligand","receptor","from","to")]

  genes1=unique(c(lr_net1$ligand,lr_net1$receptor))
  seurat1=object1[genes1,]
  allmarkers1 = FindAllMarkers(seurat1,min.pct = 0, logfc.threshold = 0,return.thresh = 1)
  allmarkers1$ct=paste0(allmarkers1$cluster,"_",sample1)

  receptor_markers1 <- allmarkers1[allmarkers1$cluster %in% unique(lr_net1$to),]
  receptor_markers1 = receptor_markers1[receptor_markers1$gene %in% lr_net1$receptor,]
  receptor_markers1 = receptor_markers1[,c("gene","ct","avg_log2FC")]
  vis_r1 = dcast(receptor_markers1,gene~ct)
  vis_r1[is.na(vis_r1)] <- 0
  rownames(vis_r1) <- vis_r1$gene
  vis_r_map1 <- vis_r1[,-1,drop=FALSE]

  ligand_markers1 <- allmarkers1[allmarkers1$cluster %in% unique(lr_net1$from),]
  ligand_markers1 = ligand_markers1[ligand_markers1$gene %in% lr_net1$ligand,]
  ligand_markers1 = ligand_markers1[,c("gene","ct","avg_log2FC")]
  vis_l1 = dcast(ligand_markers1,gene~ct)
  vis_l1[is.na(vis_l1)] <- 0
  rownames(vis_l1) <- vis_l1$gene
  vis_l_map1 <- vis_l1[,-1,drop=FALSE]

  genes2=unique(c(lr_net2$ligand,lr_net2$receptor))
  seurat2=object2[genes2,]
  allmarkers2 = FindAllMarkers(seurat2,min.pct = 0, logfc.threshold = 0,return.thresh = 1)
  allmarkers2$ct=paste0(allmarkers2$cluster,"_",sample2)

  receptor_markers2 <- allmarkers2[allmarkers2$cluster %in% unique(lr_net2$to),]
  receptor_markers2 = receptor_markers2[receptor_markers2$gene %in% lr_net2$receptor,]
  receptor_markers2 = receptor_markers2[,c("gene","ct","avg_log2FC")]
  vis_r2 = dcast(receptor_markers2,gene~ct)
  vis_r2[is.na(vis_r2)] <- 0
  rownames(vis_r2) <- vis_r2$gene
  vis_r_map2 <- vis_r2[,-1,drop=FALSE]

  ligand_markers2 <- allmarkers2[allmarkers2$cluster %in% unique(lr_net2$from),]
  ligand_markers2 = ligand_markers2[ligand_markers2$gene %in% lr_net2$ligand,]
  ligand_markers2 = ligand_markers2[,c("gene","ct","avg_log2FC")]
  vis_l2 = dcast(ligand_markers2,gene~ct)
  vis_l2[is.na(vis_l2)] <- 0
  rownames(vis_l2) <- vis_l2$gene
  vis_l_map2 <- vis_l2[,-1,drop=FALSE]

  vis_r_map=cbind(vis_r_map1,vis_r_map2)
  vis_l_map=cbind(vis_l_map1,vis_l_map2)

  result <- list(receptor_logfc=vis_r_map,ligand_logfc=vis_l_map,up.middle.mat=middle1,down.middle.mat=middle2)
  return(result)
}
################
