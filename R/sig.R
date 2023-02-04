#' @title 利用permutation计算配体受体对的pval,找出有显著意义的配体受体对
#'
#' @param object Seurat object
#' @param interaction data.frame,必须包含ligand,receptor,from,to
#' @param pval number,默认值为0.05
#' @param method "product" or "sum"
#' @return list,pvalue,有显著意义的配体受体对;
#'              significant,data.frame,行是ligand基因名,列是receptor基因名;
#'              receptor_logfc,receptor基因在特定细胞类型中的logfc值;
#'              ligand_logfc,ligand基因在特定细胞类型中的logfc值
#' @import reshape2
#' @importFrom stats quantile
#' @examples
#' load(url("https://figshare.com/ndownloader/files/37663929"))
#' interaction <- data.frame(ligand = c("MDK","LAMB3","ADM"),
#'                           receptor = c("NCL","CD44","CALCRL"),
#'                           from = "primary",to = "Macrophage")
#' sig(object=hnscc_primary,interaction,pval=0.05)
#'
#' @export
#'
sig <- function(object,interaction,pval=0.05,method="product"){

  ##ligand,receptor data.frame
  lt <- as.data.frame(matrix(0,ncol=length(unique(interaction$receptor)),nrow=length(unique(interaction$ligand))))
  colnames(lt) <- unique(interaction$receptor)
  rownames(lt) <- unique(interaction$ligand)
  for(i in 1:nrow(interaction)){
    lt[interaction$ligand[i],interaction$receptor[i]]=lt[interaction$ligand[i],interaction$receptor[i]]+1
  }
  #step1：求ligand和receptor基因在每个细胞类型中的均值
  mtx <- as.matrix(object@assays$RNA@data)
  mtx_interaction <- mtx[unique(c(interaction$ligand,interaction$receptor)),]
  cell.named.types=object$celltype
  mean_ct <- mean_celltype(mtx_interaction,cell.named.types)

  receiver=colnames(mean_ct) %in% unique(interaction$to)
  sender=colnames(mean_ct) %in% unique(interaction$from)
  #step2:求每一对interaction pairs在sender和receiver细胞类型之间的矩阵乘积
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
          #取与sender和receiver的细胞数目相同的随机细胞,做permutation求interaction的p value
          seu_sender <- subset(object,idents = send)
          seu_receiver <- subset(object,idents = rec)
          n_sender <- ncol(seu_sender@assays$RNA@data)
          n_receiver <- ncol(seu_receiver@assays$RNA@data)
          value_to_test <- as.numeric(prod)
          vector_for_test <- permutation(object,ligand,receptor,n_sender,n_receiver,method="product")
          p.value <- sum(vector_for_test > value_to_test)/length(vector_for_test)
          prod.pvalue[send,rec] <- p.value
        } else {
          prod.pvalue[send,rec] <- 1
          }
          prod.pvalue <- as.matrix(prod.pvalue)
      }
    }
    pvalue <- reshape2::melt(prod.pvalue)
    colnames(pvalue) <- c("sender","receiver","pval")
    pvalue$ligand <- rownames(l)
    pvalue$receptor <- rownames(r)

    pvalue1 <- rbind(pvalue1,pvalue)
  }

  pvalue_order <- pvalue1[order(pvalue1$pval),]
  pvalue_uni <- pvalue_order[!duplicated(pvalue_order[,c("ligand","receptor")]),]

  for(i in 1:nrow(pvalue_uni)){
    ligand2 <- pvalue_uni$ligand[i]
    receptor2 <- pvalue_uni$receptor[i]
    lt.significant[ligand2,receptor2] <- sum(pvalue_uni[c(pvalue_uni$ligand %in% ligand2 & pvalue_uni$receptor %in% receptor2),]$pval < pval)
  }

  pvalue_all <- cbind(interaction,pval=pvalue1[,"pval"])
  pvalue_sig <- pvalue_all[pvalue_all$pval <= pval,]

  sig_interaction <- list(pvalue=pvalue_sig,significant=lt.significant)
  for_map <- logfc(object,siglr=sig_interaction)
  result <- c(sig_interaction,for_map)
  return(result)
}


#######
#' @title 针对两种细胞类型中的某一ligand-receptor pairs，通过随机抽样获得100个ligand-receptor表达均值的乘积，
#' 用于后面permutation单边检验
#'
#' @param seurat Seurat object
#' @param l ligand基因名
#' @param r receptor基因名
#' @param n_sender sender_celltype的细胞数目
#' @param n_receiver receiver_celltype的细胞数目
#' @param iteration 迭代次数，默认值为100
#' @param method "product" or "sum"
#'
#' @return vector,其中的每一个值都是permutation得到的两个mean值的乘积
#' @examples
#' load(url("https://figshare.com/ndownloader/files/37663929"))
#' l <- "TIGIT"
#' r <- "PVR"
#' n_sender <- 50
#' n_receiver <- 80
#' permutation(seurat,l,r,n_sender,n_receiver,iteration=100,method="product")
#'
#' @export
#'
permutation <- function(seurat,l,r,n_sender,n_receiver,iteration=100,method=c("product","sum")){
  vector <- c()
  for(i in 1:iteration){
    l_mean <- mean(sample(seurat@assays$RNA@data[l,],n_sender, replace = FALSE))
    t_mean <- mean(sample(seurat@assays$RNA@data[r,],n_receiver, replace = FALSE))
    if(method == "product"){
    vector <- c(vector,l_mean*t_mean)
    }else{
    vector <- c(vector,l_mean+t_mean)
    }
  }
  return(vector)
}

###############################################
#' @title 计算ligand和receptor在指定细胞类型中的logfc
#'
#' @param object Seurat object
#' @param siglr data.frame,必须包含ligand,receptor,from,to，pval
#' @return list,receptor_logfc,receptor基因在特定细胞类型中的logfc值;
#'              ligand_logfc,ligand基因在特定细胞类型中的logfc值
#' @examples
#' load(url("https://figshare.com/ndownloader/files/37663929"))
#' interaction <- data.frame(ligand = c("MDK","LAMB3"),
#'                           receptor = c("NCL","CD44"),
#'                           from = "primary",to = "Macrophage",
#'                           pval = c(0,0))
#' logfc <- function(object=seurat,siglr=interaction)
#'
#' @export
#'
logfc <- function(object,siglr){

  interaction <- siglr$pvalue
  lr_net <- interaction[,c("ligand","receptor","from","to")]

  genes=unique(c(lr_net$ligand,lr_net$receptor))
  seurat=object[genes,]
  allmarkers = FindAllMarkers(seurat,min.pct = 0, logfc.threshold = 0,return.thresh = 1)

  receptor_markers <- allmarkers[allmarkers$cluster %in% unique(lr_net$to),]
  receptor_markers2 = receptor_markers[receptor_markers$gene %in% lr_net$receptor,]
  receptor_markers2 = receptor_markers2[,c("gene","cluster","avg_log2FC")]
  vis_r = dcast(receptor_markers2,gene~cluster)
  vis_r[is.na(vis_r)] <- 0
  rownames(vis_r) <- vis_r$gene
  vis_r_map <- vis_r[,-1,drop=FALSE]

  ligand_markers <- allmarkers[allmarkers$cluster %in% unique(lr_net$from),]
  ligand_markers2 = ligand_markers[ligand_markers$gene %in% lr_net$ligand,]
  ligand_markers2 = ligand_markers2[,c("gene","cluster","avg_log2FC")]
  vis_l = reshape2::dcast(ligand_markers2,gene~cluster)
  vis_l[is.na(vis_l)] <- 0
  rownames(vis_l) <- vis_l$gene
  vis_l_map <- vis_l[,-1,drop=FALSE]
  result <- list(receptor_logfc=vis_r_map,ligand_logfc=vis_l_map)
  return(result)
}
