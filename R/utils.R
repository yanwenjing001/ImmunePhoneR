#' @title Compute the geometric mean
#'
#' @param x	 a numeric vector
#' @param na.rm	 whether remove na
#' @return number,geometric mean
#' @examples
#' x <- 1:10
#' geometric.mean(x)
#'
#' @export
geometric.mean <- function (x, na.rm = TRUE)
{
  if (is.null(nrow(x))) {
    exp(mean(log(x), na.rm = na.rm))
  }
  else {
    exp(apply(log(x), 2, mean, na.rm = na.rm))
  }
}
######
#' @title 计算复合体的表达值
#'
#' @description 用复合体的多亚基的几何平均数来表示复合体表达值
#'
#' @details 为了推断细胞间的相互作用，考虑到配体受体的复合体，需要计算复合体的表达值，便于筛选出特异的相互作用，
#' 将复合体表达值附加到原表达矩阵之后。
#'
#' @param mtx 表达矩阵，行代表基因，列代表细胞
#' @param gene_name_list a list，其中每个元素是亚基的集合
#' @return mtx 将复合体表达值附加到原表达矩阵后的表达矩阵
#' @examples
#' mtx <- matrix(abs(rnorm(200)), 20, 10)
#' colnames(mtx) = paste("Cell", 1:10, sep = "")
#' rownames(mtx) = paste("Gene", 1:20, sep = "")
#' gene_name_list <- list("Gene1",c("Gene2","Gene3"),"Gene4",c("Gene5","Gene6"))
#' complex(mtx,gene_name_list)
#'
#' @export
complex <- function(mtx,gene_name_list){
  for(j in 1:length(gene_name_list)){
    print(j)
    genes=gene_name_list[[j]]
    sub_genes=intersect(genes,rownames(mtx))
    if(length(genes)> 1){
      if(length(sub_genes) > 1){
        mtx_complex=mtx[sub_genes,]
        mtx_complex=geometric.mean(mtx_complex,na.rm = FALSE)
        mtx_complex2=matrix(mtx_complex,nrow=1)
      }else if(length(sub_genes) == 1){
        mtx_complex2=matrix(mtx[sub_genes,],nrow=1)
      }else {next}
      rownames(mtx_complex2)=paste(genes,collapse = "--")
      colnames(mtx_complex2)=colnames(mtx_complex)
      mtx=rbind(mtx,mtx_complex2)
    }
  }
  return(mtx)
}
#####
#' @title Determine expressed genes of a cell type from a Seurat object single-cell RNA seq dataset
#'
#' @description \code{get_expressed_genes} Return the genes that are expressed in a given cell cluster based on
#' the fraction of cells in that cluster that should express the gene.
#' @details get_expressed_genes(ident, seurat_obj, pct = 0.10, assay_oi = NULL)
#'
#' @param ident Name of cluster identity/identities of cells
#' @param seurat_obj Single-cell expression dataset as Seurat object
#' @param pct We consider genes expressed if they are expressed in at least a specific fraction of cells of a cluster.
#' This number indicates this fraction. Default: 0.10. Choice of this parameter is important and depends largely on
#' the used sequencing platform. We recommend to require a lower fraction (like the default 0.10)
#' for 10X data than for e.g. Smart-seq2 data.
#'
#' @return A character vector with the gene symbols of the expressed genes
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' seurat_obj <- pbmc
#' get_expressed_genes(ident = "CD8 T", seurat_obj = seuratObj, pct = 0.10)
#'
#' @export
#'
get_expressed_genes <- function (seurat_obj, ident, pct = 0.1) {
  requireNamespace("Seurat")
  requireNamespace("dplyr")

  if (sum(ident %in% unique(Idents(seurat_obj))) != length(ident)) {
    stop("One or more provided cell clusters is not part of the 'Idents' of your Seurat object")
  }
  cells_oi = Idents(seurat_obj) %>% .[Idents(seurat_obj) %in% ident] %>% names()
  ident.mtx <- seurat_obj@assays$RNA@data[,cells_oi]
  number.cells=ncol(ident.mtx)
  ifpercentage <- function(x){
    sum(x>0)/number.cells >= pct
  }
  genes <- rownames(ident.mtx)[apply(ident.mtx,1,ifpercentage)]
  return(genes)
}
#####矩阵乘积
#' @title 计算矩阵乘积
#' @param vector1 某一ligand在sender_celltype的表达均值
#' @param vector2 某一receptor在receiver的表达均值
#' @return vector1和vector1的乘积
#'
#' @export
production<-function(vector1,vector2){
  #vector1：某一ligand在sender_celltype的表达均值
  #vector2：某一receptor在receiver的表达均值
  do.call("rbind",lapply(vector1, function(x){x*vector2}))
}
######
#' @title 计算细胞类型的基因表达均值
#'
#' @param mtx 表达矩阵，行代表基因，列代表细胞
#' @param cell.named.types 由细胞类型组成的向量，向量名为细胞名
#' @return mean_lr 各细胞类型的基因表达均值
#' @examples
#' mtx <- matrix(abs(rnorm(240)), 20, 12)
#' colnames(mtx) = paste("Cell", 1:12, sep = "")
#' rownames(mtx) = paste("Gene", 1:20, sep = "")
#' cell.named.types = rep(c("Bcell","CAF","Tcell"),4)
#' names(cell.named.types) = colnames(mtx)
#' mean_celltype(mtx,cell.named.types)
#'
#' @export
mean_celltype <- function(mtx,cell.named.types){
  res=data.frame()
  for(i in unique(cell.named.types)){
    print(i)
    cells=names(cell.named.types)[cell.named.types==i]
    mean_lr=apply(mtx[,colnames(mtx) %in% cells],1,mean)
    res=rbind(res,mean_lr)
  }
  mean_lr=as.data.frame(t(res))
  colnames(mean_lr)=unique(cell.named.types)
  rownames(mean_lr)=rownames(mtx)
  return(mean_lr)
}
#####
#' @title 取出sender和receiver对应的配体受体对的表达值
#'
#' @param mean_lr 各细胞类型的基因表达值
#' @param LR 先验的配体受体对
#' @param sender_celltypes character,发送ligand的细胞类型
#' @param receiver_celltypes character,产生receptor的细胞类型
#' @return list_mean_lr a list,sender的配体表达值，receiver的受体表达值，顺序与LR一致
#' @examples
#' mtx <- matrix(abs(rnorm(60)), 20, 3)
#' colnames(mtx) = c("Bcell","CAF","Tcell")
#' rownames(mtx) = paste("Gene", 1:20, sep = "")
#' LR = data.frame(ligand=paste("Gene", 1:10, sep = ""),receptor=paste("Gene", 11:20, sep = ""))
#' sender_celltypes = "Bcell"
#' receiver_celltypes = c("CAF","Tcell")
#' list_mean_lr(mean_lr,LR,sender_celltypes,receiver_celltypes)
#'
#' @export
#'
list_mean_lr <- function(mean_lr,LR,sender_celltypes,receiver_celltypes){
  mean_ligand=list()
  mean_receptor=list()
  for(col in colnames(mean_lr)){
    col_expr=mean_lr[,col]
    names(col_expr)=rownames(mean_lr)
    mean_ligand[[col]]=col_expr[LR$ligand]
    mean_receptor[[col]]=col_expr[LR$receptor]
  }
  mean.ligand=do.call(cbind,mean_ligand)
  mean.receptor=do.call(cbind,mean_receptor)
  LR_mean<- list(ligand=mean.ligand,receptor=mean.receptor)

  mean.ligand=LR_mean[["ligand"]]
  mean.ligand=mean.ligand[,sender_celltypes,drop=FALSE]
  mean.receptor=LR_mean[["receptor"]]
  mean.receptor=mean.receptor[,receiver_celltypes,drop=FALSE]
  return(list_mean=list(ligand=mean.ligand,receptor=mean.receptor))
}
######
utils::globalVariables(
  c(".","ligand","receptor")
)





