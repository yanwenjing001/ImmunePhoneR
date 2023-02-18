 
# 通过配体受体推断细胞类型之间的相互作用

在此说明文档中，您可以利用本工具进行推断、分析和可视化细胞-细胞之间的相互作用。本工具需要配体、受体相互作用的先验知识，对细胞与细胞之间的通信概率进行建模。在推断细胞间通信网络的基础上，为进一步的数据探索、分析和可视化提供了功能。 
具体来说，根据先验的配体受体库，可以预测1)使用单细胞基因表达数据，推断来自一个细胞群体(“发送者”)的哪些配体可能通过相互作用影响另一个细胞群体(“接收者”)中的哪些受体表达，
2)通过cellmarker来推断两个特定细胞类型之间的互作。   

本文档将详细指导您完成所有这些步骤。 作为单细胞表达数据的例子，我们将使用Puram等人的数据来探索头颈部鳞状细胞癌(HNSCC)肿瘤微环境中的细胞间通讯(见Puram et al. 2017)。 更具体地说，我们将研究肿瘤细胞与其他细胞类型之间的互作关系。 

## 1.加载所需的包
``` r
library(tidyverse)
library(reshape2)
library(grid)
library(gtable)
library(Seurat) 
library(iphonedev2) 

```
## 2.数据输入、选择先验的配体受体库

#### ①需要的用户输入:作为Seurat对象的单细胞数据。  

``` r
hnscc_primary <- readRDS("hnscc_primary.rds")
hnscc_primary@meta.data %>% head()
##                       orig.ident nCount_RNA nFeature_RNA      celltype
## HN28_P6_G05_S173_comb       HN28   8215.891         2723 myofibroblast
## HN25_P5_E10_S58_comb        HN25  14149.138         5064 myofibroblast
## HN26_P6_B06_S18_comb        HN26   6383.145         2595 myofibroblast
## HN25_P5_C05_S29_comb        HN25  14954.354         6423           CAF
## HN26_P5_H01_S85_comb        HN26  13097.883         5486       primary
## HN26_P5_E06_S54_comb        HN26  16132.502         7038       primary
hnscc_primary$celltype %>% table()
## CAF   Endothelial    Macrophage myofibroblast       primary        T cell 
## 218            27            63           304           691           329 

```

#### ②我们整合了八个配体受体库资源，用户可以自由选择其中的配体受体库

![](https://github.com/yanwenjing001/Iphone/blob/main/vignettes/database1.png)

![](https://github.com/yanwenjing001/Iphone/blob/main/vignettes/database2.png)

#### 下面以cellchatDB为例：

``` r
lrDB <- read.delim("C:/Users/YWJ/Desktop/库/CellchatDB/cellchatDB_human.txt", header=TRUE)

```

## 3.使用三种算法推断指定细胞类型间的互作

### 方法1: 用表达阈值法推断肿瘤与TME细胞的相互作用
``` r
interaction1_1 <- find_interactions(object=hnscc,DB=lrDB,method="expression_threshold",
                                  receiver_celltypes=c("T cell","CAF","Endothelial","Macrophage","myofibroblast"),
                                  sender_celltypes="primary",
                                  threshold=3)

interaction1_2 <- find_interactions(object=hnscc,DB=lrDB,method="expression_threshold",
                                  receiver_celltypes="primary",
                                  sender_celltypes=c("T cell","CAF","Endothelial","Macrophage","myofibroblast"),
                                  threshold=3)
```    

### 方法2: 用表达乘积法推断肿瘤与TME细胞的相互作用
``` r
interaction2_1 <- find_interactions(object=hnscc,DB=lrDB,method="expression_product",
                                  receiver_celltypes=c("T cell","CAF","Endothelial","Macrophage","myofibroblast"),
                                  sender_celltypes="primary")

interaction2_2 <- find_interactions(object=hnscc,DB=lrDB,method="expression_product",
                                  receiver_celltypes="primary",
                                  sender_celltypes=c("T cell","CAF","Endothelial","Macrophage","myofibroblast"))
```

### 方法3: 用差异组合法推断肿瘤与TME细胞的相互作用
``` r
interaction3_1 <- find_interactions(object=hnscc,DB=lrDB,method="differential_combinations",
                                  receiver_celltypes=c("T cell","CAF","Endothelial","Macrophage","myofibroblast"),
                                  sender_celltypes="primary",
                                  marker_filter=TRUE,avg_log2FC=1,padj=0.05)
                                  
interaction3_2 <- find_interactions(object=hnscc,DB=lrDB,method="differential_combinations",
                                  receiver_celltypes="primary",
                                  sender_celltypes=c("T cell","CAF","Endothelial","Macrophage","myofibroblast"),
                                  marker_filter=TRUE,avg_log2FC=1,padj=0.05)
```

## 4.找出有意义的配体受体对及可视化

#### 表达阈值法
``` r
siglr1 <- sig(object=hnscc_primary,DB=lrDB,lr=interaction1_1,pval=0.05)
lrmap1 <- map(res=siglr1)  
siglr2 <- sig(object=hnscc_primary,DB=lrDB,lr=interaction1_2,pval=0.05)
lrmap2 <- map(res=siglr2)  
```
![](https://github.com/yanwenjing001/Iphone/blob/main/vignettes/expression_threshold.png)

#### 表达乘积法
``` r
siglr1 <- sig(object=hnscc_primary,DB=lrDB,lr=interaction2_1,pval=0.05)
lrmap1 <- map(res=siglr1)  
siglr2 <- sig(object=hnscc_primary,DB=lrDB,lr=interaction2_2,pval=0.05)
lrmap2 <- map(res=siglr2)  
```
![](https://github.com/yanwenjing001/Iphone/blob/main/vignettes/expression_product.png)

#### 差异组合法
``` r
siglr1 <- sig(object=hnscc_primary,DB=lrDB,lr=interaction3_1,pval=0.05)
lrmap1 <- map(res=siglr1)  
siglr2 <- sig(object=hnscc_primary,DB=lrDB,lr=interaction3_2,pval=0.05)
lrmap2 <- map(res=siglr2)  
```
![](https://github.com/yanwenjing001/Iphone/blob/main/vignettes/differential_combination.png)

## 5.作者
闫文婧 9874703396@qq.com
徐嘉潞 xujialu13009484417@163.com
吕德康 dekanglv@126.com
         
         
         

















