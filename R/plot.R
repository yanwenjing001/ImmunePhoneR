#' @title Hierarchical cluster analysis
#'
#' @param mat matrix
#' @param row boolean values determining if rows should be clustered or hclust object
#' @param col boolean values determining if columns should be clustered or hclust object
#' @return matrix of which the order of row and column changed according to clustering
#' @examples
#' test = matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#' clustering(test)
#' @export
#'
clustering <- function(mat,row=TRUE,col=TRUE){
  if(row){
    hclust_mat = hclust(dist(mat),method = "ward.D2")
    order_mat_row = hclust_mat$labels[hclust_mat$order]
  }else{
    order_mat_row=rownames(mat)}
  if(col){
    hclust_mat = hclust(dist(t(mat)),method = "ward.D2")
    order_mat_col = hclust_mat$labels[hclust_mat$order]
  }else{
    order_mat_col=colnames(mat)}
  return(mat[order_mat_row,order_mat_col,drop=FALSE])
}
#####
#' @title 根据一组文本字符串长度，设置其水平或竖直时的宽度或高度
#'
#' @param texts vector文本
#' @param hw "height" or "width"
#' @param angle_col 角度
#' @param addtion 空白宽度或高度
#' @param fontsize 字体大小
#' @return 最长的字符串的高度或宽度
#' @import grid
#' @importFrom graphics strwidth
#' @export
#'
set_height_width <- function(texts,hw,angle_col=NULL,addtion=5,fontsize=7){
  #筛选出最上面的列名(left、bottom)中最长的字符串
  tw = strwidth(texts, units = 'in')
  longest.text = texts[which.max(tw)]
  #计算最长字符旋转后的高度，加上10bigpts后设置成最上面列名的高度
  gp = list(fontsize = fontsize)
  if(hw=="height"){
    if(is.null(angle_col)){angle_col=90}
    longest.text.grob = textGrob(longest.text, rot = angle_col, gp = do.call(gpar, gp))
    top.height.width = unit(1, units = "grobheight", data = longest.text.grob) + unit(addtion, "bigpts")
  }
  if(hw=="width"){
    if(is.null(angle_col)){angle_col=0}
    longest.text.grob = textGrob(longest.text, rot = angle_col, gp = do.call(gpar, gp))
    top.height.width = unit(1, units = "grobwidth", data = longest.text.grob) + unit(addtion, "bigpts")
  }
  return(top.height.width)
}
#####
#' @title 数据变换，把极值转换到特定范围内
#'
#' @param mat data.frame or matrix
#' @param range 转换的范围
#' @return 把极值转换到特定范围内的数据
#'
#' @export
#'
transformation=function(mat,range=c(-1,1)){
  mat=as.matrix(mat)
  group=matrix(as.numeric(cut(mat,c(-Inf,range,Inf))),ncol=ncol(mat),nrow=nrow(mat))
  mat[group==1]=range[1]
  mat[group==3]=range[2]
  return(mat)
}
#####
#' @title 绘制细胞互作热图
#'
#' @param res list,细胞互作分析结果，包含pvalue,significant,receptor_logfc,ligand_logfc
#' @param levels.mat character matrix的level
#' @param transformation 数据变换，把极值转换到特定范围内
#' @return 细胞互作热图
#'
#' @export
#'
map <- function(res,levels.mat=c("0","1"),transformation=T){

  significant_res <- res$significant

  #ligand,receptor聚类，如果只有一列，则只对行进行聚类
  if(ncol(res$ligand_logfc)==1){
    clust_l=clustering(res$ligand_logfc,row=T,col=F)
  } else{
    clust_l=clustering(res$ligand_logfc,row=T,col=T)
  }

  if(ncol(res$receptor_logfc)==1){
    clust_r=clustering(res$receptor_logfc,row=T,col=F)
  } else{
    clust_r=clustering(res$receptor_logfc,row=T,col=T)
  }
  #中间矩阵的顺序按聚类顺序调整
  left.mat=clust_l
  middle.mat=significant_res[rownames(clust_l),rownames(clust_r)]
  bottom.mat=t(clust_r)

  #数据变换，把极值转换到特定范围内
  if(transformation==T){
    left.mat=transformation(left.mat,range=c(-1,1))
    bottom.mat=transformation(bottom.mat,range=c(-1,1))
  }
  #
  left.mat=as.matrix(left.mat)
  bottom.mat=as.matrix(bottom.mat)
  middle.mat=as.matrix(middle.mat)
  mode(middle.mat)="character"

  lr2map(left.mat=left.mat,bottom.mat=bottom.mat,middle.mat=middle.mat,fontsize=7,levels.mat=c("0","1"))
}
######
#' @title 绘制热图
#'
#' @param left.mat left矩阵
#' @param middle.mat middle矩阵
#' @param bottom.mat bottom矩阵
#' @param left.colors left矩阵颜色
#' @param middle.colors middle矩阵颜色
#' @param bottom.colors bottom矩阵颜色
#' @param levels.mat character matrix的level
#' @param ... additional arguments passed on to image
#' @importFrom grDevices colorRampPalette rgb
#' @importFrom gtable gtable_add_grob gtable
#' @importFrom scales hue_pal
#' @return 热图
#'
#' @export
#'
lr2map=function(left.mat,middle.mat,bottom.mat,left.colors=NA,middle.colors=NA,bottom.colors=NA,levels.mat=NA,...){
  #给三个矩阵赋颜色
  if(is.na(left.colors)){
    rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
    r=c(157,171,185,199,216,219,228,234,237,232,227,219,211,186,164,142,120,105,88,70,60)
    g=c(56,90,113,133,156,169,187,203,221,222,225,224,226,212,199,188,175,163,149,139,128)
    b=c(75,88,101,113,124,144,165,185,206,214,224,231,236,232,231,225,217,210,206,197,187)
    coldf=data.frame(r,g,b)

    white=rgb(245, 246, 247, maxColorValue = 255)
    red=rgb(174, 61, 141, maxColorValue = 255)

    colours=rev(apply(coldf,1,rgb2hex))
  }
  if(is.na(middle.colors)){
    rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
    r=c(157,171,185,199,216,219,228,234,237,232,227,219,211,186,164,142,120,105,88,70,60)
    g=c(56,90,113,133,156,169,187,203,221,222,225,224,226,212,199,188,175,163,149,139,128)
    b=c(75,88,101,113,124,144,165,185,206,214,224,231,236,232,231,225,217,210,206,197,187)
    coldf=data.frame(r,g,b)

    white=rgb(245, 246, 247, maxColorValue = 255)
    red=rgb(174, 61, 141, maxColorValue = 255)

    colours=rev(apply(coldf,1,rgb2hex))
    middle.colors=c(white,red)
  }
  if(is.na(bottom.colors)){
    rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
    r=c(157,171,185,199,216,219,228,234,237,232,227,219,211,186,164,142,120,105,88,70,60)
    g=c(56,90,113,133,156,169,187,203,221,222,225,224,226,212,199,188,175,163,149,139,128)
    b=c(75,88,101,113,124,144,165,185,206,214,224,231,236,232,231,225,217,210,206,197,187)
    coldf=data.frame(r,g,b)

    white=rgb(245, 246, 247, maxColorValue = 255)
    red=rgb(174, 61, 141, maxColorValue = 255)
    colours=rev(apply(coldf,1,rgb2hex))
  }
  #如果矩阵只有一行或一列，颜色设置
  if(ncol(left.mat)==1){
    left.colors=colorRampPalette(c("white","#9D384B"))(256)[100:256]
  } else if(ncol(left.mat)>1){
    left.colors=colours
  }
  if(nrow(bottom.mat)==1){
    bottom.colors= colorRampPalette(c("white","#9D384B"))(256)[100:256]
  } else if(nrow(bottom.mat)>1){
    bottom.colors=colours
  }

  ##step1,首先按照mat中的值按照对应colors赋予颜色
  left.list=values2colors(left.mat, col = left.colors)
  bottom.list=values2colors(bottom.mat, col = bottom.colors)
  middle.list=values2colors(middle.mat, col = middle.colors,levels.mat=levels.mat)
  left.mat=left.list$col.mat
  middle.mat=middle.list$col.mat
  bottom.mat=bottom.list$col.mat
  left.breaks=left.list$breaks
  middle.breaks=middle.list$breaks
  bottom.breaks=bottom.list$breaks

  ##step2,根据mat的维度等设置布局
  #grid.newpage()
  gt=gtable_layout(left.mat,middle.mat,bottom.mat,6,...)
  #gtable_show_layout(gt)

  ##step3,依次将Grob对象写入gtable
  #将列名根据布局位置依次写入gtable
  left.rown=draw_rownames(rownames(left.mat),position="left",angle_col = 0,fontface=3,...)
  gt = gtable_add_grob(gt, left.rown, t = 2, l = 1, clip = "off", name = "left.rown")
  bottom.rown=draw_rownames(rownames(bottom.mat),position="right",angle_col = 0,...)
  gt = gtable_add_grob(gt, bottom.rown, t = 4, l = 5, clip = "off", name = "bottom.rown")

  #将列名根据布局位置依次写入gtable
  left.coln=draw_colnames(colnames(left.mat),position="top",angle_col = 90,...)
  gt = gtable_add_grob(gt, left.coln, t = 1, l = 2, clip = "off", name = "left.coln")
  middle.coln=draw_colnames(colnames(middle.mat),position="top",angle_col = 90,fontface=3,...)
  gt = gtable_add_grob(gt, middle.coln, t = 1, l = 4, clip = "off", name = "middle.coln")

  bottom.coln=draw_colnames(colnames(bottom.mat),position="bottom",angle_col = 90,fontface=3,...)
  gt = gtable_add_grob(gt, bottom.coln, t = 5, l = 4, clip = "off", name = "bottom.coln")
  #grid.newpage()
  #grid.draw(gt)
  #将3热图根据布局位置依次写入gtable
  left=draw_matrix(left.mat)
  middle = draw_matrix(middle.mat)
  bottom = draw_matrix(bottom.mat)
  gt = gtable_add_grob(gt, left, t = 2, l = 2, clip = "off", name = "left")
  gt = gtable_add_grob(gt, middle, t = 2, l = 4, clip = "off", name = "middle")
  gt = gtable_add_grob(gt, bottom, t = 4, l = 4, clip = "off", name = "bottom")
  #grid.newpage()
  #grid.draw(gt)

  #将legend根据布局位置依次写入gtable
  left.legend=draw_legend(color=left.colors, breaks=left.breaks,...)
  gt = gtable_add_grob(gt, left.legend, t = 5, l = 1, clip = "off", name = "left")

  bottom.legend=draw_legend(color=bottom.colors, breaks=bottom.breaks,...)
  gt = gtable_add_grob(gt, bottom.legend, t = 5, l = 2, clip = "off", name = "bottom")

  middle.legend=draw_legend(color=middle.colors, breaks=middle.breaks,just="bottom",...)
  gt = gtable_add_grob(gt, middle.legend, t = 2, l = 5, clip = "off", name = "middle")

  grid.newpage()
  grid.draw(gt)
}
#####
#' @title 将matrix中的数值转换成颜色
#'
#' @param mat matrix
#' @param col colours
#' @param levels.mat character matrix的level
#' @param na_col 默认为"white"
#' @return list,breaks为数值区间，col.mat为颜色矩阵
#'
#' @export
#'
values2colors = function(mat, col = hue_pal()(20), na_col="white",levels.mat=NA){
  mat = as.matrix(mat)#确保是matrix
  if(is.numeric(mat)){
    breaks=generate_breaks(mat, length(col))
    col.mat=matrix(value_vec_colours(as.vector(mat),
                                     col = col,
                                     breaks = breaks,
                                     na_col = na_col),
                   nrow(mat),
                   ncol(mat),
                   dimnames = list(rownames(mat),colnames(mat)))
    return(list(breaks=breaks,col.mat=col.mat))
  }
  if(is.character(mat)){
    if(is.na(levels.mat)){levels.mat=sort(unique(as.vector(mat)))}
    names(col)=levels.mat
    col.mat=matrix(col[as.vector(mat)],
                   nrow(mat),
                   ncol(mat),
                   dimnames = list(rownames(mat),colnames(mat)))
    return(list(breaks=levels.mat,col.mat=col.mat))
  }
}
#####
#' @title 根据给定数据上下限，返回区间n个上下限，从小到大排列
#'
#' @param x data.frame,vector or matrix
#' @param n 区间数
#' @param center 使用上下限的最大绝对值替换原有的上下限
#' @return n个区间的上下限
#'
#' @export
#'
generate_breaks = function(x, n, center = F){
  if(center){
    m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res = seq(-m, m, length.out = n + 1)
  }
  else{
    res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
  }
  return(res)
}
#####
#' @title 根据颜色和数值区间，把数值向量置换为颜色向量
#'
#' @param x numeric vector
#' @param col colours
#' @param breaks 数值区间
#' @param na_col 默认值，"white"
#' @return 颜色向量
#'
#' @export
#'
value_vec_colours <- function(x, col = hue_pal()(20), breaks, na_col="white"){
  res <- col[as.numeric(cut(x, breaks = breaks, include.lowest = T))]  #########???
  res[is.na(res)] <- na_col
  return(res)
}
#####
#' @title 根据left、middle和bottom中的行列名、行列数量比例，设置整体布局
#'
#' @param left 颜色矩阵
#' @param middle 颜色矩阵
#' @param bottom 颜色矩阵
#' @param gapsize 间隔距离
#' @param fontsize 字体大小
#' @return 根据需求的布局设计
#'
#' @export
#'
gtable_layout = function(left, middle, bottom, gapsize, fontsize=7){

  #首先根据字体尺寸、最长行名或列名，设置行名固定宽度和列名的固定高度
  left.rown.width=set_height_width(rownames(left),hw="width",fontsize=fontsize)
  right.rown.width=set_height_width(rownames(bottom),hw="width",fontsize=fontsize)
  top.coln.height=set_height_width(c(colnames(left),colnames(middle)),hw="height",fontsize=fontsize)
  bottom.coln.height=set_height_width(colnames(bottom),hw="height",fontsize=fontsize)
  #根据left、middle和bottom中的行列数量比例，设置三个热图的相对宽度和高度
  proportion=ncol(left)/(ncol(left)+ncol(middle))
  gapsize=unit(gapsize,"bigpts")
  matrix.width = unit(1, "npc") - left.rown.width - right.rown.width - gapsize
  left.width = proportion * matrix.width
  bottom.width = (1-proportion) * matrix.width

  proportion=nrow(left)/(nrow(left)+nrow(bottom))
  matrix.height = unit(1, "npc") - top.coln.height - bottom.coln.height - gapsize
  left.height = proportion * matrix.height
  bottom.height = (1-proportion) * matrix.height

  # Produce gtable layout
  gp = list(fontsize=fontsize)
  gt = gtable(widths = unit.c(left.rown.width, left.width, gapsize, bottom.width, right.rown.width),
              heights = unit.c(top.coln.height, left.height, gapsize, bottom.height, bottom.coln.height),
              vp = viewport(gp = do.call(gpar, gp))
  )
  return(gt)
}
#####
#' @title 根据文字的相对位置布局行名
#'
#' @param rown rownames
#' @param position "right" or "left"
#' @param ... additional arguments passed on to image
#' @return textGrob,字符串的相对位置
#'
#' @export
#'
draw_rownames = function(rown, position, ...){
  if(position=="right"){
    hjust_row = 0
    vjust_row = 0.5
    x=unit(3, "bigpts")
  }
  if(position=="left"){
    hjust_row = 1
    vjust_row = 0.5
    x=unit(1,"npc")-unit(3, "bigpts")
  }
  n=length(rown);m = 1:n
  coord = unit(m / n, "npc");size = unit(1 / n, "npc")
  y = coord - 0.5 * size
  res = textGrob(rown, x = x, y = y,
                 vjust = vjust_row, hjust = hjust_row, gp = gpar(...))

  return(res)
}
#####
#' @title 根据文字的相对位置布局列名
#'
#' @param coln colnames
#' @param position "top" or "bottom"
#' @param angle_col 旋转角度
#' @param ... additional arguments passed on to image
#' @return textGrob,字符串的相对位置
#'
#' @export
#'
draw_colnames = function(coln,position, angle_col, ...){
  if(angle_col == 0){
    hjust_col = 0.5
    vjust_col = 1
  }
  if(angle_col == 45){
    hjust_col = 1
    vjust_col = 1
  }
  if(angle_col == 90){
    if(position=="top"){
      hjust_col = 0
      vjust_col = 0.5
      y=unit(3, "bigpts")
    }
    if(position=="bottom"){
      hjust_col = 1
      vjust_col = 0.5
      y=unit(1,"npc")-unit(3, "bigpts")
    }
  }
  if(angle_col == 270){
    hjust_col = 0
    vjust_col = 0.5
  }
  if(angle_col == 315){
    hjust_col = 0
    vjust_col = 1
  }
  n=length(coln);m = 1:n
  coord = unit(m / n, "npc");size = unit(1 / n, "npc")
  x = coord - 0.5 * size
  #x,y分别表示文字的坐标位置
  res = textGrob(coln, x = x, y = y,
                 vjust = vjust_col, hjust = hjust_col, rot = angle_col, gp = gpar(...))
  return(res)
}
#####
#' @title 根据matrix中的颜色绘制成热图
#'
#' @param matrix 颜色矩阵
#' @param border_color 默认为"black"
#' @return 热图
#'
#' @export
#'
draw_matrix = function(matrix, border_color="black"){
  #gaps*，间隙位置
  y = nrow(matrix)#y轴tile数
  x = ncol(matrix)#x轴tile数

  #计算热图中tile的列的x和行的y坐标中心位置
  n=x;m = 1:n
  coord = unit(m / n, "npc");xsize = unit(1 / n, "npc")
  x = coord - 0.5 * xsize

  n=y;m = 1:n
  coord = unit(m / n, "npc");ysize = unit(1 / n, "npc")
  y = coord - 0.5 * ysize

  #计算热图中tile的中心坐标
  coord = expand.grid(y = y, x = x)

  htmp = rectGrob(x = coord$x, y = coord$y, width = xsize, height = ysize,
                  gp = gpar(fill = matrix, col = border_color))
  return(htmp)
}
#####
#' @title 绘制legend
#'
#' @param color 所有颜色
#' @param breaks 数值区间的边界
#' @param legend 图例
#' @param just "top" or "bottom"
#' @param ... additional arguments passed on to image
#' @importFrom stats dist hclust
#' @return 绘制的legend
#'
#' @export
#'
draw_legend = function(color, breaks, legend=NA, just="top",...){
  if(is.numeric(breaks)){
    #图例中的数字
    if(is.na(legend)){
      legend = grid.pretty(range(as.vector(breaks)))
    }
    color = color[!is.infinite(color)]
    breaks = breaks[!is.infinite(breaks)]

    #设置固定高度36bigpts
    height = min(unit(1, "npc"), unit(36, "bigpts"))
    #设置bar的下限
    if(just=="top"){
      start = (unit(1, "npc") - height)
    }else if(just=="bottom"){
      start= unit(0, "npc")
    }
    #图例中的数字相对位置
    legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
    #图例中的数字位置，最低点位置+相对位置
    legend_pos = start + height * legend_pos

    breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
    breaks = start + height * breaks
    #每一个tile的长度
    break_distances = breaks[-1] - breaks[-length(breaks)]


    #绘制tile，宽度10bigpts
    tile.width=10
    tile.rect = rectGrob(x = unit(5, "bigpts"), y = breaks[-length(breaks)], width = unit(tile.width, "bigpts"), height = break_distances,
                         hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
    black.rect = rectGrob(x = unit(5, "bigpts"), y = breaks[1], width = unit(tile.width, "bigpts"), height = sum(break_distances),
                          hjust = 0, vjust = 0, gp = gpar(color="black",fill="transparent"))
    axis.ticks = segmentsGrob(x0 = unit(5+tile.width, "bigpts"), y0 = legend_pos,
                              x1 = unit(5+tile.width+0.5, "bigpts"), y1 = legend_pos,
                              gp = gpar(color="black"))
    axis.title = textGrob(legend, x = unit(5+tile.width+2, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(color="black",...))
    res = grobTree(tile.rect, black.rect,axis.ticks,axis.title)
  }

  if(is.character(breaks)){
    #绘制tile，宽度和高度10bigpts,tile间隙2bigpts
    tile.size=10
    tile.gap=2
    tile.unit=tile.size + tile.gap
    height = min(unit(1, "npc"), unit(tile.unit*length(breaks), "bigpts"))

    #设置legend的下限
    if(just=="top"){
      start = (unit(1, "npc") - height)
    }else if(just=="bottom"){
      start= unit(0, "npc")
    }
    #每一个tile的长度
    tile_len = (tile.size/tile.unit)*height/length(breaks)
    y=start + (height/length(breaks))*(0:(length(breaks)-1))
    tile.rect = rectGrob(x = unit(5, "bigpts"), y = y, width = unit(tile.size, "bigpts"), height = tile_len,
                         hjust = 0, vjust = 0, gp = gpar(fill = color, col = "black"))
    title = textGrob(breaks, x = unit(5+tile.size+2, "bigpts"), y = y+tile_len/2, hjust = 0, gp = gpar(color="black",...))
    res = grobTree(tile.rect, title)
  }
  return(res)
}
######
#' @title 绘制细胞互作热图
#'
#' @param res list,细胞互作分析结果，包含pvalue,significant,receptor_logfc,ligand_logfc
#' @param levels.mat character matrix的level
#' @param transformation 数据变换，把极值转换到特定范围内
#' @return 细胞互作热图
#'
#' @export
#'
map2 <- function(res,levels.mat=c("0","1"),transformation=T){
  
  if(ncol(res$ligand_logfc)==1){
    clust_l=clustering(res$ligand_logfc,row=T,col=F)
  } else{
    clust_l=clustering(res$ligand_logfc,row=T,col=T)
  }
  
  if(ncol(res$receptor_logfc)==1){
    clust_r=clustering(res$receptor_logfc,row=T,col=F)
  } else{
    clust_r=clustering(res$receptor_logfc,row=T,col=T)
  }
  bottom.mat=t(clust_r)
  left.mat=clust_l
  
  up.middle.mat=res$up.middle.mat[rownames(clust_l),rownames(clust_r)]
  down.middle.mat=res$down.middle.mat[rownames(clust_l),rownames(clust_r)]
  up.middle.mat[up.middle.mat< 1]=0
  up.middle.mat[up.middle.mat>=1]=1
  down.middle.mat[down.middle.mat< 1]=0
  down.middle.mat[down.middle.mat>=1]=1
  if(transformation==T){
  if(min(left.mat) < 0){
    left.mat=transformation(left.mat,range=c(-1,1))
  }}
  left.mat=as.matrix(left.mat)
  up.middle.mat=as.matrix(up.middle.mat)
  down.middle.mat=as.matrix(down.middle.mat)
  mode(up.middle.mat)="character"
  mode(down.middle.mat)="character"
  if(transformation==T){
  if(min(bottom.mat) < 0){
    bottom.mat=transformation(bottom.mat,range=c(-1,1))
  }}
  bottom.mat=as.matrix(bottom.mat)
  lr2map2(left.mat,up.middle.mat,down.middle.mat,bottom.mat,fontsize=7,levels.mat=c("0","1"))
}
######
#' @title 根据matrix中的颜色绘制成三角热图
#'
#' @param matrix 表达矩阵
#' @param border_color "grey60"
#' @param border_size 0.3
#' @return 三角热图
#' @export
#'
draw_triangle_up = function(matrix, border_color="grey60",border_size=0.3){
  #gaps*，间隙位置
  y = nrow(matrix)#y轴triangle数
  x = ncol(matrix)#x轴triangle数

  #计算热图中triangle的列的x和行的y坐标中心位置
  n=x;
  xcoord_b =xcoord_t= unit(0:(n-1) / n, "npc")
  xcoord_r = unit(1:n / n, "npc")
  x3.mat=rbind(rep(xcoord_t,each=y),rep(xcoord_r,each=y),rep(xcoord_b,each=y))

  n=y;
  ycoord_t = ycoord_r =unit(1:n / n, "npc")
  ycoord_b = unit(0:(n-1) / n, "npc")
  y3=rbind(ycoord_t,ycoord_r,ycoord_b)
  y3.mat=matrix(rep(y3,x),nrow=nrow(y3))

  htmp = polygonGrob(x=x3.mat,y=y3.mat,id.lengths = rep(3,x*y),
                     gp = gpar(fill = as.vector(matrix), col = border_color,lwd=border_size))
  return(htmp)
}
#######
#' @title 根据matrix中的颜色绘制成三角热图
#'
#' @param matrix 表达矩阵
#' @param border_color "grey60"
#' @param border_size 0.3
#' @return 三角热图
#' @export
#'
draw_triangle_down = function(matrix, border_color="grey60",border_size=0.3){
  #gaps*，间隙位置
  y = nrow(matrix)#y轴triangle数
  x = ncol(matrix)#x轴triangle数

  #计算热图中triangle的列的x和行的y坐标中心位置
  n=x;
  xcoord_l=unit(0:(n-1) / n, "npc")
  xcoord_t =xcoord_b= unit(1:n / n, "npc")
  x3.mat=rbind(rep(xcoord_t,each=y),rep(xcoord_l,each=y),rep(xcoord_b,each=y))

  n=y;
  ycoord_t = unit(1:n / n, "npc")
  ycoord_b = ycoord_l= unit(0:(n-1) / n, "npc")
  y3=rbind(ycoord_t,ycoord_l,ycoord_b)
  y3.mat=matrix(rep(y3,x),nrow=nrow(y3))

  htmp = polygonGrob(x=x3.mat,y=y3.mat,id.lengths = rep(3,x*y),
                     gp = gpar(fill = as.vector(matrix), col = border_color,lwd=border_size))
  return(htmp)
}
#######
#' @title 绘制双matrix热图
#'
#' @param leftmatrix 左侧矩阵
#' @param rightmatrix 右侧矩阵
#' @return 双matrix热图
#' @export
#'
draw_diagonalmap=function(leftmatrix,rightmatrix){
  lefttriangle=draw_triangle_up(leftmatrix)
  righttriangle=draw_triangle_down(rightmatrix)
  res = grobTree(lefttriangle,righttriangle)
  return(res)
}
###
#' @title 绘制热图
#'
#' @param left.mat left矩阵
#' @param bottom.mat bottom矩阵
#' @param up.middle.mat up.middle矩阵
#' @param down.middle.mat down.middle矩阵
#' @param left.colors left矩阵颜色
#' @param middle.colors middle矩阵颜色
#' @param bottom.colors bottom矩阵颜色
#' @param levels.mat character matrix的level
#' @param ... additional arguments passed on to image
#' @return 热图
#'
#' @export
#'
lr2map2=function(left.mat,up.middle.mat,down.middle.mat,bottom.mat,left.colors=NA,middle.colors=NA,
                 bottom.colors=NA,levels.mat=NA,...){
  #input:
  if(is.na(left.colors)){
    rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
    r=c(157,171,185,199,216,219,228,234,237,232,227,219,211,186,164,142,120,105,88,70,60)
    g=c(56,90,113,133,156,169,187,203,221,222,225,224,226,212,199,188,175,163,149,139,128)
    b=c(75,88,101,113,124,144,165,185,206,214,224,231,236,232,231,225,217,210,206,197,187)
    coldf=data.frame(r,g,b)

    white=rgb(245, 246, 247, maxColorValue = 255)
    red=rgb(174, 61, 141, maxColorValue = 255)
    colours=rev(apply(coldf,1,rgb2hex))
  }
  if(is.na(middle.colors)){
    rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
    r=c(157,171,185,199,216,219,228,234,237,232,227,219,211,186,164,142,120,105,88,70,60)
    g=c(56,90,113,133,156,169,187,203,221,222,225,224,226,212,199,188,175,163,149,139,128)
    b=c(75,88,101,113,124,144,165,185,206,214,224,231,236,232,231,225,217,210,206,197,187)
    coldf=data.frame(r,g,b)

    white=rgb(245, 246, 247, maxColorValue = 255)
    red=rgb(174, 61, 141, maxColorValue = 255)

    colours=rev(apply(coldf,1,rgb2hex))
    middle.colors=c(white,red)
  }
  if(is.na(bottom.colors)){
    rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
    r=c(157,171,185,199,216,219,228,234,237,232,227,219,211,186,164,142,120,105,88,70,60)
    g=c(56,90,113,133,156,169,187,203,221,222,225,224,226,212,199,188,175,163,149,139,128)
    b=c(75,88,101,113,124,144,165,185,206,214,224,231,236,232,231,225,217,210,206,197,187)
    coldf=data.frame(r,g,b)

    white=rgb(245, 246, 247, maxColorValue = 255)
    red=rgb(174, 61, 141, maxColorValue = 255)
    colours=rev(apply(coldf,1,rgb2hex))
  }

  if(ncol(left.mat)==1){
    left.colors=colorRampPalette(c("white","#9D384B"))(256)[100:256]
  } else if(ncol(left.mat)>1){
    left.colors=colours
  }
  if(nrow(bottom.mat)==1){
    bottom.colors= colorRampPalette(c("white","#9D384B"))(256)[100:256]
  } else if(nrow(bottom.mat)>1){
    bottom.colors=colours
  }
  ##step1,首先按照mat中的值按照对应colors赋予颜色
  left.list=values2colors(left.mat, col = left.colors)
  bottom.list=values2colors(bottom.mat, col = bottom.colors)
  up.middle.mat=
    middle.mat=rbind(up.middle.mat,down.middle.mat)
  middle.list=values2colors(middle.mat, col = middle.colors,levels.mat=levels.mat)
  left.mat=left.list$col.mat
  middle.mat=middle.list$col.mat
  bottom.mat=bottom.list$col.mat
  left.breaks=left.list$breaks
  middle.breaks=middle.list$breaks
  bottom.breaks=bottom.list$breaks
  ##step2,根据

  ##step2,根据mat的维度等设置布局
  #grid.newpage()
  gt=gtable_layout(left.mat,middle.mat,bottom.mat,6,...)
  #gtable_show_layout(gt)

  ##step3,依次将Grob对象写入gtable
  #将列名根据布局位置依次写入gtable
  left.rown=draw_rownames(rownames(left.mat),position="left",angle_col = 0,fontface=3,...)
  gt = gtable_add_grob(gt, left.rown, t = 2, l = 1, clip = "off", name = "left.rown")
  bottom.rown=draw_rownames(rownames(bottom.mat),position="right",angle_col = 0,...)
  gt = gtable_add_grob(gt, bottom.rown, t = 4, l = 5, clip = "off", name = "bottom.rown")

  #将列名根据布局位置依次写入gtable
  left.coln=draw_colnames(colnames(left.mat),position="top",angle_col = 90,...)
  gt = gtable_add_grob(gt, left.coln, t = 1, l = 2, clip = "off", name = "left.coln")
  middle.coln=draw_colnames(colnames(middle.mat),position="top",angle_col = 90,fontface=3,...)
  gt = gtable_add_grob(gt, middle.coln, t = 1, l = 4, clip = "off", name = "middle.coln")

  bottom.coln=draw_colnames(colnames(bottom.mat),position="bottom",angle_col = 90,fontface=3,...)
  gt = gtable_add_grob(gt, bottom.coln, t = 5, l = 4, clip = "off", name = "bottom.coln")
  #grid.newpage()
  #grid.draw(gt)
  #将3热图根据布局位置依次写入gtable
  left=draw_matrix(left.mat)
  middle = draw_diagonalmap(middle.mat[1:(nrow(middle.mat)/2),],middle.mat[(1+(nrow(middle.mat)/2)):nrow(middle.mat),])
  bottom = draw_matrix(bottom.mat)
  gt = gtable_add_grob(gt, left, t = 2, l = 2, clip = "off", name = "left")
  gt = gtable_add_grob(gt, middle, t = 2, l = 4, clip = "off", name = "middle")
  gt = gtable_add_grob(gt, bottom, t = 4, l = 4, clip = "off", name = "bottom")
  #grid.newpage()
  #grid.draw(gt)

  #将legend根据布局位置依次写入gtable
  left.legend=draw_legend(color=left.colors, breaks=left.breaks,...)
  gt = gtable_add_grob(gt, left.legend, t = 5, l = 1, clip = "off", name = "left")

  bottom.legend=draw_legend(color=bottom.colors, breaks=bottom.breaks,...)
  gt = gtable_add_grob(gt, bottom.legend, t = 5, l = 2, clip = "off", name = "bottom")

  middle.legend=draw_legend(color=middle.colors, breaks=middle.breaks,just="bottom",...)
  gt = gtable_add_grob(gt, middle.legend, t = 2, l = 5, clip = "off", name = "middle")

  grid.newpage()
  grid.draw(gt)
}

