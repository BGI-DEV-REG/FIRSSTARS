########## 定义需要的函数
netVisual_aggregate_value = function (object, signaling, signaling.name = NULL, 
interaction = NULL,color.use = NULL, 
    thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, 
    targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE, 
    vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL, 
    weight.scale = TRUE, edge.weight.max = NULL, edge.width.max = 8, 
    layout = c("circle", "hierarchy", "chord", "spatial"), pt.title = 12, 
    title.space = 6, vertex.label.cex = 0.8, alpha.image = 0.15, 
    point.size = 1.5, group = NULL, cell.order = NULL, small.gap = 1, 
    big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, 
    legend.pos.x = 20, legend.pos.y = 20, ...) {
    layout <- match.arg(layout)
    if (is.null(vertex.weight)) {
        vertex.weight <- as.numeric(table(object@idents))
    }
    if (is.null(vertex.size.max)) {
        if (length(unique(vertex.weight)) == 1) {
            vertex.size.max <- 5
        }
        else {
            vertex.size.max <- 15
        }
    }

    net <- object@net
    pairLR.use.name <- dimnames(net$prob)[[3]]

    pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, 
        key = "pathway_name", matching.exact = T, pair.only = T)
    if (is.null(signaling.name)) {
        signaling.name <- signaling
    }

    if(is.null(interaction)==F){
      pairLR = pairLR[interaction,]
    }

    pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
    pairLR <- pairLR[pairLR.name,]
  
    prob <- net$prob
    pval <- net$pval
    prob[pval > thresh] <- 0
    if (length(pairLR.name) > 1) {
        pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 
            3, sum) != 0]
    }
    else {
        pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 
            0]
    }
    if (length(pairLR.name.use) == 0) {
        stop(paste0("There is no significant communication of ", 
            signaling.name))
    }
    else {
        pairLR <- pairLR[pairLR.name.use, ]
    }
    nRow <- length(pairLR.name.use)
    prob <- prob[, , pairLR.name.use]
    pval <- pval[, , pairLR.name.use]
    if (length(dim(prob)) == 2) {
        prob <- replicate(1, prob, simplify = "array")
        pval <- replicate(1, pval, simplify = "array")
    }
    #prob.sum <- apply(prob, c(1, 2), sum)
    return(prob)
}

netVisual_circle_modify = function (net, color.use = NULL, title.name = NULL, sources.use = NULL, 
                                    targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, 
                                    top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
                                    vertex.size.max = NULL, vertex.label.cex = 1, vertex.label.color = "black", 
                                    edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
                                    label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
                                    edge.curved = 0.2, shape = "circle", layout = in_circle(), 
                                    margin = 0.2, vertex.size = NULL, arrow.width = 1, arrow.size = 0.2) {
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    }
    else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1 - top)
  net[net < thresh] <- 0
  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use))) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) | 
                         (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g), 
                                                                             2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g), 
                                                                                                                                       2]/coords_scale[igraph::V(g), 1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)$name]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)$name]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::V(g)$color[edge.start[, 1]], alpha.edge)
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))
  
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,  2] == edge.start[, 1]), 1]]
  }
  return(g)
}
########### 开始画图
library(CellChat)
library(Seurat)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(SuppDists)
non_protein=TRUE
CellChatDB.use <- CellChatDB.mouse
lr_meta= CellChatDB.use$interaction

key_fb = c('FB_oxidative_stress_Adult','FB_inflammatory_Adult', 'FB_fascia_1_Adult','FB_oxidative_stress_EP','FB_inflammatory_EP', 'FB_fascia_EP','FB_myo_Adult')
mnp = c("Mono_Ly6c","Mono_Spp1","Macr_Ms4a7","Macr_Mmp14","Macr_Trem2" ,"Macr_Lyve1")

###################################
#############画第一个图##############
###################################

######### 读取cellchat对象，并获取里面的netP (cellchat计算结果) 
df.list = list()
for(tpi in c('Early',"Median","Late")){
    for(age in c("Adult","P5","E")){
        load(paste0("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/08.immune/Rdata/CellChat_0521_Rdata/",age,"_",tpi,".CellChat.Rdata"))
        df.net['Age'] = age
        df.net[['TPI']] = tpi
        df.list[[paste0(age,"_",tpi)]]=df.net
    }}

nep = Reduce(rbind,df.list)
######### 选取要可视化的通路
table_dir = "/Users/zhaolab3/Documents/skin_immune/table/"
Chemokine = read.csv(file.path(table_dir,'Cytokine_mouse.csv'),row.name='X')$x
cmeta = lr_meta[lr_meta$ligand%in%Chemokine,]
write.csv(cmeta,file = file.path(table_dir,'Chemokine ligands and receptors of mouse.csv'),quote = F)
pathways.show <- unique(lr_meta[lr_meta$ligand%in%Chemokine,"pathway_name"])

######### 过滤要可视化的通路和细胞类型，整理可视化的Data.frame
nep <- nep %>% 
  select(source ,target, prob,interaction_name,interaction_name_2,pathway_name, Age,TPI)
nep = nep[nep$source%in%c(mnp,key_fb)&nep$target%in%c(mnp,key_fb),]

data.nep = nep[nep$pathway_name%in%pathways.show,]

data.nep = cbind(
    data.frame(celltype = c(data.nep$source,data.nep$target),type = c(rep("source",nrow(data.nep)),rep("target",nrow(data.nep)))),
rbind(data.nep[,2:ncol(data.nep)],data.nep[,2:ncol(data.nep)]))
data.nep$type1 <- data.nep$type
data.nep$type <- paste(data.nep$Age, data.nep$type)

sum.nep <- data.nep %>%
  group_by(type1,type,celltype,TPI) %>% 
  summarise(sum_prob = sum(prob))

sum.nep <- sum.nep%>%
  mutate(sum_prob1 = if_else(type %in% c("Adult source", "P5 source","E source"), -sum_prob, sum_prob))

sum.nep$type <- factor(sum.nep$type,levels = c("E source","P5 source","Adult source","E target","P5 target","Adult target"))  

######### 开始分别画三个时期的图
age_color = c("#26c6da","#f9320c","#FFF300","#26c6da","#f9320c","#FFF300")
names(age_color) = levels(sum.nep$type)

outpath = "/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/08.immune/fig/LR_point_bar"
dir.create(outpath)

plot.list = list()
for(tpi in c('Early',"Median","Late")){
  sub.nep.1 = sum.nep[sum.nep$TPI==tpi,]
  file_name = 'Chemo_All'
  sub.nep = sub.nep.1
  data.bar <- sub.nep %>% select(type, sum_prob) %>% group_by(type) %>% summarise(sum_prob = sum(sum_prob))

  ## 画 barplot
  data.sub <- subset(sub.nep, type1== "source")
  data.sub = arrange(data.sub,desc(sum_prob))
  idx = c(intersect(data.sub$celltype,mnp),intersect(data.sub$celltype,CommonName(key_fb)))
  idx = unique(data.sub$celltype)

  p1 <- ggplot(data.bar, aes(x = type, y = sum_prob, fill = type)) +
  geom_col(position = "dodge", colour = "black", size = 0.6, width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", sum_prob)), position = position_dodge(width = 0.7),
              vjust = -0.5, size = 4, color = "white") +
  labs(title =file_name,xlab =NULL)+theme_black()+theme(axis.title.x = element_blank(),axis.text.x = element_blank())+
  scale_fill_manual(values = age_color )+
  scale_y_continuous(expand = c(0,0),limits=c(0,max(data.bar$sum_prob)*1.2))
  
  ## 画 pointplot
  sub.nep$celltype <- factor(sub.nep$celltype,levels = rev(idx))

  p2 <- ggplot(sub.nep, aes(x = sum_prob1, y = celltype, color = type, size = 8)) +  
  geom_point(alpha = 1) +
  geom_vline(xintercept = 0, color = "white",size = 1)+
  scale_color_manual(values =age_color)+
  theme_black()+coord_cartesian(xlim = c(-ceiling(max(sub.nep$sum_prob)),ceiling(max(sub.nep$sum_prob))))
  p <- patchwork::wrap_plots(plots = list(p1,p2),ncol = 1,heights=c(1/4,3/4))
  plot.list[[paste(tpi,file_name, sep = "_")]] = p
  ggsave(file.path(outpath,paste(file_name,tpi,"point_bar.pdf",sep="_")), p, width = 8, height = 16)
}

###################################
#############画网络图###############
###################################
load("./Rdata/CellChat_0521_Rdata/P5_Median.CellChat.Rdata")
#load("./Rdata/CellChat_0521_Rdata/Adult_Median.CellChat.Rdata")
# 读取不同时期的cellchat对象，分时期画

cell_level_all =levels(cellchat@idents)

print(pathways.show)  #按照pathway选
interaction_name=c("CXCL1_ACKR1" ,"CXCL2_ACKR1" ,"CXCL3_ACKR1" ,"SPP1_CD44" , "JAG1_NOTCH2",
"JAG1_NOTCH4","IL6_IL6R_IL6ST","SPP1_ITGAV_ITGB3","SPP1_ITGAV_ITGB1", "SPP1_ITGAV_ITGB5",
 "SPP1_ITGAV_ITGB6" ,"SPP1_ITGA4_ITGB1", "SPP1_ITGA9_ITGB1", "SPP1_ITGA8_ITGB1" ,"SPP1_ITGA5_ITGB1")
 #按照LR选

out_dir = "/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/08.immune/fig"

prob.sum = netVisual_aggregate_value(cellchat, 
sources.use = mnp , 
targets.use = intersect(key_fb,cell_level_all), 
signaling=pathways.show,
remove.isolate = T,
color.use = col1,
layout =  "circle")
prob.sum <- apply(prob.sum, c(1, 2), sum)

g = netVisual_circle_modify(prob.sum*5,
sources.use = mnp ,
targets.use = intersect(c("FB_fascia_1_Adult","FB_inflammatory_Adult","FB_oxidative_stress_Adult",
"FB_fascia_EP","FB_inflammatory_EP","FB_oxidative_stress_EP" ),cell_level_all), 
remove.isolate = T,
vertex.weight.max = NULL, 
vertex.size.max =6,
vertex.weight = 20, 
weight.scale = F,
color.use = col1)

radian.rescale <- function(x, start = 0, direction = 1) {
c.rotate <- function(x) (x + start)%%(2 * pi) * direction
c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
                            direction = -1, start = 0)
vertex.weight = 20
label.dist <- vertex.weight/max(vertex.weight) + 2

pdf(file = file.path(out_dir,"P5-Late-circle-Spp1.pdf"),
    width = 10, height = 6)
plot(g,  
     layout=layout_in_circle(g),  # 使用星形布局（或其他你喜欢的布局）  
     bg="black",                # 设置背景色为黑色  
     edge.label.cex= 1,
     edge.label = round(E(g)$weight,2),    # 尝试为每条边设置标签（权重）  
     edge.label.color = "black",  # 设置边标签颜色为白色  
     edge.label.family="Helvetica",
     edge.arrow.size = 1.5,       # 设置箭头大小（如果需要的话） 
     vertex.label.cex= 1,      # 设置顶点标签的大小  
     vertex.label.dist=2,     # 设置顶点标签与顶点的距离  
     vertex.label.family="Helvetica", # 设置顶点标签的字体  
     main="Black Background with White Text", # 设置图形标题 
) 
dev.off()

pdf(file = file.path(out_dir,"Adult-Late-circle-Spp1.pdf"),
    width = 10, height = 6)
plot(g,  
     layout=layout_in_circle(g),  # 使用星形布局（或其他你喜欢的布局）  
     bg="black",                # 设置背景色为黑色  
     edge.label.cex= 1,
     edge.label = round(E(g)$weight,2),    # 尝试为每条边设置标签（权重）  
     edge.label.color = "black",  # 设置边标签颜色为白色  
     edge.arrow.size = 1.5,       # 设置箭头大小（如果需要的话） 
     edge.label.family="Helvetica",
     vertex.label.cex= 1,      # 设置顶点标签的大小  
     vertex.label.dist=2,     # 设置顶点标签与顶点的距离  
     vertex.label.family="Helvetica", # 设置顶点标签的字体  
     main="Black Background with White Text", # 设置图形标题 
)  
dev.off()