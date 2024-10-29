# /jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/EP_FB_240414/FB_subset/240519_anno/marker_heatmap


library(Seurat)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)

setwd("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/EP_FB_240414/FB_subset/240519_anno/marker_heatmap")
df = read.csv("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/EP_FB_240414/FB_subset/240519_anno/marker_heatmap/conservedmarkers.csv")

obj_avg = readRDS("avg.rds")
obj_avg_scale = cbind(t(apply(obj_avg[, endsWith(colnames(obj_avg), "EP")], 1, scale)), t(apply(obj_avg[, endsWith(colnames(obj_avg), "Adult")], 1, scale)))
obj_avg_scale = apply(obj_avg_scale[df$gene, ], 2,scale)

rownames(obj_avg_scale) = df$gene
colnames(obj_avg_scale) = colnames(obj_avg)


# celltype  <-  c("FB_DP_Adult", "FB_DS_Adult", "FB_adipogenic_Adult", "FB_basement_resident_Adult", "FB_basement_activated_Adult", "FB_papi_Adult", "FB_reti_Adult","FB_muscle_Adult", "FB_fascia_1_Adult","FB_fascia_2_Adult", "FB_inflammatory_Adult", "FB_oxidative_stress_Adult","FB_exosome_Adult","FB_apoptotic_Adult", "FB_myo_Adult","FB_DP_EP", "FB_DS_EP", "FB_adipogenic_EP", "FB_basement_resident_EP", "FB_basement_activated_EP", "FB_papi_EP", "FB_reti_EP","FB_muscle_EP", "FB_fascia_EP", "FB_inflammatory_EP",  "FB_oxidative_stress_EP","FB_exosome_EP")
celltype <- c("FB_DP_EP", "FB_DS_EP", "FB_adipogenic_EP", "FB_basement_resident_EP", "FB_basement_activated_EP", "FB_papi_EP", "FB_reti_EP","FB_muscle_EP", "FB_fascia_EP", "FB_inflammatory_EP",  "FB_oxidative_stress_EP","FB_exosome_EP","FB_DP_Adult", "FB_DS_Adult", "FB_adipogenic_Adult", "FB_basement_resident_Adult", "FB_basement_activated_Adult", "FB_papi_Adult", "FB_reti_Adult","FB_muscle_Adult", "FB_fascia_1_Adult","FB_fascia_2_Adult", "FB_inflammatory_Adult", "FB_oxidative_stress_Adult","FB_exosome_Adult","FB_apoptotic_Adult", "FB_myo_Adult")


best_color<- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941",
               "#006FA6","#A30059","#FFE4E1","#0000A6","#63FFAC",
               "#B79762","#004D43","#8FB0FF","#997D87","#5A0007",
               "#809693","#1B4400","#4FC601","#3B5DFF","#FF2F80",
               "#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9",
               "#B903AA","#DDEFFF","#7B4F4B","#A1C299","#0AA6D8",
               "#00A087FF","#4DBBD5FF","#E64B35FF","#3C5488FF","#F38400",
               "#A1CAF1", "#C2B280","#848482","#E68FAC", "#0067A5", 
               "#F99379", "#604E97","#F6A600", "#B3446C","#DCD300",
               "#882D17", "#8DB600","#654522", "#E25822", "#2B3D26",
               "#191970","#000080",
               "#6495ED","#1E90FF","#00BFFF","#00FFFF","#FF1493",
               "#FF00FF","#A020F0","#63B8FF","#008B8B","#54FF9F",
               "#00FF00","#76EE00","#FFF68F","Yellow1","Gold1",
               "DarkGoldenrod4","#FF6A6A","#FF8247","#FFA54F","#FF7F24",
               "#FF3030","#FFA500","#FF7F00","#FF7256","#FF6347",
               "#FF4500","#FF1493","#FF6EB4","#EE30A7","#8B008B")
               

go_df = subset(df, df$GO != "")
cols = best_color[as.integer(as.factor(go_df$GO))]
cols[go_df$GO == "other"] = "#FFFFFF"
gp = gpar(col = cols)
ha = rowAnnotation(foo = anno_mark(at = as.integer(rownames(go_df)), 
    labels = go_df$gene, lines_gp = gp, labels_gp=gp))

col_fun <- circlize::colorRamp2(c(-2.5, 1, 2.5), viridis(30)[c(1,15, 30)])

p_heats <- Heatmap(obj_avg_scale[df$gene,celltype],
        heatmap_legend_param=list(title="expression"),
        column_title = "clustered dotplot",
        col=col_fun,
        row_order = 1:500,
        # cluster_rows = T,
        # cluster_columns = T,
        column_order = 1:27,
        column_split = endsWith(celltype, "Adult"),
        row_names_gp = gpar(type = "none"),
        # column_names_gp = gpar(type = "none"),
        # row_names_gp = gpar(col = "white"),
        column_names_gp = gpar(col = "white"),
        border = "white",
        right_annotation = ha
        )

pdf("meta_all_heatmap_scale_markers_240930.pdf", width=15, height=20)
p_heats <- draw(p_heats, background = "black")
dev.off()