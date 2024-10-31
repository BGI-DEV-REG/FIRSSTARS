library(tidyverse)
library(viridis)
library(ggdark)

first_FB_l = c("Ccl3","Il6","Sirpa",
               "Cxcl2","Cxcl3","Il1b", 
               "Cxcl10","Cxcl16","Ccl7")

FB_color_df <- read.csv("E:/CodePrograms/Immu/Immu_plot_0718/Fig2_color_seq_20240521.csv")
FB_color <- FB_color_df$colors
names(FB_color) <- FB_color_df$X

score_list <- list()
for ( age in c("Adult", "P5", "E")) {
  A_score <- read.csv(paste0("E:/CodePrograms/Immu/Immu_0621/",age,"_Early_meta.csv"))
  
  if (age == "E") {
    A_score$Age <- "E"
  }else{
    NULL
  }
  
  A_score$celltype <- sub("_[^_]+$", "", A_score$celltype)
  A_score$celltype_phase <- str_c(A_score$celltype, "_", A_score$Age)
  
  A_score_aver <- A_score %>% 
    group_by(celltype_phase) %>% 
    mutate(Overall_Score = sum(GOBP_MONOCYTE_Score1),
           Average_Mono_Score = Overall_Score / length(celltype_phase) )
  
  A_df <- as.data.frame(table(A_score_aver$Average_Mono_Score, A_score_aver$celltype_phase)) %>% 
    filter(Freq != 0)
  
  names(A_df) <- c("Score", "celltype", "cellnum") 
  
  join_A_df <- left_join(A_score, A_df, by = c("celltype_phase" = "celltype"))
  join_A_df$S_score <- as.numeric(as.character(join_A_df$Score))
  
  score_list[[age]] <- join_A_df
}

score_df <- Reduce(rbind, score_list)
FB_celltype <- grep("^FB", score_df$celltype, value = T)
FB_score_df <- score_df %>% 
  filter(celltype %in% FB_celltype) %>% 
  filter(!(celltype_phase %in% c("FB_adipogenic_Adult","FB_basement_activated_Adult","FB_basement_resident_Adult","FB_DS_Adult","FB_papi_Adult")))

FB_score_df$celltype_phase[FB_score_df$celltype_phase == "FB_fascia_1_Adult"] <- "FB_fascia_Adult"

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

FB_normalize_df <- FB_score_df %>% 
  mutate(Age = factor(Age, levels = c("E", "P5", "Adult"))) %>% 
  arrange(Age, Score) %>% 
  filter(celltype_phase %in% celltype_phase_level) # %>% 
# mutate(normalize_GOBP_MONOCYTE_Score1 = normalize(GOBP_MONOCYTE_Score1))

FB_normalize_df$celltype_phase <- factor(FB_normalize_df$celltype_phase, 
                                         rev(unique(FB_normalize_df$celltype_phase)))

color_Score <- plasma(length(unique(FB_normalize_df$S_score)))  
names(color_Score)  <- sort(as.numeric(as.character(unique(FB_normalize_df$Score))))

p_v <-  ggplot(FB_normalize_df, aes(x = celltype_phase, y = GOBP_MONOCYTE_Score1, fill = Score)) +
  coord_flip() +
  dark_theme_classic() +
  geom_violin(aes(fill = Score)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, size = 10, hjust = 1)) +
  scale_fill_manual(values = color_Score) 

ggsave("violin_plot.pdf", p_v, w = 2.5, h = 6)