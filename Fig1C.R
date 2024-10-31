library(tidyverse)

colorlist = c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941",
              "#006FA6", "#A30059", "#FFE4E1", "#0000A6", "#63FFAC",
              "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007",
              "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#FF2F80",
              "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
              "#B903AA", "#DDEFFF", "#7B4F4B", "#A1C299", "#0AA6D8",
              "#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF", "#F38400",
              "#A1CAF1", "#C2B280", "#848482", "#E68FAC", "#0067A5",
              "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
              "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
              "#191970", "#000080", "#6495ED", "#1E90FF", "#00BFFF", 
              "#00FFFF", "#FF1493", "#FF00FF", "#A020F0", "#63B8FF", 
              "#008B8B", "#54FF9F", "#00FF00", "#76EE00", "#FFF68F")

df <- read.csv("E:/CodePrograms/All_Rat/All_RAT_0910/allcell_last_obs_umap.csv")
Age_TPI <- c("E16.5_Nor", "E16.5_6H", "E16.5_1D", "E16.5_2D", "E16.5_3D", "E16.5_Born",
             "E17.5_Nor", "E17.5_6H", "E17.5_1D", "E17.5_2D", "E17.5_3D", "E17.5_Born",
             "E18.5_Nor", "E18.5_6H", "E18.5_1D", "E18.5_2D", "E18.5_3D", "E18.5_Born",
             "P0_Nor", "P5_Nor", "P5_6H", "P5_1D", "P5_2D", "P5_3D", "P5_6D",
             "Adult_Nor", "Adult_1D", "Adult_2D", "Adult_5D", "Adult_10D", "Adult_19D")
df$Age_TimePostInjury <- factor(df$Age_TimePostInjury, rev(Age_TPI))
AgeTPI_num <- as.data.frame(table(df$celltype, df$Age_TimePostInjury))
names(AgeTPI_num) <- c("celltype", "Age_TPI", "number")

AgeTPI_ratio <- AgeTPI_num %>%
  group_by(Age_TPI) %>%
  mutate(prop_number = number / sum(number))

celltype  <-  c("Melanocyte","Muscle","Endothelial","Pericyte","Schwann","Erythrocyte",
                "KC_TAC","KC_IRS", "KC_ORS", "KC_HFSC","KC_HF","KC_SG","KC_Suprabasal","KC_Basal",
                "FB_DS", "FB_DP", "FB_Reticular","FB_Papillary","FB_Fascia", "FB_Injury",  
                "Immu_Neu","Immu_pDC_B","Immu_T_NK","Immu_Myeloid")
AgeTPI_ratio$celltype <- factor(AgeTPI_ratio$celltype, celltype)

p1 <-ggplot(AgeTPI_ratio, aes(x = Age_TPI, y = prop_number, fill = celltype)) +
  coord_flip() +
  geom_col(width = 0.8) +
  dark_theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "white"),
    axis.text.y = element_text(size = 12, color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) + 
  guides(fill = FALSE) +
  scale_fill_manual(values = all_color) +
  scale_y_continuous(expand = c(0,0)) 
ggsave("all_col_stacked_SecondTest.pdf", p1, height =15 , width = 10)

celltype_num <- as.data.frame(table(df$Age_TimePostInjury))
names(celltype_num) <- c("Age_TimePostInjury", "number")
celltype_num $Age_TimePostInjury <- factor(celltype_num$Age_TimePostInjury, rev(Age_TPI))

p2 <- ggplot(celltype_num, aes(y = log10(number), x = Age_TimePostInjury, fill = Age_TimePostInjury)) +
  geom_col() +
  dark_theme_classic() +
  coord_flip() +
  scale_y_continuous(expand = c(0.01,0)) +
  theme(axis.text.x = element_text(size = 12, color = "white"),
        axis.text.y = element_text(size= 12, color = "white"),
        axis.ticks.y = element_line(size = 0.5, color = "white"),
        axis.ticks.x = element_line(size = 0.5, color = "white"), 
        plot.background = element_rect(fill = "black"),
        panel.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
  ) + 
  labs(x = "Cell Type", y = "Number", title = "Cell Type Distribution") +
  guides(fill = FALSE) +
  scale_fill_manual(values = rev(colorlist[1:31]))
ggsave("log10_bar_Age_TPIt_0910.pdf", p2, height = 15, width =8)