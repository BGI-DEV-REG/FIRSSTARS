library(tidyverse)
library(ggdark)
.libPaths(c("E:/CodePrograms/SeuratV5", .libPaths()))
library(Seurat)
packageVersion("Matrix")
packageVersion("Seurat")

first_FB_l = c("Ccl3","Il6","Sirpa",
               "Cxcl2","Cxcl3","Il1b",  
               "Cxcl10","Cxcl16","Ccl7")

first_Mono_r = c("Ccr5","Il6ra","Il6st",'Cd47',
                 "Ccr1l1","Ccr2","Ccr3","Ackr1","Ackr2",
                 "Cxcr3","Cxcr6","Cxcr1","Cxcr2","Il1r1","Il1rap")

lr_gene <- rev(c(first_FB_l, first_Mono_r))

cell_less <- c("FB_adipogenic_Adult","FB_basement_activated_Adult","FB_basement_resident_Adult","FB_DS_Adult","FB_papi_Adult")

APE_Early <- readRDS("E:/CodePrograms/Immu/Immu_0708/Immu&FB_APE_Early.rds")

celltype_phase_level <- c("FB_papi_E","FB_basement_activated_E","FB_basement_resident_E","FB_DP_E","FB_DS_E","FB_exosome_E",
                          "FB_adipogenic_E","FB_fascia_E","FB_inflammatory_E","FB_muscle_E","FB_oxidative_stress_E","FB_reti_E",
                          "FB_papi_P5","FB_basement_activated_P5","FB_basement_resident_P5","FB_DP_P5","FB_DS_P5",
                          "FB_adipogenic_P5","FB_fascia_P5","FB_inflammatory_P5","FB_muscle_P5","FB_oxidative_stress_P5","FB_reti_P5",
                          "FB_papi_Adult","FB_basement_activated_Adult","FB_basement_resident_Adult","FB_DS_Adult","FB_reti_Adult",
                          "FB_apoptotic_Adult","FB_adipogenic_Adult","FB_exosome_Adult","FB_muscle_Adult",
                          "FB_fascia_Adult","FB_fascia_2_Adult","FB_inflammatory_Adult","FB_myo_Adult","FB_oxidative_stress_Adult",
                          "Macr_resident_E","Mono_classical_E","Mono_intermediate_E","Mono_noclassical_E","Macr_neovascularization_E","Macr_remodeling_E",
                          "Macr_resident_P5","Mono_classical_P5","Mono_intermediate_P5","Mono_noclassical_P5","Macr_neovascularization_P5","Macr_remodeling_P5",
                          "Macr_resident_Adult","Mono_classical_Adult","Mono_intermediate_Adult","Mono_noclassical_Adult","Macr_neovascularization_Adult","Macr_remodeling_Adult")

FB_Early_names <-  unique(grep("^FB" ,APE_Early$celltype, value = T))
APE_FB <- subset(APE_Early, celltype %in% FB_Early_names)
APE_FB <- subset(APE_FB, celltype_phase %in% cell_less, invert = TRUE)

APE_FB$celltype_phase <- factor(APE_FB$celltype_phase, rev(celltype_phase_level))

p_first_bubble_FB  <- DotPlot(APE_FB,
                              features = first_FB_l, 
                              group.by = "celltype_phase", 
                              assay = "RNA") +
  theme(
    axis.text.x = element_text( hjust = 1, angle = 60)
  ) +
  scale_color_gradientn(colors = turbo(10)) +
  dark_theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, color = 'white', size = 10),
        axis.text.y = element_text(color = "white", size = 10),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 6, color = "white"),
        legend.title = element_text(size = 8, color = "white"),
        legend.key.size = unit(0.3, "cm"))
ggsave("fig5_first_FB_bubble_0715.pdf", p_first_bubble_FB, w = 5, h = 8)


Mono_Early_names <-  unique(grep("^M" ,APE_Early$celltype, value = T))
APE_Mono <- subset(APE_Early, celltype %in% Mono_Early_names)

p_first_bubble_Mono  <- DotPlot(APE_Mono,
                                features = first_Mono_r, 
                                group.by = "celltype_phase", 
                                assay = "RNA") +
  coord_flip() +
  theme(
    axis.text.x = element_text( hjust = 1, angle = 60)
  ) +
  scale_color_gradientn(colors = turbo(10)) +
  dark_theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, color = 'white', size = 8),
        axis.text.y = element_text(color = "white", size = 8),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 6, color = "white"),
        legend.title = element_text(size = 8, color = "white"),
        legend.key.size = unit(0.3, "cm")) +
  scale_size_continuous(range = c(1, 5))
ggsave("fig5_first_Mono_bubble.pdf", p_first_bubble_Mono, w = 4.8, h = 4.2)