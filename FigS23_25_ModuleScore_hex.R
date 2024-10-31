library(Seurat)
library(ggplot2)
library(stringr)
library(tidyverse)

### Get the parameters
parser <- argparse::ArgumentParser()
parser$add_argument('-I','--input', help='input path')
parser$add_argument('-N','--name',help='prefix name')
parser$add_argument('-O','--out', help='Output directory')
parser$add_argument('-G','--genes', nargs='+', help='Genes')
parser$add_argument('--binwidth', type="integer", default=500, help='binwidth[default %(default)s]')

args <- parser$parse_args()
print(args)
print(Sys.time())
dir.create(args$out)
setwd(args$out)


obj <- readRDS(args$input)

obj <- AddModuleScore(obj, list("genes" = args$genes), name = "ModuleScore") 
obj$ModuleScore = obj$ModuleScore1

obj$coor_y = obj@images[[names(obj@images)]]@coordinates$row
obj$coor_x = obj@images[[names(obj@images)]]@coordinates$col
results_df <- data.frame(x = obj$coor_x, y = obj$coor_y, ModuleScore=obj$ModuleScore)
write.csv(results_df, gzfile(paste0(args$name, "_ModuleScore.csv.gz")))



p <- ggplot(results_df, aes(x=x,y=y,z=ModuleScore)) + 
    stat_summary_hex(binwidth = c(args$binwidth, args$binwidth), fun = mean) + 
    scale_fill_viridis_c() + 
    coord_fixed() + 
    labs(title=args$name) + 
    theme_void() +
    theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        text = element_text(color = "white"),
        plot.title = element_text(color = "white"))

ggsave(paste0(args$name, "_ModuleScore_hex.png"),height = 16,width = 16)

print(Sys.time())
print("Done!")
