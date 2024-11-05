library(ggplot2)
library(ggsignif)
library(reshape2)

df = read.csv("adult_ep_correlation.csv",row.names=1)

sub_df = df[grepl("Adult",rownames(df)),]
sub_df = sub_df[,!(grepl("Adult",colnames(sub_df)))]
sub_df$celltype = sub("__Adult","",rownames(sub_df))
da = data.frame(stage=colnames(sub_df)[1:length(colnames(sub_df))-1])
da$celltype=sapply(da$stage,function(x){strsplit(x,"__")[[1]][1]})
da$value = 0
for (i in seq(1,nrow(da))){
    da[i,"value"]=sub_df[sub_df$celltype == da[i,"celltype"],da[i,"stage"]]
}



co = read.csv("Fig1_color_seq_20240512.csv")
colors = co$color
names(colors)=co$celltype


## 进行统计检验
library(ggpubr)
stat_test <-compare_means(value ~ celltype,  data = da, ref.group = ".all.",method = "t.test", alternative = "less")
# 过滤掉非显著性结果
stat_test <- stat_test[stat_test$p.signif != "ns",]
p <- ggboxplot(da, x = "celltype", y = "value", fill = "celltype", color="black") +scale_fill_manual(values =  colors) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 添加显著性结果
p <- p + stat_pvalue_manual(stat_test, label = "p.signif", 
                            y.position = max(da$value)*1.05 ,size = 8,  # 调整显著性符号的字体大小
                            label.size = 8,x="group2")  # 调整 y 位置以避免重叠
    

ggsave(plot=p,"correlation_EP_split_stat.pdf",w=22,h=10) 
