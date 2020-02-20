library(tidyr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(plyr)


rs = read.csv('rice_shoot_peroxide.csv', sep=";", na.strings=c("","NA"))
rr = read.csv('rice_root_peroxide.csv', sep=";", na.strings=c("","NA"))
ws = read.csv('wheat_shoot_peroxide.csv', sep=";", na.strings=c("","NA"))
wr = read.csv('wheat_root_peroxide.csv', sep=";", na.strings=c("","NA"))

pl_names <- c('rice','rice','wheat','wheat')
org_names <- c('root','shoot','root','shoot')
df_list <- list(rr,rs,wr,ws)

for (i in 1:length(list(rr,rs,wr,ws))){
  df_list[[i]]$plant=pl_names[i]
  df_list[[i]]$organ=org_names[i]
  df_list[[i]] <- gather(df_list[[i]], condition, measurement, control:Cu, factor_key=TRUE)
  df_list[[i]]$condition <- ordered(df_list[[i]]$condition,
                                    levels = c("control", "MV", "MD", "Per", "Cu"))
  df_list[[i]]$measurement <- as.numeric(gsub(",", ".", gsub("\\.", "", df_list[[i]]$measurement)))
}
merged_df <- ldply(df_list, data.frame)
merged_df$organ = factor(merged_df$organ, levels=c('shoot', 'root'))  

anno_df<-  compare_means(measurement ~ condition,group.by = c('plant','organ'), method ='wilcox.test' , merged_df) %>%
  mutate(y_pos = 0.6)

pl_df<-  compare_means(measurement ~ plant,group.by = c('condition','organ'), method ='wilcox.test' , merged_df) %>%
  mutate(y_pos = 0.6)

merged_df <- ldply(df_list, data.frame)
merged_df$organ = factor(merged_df$organ, levels=c('shoot', 'root'))  
merged_df$plant = factor(merged_df$plant)
merged_df$condition

t <- ddply(merged_df, .(plant, organ, condition),summarise,mean=mean(measurement, na.rm = T),ms=sd(measurement, na.rm = T)/sqrt(length(measurement)))
pl_df<-  compare_means(measurement ~ plant,group.by = c('condition','organ'), method ='wilcox.test' , merged_df) %>%
  mutate(y_pos = 0.6)
write.xlsx(t, 'peroxide.xlsx', sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

p <- ggbarplot(merged_df, x = "condition", 
               y = "measurement", 
               add = "mean_se",
               fill='condition',
               facet.by = c("organ",'plant'),
               position = position_dodge(0),
               width=1)+ 
  scale_y_continuous(expand = expand_scale(mult = c(0, .1)))

facet(p + theme_bw(),facet.by = c("organ",'plant'),
      panel.labs =  list(plant = c('рис', 'пшеница'), organ = c('побег', 'корень')),
      panel.labs.background = list(color = 'black',fill = 'white'),
      panel.labs.font = list(color ='black',size=21),
      scales='free_y')+
  stat_pvalue_manual(anno_df[c(3:4,33,21:23),],y.position = c(1105,200,800,480,660,1390), label = "p.signif",remove.bracket = T)+
  scale_fill_manual(values = c("white", "grey80", 'grey50','grey27','grey0'),labels=c(" контроль"," метилвиологен"," менадион"," перекись водорода"," аскорбат+Cu(II)"))+
  ylab(bquote('Содержание '*~H[2]~O[2]*', нмоль/г. сыр. массы'))+
  xlab(NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(color='black'),
        axis.title.y=element_text(size=17),
        legend.position = c(0.28, 0.89),
        legend.text = element_text(colour="black", size=18, face="bold"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor =  element_blank())