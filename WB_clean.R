library(tidyr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(plyr)
library("xlsx")
library(Cairo)
library(multcomp)
library(multcompView)
library(stringr)


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = std(x[[col]]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
} 

setwd("C:/Users/HP/Desktop/повреждения липидов и белков в окислительном стрессе/data/WB")
rs = read.csv('rice_shoot_WB.csv', sep=";", na.strings=c("","NA"))
rr = read.csv('rice_root_WB.csv', sep=";", na.strings=c("","NA"))
ws = read.csv('wheat_shoot_WB.csv', sep=";", na.strings=c("","NA"))
wr = read.csv('wheat_root_WB.csv', sep=";", na.strings=c("","NA"))

pl_names <- c('rice','rice','wheat','wheat')
org_names <- c('root','shoot','root','shoot')
df_list <- list(rr,rs,wr,ws)

for (i in 1:length(list(rr,rs,wr,ws))){
  df_list[[i]]$plant=pl_names[i]
  df_list[[i]]$organ=org_names[i]
  df_list[[i]] <- gather(df_list[[i]], condition, measurement, control:Cu, factor_key=TRUE)
  df_list[[i]]$condition <- ordered(df_list[[i]]$condition,
                                    levels = c("control","anoxia", "X24", "MV", "MD", "Per", "Cu"))
  df_list[[i]]$measurement <- as.numeric(gsub(",", ".", gsub("\\.", "", df_list[[i]]$measurement)))
}
merged_df <- ldply(df_list, data.frame)
merged_df$organ = factor(merged_df$organ, levels=c('shoot', 'root'))  
merged_df$plant = factor(merged_df$plant)
merged_df$condition
merged_df$dummy = factor(paste(merged_df$organ, merged_df$plant, merged_df$condition, sep='_'))



t <- ddply(merged_df, .(plant, organ, condition),summarise,mean=mean(measurement, na.rm = T),ms=sd(measurement, na.rm = T)/sqrt(length(measurement)))


anno_df<-  compare_means(measurement ~ condition,group.by = c('plant','organ'), method ='wilcox.test' , merged_df) %>%
  mutate(y_pos = 0.6)
anno_df_dum<-  compare_means(measurement ~ dummy, method ='wilcox.test' , merged_df) %>%
  mutate(y_pos = 0.6)

dif_dum <-as.numeric(anno_df_dum$p.format)
names(dif_dum) <- paste(anno_df_dum$group1, anno_df_dum$group2, sep="-")

  (dif_dum.T <- multcompLetters(dif_dum,Letters=c("a", "b", 'c','d','e','f','g', 'h', 'i','j','k', 'l', 'm')))
  

dif1 <-as.numeric(anno_df[c(1:21),]$p.format)
names(dif1) <- paste(anno_df[c(1:21),]$group1, anno_df[c(1:21),]$group2, sep="-")
(diff1.T <- multcompLetters(dif1))
dif2 <-as.numeric(anno_df[c(22:42),]$p.format)
names(dif2) <- paste(anno_df[c(22:42),]$group1, anno_df[c(22:42),]$group2, sep="-")
(diff2.T <- multcompLetters(dif2))
dif3 <-as.numeric(anno_df[c(43:63),]$p.format)
names(dif3) <- paste(anno_df[c(43:63),]$group1, anno_df[c(43:63),]$group2, sep="-")
(diff3.T <- multcompLetters(dif3))
dif4 <-as.numeric(anno_df[c(44:82),]$p.format)
names(dif4) <- paste(anno_df[c(44:82),]$group1, anno_df[c(44:82),]$group2, sep="-")
(diff4.T <- multcompLetters(dif4))


std <- function(x) sd(x,na.rm = T)/sqrt(length(x))
detach(package:Rmisc)
detach(package:plyr)
library(dplyr)

ag_df <- merged_df %>% 
  group_by(plant,organ, condition) %>% 
  summarise(mean = mean(measurement,na.rm = T),
            se=std(measurement),
            y_max=mean(measurement,na.rm = T)+std(measurement))
ag_df$letters <- 'letter'
ag_df[c(1:7),7] <- diff2.T
ag_df[c(8:14),7] <- diff1.T
ag_df[c(15:21),7] <- diff4.T
ag_df[c(22:28),7] <- diff3.T

ag_df$for_dum <- paste(as.character(ag_df$organ),as.character(ag_df$plant),as.character(ag_df$condition),sep='_')

for (i in 1:length(dif_dum.T$Letters)){
  ag_df[which(ag_df$for_dum==names(dif_dum.T$Letters[i])),7] <- as.character(dif_dum.T$Letters[i])
}
ag_df$factor=paste(ag_df$plant,ag_df$organ)
merged_df$factor=paste(merged_df$plant,merged_df$organ)


pl_df<-  compare_means(measurement ~ plant,group.by = c('condition','organ'), method ='wilcox.test' , merged_df) %>%
  mutate(y_pos = 0.6)




df2 <- data_summary(merged_df, varname="measurement", 
                    groupnames=c("factor", "condition"))


p <- ggplot(df2, aes(fill=condition, y=measurement, x=factor)) + 
  geom_bar(position="dodge", stat = "summary", fun.y = "mean", color="black")+
  geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=.2,
                position=position_dodge(.9)) + theme_bw()+ 
  scale_y_continuous(expand = expand_scale(mult = c(0, .1)))


ag_df_copy=ag_df


ag_df=ag_df_copy
ag_df[11,]$y_max=ag_df[11,]$y_max-0.015
ag_df[12,]$y_max=ag_df[12,]$y_max+0.015
ag_df[13,]$y_max=ag_df[13,]$y_max-0.016
ag_df[19,]$y_max=ag_df[19,]$y_max-0.014
ag_df[20,]$y_max=ag_df[20,]$y_max+0.011
CairoPDF("WB_new_mearged.pdf",width = 13.5, height =8)
p + scale_x_discrete(labels=c("Рис корень", "Рис побег", "Пшеница корень",'Пшеница побег'))+scale_fill_manual(values = c("white",'grey90' ,"grey76", 'grey60','grey47','grey20','grey5'),labels=c(" Контроль"," Аноксия", " Реаэрация"," Метилвиологен"," Менадион"," Перекись водорода"," Аскорбат +Cu(II)"))+
  ylab(bquote('Уровень карбонилирования белков'))+
  xlab(NULL)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(color='black',size=20),
        axis.text.y=element_text(color='black',size=20),
        axis.title.y=element_text(size=24),
        legend.position = c(0.35, 0.84),
        legend.text = element_text(colour="black", size=23, face="bold"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor =  element_blank())+ 
  geom_text(data =ag_df, 
            aes(x=factor, y=y_max+0.09,label=letters), 
            size = 5,position=position_dodge(0.9))
dev.off()





