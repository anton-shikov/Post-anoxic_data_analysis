library(tidyr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(plyr)
setwd("C:/Users/HP/Desktop/перекись/All_data")

std <- function(x) sd(x,na.rm = T)/sqrt(length(x))
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


rs = read.csv('rice_shoot_peroxide_agents.csv', sep=";", na.strings=c("","NA"))
rr = read.csv('rice_root_peroxide_agents.csv', sep=";", na.strings=c("","NA"))
ws = read.csv('wheat_shoot_peroxide_agents.csv', sep=";", na.strings=c("","NA"))
wr = read.csv('wheat_root_peroxide_agents.csv', sep=";", na.strings=c("","NA"))

pl_names <- c('rice','rice','wheat','wheat')
org_names <- c('root','shoot','root','shoot')
df_list <- list(rr,rs,wr,ws)

for (i in 1:length(list(rr,rs,wr,ws))){
  df_list[[i]]$plant=pl_names[i]
  df_list[[i]]$organ=org_names[i]
  df_list[[i]] <- gather(df_list[[i]], condition, measurement, control:Per, factor_key=TRUE)
  df_list[[i]]$condition <- ordered(df_list[[i]]$condition,
                                    levels = c("control", "MV", "MD", "Per"))
  df_list[[i]]$measurement <- as.numeric(df_list[[i]]$measurement)
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
merged_df$dummy = factor(paste(merged_df$organ, merged_df$plant, merged_df$condition, sep='_'))
merged_df$no_cont_dum = factor(paste(merged_df$organ, merged_df$plant, sep='_'))
merged_df_percent = merged_df %>% filter(condition!='control')
merged_df_percent_let = merged_df

detach(package:plyr)

kruskal.test(measurement ~ condition, data = merged_df)

contr_df <- merged_df %>% filter(condition=='control') %>% group_by(no_cont_dum) %>%
  summarise(mean = mean(measurement)) %>% as.data.frame()


for (i in 1:nrow(merged_df_percent)){
  merged_df_percent[i,4] <- merged_df_percent[i,4]/contr_df[which(contr_df$no_cont_dum==merged_df_percent[i,6]),2]*100
}

for (i in 1:nrow(merged_df_percent_let)){
  merged_df_percent_let[i,4] <- merged_df_percent_let[i,4]/contr_df[which(contr_df$no_cont_dum==merged_df_percent_let[i,6]),2]*100
}

library(plyr)
t <- ddply(merged_df, .(plant, organ, condition),summarise,mean=mean(measurement, na.rm = T),ms=sd(measurement, na.rm = T)/sqrt(length(measurement)))
pl_df<-  compare_means(measurement ~ plant,group.by = c('condition','organ'), method ='wilcox.test' , merged_df) %>%
  mutate(y_pos = 0.6)



anno_df<-  compare_means(measurement ~ condition,group.by = c('plant','organ'), method ='wilcox.test' , merged_df) %>%
  mutate(y_pos = 0.6)

anno_df_dum<-  compare_means(measurement ~ dummy, method ='wilcox.test' , merged_df) %>%
  mutate(y_pos = 0.6)

anno_df_percent <-  compare_means(measurement ~ condition,group.by = c('plant','organ'),method ='wilcox.test' , merged_df_percent_let) %>%
  mutate(y_pos = 0.6)

anno_df_percent$factor <- paste(paste(anno_df_percent$plant,anno_df_percent$organ))

dif_dum <-as.numeric(anno_df_dum$p.format)
names(dif_dum) <- paste(anno_df_dum$group1, anno_df_dum$group2, sep="-")

(dif_dum.T <- multcompLetters(dif_dum))

#dif_dum_percent <-as.numeric(anno_df_percent$p.format)
#names(dif_dum_percent) <- paste(anno_df_percent$group1, anno_df_percent$group2, sep="-")

#(dif_dum_percent.T <- multcompLetters(dif_dum_percent))

dif1 <-as.numeric(anno_df[c(1:6),]$p.format)
names(dif1) <- paste(anno_df[c(1:6),]$group1, anno_df[c(1:6),]$group2, sep="-")
(diff1.T <- multcompLetters(dif1))
dif2 <-as.numeric(anno_df[c(7:12),]$p.format)
names(dif2) <- paste(anno_df[c(7:12),]$group1, anno_df[c(7:12),]$group2, sep="-")
(diff2.T <- multcompLetters(dif2))
dif3 <-as.numeric(anno_df[c(13:18),]$p.format)
names(dif3) <- paste(anno_df[c(13:18),]$group1, anno_df[c(13:18),]$group2, sep="-")
(diff3.T <- multcompLetters(dif3))
dif4 <-as.numeric(anno_df[c(19:24),]$p.format)
names(dif4) <- paste(anno_df[c(19:24),]$group1, anno_df[c(19:24),]$group2, sep="-")
(diff4.T <- multcompLetters(dif4))
?cld 

detach(package:Rmisc)
detach(package:plyr)
library(dplyr)


ag_df <- merged_df %>% 
  group_by(plant,organ, condition) %>% 
  summarise(mean = mean(measurement,na.rm = T),
            se=std(measurement),
            y_max=mean(measurement,na.rm = T)+std(measurement))


ag_df_percent <- merged_df_percent %>% 
  group_by(plant,organ, condition) %>% 
  summarise(mean = mean(measurement,na.rm = T),
            se=std(measurement),
            y_max=mean(measurement,na.rm = T)+std(measurement))

ag_df$letters <- 'letter'
ag_df_percent$letters <- 'letter'
ag_df[c(1:4),7] <- diff2.T
ag_df[c(5:8),7] <- diff1.T
ag_df[c(9:12),7] <- diff4.T
ag_df[c(13:16),7] <- diff3.T
ag_df$for_dum <- paste(as.character(ag_df$organ),as.character(ag_df$plant),as.character(ag_df$condition),sep='_')
ag_df_percent$for_dum <- paste(as.character(ag_df_percent$organ),as.character(ag_df_percent$plant),as.character(ag_df_percent$condition),sep='_')



for (i in 1:length(dif_dum.T$Letters)){
  ag_df[which(ag_df$for_dum==names(dif_dum.T$Letters[i])),7] <- as.character(dif_dum.T$Letters[i])
}


#for (i in 1:length(dif_dum_percent.T$Letters)){
#  ag_df_percent[which(ag_df_percent$for_dum==names(dif_dum_percent.T$Letters[i])),7] <- as.character(dif_dum_percent.T$Letters[i])
#}


ag_df$factor=paste(ag_df$plant,ag_df$organ)
merged_df$factor=paste(merged_df$plant,merged_df$organ)

anno_df_percent=as.data.frame(anno_df_percent)
ag_df_percent=as.data.frame(ag_df_percent)

ag_df_percent$factor=paste(ag_df_percent$plant,ag_df_percent$organ)
merged_df_percent$factor=paste(merged_df_percent$plant,merged_df_percent$organ)



for (i in 1:nrow(ag_df_percent)){
  sub1=anno_df_percent[which(anno_df_percent[,12]==as.vector(ag_df_percent[i,9] )),]
  sub2=sub1[which(sub1$group1=='control' & sub1$group2==ag_df_percent[i,3]),9]
  ag_df_percent[i,7] <- ifelse(sub2!='ns',sub2,'') 
}



df2 <- data_summary(merged_df, varname="measurement", 
                    groupnames=c("factor", "condition"))
#df2 <- data_summary(merged_df_percent, varname="measurement", 
#                   groupnames=c("factor", "condition"))



write.xlsx(t, 'peroxide.xlsx', sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)



#p <- ggbarplot(merged_df, x = "condition", 
#               y = "measurement", 
#               add = "mean_se",
#               fill='condition',
#               facet.by = c("organ",'plant'),
#               position = position_dodge(0),
#               width=1)+ 
#  scale_y_continuous(expand = expand_scale(mult = c(0, .1)))


p <- ggplot(df2, aes(fill=condition, y=measurement, x=factor)) + 
  geom_bar(position="dodge", stat = "summary", fun.y = "mean", color="black")+
  geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=.2,
                position=position_dodge(.9)) + theme_bw()+ 
  scale_y_continuous(limits=c(95,320),expand = expand_scale(mult = c(0, .1)),oob = rescale_none)

p  + scale_x_discrete(labels=c("Рис корень", "Рис побег", "Пшеница корень",'Пшеница побег'))+
scale_fill_manual(values = c("white", "grey80", 'grey50','grey27'),labels=c(" Контроль"," Метилвиологен"," Менадион"," Перекись водорода"))+
  ylab(bquote('Содержание '*~H[2]~O[2]*', нмоль/г'))+
  xlab(NULL)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(color='black',size=25),
        axis.text.y=element_text(color='black',size=25),
        axis.title.y=element_text(size=28),
        legend.position = c(0.45, 0.89),
        legend.text = element_text(colour="black", size=28, face="bold"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor =  element_blank())+ 
  geom_text(data =ag_df, 
            aes(x=factor, y=y_max+25,label=letters), 
            size = 5,position=position_dodge(0.9))

p  + scale_x_discrete(labels=c("Рис корень", "Рис побег", "Пшеница корень",'Пшеница побег'))+
  scale_fill_manual(values = c( "grey80", 'grey50','grey27'),labels=c(" Метилвиологен"," Менадион"," Перекись водорода"))+
  ylab(bquote('Содержание '*~H[2]~O[2]*', % от контроля'))+
  xlab(NULL)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(color='black',size=25),
        axis.text.y=element_text(color='black',size=25),
        axis.title.y=element_text(size=28),
        legend.position = c(0.45, 0.89),
        legend.text = element_text(colour="black", size=28, face="bold"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor =  element_blank())+ 
  geom_text(data =ag_df_percent, 
            aes(x=factor, y=y_max+3,label=letters), 
            size = 5,position=position_dodge(0.9))




library(plyr)
rrs = read.csv('rice_shoot_peroxide_reaeration.csv', sep=";", na.strings=c("","NA"))
rrr = read.csv('rice_root_peroxide_reaeration.csv', sep=";", na.strings=c("","NA"))
rws = read.csv('wheat_shoot_peroxide_reaeration.csv', sep=";", na.strings=c("","NA"))
rwr = read.csv('wheat_root_peroxide_reaeration.csv', sep=";", na.strings=c("","NA"))


rpl_names <- c('rice','rice','wheat','wheat')
rorg_names <- c('root','shoot','root','shoot')
rdf_list <- list(rrr,rrs,rwr,rws)

for (i in 1:length(list(rrr,rrs,rwr,rws))){
  rdf_list[[i]]$plant=pl_names[i]
  rdf_list[[i]]$organ=org_names[i]
  rdf_list[[i]] <- gather(rdf_list[[i]], condition, measurement, control:X24h, factor_key=TRUE)
  rdf_list[[i]]$condition <- ordered(rdf_list[[i]]$condition,
                                     levels = c("control", "A","X15min","X1h", "X3h", "X6h", "X24h"))
  rdf_list[[i]]$measurement <- as.numeric(rdf_list[[i]]$measurement)
}


rmerged_df <- ldply(rdf_list, data.frame)
rmerged_df$organ = factor(rmerged_df$organ, levels=c('shoot', 'root'))  
rmerged_df$plant = factor(lmerged_df$plant)
rmerged_df$condition
rmerged_df$dummy = factor(paste(rmerged_df$organ, rmerged_df$plant, rmerged_df$condition, sep='_'))
rmerged_df$no_cont_dum = factor(paste(rmerged_df$organ, rmerged_df$plant, sep='_'))

rmerged_df_percent = rmerged_df %>% filter(condition!='control')
rmerged_df_percent_let = rmerged_df

detach(package:Rmisc)
detach(package:plyr)
library(dplyr)

rcontr_df <- rmerged_df %>% filter(condition=='control') %>% group_by(no_cont_dum) %>%
  summarise(mean = mean(measurement,na.rm = T)) %>% as.data.frame()


for (i in 1:nrow(rmerged_df_percent)){
  rmerged_df_percent[i,4] <- rmerged_df_percent[i,4]/rcontr_df[which(rcontr_df$no_cont_dum==rmerged_df_percent[i,6]),2]*100
}

for (i in 1:nrow(rmerged_df_percent_let)){
  rmerged_df_percent_let[i,4] <- rmerged_df_percent_let[i,4]/rcontr_df[which(rcontr_df$no_cont_dum==rmerged_df_percent_let[i,6]),2]*100
}



ranno_df<-  compare_means(measurement ~ condition,group.by = c('plant','organ'), method ='wilcox.test' , rmerged_df) %>%
  mutate(y_pos = 0.6)

ranno_df_dum<-  compare_means(measurement ~ dummy, method ='wilcox.test',p.adjust.method = "bonf" , rmerged_df) %>%
  mutate(y_pos = 0.6)


ranno_df_percent <-  compare_means(measurement ~ condition,group.by = c('plant','organ'),method ='wilcox.test' , rmerged_df_percent_let) %>%
  mutate(y_pos = 0.6)

ranno_df_percent$factor <- paste(paste(ranno_df_percent$plant,ranno_df_percent$organ))


dif_dum <-as.numeric(str_replace_all(ranno_df_dum$p.format, "<", ""))
names(dif_dum) <- paste(ranno_df_dum$group1, ranno_df_dum$group2, sep="-")

(dif_dum.T <- multcompLetters(dif_dum,threshold = 0.05))


dif1 <-as.numeric(ranno_df[c(1:21),]$p.format)
names(dif1) <- paste(ranno_df[c(1:21),]$group1, ranno_df[c(1:21),]$group2, sep="-")
(diff1.T <- multcompLetters(dif1))
dif2 <-as.numeric(ranno_df[c(22:42),]$p.format)
names(dif2) <- paste(ranno_df[c(22:42),]$group1, ranno_df[c(22:42),]$group2, sep="-")
(diff2.T <- multcompLetters(dif2))
dif3 <-as.numeric(ranno_df[c(43:63),]$p.format)
names(dif3) <- paste(ranno_df[c(43:63),]$group1, ranno_df[c(43:63),]$group2, sep="-")
(diff3.T <- multcompLetters(dif3))
dif4 <-as.numeric(ranno_df[c(64:84),]$p.format)
names(dif4) <- paste(ranno_df[c(64:84),]$group1, ranno_df[c(64:84),]$group2, sep="-")
(diff4.T <- multcompLetters(dif4))
detach(package:Rmisc)
detach(package:plyr)
library(dplyr)
rag_df <- rmerged_df %>% 
  group_by(plant,organ, condition) %>% 
  summarise(mean = mean(measurement,na.rm = T),
            se=std(measurement),
            y_max=mean(measurement,na.rm = T)+std(measurement))


rag_df_percent <- rmerged_df_percent %>% 
  group_by(plant,organ, condition) %>% 
  summarise(mean = mean(measurement,na.rm = T),
            se=std(measurement),
            y_max=mean(measurement,na.rm = T)+std(measurement))

rag_df$letters <- 'letter'
rag_df_percent$letters <- 'letter'
rag_df[c(1:7),7] <- diff2.T
rag_df[c(8:14),7] <- diff1.T
rag_df[c(15:21),7] <- diff4.T
rag_df[c(22:28),7] <- diff3.T



rag_df$factor=paste(rag_df$plant,rag_df$organ)
rmerged_df$factor=paste(rmerged_df$plant,rmerged_df$organ)



ranno_df_percent=as.data.frame(ranno_df_percent)
rag_df_percent=as.data.frame(rag_df_percent)

rag_df_percent$factor=paste(rag_df_percent$plant,rag_df_percent$organ)
rmerged_df_percent$factor=paste(rmerged_df_percent$plant,rmerged_df_percent$organ)


for (i in 1:nrow(rag_df_percent)){
  sub1=ranno_df_percent[which(ranno_df_percent[,12]==as.vector(rag_df_percent[i,8] )),]
  sub2=sub1[which(sub1$group1=='control' & sub1$group2==rag_df_percent[i,3]),9]
  rag_df_percent[i,7] <- ifelse(sub2!='ns',sub2,'') 
}


df4 <- data_summary(rmerged_df, varname="measurement", 
                    groupnames=c("factor", "condition"))



rp <- ggplot(df4, aes(fill=condition, y=measurement, x=factor)) + 
  geom_bar(position="dodge", stat = "summary", fun.y = "mean", color="black")+
  geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd), width=.2,
                position=position_dodge(.9)) + theme_bw()+ 
  scale_y_continuous(limits=c(0,945),expand = expand_scale(mult = c(0, .1)),oob = rescale_none)


#rag_df[8,]$y_max <- rag_df[8,]$y_max-0.005
#rag_df[4,]$y_max <- rag_df[4,]$y_max-0.005

CairoPDF("plot3_mearged",width = 12.5, height =8)
rp +
  scale_fill_manual(values = c("white", "grey93", 'grey75','grey57','grey33','grey25', 'grey0'),labels=c(" Control"," Anoxia"," 15min", "1h"," 3h"," 6h", " 24h"))+
  ylab(bquote('Protein carbonylation level ('*~OD[366]~ '/'~OD[280]*')'))+
  xlab(NULL)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(color='black',size=20),
        axis.text.y=element_text(color='black',size=20),
        axis.title.y=element_text(size=24),
        legend.position = c(0.145, 0.86),
        legend.text = element_text(colour="black", size=23, face="bold"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor =  element_blank())+ 
  geom_text(data =rag_df, aes(x=factor, y=y_max+15,label=letters), 
            size = 5,position=position_dodge(0.9))
dev.off()


rp + scale_x_discrete(labels=c("Рис корень", "Рис побег", "Пшеница корень",'Пшеница побег'))+
  scale_fill_manual(values = c(  "grey93", 'grey75','grey57','grey33', 'grey0'),labels=c(" Аноксия"," 1ч"," 3ч"," 6ч"," 24ч"))+
  ylab(bquote('Уровень карбонилирования белков, % от контроля'))+
  xlab(NULL)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(color='black',size=25),
        axis.text.y=element_text(color='black',size=25),
        axis.title.y=element_text(size=28),
        legend.position = c(0.17, 0.89),
        legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor =  element_blank())+ 
  geom_text(data =rag_df_percent, 
            aes(x=factor, y=y_max+2,label=letters), 
            size = 5,position=position_dodge(0.9))+geom_hline(yintercept=100, linetype="dashed", color = "black",size =1)

