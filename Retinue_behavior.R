##################################################################

library(ggplot2)
library(gridExtra)
library(lme4)
library(ggsignif)
library(gtools)
library(ggeffects)
library(cowplot)
library(ggprism)
library(grid)
library(dplyr)
library(tidyr)
library(nlme)
library(emojifont)
library(Cairo)

##################################################################

retinue <- read.csv("retinue_counts_2019.csv")
retinue$Proportion <- retinue$Observed/retinue$Background
retinue$Weight <- retinue$Observed+retinue$Background

retinue[retinue$Qline=="C1","Qline"] <- 0
retinue[retinue$Qline=="C3","Qline"] <- 1
retinue[retinue$Qline=="C2","Qline"] <- 0
retinue[retinue$Qline=="C4","Qline"] <- 1
retinue$Qline <- as.numeric(retinue$Qline)
retinue$Colony <- as.factor(retinue$Colony)
retinue$Block <- as.factor(retinue$Block)
retinue$Qline <- as.factor(retinue$Qline)

retinue.logit <- glm(formula=Proportion~Block/Qline,
                     weights=Weight,
                     data=retinue,
                     family="quasibinomial")
r.logit.p.b1 <- coef(summary(retinue.logit))[,4][3]
r.logit.p.b2 <- coef(summary(retinue.logit))[,4][4]
summary(retinue.logit)

rd2019.p.b1 <- data.frame(group1="C1",
                      group2="C3",
                      y.position=0.14,
                      p=as.numeric(r.logit.p.b1))
rd2019.p.b1$p <- "NS"

rd2019.p.b2 <- data.frame(group1="C2",
                          group2="C4",
                          y.position=0.14,
                          p=as.numeric(r.logit.p.b2))
rd2019.p.b2$p <- "NS"

retinue.b1 <- retinue[retinue$Block=="1",]
levels(retinue.b1$Qline) <- c(expression("A\u2640B\u2642"),
                              expression("B\u2640A\u2642"))
retinue.b2 <- retinue[retinue$Block=="2",]
levels(retinue.b2$Qline) <- c(expression("A\u2640B\u2642"),
                              expression("B\u2640A\u2642"))



retinue.b1.plot <- data.frame(Qline=retinue.b1$Qline[1:2],
                              Proportion=c(mean(retinue.b1[c(1,3),"Proportion"]),
                                           mean(retinue.b1[c(2,4),"Proportion"])),
                              std=c(sd(retinue.b1[c(1,3),"Proportion"]),
                                    sd(retinue.b1[c(2,4),"Proportion"])))

rd2019.b1 <- ggplot(data=retinue.b1.plot,
                    aes(x=Qline,y=Proportion,fill=Qline)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) +
  geom_linerange(aes(x=Qline, ymin=Proportion-std, ymax=Proportion+std)) +
  theme_prism() + theme(text = element_text(size=12),
                        axis.text.x=element_text(color="black"),
                        axis.text.y=element_text(color="black"),
                        axis.title.y=element_text(size=12),
                        legend.position="none",
                        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("") + ylab("Proportion of bees\nresponding to QMP") +
  scale_fill_manual(values=c("white", "grey70")) +
  scale_color_manual(values=c("white", "grey70")) +
  scale_y_continuous(labels = scales::percent,limits=c(0,0.15)) +
  ggtitle("Block 1")

retinue.b2.plot <- data.frame(Qline=retinue.b2$Qline[1:2],
                              Proportion=c(mean(retinue.b2[c(1,3),"Proportion"]),
                                           mean(retinue.b2[c(2,4),"Proportion"])),
                              std=c(sd(retinue.b2[c(1,3),"Proportion"]),
                                    sd(retinue.b2[c(2,4),"Proportion"])))

rd2019.b2 <- ggplot(data=retinue.b2.plot,
                    aes(x=Qline,y=Proportion,fill=Qline)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) +
  geom_linerange(aes(x=Qline, ymin=Proportion-std, ymax=Proportion+std)) +
  theme_prism() + theme(text = element_text(size=12),
                        axis.text.x=element_text(color="black"),
                        legend.position="none",
                        axis.text.y=element_text(color="white"),
                        axis.ticks.y=element_blank(),
                        axis.line.y=element_blank(),
                        plot.margin = unit(c(0,0,0,0), "cm")) + 
  xlab("") + ylab("") +
  scale_fill_manual(values=c("white", "grey70")) +
  scale_color_manual(values=c("white", "grey70")) +
  scale_y_continuous(labels = scales::percent, limits=c(0,0.15)) +
  ggtitle("Block 2")

##################################################################

retinue2 <- read.csv("retinue_counts_2021.csv")
retinue2$Proportion <- retinue2$Observed/retinue2$Background
retinue2$Weight <- retinue2$Observed+retinue2$Background

retinue2.A <- retinue2[retinue2$Block=="A",]
retinue2.A[retinue2.A$Qline=="C5","Qline"] <- 0
retinue2.A[retinue2.A$Qline=="C9","Qline"] <- 1
retinue2.A$Qline <- as.numeric(retinue2.A$Qline)

retinue2.B <- retinue2[retinue2$Block=="B",]
retinue2.B[retinue2.B$Qline=="C6","Qline"] <- 0
retinue2.B[retinue2.B$Qline=="C11","Qline"] <- 1
retinue2.B$Qline <- as.numeric(retinue2.B$Qline)

retinue2 <- rbind(retinue2.A,retinue2.B)
retinue2$Colony <- as.factor(retinue2$Colony)
retinue2$Block <- as.factor(retinue2$Block)
levels(retinue2$Block) <- c(1,2)
retinue2$Q.ID <- NULL
retinue2$Qline <- as.factor(retinue2$Qline)
retinue2$Trial <- NULL

retinue2.logit <- glm(formula=Proportion~Block/Qline,
                        weights=Weight,
                        data=retinue2,
                        family="quasibinomial")
summary(retinue2.logit)
r.logit.p1 <- coef(summary(retinue2.logit))[,4][3]
r.logit.p2 <- coef(summary(retinue2.logit))[,4][4]

r.p.b1 <- data.frame(group1="A\u2640B\u2642",
                     group2="B\u2640A\u2642",
                     y.position=0.1,
                     p=as.numeric(r.logit.p1))
r.p.b1$p <- stars.pval(r.p.b1$p)

r.p.b2 <- data.frame(group1="A\u2640B\u2642",
                     group2="B\u2640A\u2642",
                     y.position=0.13,
                     p=as.numeric(r.logit.p2))
r.p.b2$p <- stars.pval(r.p.b2$p)

retinue.b1 <- retinue2[retinue2$Block=="1",]
levels(retinue.b1$Qline) <- c(expression("A\u2640B\u2642"),
                              expression("B\u2640A\u2642"))
retinue.b2 <- retinue2[retinue2$Block=="2",]
levels(retinue.b2$Qline) <- c(expression("A\u2640B\u2642"),
                              expression("B\u2640A\u2642"))


rd2021.b1 <- ggplot(data=retinue.b1,
       aes(x=Qline,y=Proportion,fill=Qline)) +
  geom_boxplot(width=0.5,aes(fill=Qline)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize=.8) +
  theme_prism() + theme(text = element_text(size=12),
                        axis.text.x=element_text(color="black"),
                        axis.ticks.y=element_blank(),
                        axis.text.y=element_text(color="white"),
                        axis.line.y=element_blank(),
                        legend.position="none",
                        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("") + ylab("") +
  scale_color_manual(values=c("white", "grey70")) +
  scale_fill_manual(values=c("white", "grey70")) +
  scale_y_continuous(labels = scales::percent,limits=c(0,0.15)) +
  add_pvalue(r.p.b1,label.size=5,inherit.aes = FALSE) +
  ggtitle("Block 3")

rd2021.b2 <- ggplot(data=retinue.b2,
       aes(x=Qline,y=Proportion,fill=Qline)) +
  geom_boxplot(width=0.5,aes(fill=Qline)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize=.8) +
  theme_prism() + theme(text = element_text(size=12),
                        axis.text.x=element_text(color="black"),
                        legend.position="none",
                        axis.text.y=element_text(color="white"),
                        axis.ticks.y=element_blank(),
                        axis.line.y=element_blank(),
                        plot.margin = unit(c(0,0,0,0), "cm")) + 
  xlab("") + ylab("") +
  scale_fill_manual(values=c("white", "grey70")) +
  scale_color_manual(values=c("white", "grey70")) +
  scale_y_continuous(labels = scales::percent, limits=c(0,0.15)) +
  add_pvalue(r.p.b2,label.size=5,inherit.aes = FALSE) +
  ggtitle("Block 4")


rd <- plot_grid(rd2019.b1,rd2019.b2,rd2021.b1,rd2021.b2,ncol=4)
#  draw_label("Proportion of bees responding to QMP", 
#             x=0, y=0.5, vjust= 1.5, angle=90,fontface="bold",size=12)

cairo_pdf("retinue_2019_2021.pdf",width=9,height=3,
          fallback_resolution=1200)
grid.draw(rd)
dev.off()

##################################################################

data.temp <- read.csv("2022_IGC_behavior_analysis/Cordovan4_Caucasian21_retinue.csv")
data.temp[is.na(data.temp)] <- 0
data.temp$Bee <- paste(data.temp$Cross,
                       data.temp$Bee,sep=".")
data <- data.temp[,c(2:10)]
data1 <- data

data.long <- gather(data1, Trial, Score, Trial.1:Trial.3)
group_cols <- c("Lineage","ReciprocalCross","Cage","Cross","Bee","Trial")
data.AOV <- data.long[,c(1,2,3,4,6,7,8)] %>%
  group_by(across(all_of(group_cols))) %>%
  summarise_at(vars(Score),mean)
data.AOV <- data.frame(data.AOV)
group_cols <- c("Lineage","ReciprocalCross","Cross","Cage","Bee","Trial")
data.AOV[group_cols] <- lapply(data.AOV[group_cols], factor)
data.AOV$Cage.Test <- factor(paste(as.character(data.AOV$ReciprocalCross),
                                   as.character(data.AOV$Cage),sep="."))
levels(data.AOV$Lineage) <- c(0,1)
data.AOV1 <- data.AOV


data.temp <- read.csv("2022_IGC_behavior_analysis/Cordovan2_Weaver43_retinue.csv")
data.temp[is.na(data.temp)] <- 0
data.temp$Bee <- paste(data.temp$Cross,
                       data.temp$Bee,sep=".")
data <- data.temp[,c(2:10)]
data2 <- data

data.long <- gather(data2, Trial, Score, Trial.1:Trial.3)
group_cols <- c("Lineage","ReciprocalCross","Cage","Cross","Bee","Trial")
data.AOV <- data.long[,c(1,2,3,4,6,7,8)] %>%
  group_by(across(all_of(group_cols))) %>%
  summarise_at(vars(Score),mean)
data.AOV <- data.frame(data.AOV)
group_cols <- c("Lineage","ReciprocalCross","Cross","Cage","Bee","Trial")
data.AOV[group_cols] <- lapply(data.AOV[group_cols], factor)
data.AOV$Cage.Test <- factor(paste(as.character(data.AOV$ReciprocalCross),
                                   as.character(data.AOV$Cage),sep="."))
levels(data.AOV$Lineage) <- c(0,1)
data.AOV2 <- data.AOV

data.AOV <- rbind(data.AOV1,data.AOV2)
data.AOV$Block <- NA
data.AOV[data.AOV$ReciprocalCross%in%c("B1xB94","B3xB95","B5xB96"),"Block"] <- 0
data.AOV[data.AOV$ReciprocalCross%in%c("G50xY37","G52xY38","G53xY39"),"Block"] <- 1
data.AOV$Block <- factor(data.AOV$Block)

lme <- lme(Score~Block/Lineage,
           random=~1|Cage.Test,
           data=data.AOV)
s.lme <- summary(lme)
s.lme


data <- data1
metadata <- data[,c(1,2,3,4,6)]
metadata <- metadata[!duplicated(metadata),]
data.mean <- data.frame(matrix(ncol=3,nrow=0))
names(data.mean) <- c("Bee","Trial","mean")
for(i in 1:length(metadata$Bee)){
  if(metadata$Bee[i]%in%data$Bee){
    t1 <- sum(data[data$Bee==metadata$Bee[i],"Trial.1"])
    t2 <- sum(data[data$Bee==metadata$Bee[i],"Trial.2"])
    t3 <- sum(data[data$Bee==metadata$Bee[i],"Trial.3"])
    data.mean <- rbind(data.mean,data.frame(Bee=rep.int(metadata$Bee[i],3),
                                            Trial=c(1,2,3),
                                            mean=c(t1/5,t2/5,t3/5)))
  }else{
    data.mean <- rbind(data.mean,data.frame(Bee=rep.int(metadata$Bee[i],3),
                                            Trial=c(1,2,3),
                                            mean=rep.int(0,3)))
  }
}

se <- function(x){
  sd(x)/sqrt(length((x)))
}

data.plot <- data.frame(matrix(ncol=3,nrow=0))
names(data.plot) <- c("Bee","mean","se")
for(i in 1:length(metadata$Bee)){
  data.plot <- rbind(data.plot,data.frame(Bee=metadata$Bee[i],
                                          mean=mean(data.mean[data.mean$Bee==metadata$Bee[i],"mean"]),
                                          se=se(data.mean[data.mean$Bee==metadata$Bee[i],"mean"])))
}
data.plot <- data.plot[order(data.plot$mean,data.plot$se),]
data.plot <- left_join(data.plot,metadata,by="Bee")
data.plot1 <- data.plot


data <- data2
metadata <- data[,c(1,2,3,4,6)]
metadata <- metadata[!duplicated(metadata),]
data.mean <- data.frame(matrix(ncol=3,nrow=0))
names(data.mean) <- c("Bee","Trial","mean")
for(i in 1:length(metadata$Bee)){
  if(metadata$Bee[i]%in%data$Bee){
    t1 <- sum(data[data$Bee==metadata$Bee[i],"Trial.1"])
    t2 <- sum(data[data$Bee==metadata$Bee[i],"Trial.2"])
    t3 <- sum(data[data$Bee==metadata$Bee[i],"Trial.3"])
    data.mean <- rbind(data.mean,data.frame(Bee=rep.int(metadata$Bee[i],3),
                                            Trial=c(1,2,3),
                                            mean=c(t1/5,t2/5,t3/5)))
  }else{
    data.mean <- rbind(data.mean,data.frame(Bee=rep.int(metadata$Bee[i],3),
                                            Trial=c(1,2,3),
                                            mean=rep.int(0,3)))
  }
}

se <- function(x){
  sd(x)/sqrt(length((x)))
}

data.plot <- data.frame(matrix(ncol=3,nrow=0))
names(data.plot) <- c("Bee","mean","se")
for(i in 1:length(metadata$Bee)){
  data.plot <- rbind(data.plot,data.frame(Bee=metadata$Bee[i],
                                          mean=mean(data.mean[data.mean$Bee==metadata$Bee[i],"mean"]),
                                          se=se(data.mean[data.mean$Bee==metadata$Bee[i],"mean"])))
}
data.plot <- data.plot[order(data.plot$mean,data.plot$se),]
data.plot <- left_join(data.plot,metadata,by="Bee")
data.plot2 <- data.plot


data.plot <- rbind(data.plot1,data.plot2)
data.plot[data.plot$Lineage=="Cordovan2","Lineage"] <- "C7"
data.plot[data.plot$Lineage=="Weaver43","Lineage"] <- "C10"
data.plot[data.plot$Lineage=="Cordovan4","Lineage"] <- "C8"
data.plot[data.plot$Lineage=="Caucasian21","Lineage"] <- "C12"

data.plot.b1 <- data.plot[data.plot$Lineage%in%c("C7","C10"),]
data.plot.b1$Lineage <- factor(data.plot.b1$Lineage,
                               levels=c("C7","C10"))
levels(data.plot.b1$Lineage) <- c(expression("A\u2640B\u2642"),
                                  expression("B\u2640A\u2642"))
data.plot.b2 <- data.plot[data.plot$Lineage%in%c("C8","C12"),]
data.plot.b2$Lineage <- factor(data.plot.b2$Lineage,
                               levels=c("C8","C12"))
levels(data.plot.b2$Lineage) <- c(expression("A\u2640B\u2642"),
                                  expression("B\u2640A\u2642"))

r.logit.p2 <- s.lme[["tTable"]][3,5]
r.logit.p1 <- s.lme[["tTable"]][4,5]

r.p.b1 <- data.frame(group1="A\u2640B\u2642",
                     group2="B\u2640A\u2642",
                     y.position=0.95,
                     p=as.numeric(r.logit.p1))
r.p.b1$p <- stars.pval(r.p.b1$p)

r.p.b2 <- data.frame(group1="A\u2640B\u2642",
                     group2="B\u2640A\u2642",
                     y.position=0.95,
                     p=as.numeric(r.logit.p2))
r.p.b2$p <- stars.pval(r.p.b2$p)

rd2022.b1 <- ggplot(data=data.plot.b1,
       aes(x=Lineage,y=mean,fill=Lineage)) +
  geom_boxplot(width=0.5,aes(fill=Lineage)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize=.4) +
  theme_prism() + theme(text = element_text(size=12),
                        axis.text.x=element_text(color="black"),
                        axis.text.y=element_text(color="black"),
                        legend.position="none",
                        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("") + ylab("Frequency of\nresponse to QMP") +
  scale_color_manual(values=c("white","grey70")) +
  scale_fill_manual(values=c("white","grey70")) +
  scale_y_continuous(labels = scales::percent,limits=c(0,1)) +
  add_pvalue(r.p.b1,label.size=5,inherit.aes = FALSE) +
  ggtitle("Block 5")

rd2022.b2 <- ggplot(data=data.plot.b2,
       aes(x=Lineage,y=mean,fill=Lineage)) +
  geom_boxplot(width=0.5,aes(fill=Lineage)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize=.4) +
  theme_prism() + theme(text = element_text(size=12),
                        axis.text.x=element_text(color="black"),
                        legend.position="none",
                        axis.text.y=element_text(color="white"),
                        axis.ticks.y=element_blank(),
                        axis.line.y=element_blank(),
                        plot.margin = unit(c(0,0,0,0), "cm")) + 
  xlab("") + ylab("") +
  scale_fill_manual(values=c("white","grey70")) +
  scale_color_manual(values=c("white","grey70")) +
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  add_pvalue(r.p.b2,label.size=5,inherit.aes = FALSE) +
  ggtitle("Block 6")

rd <- plot_grid(rd2022.b1,rd2022.b2)

cairo_pdf("retinue_2022.pdf",width=9,height=3,
          fallback_resolution=1200)
grid.draw(rd)
dev.off()



p1 <- ggdraw() + draw_image(magick::image_read_pdf("retinue_2019_2021.pdf",density=1200),scale=0.9)
p2 <- ggdraw() + draw_image(magick::image_read_pdf("retinue_2022.pdf",density=1200),scale=0.9)
p.out <- plot_grid(p1, p2,labels=c("A","B"),ncol=1,nrow=2,label_size=200)
ggsave("retinue_AB.png",p.out,width=9,height=6,bg="white",dpi=1200)





# EXTRA


data.plot.b1.cage <- data.plot.b1
data.plot.b1.cage$Cage <- factor(data.plot.b1.cage$Cage,
                                 levels=c(1:6))
rd2022.b1.cage <- ggplot(data=data.plot.b1.cage,
                    aes(x=Cage,y=mean,fill=Lineage)) +
  geom_boxplot(width=0.5,
               position = position_dodge(width = 0.9)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize=.4,
               position = position_dodge(width = 0.9)) +
  theme_prism() + theme(text = element_text(size=12),
                        axis.text.x=element_text(color="black"),
                        axis.text.y=element_text(color="black"),
                        legend.position="none",
                        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("") + ylab("Frequency of response to QMP") +
  scale_color_manual(values=c("white","grey70")) +
  scale_fill_manual(values=c("white","grey70")) +
  scale_y_continuous(labels = scales::percent,limits=c(0,1)) +
  ggtitle("Block 5")


data.plot.b2.cage <- data.plot.b2
data.plot.b2.cage$Cage <- factor(data.plot.b2.cage$Cage,
                                 levels=c(1:6))
rd2022.b2.cage <- ggplot(data=data.plot.b2.cage,
                    aes(x=Cage,y=mean,fill=Lineage)) +
  geom_boxplot(width=0.5,
               position = position_dodge(width = 0.9)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize=.4,
               position = position_dodge(width = 0.9)) +
  theme_prism() + theme(text = element_text(size=12),
                        axis.text.x=element_text(color="black"),
                        legend.position="none",
                        axis.text.y=element_text(color="white"),
                        axis.ticks.y=element_blank(),
                        axis.line.y=element_blank(),
                        plot.margin = unit(c(0,0,0,0), "cm")) + 
  xlab("") + ylab("") +
  scale_fill_manual(values=c("white","grey70")) +
  scale_color_manual(values=c("white","grey70")) +
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  ggtitle("Block 6")

rd <- plot_grid(rd2022.b1.cage,rd2022.b2.cage)

cairo_pdf("figS1.pdf",width=10,height=5,
          fallback_resolution=1200)
grid.draw(rd)
dev.off()





data.plot.b1 <- data.plot.b1[order(data.plot.b1$mean),]
data.plot.b1$ID <- NA
data.plot.b1[data.plot.b1$Lineage=="A\u2640B\u2642","ID"] <- seq.int(1,180,2)
data.plot.b1[data.plot.b1$Lineage=="B\u2640A\u2642","ID"] <- seq.int(2,180,2)
data.plot.b1 <- data.plot.b1[order(data.plot.b1$ID),]
data.plot.b1$ID <- factor(data.plot.b1$ID,
                          levels=seq.int(1,180,1))

rd2022.b1.ID <- ggplot(data.plot.b1, aes(x=ID, y=mean, group=Lineage, color=Lineage)) + 
  geom_bar(color="black",stat="identity",aes(fill=Lineage)) +
  # geom_linerange(aes(x=ID, ymax=mean+se,ymin=mean)) +
  theme_prism() + theme(axis.text.y=element_text(size=20,face="bold"),
                        axis.title.y=element_blank(),
                        axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        legend.position="none",
                        plot.margin = unit(c(.4,0,1,1), "cm"),
                        plot.title = element_text(hjust = .1,
                                                  size=28)) + 
  scale_y_continuous(labels = scales::percent,limits=c(0,1),
                     expand = c(0, 0)) +
  scale_color_manual(values=c("white","grey70")) +
  scale_fill_manual(values=c("white","grey70")) +
  ggtitle("Block 5")


data.plot.b2 <- data.plot.b2[order(data.plot.b2$mean),]
data.plot.b2$ID <- NA
data.plot.b2[data.plot.b2$Lineage=="A\u2640B\u2642","ID"] <- seq.int(1,180,2)
data.plot.b2[data.plot.b2$Lineage=="B\u2640A\u2642","ID"] <- seq.int(2,180,2)
data.plot.b2 <- data.plot.b2[order(data.plot.b2$ID),]
data.plot.b2$ID <- factor(data.plot.b2$ID,
                          levels=seq.int(1,180,1))

rd2022.b2.ID <- ggplot(data.plot.b2, aes(x=ID, y=mean, group=Lineage, color=Lineage)) + 
  geom_bar(color="black",stat="identity",aes(fill=Lineage)) +
  # geom_linerange(aes(x=ID, ymax=mean+se,ymin=mean)) +
  theme_prism() + theme(axis.text.y=element_text(size=20,face="bold"),
                        axis.title.x=element_text(size=20,face="bold"),
                        axis.title.y=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        legend.position="none",
                        plot.margin = unit(c(.4,0,.2,1), "cm"),
                        plot.title = element_text(hjust = .1,
                                                  size=28)) + 
  xlab("Individual") +
  scale_y_continuous(labels = scales::percent,limits=c(0,1),
                     expand = c(0, 0)) +
  scale_color_manual(values=c("white","grey70")) +
  scale_fill_manual(values=c("white","grey70")) +
  ggtitle("Block 6")

rID <- plot_grid(rd2022.b1.ID,rd2022.b2.ID,ncol=1,nrow=2) +
  draw_label("Frequency of response to QMP", 
             x=0, y=0.5, vjust= 1.2, angle=90,fontface="bold",size=20)

cairo_pdf("retinue_2022_ID.pdf",width=10,height=6,
          fallback_resolution=1200)
grid.draw(rID)
dev.off()



p1 <- ggdraw() + draw_image(magick::image_read_pdf("retinue_2019_2021.pdf",density=1200),scale=0.9)
p2 <- ggdraw() + draw_image(magick::image_read_pdf("retinue_2022.pdf",density=1200),scale=0.9)
p3 <- ggdraw() + draw_image(magick::image_read_pdf("retinue_2022_ID.pdf",density=1200),scale=0.9)
p.bottom <- plot_grid(p2, p3,labels=c("B","C"),ncol=2,nrow=1,label_size=200)
p.out <- plot_grid(p1,p.bottom,nrow=2,ncol=1,labels=c("A"),label_size=200)
ggsave("retinue_all.png",p.out,width=9,height=6,bg="white",dpi=1200)


