library(ggplot2)
library(ggprism)
library(cowplot)
library(grid)
library(gridExtra)

data <- read.csv("ovary_activation.csv")
data$Block <- factor(data$Block)
data$Maternal_Lineage <- factor(data$Maternal_Lineage)
data$Maternal_Lineage <- factor(data$Maternal_Lineage, 
                                levels = c("C5","C9","C6","C11","C1","C3"))
data$Behavior <- factor(data$Behavior)
data$ReciprocalCross <- factor(data$ReciprocalCross)


data.AOV <- data
levels(data.AOV$Block) <- c(1,2,3)

model <- lm(Ovarioles~Behavior*Maternal_Lineage/Block,data=data.AOV)
ANOVA <- anova(model)
ANOVA

model <- lm(Ovarioles~Behavior+Maternal_Lineage/Behavior,
            data=data.AOV[data.AOV$Block==1,])
ANOVA <- anova(model)
ANOVA

model <- lm(Ovarioles~Behavior+Maternal_Lineage/Behavior,
            data=data.AOV[data.AOV$Block==2,])
ANOVA <- anova(model)
ANOVA

model <- lm(Ovarioles~Behavior+Maternal_Lineage/Behavior,
            data=data.AOV[data.AOV$Block==3,])
ANOVA <- anova(model)
ANOVA

data$label <- data$Maternal_Lineage

levels(data$label) <- c("A\u2640B\u2642",
                        "B\u2640A\u2642",
                        "A\u2640B\u2642",
                        "B\u2640A\u2642",
                        "A\u2640B\u2642",
                        "B\u2640A\u2642")

g1 <- ggplot(data[data$Block=="C1xC3",], 
       aes(x=Behavior, y=Ovarioles,fill=label)) +
  geom_boxplot(aes(fill=label), width=0.5,
               position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y',stackdir='center',dotsize=.8,
               position=position_dodge(0.8)) +
  theme_prism() + theme(text = element_text(size=12),
                        legend.text=element_text(size=12,face="bold",color="white"),
                        axis.text.x=element_text(color="black",vjust=-5),
                        axis.text.y=element_text(color="black"),
                        legend.position="top",
                        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("") +
  scale_fill_manual(values=c("white", "grey70"),
                    guide = guide_legend(override.aes = list(fill = "white",
                                                             color = "white"))) +
  scale_color_manual(values=c("white", "grey70"),
                     guide = guide_legend(override.aes = list(color = "white")))+
  ylim(2,10) + ggtitle("Block 1") +
  guides(color = guide_legend(override.aes = list(size=.1)))
g1

g2 <- ggplot(data[data$Block=="C5xC9",], 
             aes(x=Behavior, y=Ovarioles,fill=label)) +
  geom_boxplot(aes(fill=label), width=0.5,
               position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y',stackdir='center',dotsize=.8,
               position=position_dodge(0.8)) +
  theme_prism() + theme(legend.title=element_blank(),
                          legend.text=element_text(size=12,face="bold"),
                          strip.text.x=element_blank(),
                          text=element_text(size=12),
                          axis.text.x=element_text(color="black",vjust=-5),
                          legend.position="top",
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          axis.line.y = element_blank(),
                        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("") + ylab("") +
  scale_fill_manual(values=c("white","grey70")) +
  scale_color_manual(values=c("white","grey70")) +
  ylim(2,10) + ggtitle("Block 3") +
  guides(color = guide_legend(override.aes = list(size=.1)))
g2

g3 <- ggplot(data[data$Block=="C6xC11",], 
             aes(x=Behavior, y=Ovarioles,fill=label)) +
  geom_boxplot(aes(fill=label), width=0.5,
               position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y',stackdir='center',dotsize=.8,
               position=position_dodge(0.8)) +
  theme_prism() + theme(legend.title=element_blank(),
                        legend.text=element_text(size=12,face="bold",color="white"),
                        strip.text.x=element_blank(),
                        axis.text.x=element_text(color="black",vjust=-5),
                        text=element_text(size=12),
                        legend.position="top",
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        axis.line.y = element_blank(),
                        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("") + ylab("") +
  scale_fill_manual(values=c("white", "grey70"),
                    guide = guide_legend(override.aes = list(fill = "white",
                                                             color = "white"))) +
  scale_color_manual(values=c("white", "grey70"),
                     guide = guide_legend(override.aes = list(color = "white")))+
  ylim(2,10) + ggtitle("Block 4") +
  guides(color = guide_legend(override.aes = list(size=.1)))
g3


OCplot <- plot_grid(g1, g2, g3,nrow=1,ncol=3)

cairo_pdf("ovariole_counts.pdf",width=9.5,height=3,
          fallback_resolution=1200)
grid.draw(OCplot)
dev.off()




data.AOV <- data
levels(data.AOV$Block) <- c(1,2,3)

model <- lm(as.integer(Activation_Score)~Behavior*Maternal_Lineage/Block,data=data.AOV)
ANOVA <- anova(model)
ANOVA

model <- lm(as.integer(Activation_Score)~Behavior+Maternal_Lineage/Behavior,
            data=data.AOV[data.AOV$Block==1,])
ANOVA <- anova(model)
ANOVA

model <- lm(as.integer(Activation_Score)~Behavior+Maternal_Lineage/Behavior,
            data=data.AOV[data.AOV$Block==2,])
ANOVA <- anova(model)
ANOVA

model <- lm(as.integer(Activation_Score)~Behavior+Maternal_Lineage/Behavior,
            data=data.AOV[data.AOV$Block==3,])
ANOVA <- anova(model)
ANOVA


data.AS <- data[,c(2,3)]
data.AS <- data.AS[!duplicated(data.AS),]
data.AS <- rbind(data.AS,data.AS)
data.AS <- data.AS[order(data.AS$Behavior),]
data.AS$Activation_Score <- c(rep.int(1,6),rep.int(0,6),
                              rep.int(1,6),rep.int(0,6))

mean_sd <- function(data.sd,data){
  mean.list <- list()
  sd.list <- list()
  for(i in 1:length(row.names(data.sd))){
    data.sub <- data[data$Behavior==data.sd[i,"Behavior"]&
                       data$Maternal_Lineage==data.sd[i,"Maternal_Lineage"]&
                       data$Activation_Score==data.sd[i,"Activation_Score"],]
    crosses <- unique(as.character(data.sub$ReciprocalCross))
    crosses.count <- list()
    for(k in 1:length(crosses)){
      data.sub.sub <- data.sub[data.sub$ReciprocalCross==crosses[k],]
      crosses.count[k] <- length(row.names(data.sub.sub))
    }
    mean.list[i] <- mean(unlist(crosses.count))
    sd.list[i] <- sd(unlist(crosses.count))
  }
  data.sd$mean <- unlist(mean.list)
  data.sd$sd <- unlist(sd.list)
  return(data.sd)
}

data.AS <- mean_sd(data.AS,data)
data.AS[is.na(data.AS$sd),"sd"] <- 0
data.AS[data.AS$mean==0&data.AS$sd==0,c("mean","sd")] <- NA
data.AS$Activation_Score <- as.factor(data.AS$Activation_Score)
data.AS$Maternal_Lineage <- factor(data.AS$Maternal_Lineage, 
                                levels = c("C5","C9","C6","C11","C1","C3"))

g1 <- ggplot(data.AS[data.AS$Maternal_Lineage%in%c("C1","C3"),], 
             aes(x=Activation_Score, y=mean,
                 fill=Maternal_Lineage)) + 
  geom_bar(stat="identity",position=position_dodge(),
           aes(fill=Maternal_Lineage),color="black") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme_prism() +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=12,face="bold"),
        strip.text.x=element_blank(),
        text=element_text(size=12),
        legend.position="none",
        plot.margin = margin(b = 0)) + xlab("") +
  scale_fill_manual(values=c("white","grey70")) +
  ylab("Individuals") +
  scale_y_continuous(limits = c(0, 10.5), breaks = c(0,2,4,6,8,10)) +
  facet_grid(~Behavior)
g1

g2 <- ggplot(data.AS[data.AS$Maternal_Lineage%in%c("C5","C9"),], 
             aes(x=Activation_Score, y=mean,
                 fill=Maternal_Lineage)) + 
  geom_bar(stat="identity",position=position_dodge(),
           aes(fill=Maternal_Lineage),color="black") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme_prism() +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=12,face="bold"),
        strip.text.x=element_blank(),
        text=element_text(size=12),
        legend.position="none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        plot.margin = margin(b = 0)) + xlab("") + ylab("") +
  scale_fill_manual(values=c("white","grey70")) +
  ylim(0,10.5) +
  facet_grid(~Behavior)
g2

g3 <- ggplot(data.AS[data.AS$Maternal_Lineage%in%c("C6","C11"),], 
             aes(x=Activation_Score, y=mean,
                 fill=Maternal_Lineage)) + 
  geom_bar(stat="identity",position=position_dodge(),
           aes(fill=Maternal_Lineage),color="black") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme_prism() +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=12,face="bold"),
        strip.text.x=element_blank(),
        text=element_text(size=12),
        legend.position="none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        plot.margin = margin(b = 0)) + xlab("") + ylab("") +
  scale_fill_manual(values=c("white","grey70")) +
  ylim(0,10.5) +
  facet_grid(~Behavior)
g3

OAplot <- plot_grid(g1, g2, g3,nrow=1,ncol=3)
x.grob <- textGrob("Activation score", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))
OAplot.out <- grid.arrange(arrangeGrob(OAplot, bottom = x.grob))
ggsave("ovary_activation.png",OAplot.out,width=9,height=2,dpi=1200)




p1 <- ggdraw() + draw_image(magick::image_read_pdf("ovariole_counts.pdf",density=1200),scale=0.95)
p2 <- ggdraw() + draw_image("ovary_activation.png",scale=0.85)
p.out <- plot_grid(p1, p2,labels=c("A","B"),ncol=1,nrow=2,label_size=200)
ggsave("ovary_figure.png",p.out,width=9,height=5,bg="white",dpi=1200)


