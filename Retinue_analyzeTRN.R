WGCNA_mods <- read.csv("WGCNA_modules_locked.csv")
SVM_genes <- read.csv("SVM_GeneIDs.csv")
SVM_genes$GeneID <- paste0("LOC",SVM_genes$GeneID)
TRN <- read.csv("retinue_brain_TRN.csv")
datExpr <- read.csv("datExpr_SVM.csv",row.names=1)
background <- paste0("LOC",row.names(datExpr))
PSGEs <- read.csv("PSGE_union.csv")
modules <- data.frame(geneID=TRN$gene,module=TRN$TF)



A <- length(intersect(SVM_genes$GeneID,PSGEs$gene))
B <- length(setdiff(PSGEs$gene,SVM_genes$GeneID))
C <- length(setdiff(SVM_genes$GeneID,PSGEs$gene))
D <- length(setdiff(background,union(SVM_genes$GeneID,PSGEs$gene)))
dTable <- matrix(c(A,B,C,D),nrow=2)
fisher.test(dTable,alternative="greater")$p.value



modEnrich <- data.frame(matrix(ncol=length(unique(WGCNA_mods$module))+1,
                               nrow=length(unique(modules$module))))
names(modEnrich) <- c("TRNmodule",unique(WGCNA_mods$module))
modEnrich$TRNmodule <- unique(modules$module)
modEnrich <- modEnrich[modEnrich$TRNmodule%in%PSGEs$gene,]
modEnrich$SVM <- NA

enrichMod <- function(modules,module,genelist,background){
  mod.sub <- modules[modules$module==module,]
  A <- length(intersect(mod.sub$geneID,genelist))
  B <- length(setdiff(genelist,mod.sub$geneID))
  C <- length(setdiff(mod.sub$geneID,genelist))
  D <- length(setdiff(background,union(mod.sub$geneID,genelist)))
  dTable <- matrix(c(A,B,C,D),nrow=2)
  fisher.test(dTable,alternative="greater")$p.value
}

for(i in 1:length(modEnrich$TRNmodule)){
  for(k in 1:length(unique(WGCNA_mods$module))){
    WGCNAmod <- unique(WGCNA_mods$module)[k]
    modEnrich[i,names(modEnrich)==WGCNAmod] <- enrichMod(modules,modEnrich$TRNmodule[i],
                                WGCNA_mods[WGCNA_mods$module==WGCNAmod,"geneID"],
                                background)
  }
  modEnrich[i,"SVM"] <- enrichMod(modules,modEnrich$TRNmodule[i],
                                  SVM_genes$GeneID,background)
}

write.csv(modEnrich,"modEnrich.csv")



modules <- unique(WGCNA_mods$module)
modEnrich_WGCNA <- data.frame(module=modules)
modEnrich_WGCNA$genes <- NA
modEnrich_WGCNA$PSGEs <- NA
modEnrich_WGCNA$PSGE.p <- NA
modEnrich_WGCNA$SVM <- NA
modEnrich_WGCNA$SVM.p <- NA

enrichMod <- function(modules,module,genelist){
  mod.sub <- modules[modules$module==module,]
  A <- length(intersect(mod.sub$geneID,genelist))
  B <- length(setdiff(genelist,mod.sub$geneID))
  C <- length(setdiff(mod.sub$geneID,genelist))
  D <- length(setdiff(modules$geneID,union(mod.sub$geneID,genelist)))
  dTable <- matrix(c(A,B,C,D),nrow=2)
  fisher.test(dTable,alternative="greater")$p.value
}

for(i in 1:length(modules)){
  sub <- WGCNA_mods[WGCNA_mods$module==modules[i],]
  modEnrich_WGCNA[i,2] <- length(sub$geneID)
  modEnrich_WGCNA[i,3] <- length(intersect(sub$geneID,
                                           PSGEs$gene))
  modEnrich_WGCNA[i,4] <- enrichMod(WGCNA_mods,modules[i],PSGEs$gene)
  modEnrich_WGCNA[i,5] <- length(intersect(sub$geneID,
                                           SVM_genes$GeneID))
  modEnrich_WGCNA[i,6] <- enrichMod(WGCNA_mods,modules[i],SVM_genes$GeneID)
}

write.csv(modEnrich_WGCNA,"modEnrich_WGCNA.csv",row.names=F)



TFPSGE_targets <- TRN[TRN$TF%in%modEnrich$TRNmodule,"gene"]

WGCNA_TFPSGE <- data.frame(module=modules)
for(i in 1:length(modules)){
  sub <- WGCNA_mods[WGCNA_mods$module==modules[i],]
  WGCNA_TFPSGE[i,2] <- length(sub$geneID)
  WGCNA_TFPSGE[i,3] <- length(intersect(sub$geneID,
                                           TFPSGE_targets))
  WGCNA_TFPSGE[i,4] <- enrichMod(WGCNA_mods,modules[i],TFPSGE_targets)
}
names(WGCNA_TFPSGE) <- c("module","genes","TFPSGE_target","p")

write.csv(WGCNA_TFPSGE,"modEnrich_TFPSGE.csv",row.names=F)




figs5.dat <- read.csv("figs5_panelA.csv")
names(figs5.dat)[2] <- "value"
figs5.dat.new <- data.frame(module=c(modEnrich_WGCNA$module,
                                     modEnrich_WGCNA$module,
                                     modEnrich_WGCNA$module),
                            value=c(modEnrich_WGCNA$PSGE.p,
                                    modEnrich_WGCNA$SVM.p,
                                    WGCNA_TFPSGE$p),
                            label=c(paste(paste(modEnrich_WGCNA$PSGEs,
                                                modEnrich_WGCNA$genes,sep="/"),
                                          paste0("(",format(modEnrich_WGCNA$PSGE.p,scientific = TRUE, 
                                                            digits = 3),")"),sep="\n"),
                                    paste(paste(modEnrich_WGCNA$SVM,
                                                modEnrich_WGCNA$genes,sep="/"),
                                          paste0("(",format(modEnrich_WGCNA$SVM.p,scientific = TRUE, 
                                                            digits = 3),")"),sep="\n"),
                                    paste(paste(WGCNA_TFPSGE$TFPSGE_target,
                                                modEnrich_WGCNA$genes,sep="/"),
                                          paste0("(",format(WGCNA_TFPSGE$p,scientific = TRUE, 
                                                            digits = 3),")"),sep="\n")),
                            panel=c(rep.int("PSGE.p",length(modEnrich_WGCNA$module)),
                                    rep.int("SVM.p1",length(modEnrich_WGCNA$module)),
                                    rep.int("SVM.p2",length(modEnrich_WGCNA$module))))
figs5.dat <- rbind(figs5.dat,figs5.dat.new)
rm(figs5.dat.new)
figs5.dat$x <- "Behavioral State"
figs5.dat$module <- factor(figs5.dat$module,
                           levels=c("saddlebrown","skyblue","blue","darkturquoise",
                                    "skyblue3","darkolivegreen","violet","orange",
                                    "paleturquoise","black","sienna3","plum1",
                                    "darkred","grey60","darkmagenta","white",
                                    "darkgrey","darkorange","steelblue","lightyellow",
                                    "cyan","lightgreen"))
levels(figs5.dat$module) <- paste0("Module ",1:22)
figs5.dat$module <- factor(figs5.dat$module,
                           levels=rev(levels(figs5.dat$module)))
figs5.dat$fill <- "white"
figs5.dat[figs5.dat$value<0.05,"fill"] <- "red"

library(ggplot2)
library(ggprism)

g1 <- ggplot(figs5.dat[figs5.dat$panel=="ME.t.corr",], 
             aes(x, module,fill=fill)) +
  geom_tile(aes(fill=fill),color="black") +
  geom_text(aes(label = label),fontface="bold",
            color="black",lineheight = .75) +
  scale_fill_manual(values=c("red","grey95")) +
  theme_prism() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        plot.margin=unit(c(.8,0,0,0), "cm")) +
  ggtitle("Correlation with\nbehavioral state")
g1

g2 <- ggplot(figs5.dat[figs5.dat$panel=="PSGE.p",], 
             aes(x, module,fill=fill)) +
  geom_tile(aes(fill=fill),color="black") +
  geom_text(aes(label = label),fontface="bold",
            color="black",lineheight = .75) +
  scale_fill_manual(values=c("red","grey95")) +
  theme_prism() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.text.x=element_blank(),
        plot.margin=unit(c(.8,0,0,0), "cm")) +
  ggtitle("Enrichment for\nparent-biased genes")
g2

g3 <- ggplot(figs5.dat[figs5.dat$panel=="SVM.p1",], 
             aes(x, module,fill=fill)) +
  geom_tile(aes(fill=fill),color="black") +
  geom_text(aes(label = label),fontface="bold",
            color="black",lineheight = .75) +
  scale_fill_manual(values=c("red","grey95")) +
  theme_prism() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        plot.margin=unit(c(.8,0,0,0), "cm")) +
  ggtitle("Enrichment for\nSVM genes")
g3

g4 <- ggplot(figs5.dat[figs5.dat$panel=="SVM.p2",], 
             aes(x, module,fill=fill)) +
  geom_tile(aes(fill=fill),color="black") +
  geom_text(aes(label = label),fontface="bold",
            color="black",lineheight = .75) +
  scale_fill_manual(values=c("red","grey95")) +
  theme_prism() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        plot.margin=unit(c(.8,0,0,0), "cm")) +
  ggtitle("Enrichment for\nparent-biased TF-targets")
g4



library(cowplot)
p.out <- plot_grid(g1,g2,g3,g4,nrow=1,ncol=4,labels=c("A","B","C","D"))
ggsave("figS5.png",p.out,width=12,height=8.5,dpi=300)

