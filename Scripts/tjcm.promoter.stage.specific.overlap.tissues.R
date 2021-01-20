#qsub -t 30 "Rscript ~/scripts/tjcm.promoter.stage.specific.overlap.tissues.R"
# prewd="/data/jtan/Mammalian2020Promoter/CM2019RNAseqHuman"
# stage.order = c("b4week","b5week","b6week","b7week","b8week","b9week","b10week","b11week","b12week","b13week","b16week","b18week",
#                 "b19week","b20week","p0month","p6month","p2year","p5year","p15year","p25year","p35year","p45year","p55year")
prewd="/data/jtan/Mammalian2020Promoter/CM2019RNAseqMouse"
stage.order = c("E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5", "E17.5", "E18.5", "P0", "P03", "P14", "P28", "P63")

sample.info = read.delim(file.path(prewd,"ProcessedData","PromoterObjects","sample.sheet.txt"),header=TRUE, row.names = 1, stringsAsFactors = FALSE)
rownames(sample.info) = sample.info$newnamegender

library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(GenomicFeatures)
library(ensembldb)
library(GenomicAlignments)
library(AnnotationDbi)
library(dplyr)
library(proActiv)
library(ggplot2)
library(reshape2)
library(maSigPro)
library(doParallel)
library(foreach)
library(pheatmap)
library(webshot)
library(plotly)
library(qdapTools)
library(UpSetR)

method="Proactiv"
# method="Dexseq"
norm.method = "edger"
# norm.method = "deseq2"
gender = "Male"
# gender = "Female"
# gender = "ale"
category.order = c("Down_Down","Down_Flat","Flat_Down","Down_Up","Up_Down","Flat_Up","Up_Flat","Up_Up")

source("~/scripts/tjcm.promoter.stage.specific.function.R")

if(gender == "Male"){
  tissue.order = c("Brain", "Cerebellum", "Heart", "Kidney" , "Liver", "Testis")
} else if (gender == "Female") {
  tissue.order = c("Brain", "Cerebellum", "Heart", "Kidney" , "Liver", "Ovary")
} else {
  tissue.order = c("Brain", "Cerebellum", "Heart", "Kidney" , "Liver", "Ovary", "Testis")
}
mc.cores= 4
flat.cluster = c(0)
up.cluster = c(2) #2,3,6
down.cluster = c(1) #1,4,5

basic.objects=c("basic.objects",ls(),"name")

#############################################################
##This script is mainly for compare different tissues
##Check the overlap between different tissues
##Part 1: promoter number statistics in Up/Down cluster in each tissue
##Part 2: up cluster overlap between all tissues, so for down and flat cluster
##Part 3: table categories of overlapped genes and tissue specific genes
##Part 3: GO annotation of overlapped genes and tissue specific genes
############Three important table for all tissues overlap #####################
##cluster category (up,down,flat) info: all.tissues.promoter.cluster.category.list
##table category info: all.tissues.promoter.table.category.list ("Down_Down","Down_Flat","Flat_Down","Down_Up","Up_Down","Flat_Up","Up_Flat","Up_Up")
##overlap between tissues: all.tissues.promoter.up/down/flat.cluster.overlap
############Two important table for each tissueafterwards#####################
##in each tissue, find DE promoters common with others: all.tissues.promoter.table.category.common
##in each tissue, find DE promoters specific to this tissue: all.tissues.promoter.table.category.specific

############################################################
##Prepare cluster info and table category info from all tissues 
############################################################
if (TRUE){
    
  dir.create(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration"))
  dir.create(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender))
  dir.create(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration"))
  dir.create(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender))
  
  all.tissues.promoter.table.category.list = list()
  all.tissues.promoter.cluster.category.list = list()
  
  for (name in tissue.order){
    message(name)
    # name="Testis"
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.see.cluster0.RData")))
    load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("significant.isoform.tissue.table.category.RData")))
    
    ##merge cluster info and table category info from different tissues into one list
    all.tissues.promoter.table.category.list[[name]] = significant.isoform.tissue.table.category
    all.tissues.promoter.cluster.category.list[[name]] = rbind(significant.isoform.tissue.see.cluster, significant.isoform.tissue.see.cluster0)
  }
  
  
  save(all.tissues.promoter.table.category.list, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.table.category.list.RData"))
  save(all.tissues.promoter.cluster.category.list, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.cluster.category.list.RData"))
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
  
}

############################################################################
######promoter number statistics in Up/Down cluster in each tissue
############################################################################
if(TRUE){
  # name = "Brain" "Testis"
  
  ##load the cluster.category and table.category info
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender,"all.tissues.promoter.cluster.category.list.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender,"all.tissues.promoter.table.category.list.RData"))
  tissue.number.order = c( "Heart","Kidney","Cerebellum","Liver","Brain", "Testis")
  
  ##number statistics for up down cluster in all tissues
  all.tissue.cluster.category = sapply(all.tissues.promoter.cluster.category.list, function(x){return(table(x$cluster))})
  cluster.name = setNames(c("Flat","Up","Down"),c(flat.cluster, up.cluster, down.cluster))
  rownames(all.tissue.cluster.category) = cluster.name[rownames(all.tissue.cluster.category)]
  all.tissue.cluster.category = all.tissue.cluster.category[-1,]
  all.tissue.cluster.category.plot = melt(all.tissue.cluster.category)
  all.tissue.cluster.category.plot$Var1 = factor(all.tissue.cluster.category.plot$Var1, levels=c("Up","Down"))
  all.tissue.cluster.category.plot$Var2 = factor(as.character(all.tissue.cluster.category.plot$Var2), levels=tissue.number.order)
  colnames(all.tissue.cluster.category.plot) = c("Type","tissue","value")
  
  p1=ggplot(all.tissue.cluster.category.plot,aes(fill=Type, x=tissue,y=value)) +
    geom_bar(stat="identity", alpha=0.8,width=0.7) + #,position="fill"
    geom_text(aes(label=value),position = position_stack(reverse = FALSE, vjust=0.5), size=3.5) +
    scale_fill_manual(values=c("#f76a8c", "#77d8d8")) + #guides(fill=FALSE) +
    xlab("") +
    ylab("Promoter Number") +
    coord_flip() +
    theme_classic() 
  print(p1)
  ggsave(p1,file=file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.cluster.category.pdf")), width=4, height=4)
  
  ##number statistics for up down cluster  in all tissues
  table.category.statistics = sapply(all.tissues.promoter.table.category.list, function(x){return(table(x$category))})
  table.category.statistics.plot = melt(table.category.statistics)
  table.category.statistics.plot$Var1 = factor(table.category.statistics.plot$Var1, levels = rev(category.order))
  table.category.statistics.plot$Var2 = factor(table.category.statistics.plot$Var2, levels = tissue.number.order)
  colnames(table.category.statistics.plot) = c("Type","tissue","value")
  
  p2=ggplot(table.category.statistics.plot,aes(fill=Type, x=tissue,y=value)) +
    geom_bar(stat="identity",position="fill", alpha=0.8,width=0.7) + #,position="fill"
    geom_text(aes(label=value),position = position_fill(reverse = FALSE, vjust=0.5), size=2) +
    scale_fill_manual(values=c( "#f32053", "#f76a8c", "#fcc7d4", "#ff9206", "#ffce8f","#216353","#8cba51","#deff8b") ) + #guides(fill=FALSE) +
    scale_y_continuous(labels = scales::percent) +
    xlab("") +
    ylab("") +
    coord_polar(theta="y") +
    theme_classic() + theme(axis.text.x=element_blank())
  print(p2)
  ggsave(p2,file=file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.table.category.pdf")), width=5, height=5)
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}


##############################################################
######up cluster overlap between all tissues, also for down and flat cluster
##############################################################
if(TRUE){
  # name = "Brain" "Testis"
  
  ##load the cluster.category and table.category info
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender,"all.tissues.promoter.cluster.category.list.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender,"all.tissues.promoter.table.category.list.RData"))
  
  ##selected cluster
  build.cluster.plot = function(selected.cluster=up.cluster){
    cluster.selected.list = lapply(all.tissues.promoter.cluster.category.list, function(x){return(rownames(x[x$cluster==selected.cluster,,drop=FALSE]))} )
    cluster.selected.list.binary = t(mtabulate(cluster.selected.list)!=0)
    cluster.selected.list.plot = as.data.frame(apply(cluster.selected.list.binary, 2, as.numeric))
    rownames(cluster.selected.list.plot) = rownames(cluster.selected.list.binary)
    return(cluster.selected.list.plot)
  }   
  
  ##up cluster
  cluster.name= "up"
  all.tissues.promoter.up.cluster.overlap = build.cluster.plot(selected.cluster=up.cluster)
  p1 = upset(all.tissues.promoter.up.cluster.overlap, sets.bar.color = c("#fa7d09","#0779e4","#56B4E9","#fc5185","#ffde7d","#1fab89"),order.by = "freq", 
             main.bar.color = c("#fa7d09","#ffde7d","#0779e4","#1fab89","#fc5185","#56B4E9","black", rep("darkgrey",56)),
             nsets=ncol(all.tissues.promoter.up.cluster.overlap) ,keep.order = T, empty.intersections = "on")
  pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.overlap.", cluster.name,".pdf")), width=3, height=3)
  print(p1)
  dev.off()
  save(all.tissues.promoter.up.cluster.overlap, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.up.cluster.overlap.RData"))
  
  ##down cluster
  cluster.name= "down"
  all.tissues.promoter.down.cluster.overlap = build.cluster.plot(selected.cluster=down.cluster)
  p2 = upset(all.tissues.promoter.down.cluster.overlap, sets.bar.color = c("#fa7d09","#1fab89","#0779e4","#ffde7d","#56B4E9","#fc5185"),order.by = "freq", 
             main.bar.color = c("#fa7d09","#1fab89","#ffde7d","#0779e4","darkgrey", "#56B4E9","#fc5185","black", rep("darkgrey",55)),
             nsets=ncol(all.tissues.promoter.down.cluster.overlap) ,keep.order = T, empty.intersections = "on")
  pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.overlap.", cluster.name,".pdf")), width=3, height=3)
  print(p2)
  dev.off()
  save(all.tissues.promoter.down.cluster.overlap, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.down.cluster.overlap.RData"))
  
  ##flat cluster
  cluster.name= "flat"
  all.tissues.promoter.flat.cluster.overlap = build.cluster.plot(selected.cluster=flat.cluster)
  p3 = upset(all.tissues.promoter.flat.cluster.overlap, sets.bar.color = c("#fa7d09","#0779e4","#56B4E9","#fc5185","#ffde7d","#1fab89"),order.by = "freq", 
             main.bar.color = c("#fa7d09","#ffde7d","#0779e4","#1fab89","#fc5185","#56B4E9","black", rep("darkgrey",56)),
             nsets=ncol(all.tissues.promoter.flat.cluster.overlap) ,keep.order = T, empty.intersections = "on")
  pdf(file.path(prewd,"Output","StageSpecific",method,norm.method,"Integration",gender,paste0("all.tissues.overlap.", cluster.name,".pdf")), width=3, height=3)
  print(p3)
  dev.off()
  save(all.tissues.promoter.flat.cluster.overlap, file = file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.flat.cluster.overlap.RData"))
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

################################################
##table categories of overlapped genes and tissue specific genes
##Brain and Cerellum has special treatment
################################################
for (name in tissue.order){
  # name = "Testis" #"Brain" 
  
  ##load the table.category and up/down/flat.cluster.overlap info
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender,"all.tissues.promoter.cluster.category.list.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender,"all.tissues.promoter.table.category.list.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.up.cluster.overlap.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.down.cluster.overlap.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,"Integration",gender, "all.tissues.promoter.flat.cluster.overlap.RData"))
  
  all.tissues.promoter.cluster.category.selected = all.tissues.promoter.cluster.category.list[[name]]
  all.tissues.promoter.table.category.selected = all.tissues.promoter.table.category.list[[name]]
  table.category.statistics = data.frame(table(all.tissues.promoter.table.category.selected$category))
  
  ##common with other tissues
  all.tissues.promoter.up.cluster.common = all.tissues.promoter.up.cluster.overlap[all.tissues.promoter.up.cluster.overlap[,name]==1 & rowSums(all.tissues.promoter.up.cluster.overlap) > 1, ]
  if(name == "Brain"){all.tissues.promoter.up.cluster.2common = all.tissues.promoter.up.cluster.common[all.tissues.promoter.up.cluster.common$Cerebellum==1,]; all.tissues.promoter.up.cluster.common = all.tissues.promoter.up.cluster.common[all.tissues.promoter.up.cluster.common$Cerebellum!=1,]}
  if(name == "Cerebellum"){all.tissues.promoter.up.cluster.2common = all.tissues.promoter.up.cluster.common[all.tissues.promoter.up.cluster.common$Brain==1,]; all.tissues.promoter.up.cluster.common = all.tissues.promoter.up.cluster.common[all.tissues.promoter.up.cluster.common$Brain!=1,]}
  all.tissues.promoter.down.cluster.common = all.tissues.promoter.down.cluster.overlap[all.tissues.promoter.down.cluster.overlap[,name]==1 & rowSums(all.tissues.promoter.down.cluster.overlap) > 1, ]
  if(name == "Brain"){all.tissues.promoter.down.cluster.2common = all.tissues.promoter.down.cluster.common[all.tissues.promoter.down.cluster.common$Cerebellum==1,]; all.tissues.promoter.down.cluster.common = all.tissues.promoter.down.cluster.common[all.tissues.promoter.down.cluster.common$Cerebellum!=1,]}
  if(name == "Cerebellum"){all.tissues.promoter.down.cluster.2common = all.tissues.promoter.down.cluster.common[all.tissues.promoter.down.cluster.common$Brain==1,]; all.tissues.promoter.down.cluster.common = all.tissues.promoter.down.cluster.common[all.tissues.promoter.down.cluster.common$Brain!=1,]}
  
  all.tissues.promoter.cluster.common = c(rownames(all.tissues.promoter.up.cluster.common), rownames(all.tissues.promoter.down.cluster.common))
  all.tissues.promoter.cluster.common.geneid = unique(all.tissues.promoter.cluster.category.selected[rownames(all.tissues.promoter.cluster.category.selected) %in% all.tissues.promoter.cluster.common, "geneid"])
  all.tissues.promoter.table.category.common = all.tissues.promoter.table.category.selected[rownames(all.tissues.promoter.table.category.selected) %in% all.tissues.promoter.cluster.common.geneid, ]
  all.tissues.promoter.table.category.common.plot = data.frame(table(all.tissues.promoter.table.category.common$category))
  save(all.tissues.promoter.table.category.common, file=file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.common.RData")))
  
  all.tissues.promoter.table.category.common.plot = merge(data.frame(table(all.tissues.promoter.table.category.common$category)),table.category.statistics, by.x="Var1", by.y="Var1")
  all.tissues.promoter.table.category.common.plot$percentage = round(all.tissues.promoter.table.category.common.plot$Freq.x/all.tissues.promoter.table.category.common.plot$Freq.y * 100, 2)
  all.tissues.promoter.table.category.common.plot$Var1 = factor(all.tissues.promoter.table.category.common.plot$Var1, category.order)
  p2=ggplot(all.tissues.promoter.table.category.common.plot,aes(x=Var1,y=percentage,fill=Var1)) +
    geom_bar(stat="identity", alpha=0.8,width=0.7) + #,position="fill"
    geom_text(aes(label=percentage),position = position_stack(reverse = FALSE, vjust=0.5), size=3.5) +
    scale_fill_manual(values=c(rep("#77d8d8", 3), rep("#ffac41", 2), rep("#f76a8c", 3)) ) + guides(fill=FALSE) +
    xlab("") +
    ylab("Cateogry percentage (%)") +
    coord_flip() +
    ggtitle(paste(name,"common promoters")) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5))
  print(p2)
  ggsave(p2,file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.common.barplot.pdf")), width=3, height=3)
  
  ##common in 2 tissues: 
  if(name %in% c("Brain","Cerebellum")){
    all.tissues.promoter.cluster.2common = c(rownames(all.tissues.promoter.up.cluster.2common), rownames(all.tissues.promoter.down.cluster.2common))
    all.tissues.promoter.cluster.2common.geneid = unique(all.tissues.promoter.cluster.category.selected[rownames(all.tissues.promoter.cluster.category.selected) %in% all.tissues.promoter.cluster.2common, "geneid"])
    all.tissues.promoter.table.category.2common = all.tissues.promoter.table.category.selected[rownames(all.tissues.promoter.table.category.selected) %in% all.tissues.promoter.cluster.2common.geneid, ]
    all.tissues.promoter.table.category.2common.plot = data.frame(table(all.tissues.promoter.table.category.2common$category))
    save(all.tissues.promoter.table.category.2common, file=file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.2common.RData")))
    
    all.tissues.promoter.table.category.2common.plot = merge(data.frame(table(all.tissues.promoter.table.category.2common$category)),table.category.statistics, by.x="Var1", by.y="Var1")
    all.tissues.promoter.table.category.2common.plot$percentage = round(all.tissues.promoter.table.category.2common.plot$Freq.x/all.tissues.promoter.table.category.2common.plot$Freq.y * 100, 2)
    all.tissues.promoter.table.category.2common.plot$Var1 = factor(all.tissues.promoter.table.category.2common.plot$Var1, category.order)
    p2=ggplot(all.tissues.promoter.table.category.2common.plot,aes(x=Var1,y=percentage,fill=Var1)) +
      geom_bar(stat="identity", alpha=0.8,width=0.7) + #,position="fill"
      geom_text(aes(label=percentage),position = position_stack(reverse = FALSE, vjust=0.5), size=3.5) +
      scale_fill_manual(values=c(rep("#77d8d8", 3), rep("#ffac41", 2), rep("#f76a8c", 3)) ) + guides(fill=FALSE) +
      xlab("") +
      ylab("Cateogry percentage (%)") +
      coord_flip() +
      ggtitle(paste("Brain Cerebellum common promoters")) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5))
    print(p2)
    ggsave(p2,file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.2common.barplot.pdf")), width=3, height=3)
  }
  
  ##specific to one tissue
  all.tissues.promoter.up.cluster.specific = all.tissues.promoter.up.cluster.overlap[all.tissues.promoter.up.cluster.overlap[,name]==1 & rowSums(all.tissues.promoter.up.cluster.overlap) == 1, ]
  all.tissues.promoter.down.cluster.specific = all.tissues.promoter.down.cluster.overlap[all.tissues.promoter.down.cluster.overlap[,name]==1 & rowSums(all.tissues.promoter.down.cluster.overlap) == 1, ]
  all.tissues.promoter.cluster.specific = c(rownames(all.tissues.promoter.up.cluster.specific), rownames(all.tissues.promoter.down.cluster.specific))
  all.tissues.promoter.cluster.specific.geneid = unique(all.tissues.promoter.cluster.category.selected[rownames(all.tissues.promoter.cluster.category.selected) %in% all.tissues.promoter.cluster.specific, "geneid"])
  all.tissues.promoter.cluster.specific.geneid = all.tissues.promoter.cluster.specific.geneid[!all.tissues.promoter.cluster.specific.geneid %in% all.tissues.promoter.cluster.common.geneid]
  if(name %in% c("Brain","Cerebellum")){all.tissues.promoter.cluster.specific.geneid = all.tissues.promoter.cluster.specific.geneid[!all.tissues.promoter.cluster.specific.geneid %in% all.tissues.promoter.cluster.2common.geneid]}
  all.tissues.promoter.table.category.specific = all.tissues.promoter.table.category.selected[rownames(all.tissues.promoter.table.category.selected) %in% all.tissues.promoter.cluster.specific.geneid, ]
  save(all.tissues.promoter.table.category.specific, file=file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.specific.RData")))
  
  all.tissues.promoter.table.category.specific.plot = merge(data.frame(table(all.tissues.promoter.table.category.specific$category)),table.category.statistics, by.x="Var1", by.y="Var1")
  all.tissues.promoter.table.category.specific.plot$percentage = round(all.tissues.promoter.table.category.specific.plot$Freq.x/all.tissues.promoter.table.category.specific.plot$Freq.y * 100, 2)
  all.tissues.promoter.table.category.specific.plot$Var1 = factor(all.tissues.promoter.table.category.specific.plot$Var1, category.order)
  p1=ggplot(all.tissues.promoter.table.category.specific.plot,aes(x=Var1,y=percentage,fill=Var1)) +
    geom_bar(stat="identity", alpha=0.8,width=0.7) + #,position="fill"
    geom_text(aes(label=percentage),position = position_stack(reverse = FALSE, vjust=0.5), size=3.5) +
    scale_fill_manual(values=c(rep("#77d8d8", 3), rep("#ffac41", 2), rep("#f76a8c", 3)) ) + guides(fill=FALSE) +
    xlab("") +
    ylab("Cateogry percentage (%)") +
    coord_flip() +
    ggtitle(paste(name,"specific promoters")) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5))
  print(p1)
  ggsave(p1,file=file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.specific.barplot.pdf")), width=3, height=3)
  
  
  # all.tissues.promoter.up.cluster.overlap.brain.cere = all.tissues.promoter.up.cluster.overlap[all.tissues.promoter.up.cluster.overlap$Brain==1 & all.tissues.promoter.up.cluster.overlap$Cerebellum==1 & rowSums(all.tissues.promoter.up.cluster.overlap[,-1:-2]) == 0, ]
  # all.tissues.promoter.down.cluster.overlap.brain.cere = all.tissues.promoter.down.cluster.overlap[all.tissues.promoter.down.cluster.overlap$Brain==1 & all.tissues.promoter.down.cluster.overlap$Cerebellum==1 & rowSums(all.tissues.promoter.down.cluster.overlap[,-1:-2]) == 0, ]
  # all.tissues.promoter.up.cluster.overlap.common = all.tissues.promoter.up.cluster.overlap[rowSums(all.tissues.promoter.up.cluster.overlap) >=2 & !(rownames(all.tissues.promoter.up.cluster.overlap) %in% rownames(all.tissues.promoter.up.cluster.overlap.brain.cere)) , ]
  # all.tissues.promoter.down.cluster.overlap.common = all.tissues.promoter.down.cluster.overlap[rowSums(all.tissues.promoter.down.cluster.overlap) >=2 & 
  #                                                                                             !(rownames(all.tissues.promoter.down.cluster.overlap) %in% rownames(all.tissues.promoter.down.cluster.overlap.brain.cere)) &
  #                                                                                               !(rownames(all.tissues.promoter.down.cluster.overlap) %in% rownames(all.tissues.promoter.down.cluster.overlap.testis.liver)), ]
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}



################################################
##GO annotation of overlapped genes and tissue specific genes
################################################
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(DOSE)

for (name in tissue.order){
  # name = "Testis" #"Brain" 
  Pcutoff=0.05
  Qcutoff=0.05
  
  #load common and specific DE promoter list
  load(file.path(prewd,"ProcessedData","PromoterObjects",method,norm.method,"gene.expression.star.mean.RData"))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.common.RData")))
  load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.specific.RData")))
  if(name %in% c("Brain","Cerebellum")){load(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.2common.RData")))}
  ###covert GeneID to entrez ID
  database.name = ifelse(length(grep("Mouse",prewd))!=0, "org.Mm.eg.db", "org.Hs.eg.db")
  species.name = ifelse(length(grep("Mouse",prewd))!=0, "Mus musculus", "Homo sapiens")
  expressed.genes.entrezid = bitr(row.names(gene.expression.star.mean[rowSums(gene.expression.star.mean)>10,]), fromType="ENSEMBL", toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb=database.name)
  all.genes.entrezid = bitr(row.names(gene.expression.star.mean), fromType="ENSEMBL", toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb=database.name)
  
  ##convert common and specific DE promoter list ###
  DEgenes.name.list = list(common = all.genes.entrezid[all.genes.entrezid$ENSEMBL %in% rownames(all.tissues.promoter.table.category.common),"ENTREZID"], 
                           specific = all.genes.entrezid[all.genes.entrezid$ENSEMBL %in% rownames(all.tissues.promoter.table.category.specific),"ENTREZID"])
  if(name %in% c("Brain","Cerebellum")){
    DEgenes.name.list = list(common = all.genes.entrezid[all.genes.entrezid$ENSEMBL %in% rownames(all.tissues.promoter.table.category.common),"ENTREZID"],
                             specific = all.genes.entrezid[all.genes.entrezid$ENSEMBL %in% rownames(all.tissues.promoter.table.category.specific),"ENTREZID"], 
                             common2 = all.genes.entrezid[all.genes.entrezid$ENSEMBL %in% rownames(all.tissues.promoter.table.category.2common),"ENTREZID"] )
  }
  #####################
  ##GO Analysis
  #####################
  all.tissues.promoter.table.category.go = try(compareCluster(DEgenes.name.list, fun = "enrichGO", OrgDb = database.name, ont="BP", pvalueCutoff = Pcutoff, qvalueCutoff =  Qcutoff, universe = expressed.genes.entrezid$ENTREZID, readable = TRUE),silent = TRUE)
  save(all.tissues.promoter.table.category.go, file=file.path(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.go.RData"))))
  #####dotPlot of enrich annotations
  pdf( file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.go.pdf") ))
  print(dotplot(all.tissues.promoter.table.category.go, showCategory=14))
  dev.off()
  
  #####################
  ##MSigDb C5
  #####################
  #####common DE promoter
  all.tissues.promoter.table.category.common.msigc5 = try(enricher(DEgenes.name.list[[1]], TERM2GENE=msigdbr(species = species.name , category = "C5") %>% dplyr::select(gs_name, entrez_gene), pvalueCutoff = Pcutoff, qvalueCutoff =  Qcutoff, universe = expressed.genes.entrezid$ENTREZID),silent = TRUE)
  save(all.tissues.promoter.table.category.common.msigc5, file=file.path(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.common.msigc5.RData"))))
  #####specific DE promoter
  all.tissues.promoter.table.category.specific.msigc5 = try(enricher(DEgenes.name.list[[2]], TERM2GENE=msigdbr(species = species.name , category = "C5") %>% dplyr::select(gs_name, entrez_gene), pvalueCutoff = Pcutoff, qvalueCutoff =  Qcutoff, universe = expressed.genes.entrezid$ENTREZID),silent = TRUE)
  save(all.tissues.promoter.table.category.specific.msigc5, file=file.path(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.specific.msigc5.RData"))))
  #####DE promoter that are common in "Brain","Cerebellum" 
  if(name %in% c("Brain","Cerebellum")){
    all.tissues.promoter.table.category.2common.msigc5 = try(enricher(DEgenes.name.list[[3]], TERM2GENE=msigdbr(species = species.name , category = "C5") %>% dplyr::select(gs_name, entrez_gene), pvalueCutoff = Pcutoff, qvalueCutoff =  Qcutoff, universe = expressed.genes.entrezid$ENTREZID),silent = TRUE)
    save(all.tissues.promoter.table.category.2common.msigc5, file=file.path(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.2common.msigc5.RData"))))
  }
  #####dotPlot of enrich annotations
  pdf( file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.msigc5.pdf") ))
  print(dotplot(all.tissues.promoter.table.category.common.msigc5, showCategory=14))
  if(name %in% c("Brain","Cerebellum")){print(dotplot(all.tissues.promoter.table.category.2common.msigc5, showCategory=14))}
  print(dotplot(all.tissues.promoter.table.category.specific.msigc5, showCategory=14))
  dev.off()
  
  if(length(grep("Mouse",prewd))==0){
    #####################
    ### Disease Analysis ###
    #####################
    ##DO
    all.tissues.promoter.table.category.do = try(compareCluster(DEgenes.name.list, fun = "enrichDO", ont="DO", pvalueCutoff = Pcutoff, qvalueCutoff =  Qcutoff, universe = expressed.genes.entrezid$ENTREZID, readable = TRUE),silent = TRUE)
    save(all.tissues.promoter.table.category.do, file=file.path(file.path(prewd,"ProcessedData","StageObjects",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.do.RData"))))
    #####dotPlot of enrich annotations
    pdf( file.path(prewd,"Output","StageSpecific",method,norm.method,name,gender,paste0("all.tissues.promoter.table.category.do.pdf") ))
    print(dotplot(all.tissues.promoter.table.category.do, showCategory=14))
    dev.off()
  }
  
  
  rm(list=setdiff(ls(),basic.objects))
  gc()
}

