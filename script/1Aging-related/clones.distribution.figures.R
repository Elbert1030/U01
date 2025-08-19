library(ggpubr)
library(RColorBrewer)
library(EnhancedVolcano)
library(DescTools)
library(reshape2)
library(ggalluvial)
library(dplyr)
# install.packages("treemapify")
library(treemapify)
library(Seurat)
library(patchwork)  # recommended


lnWD<-"/project2/sli68423_1316/users-sync/lamis/U01_Projects/AllExpDownStream/"
setwd(lnWD)
seurat.vivo1<-readRDS("hspc.combined.vivo.ann.exp1.RDS")
colnames(seurat.vivo1@meta.data)
seurat.HSC.vivo1<-subset(seurat.vivo1,subset=celltype%in%c("LT-HSC","ST-HSC","HSC"))
#rm(seurat.vivo1)

seurat.vivo2<-readRDS("hspc.combined.vivo.ann.exp2.RDS")
seurat.HSC.vivo2<-subset(seurat.vivo2,subset=celltype%in%c("LT-HSC","ST-HSC","HSC"))
#rm(seurat.vivo2)

seurat.vivo.cross<-readRDS("hspc.combined.vivo.ann.cross.RDS")
colnames(seurat.vivo.cross@meta.data)
head(seurat.vivo.cross@meta.data$Cells)
seurat.HSC.vivo.cross<-subset(seurat.vivo.cross,subset=celltype%in%c("LT-HSC","ST-HSC","HSC"))
#rm(seurat.vivo.cross)

seurat.vivo.supp<-readRDS("hspc.combined.vivo.ann.supp.RDS")
seurat.HSC.vivo.supp<-subset(seurat.vivo.supp,subset=celltype%in%c("LT-HSC","ST-HSC","HSC"))
#rm(seurat.vivo.supp)

seurat.vitro1<-readRDS("hspc.combined.vitro.ann.exp1.RDS")
seurat.HSC.vitro1<-subset(seurat.vitro1,subset=celltype%in%c("LT-HSC","ST-HSC","HSC"))
#rm(seurat.vitro1)

seurat.vitro2<-readRDS("hspc.combined.vitro.ann.exp2.RDS")
seurat.HSC.vitro2<-subset(seurat.vitro2,subset=celltype%in%c("LT-HSC","ST-HSC","HSC"))
#rm(seurat.vitro2)

seurat.vitro.cross<-readRDS("hspc.combined.vitro.ann.cross.RDS")
colnames(seurat.vitro.cross@meta.data)
seurat.HSC.vitro.cross<-subset(seurat.vitro.cross,subset=celltype%in%c("LT-HSC","ST-HSC","HSC"))
#rm(seurat.vitro.cross)

seurat.vitro.supp<-readRDS("hspc.combined.vitro.ann.supp.RDS")
seurat.HSC.vitro.supp<-subset(seurat.vitro.supp,subset=celltype%in%c("LT-HSC","ST-HSC","HSC"))
#rm(seurat.vitro.supp)

Sample.vivo1=c("1_YOUNG","2_Old")
Sample.vivo2=c("Y1B","Y2B","Y3B","O1B","O2B","O3B")
Sample.vivo.cross=c("OO","OY","YY1","YY2")
Sample.vivo.supp=c("SC2400270_U01s-O-Veh-1","SC2400279_U01s-Y-Veh-2","SC2400282_U01s-Y-Veh-3",
                   "SC2400272_U01s-O-DQ-1","SC2400280_U01s-O-Veh-2","SC2400283_U01s-O-Veh-3"
                   ,"SC2400268_U01s-Y-Veh-1","SC2400281_U01s-O-DQ-2","SC2400284_U01s-O-DQ-3")

Sample.vitro1<-c("1_Y2","2_Y3","3_Y4","4_O2","5_O3","6_O4")
Sample.vitro2<-c("1_Ya","2_Ya","3_Ya","4_Oa","5_Oa","6_Oa")
Sample.vitro.cross<-c("O_vitro","Y_vitro")
Sample.vitro.supp<-"SC2400110_U01s-pre-HSCT"

Sample.Trace1<-c("Y_Traced","O_Traced")
Sample.Trace2<-c("1_Y1C","2_Y2C","3_Y3C","4_O1C","5_O2C","6_O3C")
Sample.Trace.cross<-c("OO_C","OY_C","YO_C","YY_C") #### switching YO and YY --> fix
Sample.Trace.supp<-c("O-Veh-T-1","Y-Veh-T-2","Y-Veh-T-3","O-DQ-T-1","O-Veh-T-2","O-Veh-T-3","Y-Veh-T-1","O-DQ-T-2","O-DQ-T-3")


Rep.vitro.cross<-c("O_Vitro_Rep1","O_Vitro_Rep2","O_Vitro_Rep3" , "Y_Vitro_Rep1","Y_Vitro_Rep2","Y_Vitro_Rep3")
Rep.vitro.cross<-c("OO_Rep1","OO_Rep2" , "OO_Rep3" , "OY_Rep1" , "OY_Rep2",  "OY_Rep3" , "YO_Rep1",
                      "YO_Rep2",  "YO_Rep3" , "YY_Rep1","YY_Rep2" , "YY_Rep3")
##########
######-----------------------------
Clones<-function(seurat,unique,start){
  clones<-type<-NULL
# x<-seurat$CloneID
  x <- if (unique == 1) seurat$uniqueClonesTraced else seurat$CloneID
  type[which(seurat$celltype%in%c("HSC","LT-HSC","ST-HSC"))]<-"HSC"
  type[which(!seurat$celltype %in% c("HSC", "LT-HSC", "ST-HSC"))]<-"P"
  if(unique==1){ 
    x<-seurat$uniqueClonesTraced}
  
  type<-type[-which(x=="0")]
  x<-x[-which(x=="0")]
  names(x)<-type
  df <- data.frame(
    CloneID = x,
    CellType = type,
    stringsAsFactors = FALSE)
  clones <- df %>%
    group_by( CloneID, CellType) %>%
    summarise(clone.freq = n(), .groups = "drop")
  clones$Sample<-sub("^[^_]+_", "", clones$CloneID)
  if(start==1){clones$Sample<-sub("_[0-9]+$", "", clones$CloneID)}
  return(clones)
}


Clones2<-function(seurat,unique,start){
  clones<-type<-NULL
#  x<-seurat$CloneID
  x <- if(unique == 1) seurat$uniqueClonesTraced else seurat$CloneID
    if(unique==1){ 
    x<-seurat$uniqueClonesTraced}
  x<-x[-which(x=="0")]
  df <- data.frame(
    CloneID = x,
    stringsAsFactors = FALSE)
  clones <- df %>%
    group_by( CloneID) %>%
    summarise(clone.freq = n(), .groups = "drop")
  clones$Sample<-sub("^[^_]+_", "", clones$CloneID)
  if(start==1){clones$Sample<-sub("_[0-9]+$", "", clones$CloneID)}
  return(clones)
}

PieChart2<-function(clones,samples){
  # Create Data
  clones=clones[which(clones$Sample%in%samples),]
  data <- data.frame(
    group=clones$CloneID,
    value=as.numeric(clones$clone.freq) )
  
  colourCount = length(unique(data$group))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  clones.num<-length(clones$CloneID)
  cell.num<-sum(as.numeric(clones$clone.freq) )
  title_text <- paste0(samples,"\nClones: ", clones.num, "\nCells: ", cell.num)
  # Basic piechart
  P<-ggplot(data, aes(x="", y=value, fill=getPalette(colourCount))) +
    # geom_bar(stat="identity", color="white", width =9, show.legend = FALSE) + scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colourCount))+
    geom_bar(stat = "identity", width = 1, show.legend = FALSE)+ scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colourCount))+
    ggtitle(title_text) +
    coord_polar("y", start=0)+theme_void()+theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  return(P)
}

avg_proliferation_expression_by_clone <- function(seu, clone_id, Sample, clone_col , assay = "RNA") {
  genes <- c("Mki67", "Top2a")
  seu<-subset(seu,orig.ident%in%Sample)
  seu.data<-GetAssayData(seu,slot="data",assay="RNA")
  genes<-intersect(genes,rownames(seu.data))
  seu.data.genes<-seu.data[genes,]
  seu.data.genes<-as.matrix(seu.data.genes)
  
  Idents(seu) <- clone_col
  tmp_cells<-subset(seu, idents=clone_id)
  tmp_cells<-colnames(tmp_cells)
  tmp_data<-seu.data.genes[,tmp_cells]
  if(length(tmp_cells)>1){tmp_mean<-apply(tmp_data, 1, mean)}else {tmp_mean<-tmp_data}
  return(mean(tmp_mean))
}

#####----------
BarplotClones <- function(clones,clones2 ,seu,samples,col,Title,Sample,clone_col) {
  clones <- clones[clones$Sample %in% samples, ]
  clones2<-clones2[clones2$Sample %in% samples, ]
  data <- data.frame(
    group = clones$CloneID,
    value = as.numeric(clones$clone.freq),
    celltype = clones$CellType
  )
  # Optional: filter out very small clones
  clones2<- clones2[clones2$clone.freq > 0.02 * sum(clones2$clone.freq ), ]
  clones2 <- clones2[order(clones2$clone.freq, decreasing = TRUE), ]
  
  data <- data [which(data$group%in% clones2$CloneID),]
  # Set factor levels for ordering on x-axis
  data$group <- factor(data$group, levels = clones2$CloneID)
  
  for (j in 1:length(data$group)){
    data$proliferation.exp[j]<-avg_proliferation_expression_by_clone(seu,data$group[j], Sample, clone_col, assay = "RNA")}

  ggplot(data, aes(x = group,y = value,
    fill = ifelse(celltype == "P", proliferation.exp, NA), color = proliferation.exp )) +
    geom_bar(
      stat = "identity",
      position = "stack",
      width = 0.7
    ) +
    scale_fill_gradient(
      low = "grey90", high = col, name = "Proliferation", na.value = "white"
    ) +
    scale_color_gradient(
      low = "grey90", high = col, name = "Proliferation"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.5, "lines")
    ) +
    ggtitle(Title) +
    xlab("Top clones") +
    ylab("Frequency")
  
}

#######
DotplotClones <- function(clones,clones2 ,seu,samples,col,Title,Sample,clone_col ) {
clones <- clones[clones$Sample %in% samples, ]
clones2 <- clones2[clones2$Sample %in% samples, ]
data <- data.frame(
  group = clones$CloneID,
  value = as.numeric(clones$clone.freq),
  celltype = clones$CellType
)

# Filter out very small clones
clones2 <- clones2[clones2$clone.freq > 0.02 * sum(clones2$clone.freq), ]
clones2 <- clones2[order(clones2$clone.freq, decreasing = TRUE), ]

# Subset data to top clones
data <- data[data$group %in% clones2$CloneID,]
data$group <- factor(data$group, levels = clones2$CloneID)
data$proliferation.exp<-NULL
for (j in 1:length(data$group)){
  data$proliferation.exp[j]<-avg_proliferation_expression_by_clone(seu,data$group[j], Sample, clone_col, assay = "RNA")}

ggplot(data, aes(x = group, y = value.vitro, color = proliferation.exp, shape = celltype)) +
    geom_point(size = 4, stroke = 1.2) +
    scale_color_gradient(low = "grey90", high = col, name = "Proliferation") +
    scale_shape_manual(values = c("P" = 16, "HSC" = 1)) +  # solid vs hollow
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
      # ,
      # legend.position = "none"
    )+
    ggtitle(Title) +
    xlab("Top clones") +
    ylab("Frequency")}


DotplotClones.vitro <- function(clones,clones2 ,seu,seu.vitro,samples,col,Title,Sample.vivo,Sample.vitro,clone_col ) {
  clones <- clones[clones$Sample %in% samples, ]
  clones2 <- clones2[clones2$Sample %in% samples, ]
  
  data <- data.frame(
    group = clones$CloneID,
    value = as.numeric(clones$clone.freq),
    celltype = clones$CellType
  )
  
  # Filter out very small clones
  clones2 <- clones2[clones2$clone.freq > 0.02 * sum(clones2$clone.freq), ]
  clones2 <- clones2[order(clones2$clone.freq, decreasing = TRUE), ]
  
  # Subset data to top clones
  data <- data[data$group %in% clones2$CloneID,]
  data$group <- factor(data$group, levels = clones2$CloneID)
  data$proliferation.exp<-NULL
  for (j in 1:length(data$group)){
    data$proliferation.exp[j]<-avg_proliferation_expression_by_clone(seu,data$group[j], Sample.vivo, clone_col, assay = "RNA")}

  data.vitro <- data.frame(
    group = clones.vitro$CloneID,
    value.vitro = as.numeric(clones.vitro$clone.freq),
    celltype = clones.vitro$CellType)
  for (j in 1:length(data.vitro$group)){
    data.vitro$proliferation.exp.vitro[j]<-avg_proliferation_expression_by_clone(seurat.vitro1,data.vitro$group[j], Sample.vitro, clone_col, assay = "RNA")}
  
  clones.vitro<-clones.vitro[which(clones.vitro$CloneID%in%data$group),]
  data<-merge(data,data.vitro,by=c("group","celltype"),all.x=T)
  data$value.vitro[is.na(data$value.vitro)] <- 0
  # Dot plot
  ggplot(data, aes(x = group, y = value.vitro, color = proliferation.exp, shape = celltype)) +
    geom_point(size = 4, stroke = 1.2) +
    scale_color_gradient(low = "grey90", high = col, name = "Proliferation") +
    scale_shape_manual(values = c("P" = 16, "HSC" = 1)) +  # solid vs hollow
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
      # ,
      # legend.position = "none"
    )+
    ggtitle(Title) +
    xlab("Top clones") +
    ylab("Frequency")
}
colnames(seurat.vitro1@meta.data)
colnames(seurat.vitro2@meta.data)
colnames(seurat.vitro.cross@meta.data)
colnames(seurat.vitro.supp@meta.data)
colnames(seurat.vivo1@meta.data)
colnames(seurat.vivo2@meta.data)
colnames(seurat.vivo.cross@meta.data)
colnames(seurat.vivo.supp@meta.data)
colnames(seurat.HSC.vitro1@meta.data)
colnames(seurat.HSC.vitro2@meta.data)
colnames(seurat.HSC.vitro.cross@meta.data)
colnames(seurat.HSC.vitro.supp@meta.data)
colnames(seurat.HSC.vivo1@meta.data)
colnames(seurat.HSC.vivo2@meta.data)
colnames(seurat.HSC.vivo.cross@meta.data)
colnames(seurat.HSC.vivo.supp@meta.data)
head(seurat.vivo1$CloneID)
head(seurat.vivo.cross@meta.data$uniqueClonesTraced)
head(seurat.vivo.supp@meta.data$uniqueClonesTraced)

clones.vivo1<-Clones(seurat.vivo1,1,1)
clones.vivo2<-Clones(seurat.vivo2,0,0)
clones.vivo.cross<-Clones(seurat.vivo.cross,1,1)
clones.vivo.supp<-Clones(seurat.vivo.supp,0,0)

clones.vitro1<-Clones(seurat.vitro1,1,1) 
clones.vitro2<-Clones(seurat.vitro2,0,0)
clones.vitro.cross<-Clones(seurat.vitro.cross,1,1)
clones.vitro.supp<-Clones(seurat.vitro.supp,0,0)

clones.vivo1.total<-Clones2(seurat.vivo1,1,1)
clones.vivo2.total<-Clones2(seurat.vivo2,0,0)
clones.vivo.cross.total<-Clones2(seurat.vivo.cross,1,1)
clones.vivo.supp.total<-Clones2(seurat.vivo.supp,0,0)

clones.vitro1.total<-Clones2(seurat.vitro1,1,1) 
clones.vitro2.total<-Clones2(seurat.vitro2,0,0)
clones.vitro.cross.total<-Clones2(seurat.vitro.cross,1,1)
clones.vitro.supp.total<-Clones2(seurat.vitro.supp,0,0)


clones.HSC.vivo1<-Clones(seurat.HSC.vivo1,1,1)
clones.HSC.vivo2<-Clones(seurat.HSC.vivo2,0,0)
clones.HSC.vivo.cross<-Clones(seurat.HSC.vivo.cross,1,1)
clones.HSC.vivo.supp<-Clones(seurat.HSC.vivo.supp,0,0)

clones.HSC.vitro1<-Clones(seurat.HSC.vitro1,1,1) 
clones.HSC.vitro2<-Clones(seurat.HSC.vitro2,0,0)
clones.HSC.vitro.cross<-Clones(seurat.HSC.vitro.cross,1,1)
clones.HSC.vitro.supp<-Clones(seurat.HSC.vitro.supp,0,0)

clones.HSC.vivo1.total<-Clones2(seurat.HSC.vivo1,1,1)
clones.HSC.vivo2.total<-Clones2(seurat.HSC.vivo2,0,0)
clones.HSC.vivo.cross.total<-Clones2(seurat.HSC.vivo.cross,1,1)
clones.HSC.vivo.supp.total<-Clones2(seurat.HSC.vivo.supp,0,0)

clones.HSC.vitro1.total<-Clones2(seurat.HSC.vitro1,1,1) 
clones.HSC.vitro2.total<-Clones2(seurat.HSC.vitro2,0,0)
clones.HSC.vitro.cross.total<-Clones2(seurat.HSC.vitro.cross,1,1)
clones.HSC.vitro.supp.total<-Clones2(seurat.HSC.vitro.supp,0,0)



col<-c("grey30","chocolate")
for (i in 1:length(Sample.Trace1)){
P<-BarplotClones(clones.vivo1,clones.vivo1.total,seurat.vivo1,Sample.Trace1[i],col[i],Sample.vivo1[i],Sample.vivo1[i],"uniqueClonesTraced")
assign(paste0("P.vivo1","_",i),P)
#   P<-DotplotClones(clones.vivo1,clones.vivo1.total,seurat.vivo1,Sample.Trace1[i],col[i],Sample.vivo1[i],Sample.vivo1[i],"uniqueClonesTraced")
# assign(paste0("P.vivo1.dot","_",i),P)
}

# i=1
# P<-DotplotClones.vitro(clones.vivo1,clones.vivo1.total,seurat.vivo1,seurat.vitro1,Sample.Trace1[i],col[i],Sample.Trace1[i],Sample.vivo1[i],Sample.vitro1[1:3],"uniqueClonesTraced")
# assign(paste0("P.vitro1.dot","_",i),P)
# 
# i=2
# P<-DotplotClones.vitro(clones.vivo1,clones.vivo1.total,seurat.vivo1,seurat.vitro1,Sample.Trace1[i],col[i],Sample.Trace1[i],Sample.vivo1[i],Sample.vitro1[4:6],"uniqueClonesTraced")
# assign(paste0("P.vitro1.dot","_",i),P)
# 

col<-c(rep("grey30",3),rep("chocolate",3))
for (i in 1:length(Sample.Trace2)){
  P<-BarplotClones(clones.vivo2,clones.vivo2.total,seurat.vivo2,Sample.Trace2[i],col[i],Sample.vivo2[i],Sample.vivo2[i],"CloneID")
  assign(paste0("P.vivo2","_",i),P)
  # P<-DotplotClones(clones.vivo2,clones.vivo2.total,seurat.vivo2,Sample.Trace2[i],col[i],Sample.vivo2[i],Sample.vivo2[i],"CloneID")
  # assign(paste0("P.vivo2.dot","_",i),P)
  }

col<-c("chocolate","#65A593","grey30","dodgerblue4")
for (i in 1:length(Sample.Trace.cross)){
  P<-BarplotClones(clones.vivo.cross,clones.vivo.cross.total,seurat.vivo.cross,Sample.Trace.cross[i],col[i],Sample.vivo.cross[i],Sample.vivo.cross[i],"uniqueClonesTraced")
  assign(paste0("P.cross","_",i),P)
  # P<-DotplotClones(clones.vivo.cross,clones.vivo.cross.total,seurat.vivo.cross,Sample.Trace.cross[i],col[i],Sample.vivo.cross[i],Sample.vivo.cross[i],"uniqueClonesTraced")
  # assign(paste0("P.cross.dot","_",i),P)
  }



col<-c("chocolate","#65A593","#65A593","red4","chocolate","chocolate","#65A593","red4","red4")
for (i in 1:length(Sample.Trace.supp)){
  P<-BarplotClones(clones.vivo.supp,clones.vivo.supp.total,seurat.vivo.supp,Sample.Trace.supp[i],col[i],Sample.Trace.supp[i],Sample.vivo.supp[i],"CloneID")
  assign(paste0("P.supp","_",i),P)
  # P<-DotplotClones(clones.vivo.supp,clones.vivo.supp.total,seurat.vivo.supp,Sample.Trace.supp[i],col[i],Sample.vivo.supp[i],Sample.vivo.supp[i],"CloneID")
  # assign(paste0("P.supp.dot","_",i),P)
  }



for (i in 1:length(Sample.Trace1)){
  P<-PieChart2(clones.vivo1.total,Sample.Trace1[i])
  assign(paste0("Pie.vivo1","_",i),P)}

for (i in 1:length(Sample.Trace2)){
  P<-PieChart2(clones.vivo2.total,Sample.Trace2[i])
  assign(paste0("Pie.vivo2","_",i),P)}

for (i in 1:length(Sample.Trace.cross)){
  P<-PieChart2(clones.vivo.cross.total,Sample.Trace.cross[i])
  assign(paste0("Pie.vivo.cross","_",i),P)}

for (i in 1:length(Sample.Trace.supp)){
  P<-PieChart2(clones.vivo.supp.total,Sample.Trace.supp[i])
  assign(paste0("Pie.vivo.supp","_",i),P)}


for (i in 1:length(Sample.Trace1)){
  P<-PieChart2(clones.vitro1.total,Sample.Trace1[i])
  assign(paste0("Pie.vitro1","_",i),P)}

for (i in 1:length(Sample.Trace2)){
  P<-PieChart2(clones.vitro2.total,Sample.Trace2[i])
  assign(paste0("Pie.vitro2","_",i),P)}

for (i in 1:length(Sample.Trace.cross)){
  P<-PieChart2(clones.vitro.cross.total,Sample.Trace.cross[i])
  assign(paste0("Pie.vitro.cross","_",i),P)}

for (i in 1:length(Sample.Trace.supp)){
  P<-PieChart2(clones.vitro.supp.total,Sample.Trace.supp[i])
  assign(paste0("Pie.vitro.supp","_",i),P)}

###### pie for HSC only:
for (i in 1:length(Sample.Trace1)){
  P<-PieChart2(clones.HSC.vivo1,Sample.Trace1[i])
  assign(paste0("Pie.HSC.vivo1","_",i),P)}

for (i in 1:length(Sample.Trace2)){
  P<-PieChart2(clones.HSC.vivo2.total,Sample.Trace2[i])
  assign(paste0("Pie.HSC.vivo2","_",i),P)}

for (i in 1:length(Sample.Trace.cross)){
  P<-PieChart2(clones.HSC.vivo.cross.total,Sample.Trace.cross[i])
  assign(paste0("Pie.HSC.vivo.cross","_",i),P)}

for (i in 1:length(Sample.Trace.supp)){
  P<-PieChart2(clones.HSC.vivo.supp.total,Sample.Trace.supp[i])
  assign(paste0("Pie.HSC.vivo.supp","_",i),P)}


for (i in 1:length(Sample.Trace1)){
  P<-PieChart2(clones.HSC.vitro1.total,Sample.Trace1[i])
  assign(paste0("Pie.HSC.vitro1","_",i),P)}

for (i in 1:length(Sample.Trace2)){
  P<-PieChart2(clones.HSC.vitro2.total,Sample.Trace2[i])
  assign(paste0("Pie.HSC.vitro2","_",i),P)}

for (i in 1:length(Sample.Trace.cross)){
  P<-PieChart2(clones.HSC.vitro.cross.total,Sample.Trace.cross[i])
  assign(paste0("Pie.HSC.vitro.cross","_",i),P)}

for (i in 1:length(Sample.Trace.supp)){
  P<-PieChart2(clones.HSC.vitro.supp.total,Sample.Trace.supp[i])
  assign(paste0("Pie.HSC.vitro.supp","_",i),P)}

###############
setwd("/project2/sli68423_1316/users/Kailiang/Test_Rcode/U1 test code/1Age/output")
dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
out_dir <- "/project2/sli68423_1316/users/Kailiang/Test_Rcode/U1 test code/1Age/output"
  # 你希望生成 Figures 的路径

graphics.off()
pdf("./Figures/TopClonesVivo1.pdf", width = 5,height = 1.5)
ggarrange(P.vivo1.dot_1,P.vivo1.dot_2, font.label = list(size = 25),hjust=-0.2, ncol = 2, nrow = 1)
dev.off()

graphics.off()
pdf("./Figures/TopClonesVivo1Bar.pdf", width = 6,height = 1.5)
ggarrange(P.vivo1_1,P.vivo1_2, font.label = list(size = 25),hjust=-0.2, ncol = 2, nrow = 1)
dev.off()

graphics.off()
pdf("./Figures/TopClonesVivo2Bar.pdf", width = 7,height = 3)
ggarrange(P.vivo2_1,P.vivo2_2,P.vivo2_3,P.vivo2_4,P.vivo2_5,P.vivo2_6, font.label = list(size = 25),hjust=-0.2, ncol = 3, nrow = 2)
dev.off()

graphics.off()
pdf("./Figures/TopClonesVivoSuppBar.pdf", width = 7,height = 4)
ggarrange(P.supp_1,P.supp_5,P.supp_6, P.supp_2,P.supp_3,P.supp_7,P.supp_4,P.supp_8,P.supp_9,font.label = list(size = 25),hjust=-0.2, ncol = 3, nrow = 3)
dev.off()

#ggarrange(P.vivo2.dot_1,P.vivo2.dot_2,P.vivo2.dot_3,P.vivo2.dot_4,P.vivo2.dot_5,P.vivo2.dot_6, font.label = list(size = 25),hjust=-0.2, ncol = 6, nrow = 1)

graphics.off()
pdf("./Figures/TopClonesVivoCrossBar.pdf", width = 8,height = 1.5)
ggarrange(P.cross_3,P.cross_4,P.cross_1,P.cross_2, font.label = list(size = 25),hjust=-0.2, ncol = 4, nrow = 1)
dev.off()


graphics.off()
pdf("./Figures/PieExp1.pdf", width = 8,height = 8)
ggarrange(Pie.vitro1_1,Pie.vitro1_2,Pie.vivo1_1,Pie.vivo1_2, font.label = list(size = 25),hjust=-0.2, ncol = 2, nrow = 2)
dev.off()

graphics.off()
pdf("./Figures/PieExp2.pdf", width = 24,height = 8)
ggarrange(Pie.vitro2_1,Pie.vitro2_2,Pie.vitro2_3,Pie.vitro2_4,Pie.vitro2_5,Pie.vitro2_6,
          Pie.vivo2_1,Pie.vivo2_2, Pie.vivo2_3,Pie.vivo2_4, Pie.vivo2_5,Pie.vivo2_6, 
          font.label = list(size = 25),hjust=-0.2, ncol = 6, nrow = 2)
dev.off()

graphics.off()
pdf("./Figures/PieExp.cross.pdf", width = 16,height = 8)
ggarrange(Pie.vitro.cross_1,Pie.vitro.cross_2,Pie.vitro.cross_3,Pie.vitro.cross_4,
          Pie.vivo.cross_1,Pie.vivo.cross_2, Pie.vivo.cross_3,Pie.vivo.cross_4,  
          font.label = list(size = 25),hjust=-0.2, ncol = 4, nrow = 2)
dev.off()



graphics.off()
pdf("./Figures/PieExp.suppVitro.pdf", width = 12,height = 12)
ggarrange(Pie.vitro.supp_1,Pie.vitro.supp_2,Pie.vitro.supp_3,Pie.vitro.supp_4,
          Pie.vitro.supp_5,Pie.vitro.supp_6,Pie.vitro.supp_7,Pie.vitro.supp_8,Pie.vitro.supp_9,
          font.label = list(size = 25),hjust=-0.2, ncol = 3, nrow = 3)
dev.off()

graphics.off()
pdf("./Figures/PieExp.suppVivoo.pdf", width = 12,height = 12)
ggarrange(Pie.vivo.supp_1,Pie.vivo.supp_2,Pie.vivo.supp_3,Pie.vivo.supp_4,Pie.vivo.supp_5,
          Pie.vivo.supp_6,Pie.vivo.supp_7,Pie.vivo.supp_8,Pie.vivo.supp_9,
          font.label = list(size = 25),hjust=-0.2, ncol = 3, nrow = 3)
dev.off()



########################
graphics.off()
pdf("./Figures/Pie.HSCExp1.pdf", width = 8,height = 8)
ggarrange(Pie.HSC.vitro1_1,Pie.HSC.vitro1_2,Pie.HSC.vivo1_1,Pie.HSC.vivo1_2, font.label = list(size = 25),hjust=-0.2, ncol = 2, nrow = 2)
dev.off()

graphics.off()
pdf("./Figures/Pie.HSCExp2.pdf", width = 24,height = 8)
ggarrange(Pie.HSC.vitro2_1,Pie.HSC.vitro2_2,Pie.HSC.vitro2_3,Pie.HSC.vitro2_4,Pie.HSC.vitro2_5,Pie.HSC.vitro2_6,
          Pie.HSC.vivo2_1,Pie.HSC.vivo2_2, Pie.HSC.vivo2_3,Pie.HSC.vivo2_4, Pie.HSC.vivo2_5,Pie.HSC.vivo2_6, 
          font.label = list(size = 25),hjust=-0.2, ncol = 6, nrow = 2)
dev.off()

graphics.off()
pdf("./Figures/Pie.HSCExp.cross.pdf", width = 16,height = 8)
ggarrange(Pie.HSC.vitro.cross_1,Pie.HSC.vitro.cross_2,Pie.HSC.vitro.cross_3,Pie.HSC.vitro.cross_4,
          Pie.HSC.vivo.cross_1,Pie.HSC.vivo.cross_2, Pie.HSC.vivo.cross_3,Pie.HSC.vivo.cross_4,  
          font.label = list(size = 25),hjust=-0.2, ncol = 4, nrow = 2)
dev.off()



graphics.off()
pdf("./Figures/Pie.HSCExp.suppVitro.pdf", width = 12,height = 12)
ggarrange(Pie.HSC.vitro.supp_1,Pie.HSC.vitro.supp_2,Pie.HSC.vitro.supp_3,Pie.HSC.vitro.supp_4,
          Pie.HSC.vitro.supp_5,Pie.HSC.vitro.supp_6,Pie.HSC.vitro.supp_7,Pie.HSC.vitro.supp_8,Pie.HSC.vitro.supp_9,
          font.label = list(size = 25),hjust=-0.2, ncol = 3, nrow = 3)
dev.off()

graphics.off()
pdf("./Figures/Pie.HSCExp.suppVivoo.pdf", width = 12,height = 12)
ggarrange(Pie.HSC.vivo.supp_1,Pie.HSC.vivo.supp_2,Pie.HSC.vivo.supp_3,Pie.HSC.vivo.supp_4,Pie.HSC.vivo.supp_5,
          Pie.HSC.vivo.supp_6,Pie.HSC.vivo.supp_7,Pie.HSC.vivo.supp_8,Pie.HSC.vivo.supp_9,
          font.label = list(size = 25),hjust=-0.2, ncol = 3, nrow = 3)
dev.off()




x<-seurat.vivo2$CloneID
x<-x[-which(x=="0")]
hist(data$freq)





#clones.vivo1.total<-
clones.vivo2.total<-Clones2(seurat.vivo2,0,0)
clones.vivo.cross.total<-clones.vivo.cross.total[order(clones.vivo.cross.total$clone.freq,decreasing = T),]
clones.vivo.supp.total<-clones.vivo.supp.total[order(clones.vivo.supp.total$clone.freq,decreasing = T),]


seu<-seurat.vivo.supp
seu<-seurat.vivo2

library(Seurat)
library(ggplot2)
library(ggpubr)

# Specify your clone of interest
clone_of_interest <-"2_O-Veh-T-1"
clone_of_interest <- "1_3_Y3C"

# Add logical column for plotting
seu$IsClone <- ifelse(seu$CloneID== clone_of_interest, clone_of_interest, "Other")

# Basic highlighting using DimPlot
p <- DimPlot(seu, group.by = "IsClone", cols = c("#65A593","grey90")) +
  ggtitle(paste("UMAP: Clone", clone_of_interest)) +
  theme_classic()
p


