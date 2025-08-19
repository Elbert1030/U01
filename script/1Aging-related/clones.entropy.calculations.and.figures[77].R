library(ggpubr)
library(RColorBrewer)
library(EnhancedVolcano)
library(DescTools)
library(reshape2)
library(ggalluvial)
library(dplyr)

#clones entropy
# Pie_chart of clone distribution
# histogram of clone distribution
# connected barlpoy


lnWD<-"/project2/sli68423_1316/users-sync/lamis/U01_Projects/AllExpDownStream/"
setwd(lnWD)
seurat.vivo1<-readRDS("hspc.combined.vivo.ann.exp1.RDS")
seurat.vivo2<-readRDS("hspc.combined.vivo.ann.exp2.RDS")
seurat.vivo.cross<-readRDS("hspc.combined.vivo.ann.cross.RDS")
seurat.vivo.supp<-readRDS("hspc.combined.vivo.ann.supp.RDS")

seurat.vitro1<-readRDS("hspc.combined.vitro.ann.exp1.RDS")
seurat.vitro2<-readRDS("hspc.combined.vitro.ann.exp2.RDS")
seurat.vitro.cross<-readRDS("hspc.combined.vitro.ann.cross.RDS")
seurat.vitro.supp<-readRDS("hspc.combined.vitro.ann.supp.RDS")

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
Sample.Trace.cross<-c("OO_C","OY_C","YY_C","YO_C")
Sample.Trace.supp<-c("O-Veh-T-1","Y-Veh-T-2","Y-Veh-T-3","O-DQ-T-1","O-Veh-T-2","O-Veh-T-3","Y-Veh-T-1","O-DQ-T-2","O-DQ-T-3")

##########
seurat<-seurat.vitro.cross
sample<-rep.vitro[i]
Hasthtag=unique=1

Entropy<-function(seurat,sample,Hasthtag,unique){
  if(Hasthtag==1){
  seurat<-subset(seurat,RepExtended==sample)}
  if(Hasthtag==0){
    seurat<-subset(seurat,orig.ident==sample)}
  if (unique==1){
  clone.id<-names(table(seurat$uniqueClonesTraced))
  freq<-as.numeric(table(seurat$uniqueClonesTraced))}
  if (unique==0){
    clone.id<-names(table(seurat$CloneID))
    freq<-as.numeric(table(seurat$CloneID))}
  i<-which(clone.id=="0")
  clone.id<-clone.id[-i]
  freq<-freq[-i]
  p<-freq/sum( freq)
  S<-sum(-p*log(p))/log(length(clone.id))
  return(S)
}

Entropy.vitro1<-Entropy.vitro2<-Entropy.vitro.cross<-Entropy.vitro.supp<-Entropy.vivo1<-Entropy.vivo2<-Entropy.vivo.cross<-Entropy.vivo.supp<-NULL
for (i in 1:length(Sample.vitro1)){
  Entropy.vitro1[i]<-Entropy(seurat.vitro1,Sample.vitro1[i],0,1)}

for (i in 1:length(Sample.vitro2)){
  Entropy.vitro2[i]<-Entropy(seurat.vitro2,Sample.vitro2[i],0,0)}

rep.vitro<-c("O_Vitro_Rep1","O_Vitro_Rep2","O_Vitro_Rep3","Y_Vitro_Rep1","Y_Vitro_Rep2","Y_Vitro_Rep3")
for (i in 1:length(rep.vitro)){ 
  Entropy.vitro.cross[i]<-Entropy(seurat.vitro.cross,rep.vitro[i],1,1)}

for (i in 1:length(Sample.vitro.supp)){
  Entropy.vitro.supp[i]<-Entropy(seurat.vitro.supp,Sample.vitro.supp[i],0,0)}

for (i in 1:length(Sample.vivo1)){
  Entropy.vivo1[i]<-Entropy(seurat.vivo1,Sample.vivo1[i],0,1)}

for (i in 1:length(Sample.vivo2)){
  Entropy.vivo2[i]<-Entropy(seurat.vivo2,Sample.vivo2[i],0,0)}

rep.vivo<-names(table(seurat.vivo.cross$RepExtended))
rep.vivo<-rep.vivo[-which(rep.vivo=="UnMapped")]
for (i in 1:length(rep.vivo)){ 
  Entropy.vivo.cross[i]<-Entropy(seurat.vivo.cross,rep.vivo[i],1,1)}

for (i in 1:length(Sample.vivo.supp)){
  Entropy.vivo.supp[i]<-Entropy(seurat.vivo.supp,Sample.vivo.supp[i],0,0)}
 ######################

Sample.vivo22<-c("1Y1B","2Y2B" ,"3Y3B", "4O1B" ,"5O2B", "6O3B")
rep.vivo1<-c("2OO_Rep1", "2OO_Rep2" ,"2OO_Rep3" ,"4OY_Rep1" ,"4OY_Rep2" ,"4OY_Rep3", "3YO_Rep1", "3YO_Rep2" ,"3YO_Rep3",
             "1YY_Rep1", "1YY_Rep2", "1YY_Rep3")
rep.vitro1<-c("2O_Vitro_Rep1", "2O_Vitro_Rep2", "2O_Vitro_Rep3", "1Y_Vitro_Rep1", "1Y_Vitro_Rep2" , "1Y_Vitro_Rep3")
Sample.vivo.supp1<-c("2O-Veh-1" ,"1Y-Veh-2", "1Y-Veh-3", "3O-DQ-1" ,
                    "2O-Veh-2", "2O-Veh-3" ,"1Y-Veh-1" ,"3O-DQ-2" , "3O-DQ-3" )
data<-data.frame(Age=c(rep("1Young",3),rep("2Old",3),"1Young","2Old",rep("1Young",3),rep("2Old",3),
                       rep("1Young",3),rep("2Old",3),rep("2Old",3),rep("1Young",3),rep("2Old",3),
                       rep("4Old-Young",3),rep("3Young-Old",3),rep("1Young",3),"2Old","2Old",rep("4Old-Young",2),"5Old-DQ",
                       rep("2Old",2),"4Old-Young",rep("5Old-DQ",2)),
                 Entropy=c(Entropy.vitro1,Entropy.vivo1,Entropy.vitro2,Entropy.vivo2,Entropy.vitro.cross,
                           Entropy.vivo.cross,Entropy.vitro.supp,Entropy.vivo.supp),
                 Sample=c(Sample.vitro1,Sample.vivo1,Sample.vitro2,Sample.vivo22,rep.vitro1,rep.vivo1,Sample.vitro.supp,Sample.vivo.supp1),
                 Stage=c(rep("1pre.exp1",6),rep("2post.exp1",2),
                         rep("3pre.exp2",6),rep("4post.exp2",6),
                         rep("5pre.cross",6),rep("6post.cross",12),
                         "7pre.supp",rep("8post.supp",9)))

p_pre1  <- t.test(Entropy ~ Age, data = subset(data, Stage == "1pre.exp1"),alternative = "less")$p.value
# p_post <- t.test(Entropy ~ Age, data = subset(data, Stage == "2post.exp1"),alternative = "less")$p.value
p_pre2  <- t.test(Entropy ~ Age, data = subset(data, Stage == "3pre.exp2"),alternative = "less")$p.value
p_post2 <- t.test(Entropy ~ Age, data = subset(data, Stage == "4post.exp2"),alternative = "less")$p.value
p_pre3 <- t.test(Entropy ~ Age, data = subset(data, Stage == "5pre.cross"),alternative = "less")$p.value
p_post3 <- t.test(Entropy ~ Age, data = subset(subset(data, Stage == "6post.cross"),Age%in%c("1Young","2Old")),alternative = "less")$p.value
p_post.cross <- t.test(Entropy ~ Age, data = subset(subset(data, Stage == "6post.cross"),Age%in%c("3Young-Old","4Old-Young")),alternative = "less")$p.value
p_post.supp <- t.test(Entropy ~ Age, data = subset(subset(data, Stage == "8post.supp"),Age%in%c("2Old","4Old-Young")))$p.value
p_post.supp2 <- t.test(Entropy ~ Age, data = subset(subset(data, Stage == "8post.supp"),Age%in%c("2Old","5Old-DQ")))$p.value

p<-ggplot(data) +
  geom_col(
    aes(x = Stage, y = Entropy, color = Age, group = Sample), 
    fill = "white",
    position = position_dodge(width = 0.8),  # wider dodge = more space between groups
    width = 0.6,  # narrower bars = more space between bars
    size = 1
  ) +
  scale_color_manual(values = c("1Young" = "grey50", "2Old" = "chocolate","3Young-Old"="dodgerblue4","4Old-Young"="#65A593","5Old-DQ"="red4")) +
  annotate("text", x = 1, y = max(data$Entropy) * 1.1,
           label = paste0("P = ", signif(p_pre1, 1)), size = 4) +
  annotate("text", x = 3, y = max(data$Entropy) * 1.1,
           label = paste0("P = ", signif(p_pre2, 1)), size = 4) +
  annotate("text", x = 4, y = max(data$Entropy) * 1.1,
           label = paste0("P = ", signif(p_post2, 1)), size = 4) +
  annotate("text", x = 5, y = max(data$Entropy) * 1.1,
           label = paste0("P = ", signif(p_pre3, 1)), size = 4) +
  
  annotate("text", x = 5.75, y = max(data$Entropy) * 1.1-0.05,
           label = paste0("P = ", signif(p_post3, 1)), size = 4) +
  annotate("text", x = 6.25, y = max(data$Entropy) * 1.1+0.05,
           label = paste0("P = ", signif(p_post.cross, 1)), size = 4) +
  
  annotate("text", x = 7.8, y = max(data$Entropy) * 1.1-0.05,
           label = paste0("P = ", signif(p_post.supp, 1)), size = 4) +
  annotate("text", x = 8.2, y = max(data$Entropy) * 1.1+0.05,
           label = paste0("P = ", signif(p_post.supp2, 1)), size = 4) +
  
  theme_classic(base_size = 17) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

out_dir <- "/project2/sli68423_1316/users/Kailiang/Test_Rcode/U1 test code/1Age/output"
dir.create("Figures")
graphics.off()

outfn <- file.path(out_dir, "Figures/cloneDistributionEntropy.pdf")
pdf(outfn,width = 12,height = 4)
print(p)
dev.off()
######-----------------------------