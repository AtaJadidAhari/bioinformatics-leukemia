WD <- getwd()
setwd(WD)

library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(plyr)


# some comments are present in source.R.
#in this code, we group datas into to groups. Normal and Leukimeia



## make a new folder named "Data" in the current working directory. Put GSE48558_series_matrix.txt.gz and GPL6244.annot.gz in Data folder.
#loading the data:
gset <- getGEO( GSEMatrix =TRUE, filename = "Data/GSE48558_series_matrix.txt.gz", AnnotGPL=TRUE, destdir =  "Data/")

#zeros denote leukimia and 1's are normal. the X's are excluded from gset
gsms <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXX",
               "XXXXXXXXXXXXXXXXXX1X1XXXXX1111X1XX11XX11X1X1X1X1X1",
               "XXX1XXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111001000",
               "11111111111111111111")
AMLPatient <- "AML.Patient"
normal <- "Normal"
gr <- c(rep(AMLPatient, 13), rep(normal, 27), 
         rep(AMLPatient, 2), normal,   rep(AMLPatient, 3),rep(normal, 20)  )
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }


# eliminate samples marked as "X"
sel <- which(sml != "X")

gset <- gset[ ,sel]

ex <-exprs(gset)



pdf("Results/boxplot2.pdf", width=50)
boxplot(ex)
dev.off()

#ex <- normalizeQuantiles(ex)
#exprs(gset) <- ex

#### Corelation heatmap of 2 groups
pdf("Results/CoreHeatmapend_of_two_groups.pdf", width =  50, height = 50)
pheatmap(cor(ex),  labels_row = gr, labels_col = gr)
dev.off()





#### Principal Components

pc <- prcomp(ex)

pdf("Results/PC2.pdf")
plot(pc)
plot(pc$x[, 1:2])
dev.off()


ex_scale <- t(scale(t(ex), scale = FALSE)) #miangine gene ha sefr mishavand
pc <- prcomp(ex_scale)
pdf("Results/PC_scaled2.pdf")
plot(pc)
plot(pc$x[, 1:2])
dev.off()


############################################################ pc_scaled

pcr <- data.frame(pc$rotation[, 1:3], Group=gr)

pdf("Results/PCA_samples2.pdf")
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()



#### differenctial Expressional Analisys

gr <- factor(gr)
gset$description <- gr

design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(gr)

fit <- lmFit(gset, design )
cont.matrix <- makeContrasts(Normal - AML.Patient, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01 )

tT <- topTable(fit2, adjust="fdr", sort.by = "B", number = Inf)
tT <- subset(tT, select = c("Gene.symbol", "Gene.ID" ,"adj.P.Val", "logFC" ))
write.table(tT, "Results/AML.Patient-Normal.txt", row.names = F, sep = "\t", quote = F)

aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(aml.up$Gene.symbol)

aml.up.genes <- unique(as.character(strsplit2(aml.up.genes, "///")))
write.table(aml.up.genes, file="Results/AML.Patient-Normal_Up.txt", quote = F, row.names = F, col.names = F)


aml.down <- subset(tT, logFC < 1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, file="Results/AML.Patient-Normal_Down.txt", quote = F, row.names = F, col.names = F)
