
### Documentation
### Github/Shark



import os
import sys
import rpy2
import rpy2.robjects as robjects
import rpy2.rpy_classic as rpy
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
rpy.set_default_mode(0)
import itertools

import urllib

r = robjects.r
pvalue = 0.05





## Read count table and filter using R


r.library("edgeR")
r('data<-read.csv("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Project/CountTable/gene_counts.csv", row.names = 1, header = TRUE)')
r('nf <- calcNormFactors(data)')
r('y <- voom(data, plot=FALSE, lib.size=colSums(data)*nf)')
r('datacpm  <- rowSums(cpm(data) >= 2) > 38') ## Filter on Cpm where 2 in more than 50% of samples
#r('datacpm2  <- rowSums(cpm(y) >= 2) > 39') ## Filter on Cpm where 2 in more than 50% of samples
r('datar <- y[datacpm,]')

#--------------------------





## PCA on PKDVsWT


r('library("ggbiplot")')
r('library("ggplot2")')

##  Read targets file to select samples


r('targets <- read.table("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/PCA_Targets/PKD_Profile_adultOnly.txt", header = TRUE)')
r('samplesIds <- rownames(targets)')
r('samples <- targets[,2]')
r('samples2 <- targets[,4]')
r('datar_selected <- datar[, samplesIds]')
print r('colnames(datar_selected)')
r('pdf("SampleGraph.pdf",width=7,height=5)')
r('ir.pca <- prcomp(t(datar_selected$E), center = TRUE, scale. = TRUE)')
r('g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, groups = samples, ellipse = TRUE, circle = TRUE, choices = 1:2, varname.size = 0, var.axes = TRUE)')
r('print(g)')
r('dev.off()')


##Extract PC1 genes
r('loadings <- ir.pca$rotation[,1:3]')
r('x<-loadings[order(loadings[,1]),]')
r('y<-loadings[order(-loadings[,1]),]')
print r('newGenesNeg <- rownames(head(x[50:750,1:3],1000))')
print r('newGenesPos <- rownames(head(y[50:750,1:3],1000))')

r('newGenes <- rbind(newGenesNeg, newGenesPos)')

print r('newPCA <- datar_selected[newGenes,]')
r('pdf("SampleGraph2.pdf",width=7,height=5)')
r('ir.pca <- prcomp(t(newPCA$E), center = TRUE, scale. = TRUE)')
r('g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, ,   group=interaction(samples2), ellipse = TRUE, circle = TRUE, choices = 1:2, varname.size = 0, var.axes = TRUE)')
r('print(g)')
r('dev.off()')



#### Read different treatments

r('targets <- read.table("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/PCA_Targets/PKD_Profile_adultCurcumin.csv", header = TRUE)')
r('samplesIds <- rownames(targets)')
r('samples <- targets[,2]')
r('samples2 <- targets[,4]')
r('datar_selected <- datar[newGenes, samplesIds]')

##PCA

r('pdf("SampleGraphTreatment.pdf",width=7,height=5)')
r('ir.pca <- prcomp(t(datar_selected$E), center = TRUE, scale. = TRUE)')

for ii in range(2,20):
	command = 'g1 <- ggbiplot(ir.pca,  obs.scale = 1, var.scale = 1, ,   groups=samples, ellipse = TRUE, circle = TRUE, choices = c(1,' + str(ii) + '), varname.size = 0, var.axes = TRUE)'
	r(command)
	r('print(g1)')
r('dev.off()')

### Try after selection of treatment genes



r('loadingsTreatment <- ir.pca$rotation[,1:4]')
r('x<-loadingsTreatment[order(loadingsTreatment[,4]),]')
r('y<-loadingsTreatment[order(-loadingsTreatment[,4]),]')
print r('newGenesNeg <- rownames(head(x[0:5,],1000))')
print r('newGenesPos <- rownames(head(y[0:5,],1000))')
r('newGenes <- rbind(newGenesNeg, newGenesPos)')
r('datar_selected <- datar[newGenes, samplesIds]')


r('ir.pca <- prcomp(t(datar_selected$E), center = TRUE, scale. = TRUE)')
r('pdf("SampleGraphTreatment2.pdf",width=7,height=5)')

for ii in range(2,20):
	command = 'g1 <- ggbiplot(ir.pca,  obs.scale = 0.5, var.scale = 0.5 , ellipse = TRUE, circle = TRUE, groups = samples2, choices = c(1,' + str(ii) + '), varname.size = 1, labels.size = 3, var.axes = TRUE) + geom_point(aes(shape=samples))'
	r(command)
	r('print(g1)')
r('dev.off()')


#### Connect to bioMart localDB

## Read biomaRt DB into R

biomartDB = r('biomartDB <- read.csv("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/biomArt_results2.csv", header = TRUE) ')

r('newPCA$E <- cbind(ensembl_gene_id = rownames(newPCA$E), newPCA$E)')
r('rownames(newPCA$E) <- NULL')

total = r('total <- merge(newPCA$E,biomartDB,by="ensembl_gene_id")')

genes = r('genes <- (total[,31])')   ## 31 is the column that contains the gene symbols

geneSymbols = []

for ii in range(1, len(genes)+1):
	geneSymbols.append(str(genes.rx2(ii)).split("]")[1].split("\n")[0].strip(" ").upper())




### Read signature genes




	






