
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
r('datar <- y[datacpm,]')
#r('data2[,-1]<-round(data2[,-1],3)')
#r('write.table(datar, file="/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Project/CountTable/gene_counts_n2.csv", sep = "\t")')


## Read biomaRt DB into R

biomartDB = r('biomartDB <- read.csv("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/biomArt_results2.csv", header = TRUE) ')


## Read experimental design

for file in os.listdir(sys.argv[1]):
	if file.find("csv") >=0 and file.find("lock")<0:
		r.library("limma")
		#print r('colnames(datar)')
		#print r('rownames(datar)')
		targets_path = "targets <- read.table(" + '"' + sys.argv[1]  + "/" + file + '"' + ', header = TRUE)'	
		targets = r(targets_path)
		targets_sorted = r('targets.sort<-targets[order(targets[,1]),]') ## Sort by sample id
		r('rownames(targets.sort) <- targets.sort[,1]')## Change column names 
		#print r('targets.sort[,1]')
		data_selected = r('  datar_selected <- datar[, rownames(targets.sort)]') ## Select data from the complete counts table
		#print r('datar_selected')
		targetsNames =  r('colnames(targets.sort)')

		print file

		## design sheet
		x = "paste("
		for item in targetsNames:
			if item.find("sample_ID") < 0:
				x = x + "targets.sort$" + item + ", "

		x  = "TS <- " + x + 'sep="_")'
		#print x
		TS = r(x)  
		r('TS1 <- factor(TS, levels=c(unique(TS)))')
		r('design2 <- model.matrix(~0+TS1)')
		r('rownames(design2) <- rownames(targets.sort)')
		r('colnames(design2) <- levels(TS1)')
		TSUnique = r('unique(TS)')
	
		## Create contrasts
		contrastsAll = ""
		for item in itertools.combinations(TSUnique, 2):
			if item[0].find("Sham") >= 0 or  item[0].find("sham") >= 0 or  item[0].find("WT") >= 0 or  item[0].find("N") >= 0  or  item[0].find("PBS") >= 0:  ## Have WT in the beginning of the contrast for constant meaning of foldChange
				contrast =  item[1].replace(".","").replace("_", "").strip() + "Vs" + item[0].replace(".","").replace("_", "").strip() + "=" + item[1].strip() + "-" + item[0].strip()
			else:
				contrast = item[0].replace(".","").replace("_", "").strip() + "Vs" + item[1].replace(".","").replace("_", "").strip() + "="+ item[0].strip() + "-" + item[1].strip()
			contrastsAll = contrastsAll +  contrast + ", "
		contrastsCommand = "cont.matrix <- makeContrasts(" + contrastsAll +  " levels = design2)"
		r(contrastsCommand)


		## Fit model
		print r('design2')
		print r('colnames(datar_selected)')
		r('fit <- lmFit(datar_selected, design2)')
		r('fit2 <- contrasts.fit(fit, cont.matrix)')
		r('fit2e <- eBayes(fit2)')

	


		command = "mkdir " +    sys.argv[1] + "/" + file.split(".")[0]

		#path = "/home/tareq/Documents/leid/Work/PKD_Data/AutomatedData/" + type + "/" + geo + "/Results/"
		os.system(command)
		### Write the output of each contrast in its file
	 	for ii in range (0, contrastsAll.count(',')):
			#print contrastsAll.split(",")[ii].strip().split("=")[0]
		
			toptable = 'toptable <- topTable(fit2e, number=nrow(datar_selected), adjust="fdr", ' + str(ii + 1) + ")"
			contrastOut = r(toptable)
			#print contrastOut
			sum1 = 'sum(toptable$adj.P.Val<=  ' + str(pvalue) + ')'
			filter1 =' filter1 <- toptable[toptable$adj.P.Val<  ' + str(pvalue) + ", ]"   ## Filter based on adj.pvalue
			r(filter1)


			directions = ["up", "down"]
		
			### Retrieve gene id using BiomaRt local DB 
			for item in directions:
				if item =="up":
					filter2 =' filter2 <- filter1[filter1$logFC >=0 , ]'   ## Filter based on logFC sign.
					r(filter2)
					genes = r('genes <- rownames(filter2)')
					output = 'write.table(filter2,  ' + '"' + sys.argv[1] + "/" + file.split(".")[0]+ '/' + contrastsAll.split(",")[ii].strip().split("=")[0].strip() + '_Up.txt", sep="\t")'
					r(output)
					r('filter2 <- cbind(ensembl_gene_id = rownames(filter2), filter2)')
					r('rownames(filter2) <- NULL')

					if len(genes)>0:					
						r('total <- merge(filter2,biomartDB,by="ensembl_gene_id")')
						r('total <-total[ order(-total[,4]), ]')    ### Sort by adj. Pvalue
						r('rownames(total) <- NULL')
						output = 'write.table(total,  ' + '"' + sys.argv[1] + "/" + file.split(".")[0]+ '/' + contrastsAll.split(",")[ii].strip().split("=")[0].strip() + '_Up_edited.txt", sep="\t", col.names=NA)'
						r(output)
				else:
					

					filter2 =' filter2 <- filter1[filter1$logFC <0 , ]'   ## Filter based on logFC sign.
					r(filter2)
					genes = r('genes <- rownames(filter2)')
					output = 'write.table(filter2,  ' + '"' + sys.argv[1] + "/" + file.split(".")[0]+ '/' + contrastsAll.split(",")[ii].strip().split("=")[0].strip() + '_Down.txt", sep="\t")'
					r(output)
					r('filter2 <- cbind(ensembl_gene_id = rownames(filter2), filter2)')
					r('rownames(filter2) <- NULL')

					if len(genes)>0:					
						r('total <- merge(filter2,biomartDB,by="ensembl_gene_id")')
						r('total<- total[ order(-total[,4]), ]')    ### Sort by adj. Pvalue
						r('rownames(total) <- NULL')
						output = 'write.table(total,  ' + '"' + sys.argv[1] + "/" + file.split(".")[0]+ '/' + contrastsAll.split(",")[ii].strip().split("=")[0].strip() + '_Down_edited.txt", sep="\t" , col.names=NA)'
						r(output)
				

				

