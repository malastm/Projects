
### Documentation
### Github/Shark



import os
import sys
import rpy2
import rpy2.robjects as robjects
import rpy2.rpy_classic as rpy
rpy.set_default_mode(0)
import itertools

import urllib

r = robjects.r
pvalue = 0.05

### Get data from db2db

def getData(mylist, myoutput, genesDEdata):
	ids = ""
	Idtype =""
	for ii in range (0, len(mylist)): 
		Idtype = "EnsemblGeneID" 
		ids = ids + "," + mylist[ii].strip().split(".")[0]
	ids = ids[1:]
	url = "http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi?method=db2db&input=" + Idtype+ "&inputValues=" + ids + "&format=row&&outputs=homolog-humangeneid,geneinfo,keggpathwayinfo,go-biologicalprocess"
	#print url
	u = urllib.urlopen(url)
	response = u.read()
	for ii in range(0, response.count("}")):
		responsex = response.split("}")[ii][1:]
		print responsex
		homolog =  responsex.split(",")[1].split(":")[1].strip('"').split('"')[1]
		#homolog=""
		if responsex.split(",")[2].find("Gene Symbol")>=0:
			genesymbol =  responsex.split(",")[2].split("Gene Symbol")[1].split("]")[0].split(":")[1].strip()
			description =  responsex.split(",")[2].split("Description")[1].split("]")[0].split(":")[1].strip()
						#print description				
		else:				
			genesymbol = "-"
			description = "-"
		#kegg = responsex.split('KEGG Pathway Info":')[1].split(",")[0].strip()
		#BP = responsex.split('GO - Biological Process":')[1].strip().split("}")[0].strip()
		myoutput.write(genesymbol + "\t" + description + "\t" + genesDEdata[mylist[ii]])







## Read count table and filter using R


r.library("edgeR")
r('data<-read.csv("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Project/CountTable/gene_counts.csv", row.names = 1, header = TRUE)')
r('nf <- calcNormFactors(data)')
r('y <- voom(data, plot=TRUE, lib.size=colSums(data)*nf)')
r('datacpm  <- rowSums(cpm(data) >= 2) > 39') ## Filter on Cpm where 2 in more than 50% of samples
r('data2 <- y$E[datacpm,]')
r('data2[,-1]<-round(data2[,-1],3)')
r('write.table(data2, file="/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Project/CountTable/gene_counts_n2.csv", sep = "\t")')


## Read experimental design

  


for file in os.listdir(sys.argv[1]):
	if file.find("csv") >=0:
		r.library("limma")
		data = r(' datar <- read.table("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Project/CountTable/gene_counts_n2.csv", header = TRUE , row.names = 1)')
		print r('colnames(datar)')
		#print r('rownames(datar)')
		targets_path = "targets <- read.table(" + '"' + sys.argv[1]  + "/" + file + '"' + ', header = TRUE)'	
		targets = r(targets_path)
		targets_sorted = r('targets.sort<-targets[order(targets[,1]),]') ## Sort by sample id 
		data_selected = r('  datar_selected <- datar[, which(names(datar) %in% targets.sort[,1])]') ## Select data from the complete counts table

		targetsNames =  r('colnames(targets.sort)')


		print r('colnames(datar_selected)')  
		print file

		## design sheet
		x = "paste("
		for item in targetsNames:
			if item.find("sample_ID") < 0:
				x = x + "targets.sort$" + item + ", "

		x  = "TS <- " + x + 'sep="_")'
		print x
		TS = r(x)  
		r('TS1 <- factor(TS, levels=c(unique(TS)))')
		r('design2 <- model.matrix(~0+TS1)')
		r('rownames(design2) <- rownames(targets.sort)')
		r('colnames(design2) <- levels(TS1)')
		TSUnique = r('unique(TS)')
	
		## Create contrasts
		contrastsAll = ""
		for item in itertools.combinations(TSUnique, 2):
			if item[0].find("Sham") >= 0 or  item[0].find("sham") >= 0 or  item[0].find("WT") >= 0 or  item[0].find("N") >= 0 :  ## Have WT in the beginning of the contrast for constant meaning of foldChange
				contrast =  item[1].replace(".","").replace("_", "").strip() + "Vs" + item[0].replace(".","").replace("_", "").strip() + "=" + item[1].strip() + "-" + item[0].strip()
			else:
				contrast = item[0].replace(".","").replace("_", "").strip() + "Vs" + item[1].replace(".","").replace("_", "").strip() + "="+ item[0].strip() + "-" + item[1].strip()
			contrastsAll = contrastsAll +  contrast + ", "
		contrastsCommand = "cont.matrix <- makeContrasts(" + contrastsAll +  " levels = design2)"
		r(contrastsCommand)


		## Fit model
		print r('design2')
		#print r('datar_selected')
		r('fit <- lmFit(datar_selected, design2)')
		r('fit2 <- contrasts.fit(fit, cont.matrix)')
		r('fit2e <- eBayes(fit2)')

	


		command = "mkdir " +    sys.argv[1] + "/" + file.split(".")[0]

		#path = "/home/tareq/Documents/leid/Work/PKD_Data/AutomatedData/" + type + "/" + geo + "/Results/"
		os.system(command)
		### Write the output of each contrast in its file
	 	for ii in range (0, contrastsAll.count(',')):
			print contrastsAll.split(",")[ii].strip().split("=")[0]
		
			toptable = 'toptable <- topTable(fit2e, number=nrow(datar_selected), adjust="fdr", ' + str(ii + 1) + ")"
			r(toptable)
			sum1 = 'sum(toptable$adj.P.Val<=  ' + str(pvalue) + ')'
			filter1 =' filter1 <- toptable[toptable$adj.P.Val<=  ' + str(pvalue) + ", ]"   ## Filter based on adj.pvalue
			r(filter1)


			directions = ["up", "down"]
		
			### Retrieve gene id using BiomaRt 
			r('library(biomaRt)')
			r('ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")')
			for item in directions:
				if item =="up":
					filter2 =' filter2 <- filter1[filter1$logFC >=0 , ]'   ## Filter based on logFC sign.
					r(filter2)
					#print r('filter2')
					#print r(sum1)

					output = 'write.table(filter2,  ' + '"' + sys.argv[1] + "/" + file.split(".")[0]+ '/' + contrastsAll.split(",")[ii].strip().split("=")[0].strip() + '_Up.txt", sep="\t")'
					r(output)
					genes = r('genes <- rownames(filter2)')
					r('filter2 <- cbind(ensembl_gene_id = rownames(filter2), filter2)')
					r('rownames(filter2) <- NULL')
					r('ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")')

					if len(genes)>0:					
					
						r('result <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "description"), filters ="ensembl_gene_id", values = genes, mart =ensembl)')
						#print r('result')
						r('total <- merge(filter2,result,by="ensembl_gene_id")')
						output = 'write.table(total,  ' + '"' + sys.argv[1] + "/" + file.split(".")[0]+ '/' + contrastsAll.split(",")[ii].strip().split("=")[0].strip() + '_Up_edited.txt", sep="\t")'
						r(output)
				else:
					filter2 =' filter2 <- filter1[filter1$logFC <0 , ]'   ## Filter based on logFC sign.
					r(filter2)
					#print r('filter2')


					output = 'write.table(filter2,  ' + '"' + sys.argv[1] + "/" + file.split(".")[0]+ '/' + contrastsAll.split(",")[ii].strip().split("=")[0].strip() + '_Down.txt", sep="\t")'
					r(output)
					genes = r('genes <- rownames(filter2)')
					r('filter2 <- cbind(ensembl_gene_id = rownames(filter2), filter2)')
					r('rownames(filter2) <- NULL')
					r('ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")')
					if genes >0:
						r('result <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "description"), filters ="ensembl_gene_id", values = genes, mart =ensembl)')
						#print r('result')
						r('total <- merge(filter2,result,by="ensembl_gene_id")')
						output = 'write.table(total,  ' + '"' + sys.argv[1] + "/" + file.split(".")[0]+ '/' + contrastsAll.split(",")[ii].strip().split("=")[0].strip() + '_Down_edited.txt", sep="\t")'
						r(output)
				
				

				
				
