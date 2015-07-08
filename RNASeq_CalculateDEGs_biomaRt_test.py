
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
r('y <- voom(data, plot=FALSE, lib.size=colSums(data)*nf)')
r('datacpm  <- rowSums(cpm(y) >= 2) > 39') ## Filter on Cpm where 2 in more than 50% of samples
r('data2 <- y[datacpm,]')
#r('data2[,-1]<-round(data2[,-1],3)')
#r('write.table(data2, file="/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Project/CountTable/gene_counts_n2.csv", sep = "\t")')




## Read experimental design

  


for file in os.listdir(sys.argv[1]):
	if file.find("csv") >=0 and file.find("lock")<0:
		r.library("limma")
		print file
		#data = r(' datar <- read.table("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Project/CountTable/gene_counts_n2.csv", header = TRUE , row.names = 1)')
		r('datar<- data2')
		print r('colnames(datar)')
		#print r('rownames(datar)')
		targets_path = "targets <- read.table(" + '"' + sys.argv[1]  + "/" + file + '"' + ', header = TRUE)'	
		targets = r(targets_path)
		targets_sorted = r('targets.sort<-targets[order(targets[,1]),]') ## Sort by sample id 
		print r('targets.sort')
		data_selected = r('  datar_selected <- datar[,  targets.sort[,1]]') ## Select data from the complete counts table
		print r('datar_selected')
		
