
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



### Read PC1 loadings into positive and negative groups

positive = open("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/loadingsPositive.csv", 'r+')
negative = open("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/loadingsNegative.csv", 'r+')


count = 0  ## Get top 500

pos  = []
neg = []



for line in positive:
	count = count  + 1
	if count <= 55000:
		pos.append(line.split(",")[0].strip('"').strip())

count = 0 
for line in negative:
	count = count  + 1
	if count <= 55000:
		neg.append(line.split(",")[0].strip('"').strip())	
	
	
#print neg



path = "/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Single_Comparisons/Run1July/"


filenames = next(os.walk("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Single_Comparisons/Run1July/"))[1]

for item in filenames:
	print item
	for subfile in  os.listdir("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Single_Comparisons/Run1July/" +  item):
		if subfile.find("edited") >= 0:
			f2 = open("/home/tareq/Documents/leid/Work/RNASeq_TranCyst/Single_Comparisons/Run1July/" +  item + "/" + subfile, 'r+')
			genesinSubfile = []
			for line in f2:
				genesinSubfile.append(line.split("\t")[1].strip('"'))
			#print genesinSubfile
			print subfile + "\t" + str(len(genesinSubfile)) + "\t" + str(len(set(genesinSubfile).intersection(pos))) + "\t" + str(len(set(genesinSubfile).intersection(neg))) 


		
