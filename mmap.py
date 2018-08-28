#!/usr/bin/env pypy3 python3

import os
import glob
import time
import math
from itertools import groupby
from multiprocessing import Pool
import re
from collections import OrderedDict
import shutil
import sys
wd = os.getcwd()

#----------------------------------------------------------------------------------
#///////// mMap - Mouse Genetype to Phenotype Mapping scripts /////////////////////
#----------------------------------------------------------------------------------
	
def gene_file(genefile):

	global p2
	func = []
	with open(genefile, 'rU') as ge:
		for g in ge:
			g = g.rstrip().split("\t")
			p = ge.name
			p1 = p.split('.')
			p2 = p1[0].split("_")
			with open(os.path.join(wd+'/'+"Functional datafiles"+"/"+"Functional-data.txt")) as rr:
				for r in rr:
					r = r.rstrip().split("\t")
					if g[0] == r[-1] and int(g[1]) >= int(r[2]) and int(r[3]) >= int(g[1]):
							func.append('\t'.join((g[0], r[1], g[1], r[4])))

		# PPF = Protein Functional Region, MGL = Mutated Genomic Region, MAA = Muated Amino Acid, D-PFF = Description of PFF
		header = ['Gene', 'PFF', 'MAA', 'D-PFF']
		with open(p2[0]+"_functional-accessment.txt", 'a') as ii:
			ii.write('\t'.join(header) + '\n')
			for i in func:
				ii.write(i+'\n')

def downstream(func_file): 
	phen = []
	with open(func_file, 'rU') as gg:
		for g in gg:
			g = g.rstrip().split()
			p = gg.name
			p1 = p.split('.')
			p2 = p1[0].split("_")
			with open(os.path.join(wd+'/'+"Functional datafiles"+"/"+"Phenotype.txt")) as rr:  
				for r in rr:
					d =  r.rstrip().split("\t")
					if d[0] == g[0]:
						phen.append('\t'.join(d))

	with open(p2[0]+"_phenotype_accessment.txt", 'a') as ii:
		for i in phen:
			ii.write(i+'\n')

def network(functional):

	import json

	import urllib.request

	from urllib.request import urlopen

	with open(functional, 'rU') as func:
		p = func.name
		p1 = p.split('.')
		p2 = p1[0].split("_")
		enrich = []
		k = [(line.split())[0] for line in func]
		for i in k:
			enrich.append(i)

		string_api_url = "https://string-db.org/api"
		output_format = "json"
		method_enrich = "enrichment"
		species = "10090"
		my_app = "www.aweseome_app.org"
		caller_identity = 'aarslan'
		request_url = string_api_url + "/" + output_format + "/" + method_enrich + "?"
		request_url += "identifiers=" + "%0d".join(enrich[1:])
		request_url += "&" + "species=" + species
		request_url += "&" + "caller_identity=" + my_app
		response = urllib.request.urlopen(request_url)
		result = response.read()
		if result:
				data = json.loads(result.decode('utf-8'))
				for row in data:
					term = row["term"]
					preferred_names = ",".join(row["preferredNames"])
					fdr = row["fdr"]
					description = row["description"]
					with open(p2[0]+'_enrichment.txt', 'a') as en:
						en.write("\t".join([term, str(fdr), description, preferred_names])+'\n')

		output_format_figure = "image"
		method = "network"
		request_urlf = string_api_url + "/" + output_format_figure + "/" + method + "?"
		request_urlf += "identifiers=%s"
		request_urlf += "&" + "species=" + species
		request_urlf += "&" + "add_white_nodes=0"
		request_urlf += "&" + "network_flavor=actions"
		request_urlf += "&" + "caller_identity=aarslan" + my_app	
		c = "%0d".join(enrich[1:])
		urllib.request.urlretrieve(request_urlf % c, "%s.png" % "".join(p2[0]+"_Protein_Network"))


def noncoding(genefile):

	global p1
	c = []
	with open(genefile, 'rU') as mut:
        		for a in mut:
        			a = a.rstrip().split('\t')
        			p = mut.name
        			p1 = p.split('.')
        			c.append('%s %s' % (a[0], a[1]))
	for m in c:
		m = m.rstrip().split()
		vista = []
		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'Vista_enhancers_flankinGenes.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and m[1] >= v[1] and v[2] >= m[1]:
					vista.append('%s %s %s %s %s' % (v[0], m[1], v[-2], v[-1], 'Enhancer'))
					with open(p1[0]+"_regulatory-analysis.txt", 'a') as k:
						k.write('\t'.join(vista)+'\n')
		pro = []
		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'EPDnew_promoter.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and m[1] >= v[1] and v[2] >= m[1]:
					pro.append('%s %s %s %s' % (v[0], m[1], v[3], 'Promoter'))
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join(pro)+'\n')
		tss= [] 
		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'TSS_promoter.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and m[1] == v[2]:
					tss.append('%s %s %s %s %s' % (v[0], m[1], v[1], v[-1], 'TSS'))
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join(tss)+'\n')
		cpg = []	
		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'cpgIsland.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and m[1] >= v[1] and v[2] >= m[1]:
					cpg.append('%s %s %s %s' % (v[0], m[1], v[-1], 'CpG'))
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join(cpg)+'\n')
		ctcf = []
		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'insulator-CTCF-binding-sites.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and m[1] >= v[1] and v[2] >= m[1]:
					ctcf.append('%s %s %s %s %s' % (v[0], m[1], v[-2], v[-1], 'Insulator'))
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join(ctcf)+'\n')
						#print(ctcf)


def ncnetwork(functional):

	import json

	import urllib.request

	from urllib.request import urlopen

	with open(functional, 'rU') as nc:
		p = nc.name
		p1 = p.split('.')
		p2 = p1[0].split("_")
		enrich = []
		for n in nc:
			n = n.split()
			if n[-1] == 'Enhancer':
				enrich.append(n[-3] +'\t'+n[-2])
			if n[-1] == 'Promoter' or n[-1] == 'TSS':
				take1 = n[-2].split('_')
				enrich.append(take1[0])
	phen = []
	with open(os.path.join("mgi-phenotype.txt")) as rr:  
		for r in rr:
			d =  r.rstrip().split("\t")
			for g in enrich:
				g = g.split()
				if d[0] == g[0]:
					phen.append('\t'.join(d))

	with open(p1[0]+"_phenotype_accessment.txt", 'a') as ii:
		for i in phen:
			ii.write(i+'\n')
		
		i = [(line.split())[0] for line in enrich]
		string_api_url = "https://string-db.org/api"
		output_format = "json"
		method_enrich = "enrichment"
		species = "10090"
		my_app = "www.aweseome_app.org"
		caller_identity = 'aarslan'
		request_url = string_api_url + "/" + output_format + "/" + method_enrich + "?"
		request_url += "identifiers=" + "%0d".join(i[1:])
		request_url += "&" + "species=" + species
		request_url += "&" + "caller_identity=" + my_app
		response = urllib.request.urlopen(request_url)
		result = response.read()
		if result:
				data = json.loads(result.decode('utf-8'))
				for row in data:
					term = row["term"]
					preferred_names = ",".join(row["preferredNames"])
					fdr = row["fdr"]
					description = row["description"]
					#if fdr < 0:
					with open(p2[0]+'_enrichment.txt', 'a') as en:
						en.write("\t".join([term, str(fdr), description, preferred_names])+'\n')
		else:
			with open(p2[0]+'_enrichment.txt', 'a+') as en:
				en.write("No significant enrichment detected")

		output_format_figure = "image"
		method = "network"
		request_urlf = string_api_url + "/" + output_format_figure + "/" + method + "?"
		request_urlf += "identifiers=%s"
		request_urlf += "&" + "species=" + species
		request_urlf += "&" + "add_white_nodes=0"
		request_urlf += "&" + "network_flavor=actions"
		request_urlf += "&" + "caller_identity=aarslan" + my_app	
		c = "%0d".join(i[1:])
		urllib.request.urlretrieve(request_urlf % c, "%s.png" % "".join(p2[0]+"_Protein_Network"))

def nocode_mutation(file1):
	try:
		noncoding(os.path.join(file1))

		ncnetwork(os.path.join(p1[0]+"_regulatory-analysis.txt"))

	except IOError:
		pass
	try:
		os.system("mkdir "+wd+"/"+p1[0]+"_mutations")
		shutil.move(wd+"/"+p1[0]+"_regulatory-analysis.txt", wd+"/"+p1[0]+"_mutations")
		shutil.move(wd+"/"+p1[0]+"_Protein_Network.png", wd+"/"+p1[0]+"_mutations")
		shutil.move(wd+"/"+p1[0]+"_enrichment.txt", wd+"/"+p1[0]+"_mutations")
		shutil.move(wd+"/"+p1[0]+"_regulatory-analysis_phenotype_accessment.txt", wd+"/"+p1[0]+"_mutations")
	except IOError:
		pass

	return "Done"

def gene_mutation(genefile):

	try:
			gene_file(os.path.join(genefile))
			
			network(os.path.join(p2[0]+"_functional-accessment.txt"))

			downstream(os.path.join(p2[0]+"_functional-accessment.txt"))

	except IOError:
			pass

	try:
		os.system("mkdir "+wd+"/"+p2[0]+"_mutations")
		shutil.move(wd+"/"+p2[0]+"_functional-accessment.txt", wd+"/"+p2[0]+"_mutations")
		shutil.move(wd+"/"+p2[0]+"_Protein_Network.png", wd+"/"+p2[0]+"_mutations")
		shutil.move(wd+"/"+p2[0]+"_enrichment.txt", wd+"/"+p2[0]+"_mutations")
		shutil.move(wd+"/"+p2[0]+"_phenotype_accessment.txt", wd+"/"+p2[0]+"_mutations")
	except IOError:
		pass
	return "Done"


if __name__=='__main__':
	
	import argparse
	p = Pool()
	sys.setrecursionlimit(2000)

	parser = argparse.ArgumentParser()

	parser.add_argument('-g','--genes', action='append', help='performs the msMap on genes level mutation positions')
	parser.add_argument('-nc', '--noncoding', action='append', help='performs the mMap on proteins level mutation positions')

	args = parser.parse_args()
	if args.genes:
		try:
			genefile =  [f for f in glob.glob(os.path.join("*.txt")) if f.endswith(sys.argv[-1])]
			p.map(gene_mutation, genefile)
			p.close()
			p.join()
		except IOError:
			pass
	elif args.noncoding:
		try:
			gene1 =  [f for f in glob.glob(os.path.join("*.txt")) if f.endswith(sys.argv[-1])]
			p.map(nocode_mutation, gene1)
			p.close()
			p.join()
		except IOError:
			pass
	else:
		print ("to run a function seek help")
