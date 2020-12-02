#!/usr/bin/env pypy3 python3

import os
import time
import glob
import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from collections import OrderedDict
from decimal import Decimal
from scipy.stats import hypergeom
import math
import mechanize
from urllib.error import HTTPError
import numpy as np
import shutil
start_time = time.time()
wd = os.getcwd()

def gene_file(genefile):

	cons = glob.glob(wd+'/'+"Functional-datafiles"+"/"+"conversation/*")
	global p2
	func = []
	funcC = []
	with open(genefile, 'r') as ge:
		for g in ge:
			g = g.rstrip().split("\t")
			p = ge.name
			p2 = p.split('.')
			with open(os.path.join(wd+'/'+"Functional-datafiles"+"/"+"Functional-data.txt")) as rr:
				for r in rr:
					r = r.rstrip().split("\t")
					if g[0] == r[-1] and int(g[1]) >= int(r[2]) and int(r[3]) >= int(g[1]):
							func.append('\t'.join((g[0], r[1], g[1], r[4])))

		for i in func:
			i = i.rstrip().split("\t")
			for c in cons:
				c1 = c.split('/')
				c2 = c1[-1].split('.')
				if i[0] == c2[0]:
					with open(c) as cc:
						for ci in cc:
							con = ci.rstrip().split('\t')
							if int(i[2]) == int(con[0]):
								funcC.append('\t'.join((i[0], i[1], i[2], i[-1], con[-1])))

		# PPF = Protein Functional Region, MGL = Mutated Genomic Region, MAA = Muated Amino Acid, D-PFF = Description of PFF, Cons = Conservation
		header = ['Gene', 'PFF', 'MAA', 'D-PFF', "Cons"]
		with open(p2[0]+"_functional-accessment.txt", 'a') as ii:
			ii.write('\t'.join(header) + '\n')
			for fc in OrderedDict.fromkeys(funcC):
				ii.write(fc+'\n')
	
	if p2[0]+"_functional-accessment.txt":

		rat1 = pd.read_csv(p2[0]+"_functional-accessment.txt", sep='\t')
		rat = rat1.sort_values(by=['Cons'])

		e_count = []

		for index, gene in enumerate(rat['Gene']):
			subset = rat.iloc[:index + 1]
			count = len(subset[subset['Gene'] == gene])
			e_count.append(count)

		x = rat['Gene']
		y = e_count
		hue = rat['Cons']

		plt.figure(figsize=(6.5, 4.5))
		ax = sns.scatterplot(x, y, hue=hue, s=15, legend="full", palette="RdYlGn")
		ax.grid(False)  # Remove grid
		ax.get_legend().remove()  # Delete default legend
		scale_legend = plt.Normalize(hue.min() - 1, hue.max())  # Create a scale for the colormap.
		color_map = plt.cm.ScalarMappable(cmap="RdYlGn", norm=scale_legend)  # Colormap used in legend.
		color_map.set_array([])  # Dummy variable needed to create a colormap.
		cb = ax.figure.colorbar(color_map, shrink=0.40)  # Add colormap as a legend.
		cb.ax.tick_params(labelsize=3)
		plt.xticks(rotation='vertical', fontsize=3)
		plt.ylabel(" SNP count", size=7)
		plt.xlabel('Protein', fontsize=7)
		plt.yticks(fontsize=5)
		plt.gcf().text(0.85, 0.45, "Conservation*", fontsize=4, rotation=90)  # Label used for colormap.
		plt.savefig(p2[0]+'_Mutated-Proteins-at-Conserved-Regions.png', format='png', dpi=900, bbox_inches='tight')
	print(p2[0]+" conservation graph is done")

	sc = glob.glob(wd+'/'+"Functional-datafiles"+"/"+"numbers/*")
	df = rat1.drop_duplicates(subset = ["Gene"], keep='last')
	for c in sc:
		f1 = c.split('/')
		f2 = f1[-1].split('.')
		for d in df['Gene']:
			if d == f2[0]:
				with open(c) as ij:
					for i in ij:
						if not 'protein-coding gene' in i and not 'Unknown' in i:
							with open(p2[0]+"_scExpression-data.txt", 'a') as k:
								k.write(d +'\t'+ str(i))
	print(p2[0]+" sc-Expression analysis is done")

	if p2[0]+"_scExpression-data.txt":
		plt.figure(figsize=(6.5, 4.5))
		palette = np.repeat(np.array(sns.color_palette("deep")),1, axis=0)
		plt.yticks(fontsize=3)
		plt.xticks(fontsize=4)
		scx = pd.read_csv(p2[0]+"_scExpression-data.txt", sep='\t', index_col=None)
		scx.columns = ["Genes", "Cell", "Number"]
		sns.barplot(x='Number', y='Cell', data=scx, hue='Genes', orient='h',palette=palette)
		plt.ylabel(" cell-type", size=7)
		plt.xlabel("number of clusters", fontsize=7)
		sns.set_color_codes("pastel")
		plt.legend(loc='upper left', fontsize="x-small", title='Genes', title_fontsize=3, bbox_to_anchor=(1.05, 1))
		plt.savefig(p2[0]+'_scExp.png', format='png', dpi=900, bbox_inches='tight')

	print(p2[0]+" sc-Expression analysis visualisation is done")

def downstream(func_file): 
	
	global p3
	phen = []
	with open(func_file, 'r') as gg:
		for g in gg:
			g = g.rstrip().split()
			p = gg.name
			p1 = p.split('.')
			p3 = p1[0].split("_")
			with open(os.path.join(wd+'/'+"Functional-datafiles"+"/"+"Phenotype.txt")) as rr:  
				for r in rr:
					d =  r.rstrip().split("\t")
					if d[0] == g[0]:
						phen.append('\t'.join(d))

	with open(p3[0]+"_phenotype_accessment.txt", 'a') as ii:
		for i in phen:
			ii.write(i+'\n')

def testScore(f):

	global p4
	l = []
	with open(f, 'r') as ij:
		p = ij.name
		p1 = p.split('.')
		p4 = p1[0].split("_")
		for i in ij.readlines()[1:]:
			i = i.rstrip().split('\t')
			l.append('\t'.join((i[1], i[2])))
	de = {}
	lines = (line.rstrip('\t') for line in l)
	unique = OrderedDict.fromkeys( (line for line in lines if line) )
	k = [(line.split('\t'))[0] for line in unique]

	for word in k:
		de[word] = de.get(word, 0) + 1

	word_freq = []
	for key, value in de.items():
		word_freq.append(list((value, key)))

	pathway = []
	word_freq.sort(reverse=True)
	for i in word_freq:

		pv = [(	'Binding site'	,	'14'	),
		(	'Coiled coil'	,	'1337'	),
		(	'Cross-link'	,	'12'	),
		(	'DNA binding'	,	'52'	),
		(	'Disulfide bond'	,	'4117'	),
		(	'Domain'	,	'9114'	),
		(	'Glycosylation'	,	'84'	),
		(	'Helix'	,	'592'	),
		(	'Initiator methionine'	,	'2'	),
		(	'Metal binding'	,	'11'	),
		(	'Modified residue'	,	'164'	),
		(	'Motif'	,	'74'	),
		(	'Nucleotide binding'	,	'68'	),
		(	'Region'	,	'5472'	),
		(	'Repeat'	,	'2049'	),
		(	'Signal peptide'	,	'752'	),
		(	'Avtive Site'	,	'6'	),
		(	'Transit peptide '	,	'208'	),
		(	'Transmembrane'	,	'1596'	),
		(	'Turn'	,	'82'	),
		(	'Zinc finger'	,	'273'	)]
	
		for v in pv:
			if v[0] == i[1]:
				# Number of SNPs
				N = len(k)
				# functional plus regulatory SNPs
				M = 340091 
				p = f"{Decimal(hypergeom(M, int(v[1]), N).pmf(i[0])):.2E}"
				pathway.append((v[0], p))

	pathway.sort(key = lambda x: float(x[1]), reverse = False)
	
	for pa in pathway:	
		with open(p4[0]+"_Enriched-functional-regions.txt",'a') as eg:
			eg.write('\t'.join(pa)+'\n')

	if p4[0]+'_Enriched-functional-regions.txt':

			fig, ax = plt.subplots(figsize=(7, 5.5))
			x = []
			y = []
			with open(p4[0]+'_Enriched-functional-regions.txt','r') as ij:
				for i in ij.readlines():
					i = i.rstrip().split('\t')
					x.append(i[0])
					try:
						y.append(abs(math.log10(float(i[1]))))
					except ValueError:
						pass

	xt = tuple(x)
	ys, xs = zip(*sorted(zip(y, xt)))
	my_range=list(range(1,len(y)+1))
	plt.hlines(y=my_range, xmin=ys, xmax=ys[0], color='indigo', alpha=0.2, linewidth=3)
	plt.plot(ys, my_range, "o", markersize=3, color='indigo', alpha=0.6)
	plt.yticks(my_range, xs)
	ax.tick_params(axis='both', which='major', labelsize=5)
	plt.xlabel ('-log10(p-value)', size=7)
	plt.ylabel(" regions", size=7)
	fig.savefig(p4[0]+'_Enriched-functional-regions.png', format='png', dpi=900, bbox_inches='tight')

def network(functional):

	global p5

	import json

	import urllib.request

	from urllib.request import urlopen

	with open(functional, 'r') as func:
		p = func.name
		p1 = p.split('.')
		p5 = p1[0].split("_")
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
		request_url += "identifiers=" + "%0d".join(enrich[1:200])
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
					category = row["category"]
					description = row["description"]
					if fdr <= 0.01 and category == 'Process':
						with open(p5[0]+'_BiologicalProcess.txt', 'a') as biop, open(p5[0]+'_BiologicalProcess.txt', 'a') as biolp:
							biop.write("\t".join([term, str(fdr), description, preferred_names])+'\n')
							biolp.write("\t".join([term, str(fdr)])+'\n')
					if fdr <= 0.01 and category == 'RCTM' or category == 'KEGG':
						with open(p5[0]+'_Pathways.txt', 'a') as pathw:
							pathw.write("\t".join([term, str(fdr), description, preferred_names])+'\n')
					else:
						with open(p5[0]+'_enrichment.txt', 'a') as en:
							en.write("\t".join([term, str(fdr), description, preferred_names])+'\n')

		else:
			with open(p5[0]+'_enrichment.txt', 'a+') as en:
				en.write("No significant enrichment detected")					

		output_format_figure = "image"
		method = "network"
		request_urlf = string_api_url + "/" + output_format_figure + "/" + method + "?"
		request_urlf += "identifiers=%s"
		request_urlf += "&" + "species=" + species
		request_urlf += "&" + "add_white_nodes=0"
		request_urlf += "&" + "network_flavor=actions"
		request_urlf += "&" + "caller_identity=aarslan" + my_app	
		c = "%0d".join(enrich[1:200])
		urllib.request.urlretrieve(request_urlf % c, "%s.png" % "".join(p2[0]+"_Protein_Network"))

		url = "http://revigo.irb.hr/"
		mec = mechanize.Browser()
		mec.set_handle_robots(False)
		mec.open(url)
		mec.select_form(name="submitToRevigo")

		try:	

			x = []
			y = []
			
			with open(p5[0]+'_BiologicalProcess.txt', 'r') as ij:
				data = ij.read().replace('\t', ' ')
				mec["goList"] = data
				res = mec.submit()
				content = res.read()
				f = mec.retrieve('http://revigo.irb.hr/export.jsp?table=1','REVIGO.csv')[0]
				g = open(f, 'r')
				g1 = g.readlines()[1:]
				for g2 in g1:
					g2 = g2.rstrip().split(',')
					x.append(g2[1])
					y.append(abs(float(g2[6])))
					with open(p5[0]+'_BiologicalProcess_Revigo.txt', 'a') as revi:
						revi.write('\t'.join(g2)+'\n')
		except HTTPError as q:
			print(q)
		
		if not y:

			with open(p5[0]+"_BiologicalProcess.txt", 'r') as ij:
				for i in ij.readlines()[:70]:
					ii = i.rstrip().split('\t')
					if float(ii[1]) < 0.005 and len(ii) >= 3:
							y.append(abs(math.log10(float(ii[1]))))
							x.append(ii[2])

		os.remove('REVIGO.csv')
		fig, ax = plt.subplots(figsize=(7, 5.5)) #very good size
		xt = tuple(x)
		ys, xs = zip(*sorted(zip(y, xt)))
		my_range=list(range(1,len(y)+1))
		plt.hlines(y=my_range, xmin=0, xmax=ys, color='#007acc', alpha=0.2, linewidth=3)
		plt.plot(ys, my_range, "o", markersize=3, color='#007acc', alpha=0.6)
		plt.yticks(my_range, xs)
		plt.yticks(fontsize=5, rotation=0)
		plt.xlabel ('-log10(p-value)', fontsize=6)
		fig.savefig(p5[0]+'_Significantly-Altered-Biological-Process.png', format='png', dpi=900, bbox_inches='tight')

		
		if p5[0]+'_Pathways.txt':
			fig, ax = plt.subplots(figsize=(7, 5.5))
			x = []
			y = []
			with open(p5[0]+'_Pathways.txt','r') as ij:
				for i in ij.readlines()[:30]:
					i = i.rstrip().split('\t')
					x.append(i[2])
					y.append(abs(math.log10(float(i[1]))))

		xt = tuple(x)
		ys, xs = zip(*sorted(zip(y, xt)))
		my_range=list(range(1,len(y)+1))
		plt.hlines(y=my_range, xmin=ys, xmax=ys[0], color='#197213', alpha=0.2, linewidth=3)
		plt.plot(ys, my_range, "o", markersize=3, color='#197213', alpha=0.6)
		plt.yticks(my_range, xs)
		ax.tick_params(axis='both', which='major', labelsize=6)
		plt.xlabel ('-log10(p-value)', fontsize=7)
		fig.savefig(p5[0]+'_Significantly-Altered-Biological-Pathways.png', format='png', dpi=900, bbox_inches='tight')

def go():

	try:
		if p2[0]+'_enrichment.txt':

				fig, ax = plt.subplots(figsize=(9, 5.5))
				plt.tick_params(labelsize=6)
				ij = pd.read_csv(p2[0]+'_enrichment.txt', sep='\t')
				ij.columns = ['GOt', 'FDR', 'Desc', "Genes"]
				ij = ij[ij.iloc[:,0].str.contains("GO.")]
				ij = ij[(ij.iloc[:,1] <= 0.01 )]
				ij = ij.sort_values(by='FDR', ascending=True)
				ij = ij.head(50)
				df1 = ij.iloc[:,-1].str.split(',')
				fdr = round(abs(np.log10(ij['FDR'])),1)
				scatter = ax.scatter(fdr,  ij['Desc'], c=ij['FDR'], s=[i for i in range(len(ij['FDR']))])
				handles, labels = scatter.legend_elements(prop="sizes", alpha=0.3)
				legend2 = ax.legend(handles, labels, title='number of genes', fontsize=4, bbox_to_anchor=(1.15, 1))
				legend2.get_title().set_fontsize('5')
				plt.xlabel(" -log10 FDR", fontsize=7)
				fig.savefig(p2[0]+'_GO-term-Enrichment.png', format='png', dpi=900, bbox_inches='tight')
	except ValueError as q:
		print(q)


def noncoding(genefile):

	global p1
	c = []
	with open(genefile, 'r') as mut:
		for a in mut:
			a = a.rstrip().split('\t')
			p = mut.name
			p1 = p.split('.')
			c.append('%s %s' % (a[0], a[1]))
	
	for m in c:
		m = m.rstrip().split()
		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'Vista_enhancers_flankinGenes.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and int(m[1]) >= int(v[1]) and int(v[2]) >= int(m[1]):
					with open(p1[0]+"_regulatory-analysis.txt", 'a') as k:
						k.write('\t'.join((v[0], m[1], v[-2], v[-1], 'Enhancer'))+'\n')

		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'EPDnew_promoter.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and int(m[1]) >= int(v[1]) and int(v[2]) >= int(m[1]):
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join((v[0], m[1], v[3], 'Promoter'))+'\n')

		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'TSS_promoter.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and int(m[1]) == int(v[2]):
					# tss.append('%s %s %s %s %s' % (v[0], m[1], v[1], v[-1], 'TSS'))
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join((v[0], m[1], v[1], v[-1], 'TSS'))+'\n')

		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'cpgIsland.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and int(m[1]) >= int(v[1]) and int(v[2]) >= int(m[1]):
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join((v[0], m[1], v[-1], 'CpG'))+'\n')

		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'insulator-CTCF-binding-sites.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and int(m[1]) >= int(v[1]) and int(v[2]) >= int(m[1]):
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join((v[0], m[1], v[-2], v[-1], 'Insulator'))+'\n')


		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'miRNA.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and int(m[1]) == int(v[1]):
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join((v))+'\n')

		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'AltCodon.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and int(m[1]) == int(v[2]):
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join((v[0], v[2], v[1], 'AltCodon'))+'\n')

		with open(os.path.join(wd+'/'+'Regulatory datafiles'+'/'+'SPLICE_sires-snps.txt')) as vi:
			for v in vi:
				v = v.rstrip().split('\t')
				if m[0] == v[0] and int(m[1]) == int(v[1]):
					with open(p1[0]+"_regulatory-analysis.txt", 'a+') as k:
						k.write('\t'.join((v[0], m[1], v[-1], 'Splice-site'))+'\n')

def ncnetwork(functional):

	global p1

	import json

	import urllib.request

	from urllib.request import urlopen

	with open(functional, 'r') as nc:
		p = nc.name
		p1 = p.split('.')
		enrich = []
		for n in nc:
			n = n.rstrip().split('\t')
			if n[-1] == 'Enhancer':
				enrich.append(n[-3] +'\t'+n[-2])
			if n[-1] == 'Promoter' or n[-1] == 'TSS':
				take1 = n[-2].split('_')
				enrich.append(take1[0])
			if n[-1] == 'miRNA':
				enrich.append(n[-3] +'\t'+n[-2])
			if n[-1] == 'Alternative-Translation-site':
				enrich.append(n[-3] +'\t'+n[-2])
			if n[-1] == 'Splice-site':
				enrich.append(n[-3] +'\t'+n[-2])

		de = {}
		lines = (line.rstrip('\t') for line in enrich)
		unique = OrderedDict.fromkeys( (line for line in lines if line) )
		k = [(line.split())[0] for line in unique]
		for word in k:
			de[word] = de.get(word, 0) + 1

		word_freq = []
		for key, value in de.items():
			word_freq.append(list((value, key)))

		word_freq.sort(reverse=True)
		for i in word_freq:
				with open(p1[0]+"_Most-Altered-Mouse-Genes.txt", 'a') as fil:
					fil.write(i[1] +'\t'+ str(i[0]) +'\t'+ str(round(i[0] / len(k) * 100, 2))+'%'+'\n')
		
		y = [j[0] for j in word_freq]
		x = [z[1] for z in word_freq]

		fig, ax = plt.subplots(figsize=(6, 3.5))
		xt = tuple(x)
		my_range=list(range(1,len(y)+1))

		plt.plot(my_range, y, "o", markersize=2, color='#55342d', alpha=0.6)
		plt.xticks(my_range, xt , rotation=90)
		ax.tick_params(axis='both', which='major', labelsize=2)
		plt.ylabel ('Number of SNPs', size=7)
		plt.xlabel ('Protein', size=7)
		fig.savefig(p1[0]+'_Most-Altered-Mouse-Genes.png', format='png', dpi=900, bbox_inches='tight')

		sc = glob.glob(wd+'/'+"Functional-datafiles"+"/"+"numbers/*")
		for c in sc:
			f1 = c.split('/')
			f2 = f1[-1].split('.')
			k = [(line.split())[0] for line in unique]
			for di in k:
				if di == f2[0]:
					with open(c) as ij:
						for i in ij:
							if not 'protein-coding gene' in i and not 'Unknown' in i:
								with open(p1[0]+'_scExpression-data.txt', 'a') as k:
									k.write(di +'\t'+ str(i))
		print(p1[0]+" sc-Expression analysis is done")

		if p1[0]+'_scExpression-data.txt':
			plt.figure(figsize=(6, 4))
			#plt.yscale('symlog')
			palette = np.repeat(np.array(sns.color_palette("deep")),1, axis=0)
			plt.yticks(fontsize=3)
			plt.xticks(fontsize=4)
			scx = pd.read_csv(p1[0]+"_scExpression-data.txt", sep='\t', index_col=None)
			scx.columns = ["Genes", "Cell", "Number"]
			sns.barplot(x='Number', y='Cell', data=scx, hue='Genes', orient='h',palette=palette)
			plt.ylabel(" cell-type", size=7)
			plt.xlabel("number of clusters", fontsize=7)
			sns.set_color_codes("pastel")
			plt.legend(loc='upper left', fontsize=4, title='Genes', title_fontsize=4, bbox_to_anchor=(1.05, 1))
			plt.savefig(p1[0]+'_scExp.png', format='png', dpi=350, bbox_inches='tight')

		print(p1[0]+" sc-Expression analysis visualisation is done")

	
	phen = []
	with open(os.path.join(wd+'/'+"Functional-datafiles"+"/"+"Phenotype.txt")) as rr:  
		for r in rr:
			d =  r.rstrip().split("\t")
			for g in enrich:
				g = g.rstrip().split('\t')
				if d[0] == g[0]:
					phen.append('\t'.join(d))


	with open(p1[0]+"_phenotype_accessment.txt", 'a') as pi:
		for pe in OrderedDict.fromkeys(phen):
			pi.write(pe+'\n')

	ei = [(line.split())[0] for line in enrich]
	
	string_api_url = "https://string-db.org/api"
	output_format = "json"
	method_enrich = "enrichment"
	species = "10090"
	my_app = "www.aweseome_app.org"
	caller_identity = 'aarslan'
	request_url = string_api_url + "/" + output_format + "/" + method_enrich + "?"
	request_url += "identifiers=" + "%0d".join(ei[:200])
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
						category = row["category"]
						fig, ax = plt.subplots(figsize=(9, 5.5))
						
						if fdr <= 0.05 and category == 'Process':
							with open(p1[0]+'_BiologicalProcess.txt', 'a') as biop:
								biop.write("\t".join([term, str(fdr), description, preferred_names])+'\n')
						else:
							with open(p1[0]+'_enrichment.txt', 'a+') as en:
								en.write("No significant enrichment detected")

						if fdr <= 0.05 and category == 'RCTM' or category == 'KEGG':
							with open(p1[0]+'_Pathways.txt', 'a') as pathw:
								pathw.write("\t".join([term, str(fdr), description, preferred_names])+'\n')
						else:
							with open(p1[0]+'_enrichment.txt', 'a+') as en:
								en.write("No significant enrichment detected")

	else:
		with open(p1[0]+'_enrichment.txt', 'a+') as en:
			en.write("No significant enrichment detected")

	output_format_figure = "image"
	method = "network"
	request_urlf = string_api_url + "/" + output_format_figure + "/" + method + "?"
	request_urlf += "identifiers=%s"
	request_urlf += "&" + "species=" + species
	request_urlf += "&" + "add_white_nodes=0"
	request_urlf += "&" + "network_flavor=actions"
	request_urlf += "&" + "caller_identity=aarslan" + my_app	
	c = "%0d".join(ei[:200])
	urllib.request.urlretrieve(request_urlf % c, "%s.png" % "".join(p1[0]+"_Protein_Network"))

	try:
		if p1[0]+'_Pathways.txt':

				fig, ax = plt.subplots(figsize=(7.5, 5.5))
				x = []
				y = []
				with open(p1[0]+'_Pathways.txt') as ij:
					for i in ij.readlines()[:30]:
						i = i.rstrip().split('\t')
						x.append(i[2])
						y.append(abs(math.log10(float(i[1]))))

		xt = tuple(x)
		ys, xs = zip(*sorted(zip(y, xt)))
		my_range=list(range(1,len(y)+1))
		plt.hlines(y=my_range, xmin=ys, xmax=ys[0], color='#197213', alpha=0.2, linewidth=3)
		plt.plot(ys, my_range, "o", markersize=3, color='#197213', alpha=0.6)
		plt.yticks(my_range, xs)
		ax.set_xticks(ax.get_xticks()[::2])
		ax.tick_params(axis='both', which='major', labelsize=6)
		plt.xlabel ('-log10(p-value)', size=7)
		fig.savefig(p1[0]+'_Significantly-Altered-Biological-Pathways.png', format='png', dpi=900, bbox_inches='tight')
	
	except (StopIteration, FileNotFoundError) as q:
		print(q)
	
	try:
		if p1[0]+'_BiologicalProcess.txt':

			fig, ax = plt.subplots(figsize=(7.5, 5.5))
			x = []
			y = []
			with open(p1[0]+'_BiologicalProcess.txt','r') as ij:
				for i in ij.readlines()[:30]:
					i = i.rstrip().split('\t')
					x.append(i[2])
					y.append(abs(math.log10(float(i[1]))))

		xt = tuple(x)
		ys, xs = zip(*sorted(zip(y, xt)))
		my_range=list(range(1,len(y)+1))
		plt.hlines(y=my_range, xmin=ys, xmax=ys[0], color='#007acc', alpha=0.2, linewidth=3)
		plt.plot(ys, my_range, "o", markersize=3, color='#007acc', alpha=0.6)
		plt.yticks(my_range, xs)
		ax.set_xticks(ax.get_xticks()[::2])
		ax.tick_params(axis='both', which='major', labelsize=6)
		plt.xlabel ('-log10(p-value)', size=7)
		fig.savefig(p1[0]+'_Significantly-Altered_BiologicalProcess.png', format='png', dpi=900, bbox_inches='tight')
	except (StopIteration,FileNotFoundError) as q:
		print(q)


def testScorer(functional):

	global p1
	l = []
	with open(functional) as ij:
		p = ij.name
		p1 = p.split('.')
		for i in ij:
			i = i.rstrip().split('\t')
			l.append('\t'.join((i[-1], i[1])))

	de = {}
	lines = (line.rstrip('\t') for line in l)
	unique = OrderedDict.fromkeys( (line for line in lines if line) )
	k = [(line.split('\t'))[0] for line in unique]

	for word in k:
		de[word] = de.get(word, 0) + 1

	word_freq = []
	for key, value in de.items():
		word_freq.append(list((value, key)))

	pathw = []
	word_freq.sort(reverse=True)
	for i in word_freq:

		pv = [(	'Enhancer'	,	'28947'	),
		(	'Splice-site'	,	'486'	),
		(	'CpG'	,	'8938'	),
		(	'Insulator'	,	'12640'	),
		(	'Promoter'	,	'60136'	),
		(	'TSS'	,	'22'	),
		(	'AltCodon'	,	'15'	),
		(	'miRNA'	,	'425'	)]
	
		for v in pv:
			if v[0] == i[1]:
				# Number of SNPs
				N = len(k)
				# functional plus regulatory SNPs
				M = 340091 
				p = f"{Decimal(hypergeom(M, int(v[1]), N).pmf(i[0])):.2E}"
				pathw.append((v[0], p))

	pathw.sort(key = lambda x: float(x[1]), reverse = False)
	for pa in pathw:		
		with open(p1[0]+"_Enriched.txt",'a') as eg:
			eg.write('\t'.join(pa)+'\n')

	if p1[0]+'_Enriched.txt':

			fig, ax = plt.subplots(figsize=(7.5, 5.5))
			x = []
			y = []
			with open(p1[0]+'_Enriched.txt','r') as ij:
				for i in ij.readlines():
					i = i.rstrip().split('\t')
					x.append(i[0])
					try:
						y.append(abs(math.log10(float(i[1]))))
					except ValueError:
						pass

	xt = tuple(x)
	ys, xs = zip(*sorted(zip(y, xt)))
	my_range=list(range(1,len(y)+1))
	plt.hlines(y=my_range, xmin=ys, xmax=ys[0], color='indigo', alpha=0.2, linewidth=3)
	plt.plot(ys, my_range, "o", markersize=3, color='indigo', alpha=0.6)
	plt.yticks(my_range, xs)
	ax.tick_params(axis='both', which='major', labelsize=5)
	plt.xlabel ('-log10(p-value)', size=7)
	fig.savefig(p1[0]+'_Enriched-regulatory-regions.png', format='png', dpi=900, bbox_inches='tight')

	print(p1[0]+'_Enriched-regulatory-regions is ready')

def mainResultsR():

	try:
		pathway = p1[0]+"_Significantly-Altered-Biological-Pathways.png"
	except FileNotFoundError as e:
		print(e)
	try:
		process = p1[0]+"_Significantly-Altered_BiologicalProcess.png"
	except FileNotFoundError as e:
		print(e)
	try:
		genes = p1[0]+'_Most-Altered-Mouse-Genes.png'
	except FileNotFoundError as e:
		print(e)
	try:
		domain = p1[0]+'_Enriched-regulatory-regions.png'
	except FileNotFoundError as e:
		print(e)
	try:
		network = p1[0]+'_Protein_Network.png'
	except FileNotFoundError as e:
		print(e)
	try:
		singlecell = p1[0]+'_scExp.png'
	except FileNotFoundError as e:
		print(e)
	try:
		GO = p1[0]+"_GO-term-Enrichment.png"
	except FileNotFoundError as e:
		print(e)
	
	with open(p1[0]+'.html', 'a') as ht:
		
		html0 =  f"""<html> 
		<br> <body> <link rel="stylesheet" href="style.css">
		<body><img src="mmap.png" width="20%" alt="logo"  />
		<br>
		</br>
		<div class="wrapper"> <center>
				<h1> mMap results for {p1[0]} </div>
		<br>
		</br>
		<br><center>
		<div class="container">
	    	<input type="checkbox" id="zoomCheck1">
		    <label for="zoomCheck1">
		    <h2> Proteins with Alternative Alleles and Allele Conversation</h2>
		    <img src={genes}  onerror="WARNGING: no mutated genes found!" alt="WARNING: no mutated genes found!">
		    </label>
		    <br> *conservation scale = more nagative means more conserved site and vice-versa </br>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		    <input type="checkbox" id="zoomCheck2">
		    <label for="zoomCheck2">
		    <h2> Altered Protein Functional Regions Enrichment</h2>
		    <img src={domain}  onerror="WARNGING: no mutated domain found!" alt="WARNING: no mutated genes found!">
		    </label>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		    <input type="checkbox" id="zoomCheck3">
		    <label for="zoomCheck3">
		    <h2> Altered GO terms Enrichment</h2>
		    <img src={GO} onerror="WARNGING: no mutated GO term found!" alt="WARNING: no mutated genes found!" >
		    </label>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		    	<input type="checkbox" id="zoomCheck4">
		    	<label for="zoomCheck4">
		    	<h2> Altered Biological Pathways Enrichment</h2>
		    <br>
		   		<img src={pathway} onerror="WARNGING: no enrichment of pathways found!" alt="WARNING: no enrichment of pathways found!" >
		    </label>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		    	<input type="checkbox" id="zoomCheck5">
		    	<label for="zoomCheck5">
		    	<h2> Altered Biological processes Enrichment</h2>
		    	<img src={process}  onerror="WARNGING: no enrichment of pathways found!" alt="WARNING: no enrichment of pathways found!" >
		    </label>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		   		<input type="checkbox" id="zoomCheck6">
		    	<label for="zoomCheck6">
		    	<h2> Alternative Allele containing Protein(s) Interaction Network</h2>
		    	<img src={network} onerror="WARNGING: no protein network available for in the input data!" alt="WARNING: no mutated genes found!" >
		    </label>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		    	<input type="checkbox" id="zoomCheck7">
		   		<label for="zoomCheck7">
		    	<h2> Alternative Allele containing Protein(s) Single Cell RNA Expression</h2>
		    	<img src={singlecell} onerror="WARNGING: no sc data found!" alt="WARNING: no mutated genes found!" >
		    </label>
		    <br>
		    </br>
		    <br>
		    </br>
		    <div id="footer_container">
		    	<font color="white" text-align= "left">
		    		<div id="footer"> © 2020 Developed by Ahmed Arslan <aarslan@stanford.edu> @ Peltz Lab Stanford University School of Medicine
		    </div>
		</div>
		</body> </html>""" 

		ht.write(html0)
		print(p1[0]+" mMap results are ready")

def nocode_mutation(file1):
	try:
		noncoding(os.path.join(file1))
		print('Regulatory analysis is done...')
	except IOError:
		pass
	try:
		ncnetwork(os.path.join(p1[0]+"_regulatory-analysis.txt"))
		print('Network analysis is done...')
	except IOError:
		pass
	try:
		testScorer(os.path.join(p1[0]+".txt"))
		print('Enriched Regions analysis is done...')
	except IOError:
		pass
	try:
		mainResultsR()
		print('Final HTML page of mMap resutlts is ready ...')
	except IOError:
		pass
	
	try:
		mg = (time.time() - start_time)
		os.system("mkdir "+wd+"/"+p1[0]+"_mutations_"+ str(mg))
		try:
				shutil.move(wd+"/"+p1[0]+".txt", wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+"_Enriched.txt", wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+'_Enriched-regulatory-regions.png', wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+'_Most-Altered-Mouse-Genes.txt', wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+'_Most-Altered-Mouse-Genes.png', wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+"_phenotype_accessment.txt", wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+'_scExp.png', wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+'_scExpression-data.txt', wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+'_Significantly-Altered-Biological-Pathways.png', wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+'_Significantly-Altered_BiologicalProcess.png', wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+'_BiologicalProcess.txt', wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+'_enrichment.txt', wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+"_Protein_Network.png", wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
				shutil.move(wd+"/"+p1[0]+"_Pathways.txt", wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
				print(e)
		try:
			shutil.move(wd+"/"+p1[0]+".html", wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
			print(e)
		try:
			shutil.copy(wd+"/"+'style.css', wd+"/"+p1[0]+"_mutations_"+ str(mg))
			shutil.copy(wd+"/"+'mmap.png', wd+"/"+p1[0]+"_mutations_"+ str(mg))
		except FileNotFoundError as e:
			print(e)
	except IOError as e:
		print(e)


def mainResults():

	try:
		GO = p2[0]+"_GO-term-Enrichment.png"
	except FileNotFoundError as e:
		print(e)
	try:
		pathway = p2[0]+"_Significantly-Altered-Biological-Pathways.png"
	except FileNotFoundError as e:
		print(e)
	try:
		process = p2[0]+"_Significantly-Altered-Biological-Process.png"
	except FileNotFoundError as e:
		print(e)
	try:
		genes = p2[0]+'_Mutated-Proteins-at-Conserved-Regions.png'
	except FileNotFoundError as e:
		print(e)
	try:
		domain = p2[0]+"_Enriched-functional-regions.png"
	except FileNotFoundError as e:
		print(e)
	try:
		network = p2[0]+'_Protein_Network.png'
	except FileNotFoundError as e:
		print(e)
	try:
		singlecell = p2[0]+'_scExp.png'
	except FileNotFoundError as e:
		print(e)

	with open(p2[0]+'.html', 'a') as ht:

		html0 =  f"""<html> 
		<br> <body> <link rel="stylesheet" href="style.css">
		<body><img src="mmap.png" width="20%" alt="logo"  />
		<br>
		</br>
		<div class="wrapper"> <center>
				<h1> mMap results for {p2[0]} </div>
		<br>
		</br>
		<br><center>
		<div class="container">
	    	<input type="checkbox" id="zoomCheck1">
		    <label for="zoomCheck1">
		    <h2> Proteins with Alternative Alleles and Allele Conversation</h2>
		    <img src={genes}  onerror="WARNGING: no mutated genes found!" alt="WARNING: no mutated genes found!">
		    </label>
		    <br> *conservation scale = more nagative means more conserved site and vice-versa </br>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		    <input type="checkbox" id="zoomCheck2">
		    <label for="zoomCheck2">
		    <h2> Altered Protein Functional Regions Enrichment</h2>
		    <img src={domain}  onerror="WARNGING: no mutated domain found!" alt="WARNING: no mutated genes found!">
		    </label>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		    <input type="checkbox" id="zoomCheck3">
		    <label for="zoomCheck3">
		    <h2> Altered GO terms Enrichment</h2>
		    <img src={GO} onerror="WARNGING: no mutated GO term found!" alt="WARNING: no mutated genes found!" >
		    </label>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		    	<input type="checkbox" id="zoomCheck4">
		    	<label for="zoomCheck4">
		    	<h2> Altered Biological Pathways Enrichment</h2>
		    <br>
		   		<img src={pathway} onerror="WARNGING: no enrichment of pathways found!" alt="WARNING: no enrichment of pathways found!" >
		    </label>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		    	<input type="checkbox" id="zoomCheck5">
		    	<label for="zoomCheck5">
		    	<h2> Altered Biological processes Enrichment</h2>
		    	<img src={process}  onerror="WARNGING: no enrichment of pathways found!" alt="WARNING: no enrichment of pathways found!" >
		    </label>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		   		<input type="checkbox" id="zoomCheck6">
		    	<label for="zoomCheck6">
		    	<h2> Alternative Allele containing Protein(s) Interaction Network</h2>
		    	<img src={network} onerror="WARNGING: no protein network available for in the input data!" alt="WARNING: no mutated genes found!" >
		    </label>
		</div>
		<br>
		</br>
		<br>
		</br>
		<div><center>
			<div class="container">
		    	<input type="checkbox" id="zoomCheck7">
		   		<label for="zoomCheck7">
		    	<h2> Alternative Allele containing Protein(s) Single Cell RNA Expression</h2>
		    	<img src={singlecell} onerror="WARNGING: no sc data found!" alt="WARNING: no mutated genes found!" >
		    </label>
		    <br>
		    </br>
		    <br>
		    </br>
		    <div id="footer_container">
		    	<font color="white" text-align= "left">
		    		<div id="footer"> © 2020 Developed by Ahmed Arslan <aarslan@stanford.edu> @ Peltz Lab Stanford University School of Medicine
		    </div>
		</div>
		</body> </html>""" 

		ht.write(html0)
		print(p2[0]+" mMap results are ready")

def gene_mutation(genefile):

	try:
			gene_file(os.path.join(genefile))
			print(p2[0]+' conservation analysis of functional SNPs is finished')
	except Exception as e:
			print(e)
	try:
			testScore(os.path.join(p2[0]+"_functional-accessment.txt"))
			print('testScore analysis is finished')
	except Exception as e:
			print(e)
	try:
			network(os.path.join(p2[0]+"_functional-accessment.txt"))
			print('Network is generated')
	except Exception as e:
			print(e)
	try:
		downstream(os.path.join(p2[0]+"_functional-accessment.txt"))
		print('Disease-Gene NLP analysis is done...')
	except Exception as e:
		print(e)
	try:
			go()
			print('GO plot is generated')
	except Exception as e:
			print(e)
	try:
			mainResults()
			print('mMap analysis is done and Final HTML page of mMap resutlts is ready ...')
	except Exception as e:
			print(e)

	try:
		mg = (time.time() - start_time)
		os.system("mkdir "+wd+"/"+p2[0]+"_mutations_"+ str(mg))
		try:
			shutil.move(wd+"/"+p2[0]+"_functional-accessment.txt", wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:
			shutil.move(wd+"/"+p2[0]+"_Protein_Network.png", wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:
			shutil.move(wd+"/"+p2[0]+"_enrichment.txt", wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:
			shutil.move(wd+"/"+p2[0]+'_Mutated-Proteins-at-Conserved-Regions.png', wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:
			shutil.move(wd+"/"+p2[0]+'_Significantly-Altered-Biological-Pathways.png', wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:
			shutil.move(wd+"/"+p2[0]+'_Pathways.txt', wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:
			shutil.move(wd+"/"+p2[0]+'_BiologicalProcess.txt', wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:
			shutil.move(wd+"/"+p2[0]+'_Significantly-Altered-Biological-Process.png', wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:
			shutil.move(wd+"/"+ p2[0]+'_Enriched-functional-regions.txt', wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:
			shutil.move(wd+"/"+p2[0]+'_BiologicalProcess_Revigo.txt', wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:
			shutil.move(wd+"/"+ p2[0]+'_Enriched-functional-regions.png', wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass	
		try:
			shutil.move(wd+"/"+ p2[0]+'_GO-term-Enrichment.png', wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:	
			shutil.move(wd+"/"+p2[0]+"_phenotype_accessment.txt", wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:	
			shutil.move(wd+"/"+p2[0]+"_scExpression-data.txt", wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:	
			shutil.move(wd+"/"+p2[0]+'_scExp.png', wd+"/"+p2[0]+"_mutations_" + str(mg))
		except IOError:
			pass
		try:	
			shutil.move(wd+"/"+p2[0]+'.html', wd+"/"+p2[0]+"_mutations_" + str(mg))
			shutil.copy(wd+"/"+'style.css', wd+"/"+p2[0]+"_mutations_"+ str(mg))
			shutil.copy(wd+"/"+'mmap.png', wd+"/"+p2[0]+"_mutations_"+ str(mg))
		except IOError:
			pass
			print('All analyses for ' +p2[0]+' are completed!' +'\n'+ '<<< mMap curated and maintained by Ahmed Arslan <aarslan@stanford.edu> >>>')

	except IOError:
		pass
	return "Done"

if __name__=='__main__':
	
	import argparse
	from multiprocessing import Pool
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
		except (IOError, IndexError):
			pass
	elif args.noncoding:
		try:
			gene1 =  [f for f in glob.glob(os.path.join("*.txt")) if f.endswith(sys.argv[-1])]
			p.map(nocode_mutation, gene1)
			p.close()
			p.join()
		except (IOError, IndexError):
			pass
	else:
		print ("to run a mMap function seek help")
