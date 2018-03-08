#!/usr/bin/env python

import re
import os
import sys
from scipy import *
from numpy import *
import numpy as np
import scipy.stats as sp
from scipy.stats import norm
import matplotlib
import matplotlib.mlab as mlab
matplotlib.use('Agg')
import matplotlib.pyplot as plt

'''
This program requires several input files, including:
DATA
	expression_data (entrez genes and RPKM across tissues)
KEYS
	ec.ec (kegg ec pairs)
	ec.linked (ec pairs linked in fungi)
	entrez.ensembl (entrez genes and associated ensembl ids)
	entrez.genes (entrez genes and common names)
	entrez.pathway (entrez genes and associated pathways)
	entrez.ec (entrez genes and associated ec numbers)
	ec.compound (ec numbers and associated compound numbers)
	compound.names (compound numbers and associated names)
'''

# ---PART 1-----------------------------------------------------
# human-specific analysis using non-normalized data-------------

# parse Jason's EC_EC.crapml file to extract all EC pairs

ec_jason = open('ec.jason', 'r')
ec_ec = open('ec.ec', 'w')

for line in ec_jason:
    if 'parent' in line:
        ec1 = line[line.find('parent="') + len('parent="'):line.find('" ')]
        if line.find('M_0=" e') != -1:
            ecs = line[line.find('M_0=" e') + len('M_0=" '):line.find('" M_1')].split(' ')
        else:
            ecs = line[line.find('M_0="e') + len('M_0="'):line.find('" M_1')].split(' ')
        for i in ecs:
            if i:
                ec_ec.write('%s\t%s\n' %
                            (min(ec1[3::], i[3::]), max(ec1[3::], i[3::])))

ec_jason.close()
ec_ec.close()


# pull linked EC pairs sent from Kris out of Jason's file and call the rest unlinked

ec_ec = open('ec.ec', 'r')
ec_linked = open('ec.linked', 'r')
ec_unlinked = open('ec.unlinked', 'w')

ec_ec_dict = {}
for line in ec_ec:
    ec_ec_dict[line] = 1
ec_linked_dict = {}
for line in ec_linked:
    ec_linked_dict[line] = 1
for i in ec_ec_dict:
    if i in ec_linked_dict:
        continue
    else:
        ec_unlinked.write(i)

ec_ec.close()
ec_linked.close()
ec_unlinked.close()


# make linked and unlinked gene pair files

entrez_ec = open('entrez.ec', 'r')
ec_linked = open('ec.linked', 'r')
ec_unlinked = open('ec.unlinked', 'r')
gene_linked = open('gene.linked', 'w')
summary_table = open('summary.table', 'w')
gene_unlinked = open('gene.unlinked', 'w')
dict = {}
entrez_ec_dict = {}
for line in entrez_ec:
    cols = line.strip().split('\t')
    entrez_ec_dict.setdefault(cols[1], []).append(cols[0])
for line in ec_linked:
    cols = line.strip().split('\t')
    if cols[0] in entrez_ec_dict and cols[1] in entrez_ec_dict:
        genes1 = entrez_ec_dict[cols[0]]
        genes2 = entrez_ec_dict[cols[1]]
        for i in genes1:
            for j in genes2:
                ecpair = '(' + str(cols[0]) + ',' + str(cols[1]) + ')'
                genepair = '(' + str(i) + ',' + str(j) + ')'
                dict.setdefault(ecpair, []).append(genepair)
for x in dict:
    summary_table.write('%s\t%s\n' % (x, dict[x]))
for line in ec_unlinked:
    cols = line.strip().split('\t')
    if cols[0] in entrez_ec_dict and cols[1] in entrez_ec_dict:
        genes1 = entrez_ec_dict[cols[0]]
        genes2 = entrez_ec_dict[cols[1]]
        for i in genes1:
            for j in genes2:
                gene_unlinked.write('%s\t%s\n' % (min(i, j), max(i, j)))

entrez_ec.close()
ec_linked.close()
ec_unlinked.close()
gene_linked.close()
gene_unlinked.close()
summary_table.close()


# sum gene data to create EC data

entrez_ec = open('entrez.ec', 'r')
gene_expression_data = open('gene.expression_data', 'r')
ec_expression_data_all = open('ec.expression_data_all', 'w')

entrez_ec_dict = {}
for line in entrez_ec:
    cols = line.strip().split('\t')
    entrez_ec_dict.setdefault(cols[0], []).append(cols[1])
for line in gene_expression_data:
    cols = line.strip().split('\t')
    if cols[0] in entrez_ec_dict:
        ec_expression_data_all.write(line.replace(
            str(cols[0]), str(entrez_ec_dict[cols[0]])[2:-2], 1))

entrez_ec.close()
gene_expression_data.close()
ec_expression_data_all.close()

ec_expression_data_all = open('ec.expression_data_all', 'r')
ec_expression_data = open('ec.expression_data', 'w')

ecs = {}
for line in ec_expression_data_all:
    cols = line.strip().split('\t')
    if cols[0] in ecs:
        ecs[cols[0]] += np.array([float(e) for e in cols[1:]])
    else:
        ecs[cols[0]] = np.array([float(e) for e in cols[1:]])
for key in ecs:
    ec_expression_data.write(key + '\t' + str(ecs[key][0]) + '\t' + str(ecs[key][1]) + '\t' + str(ecs[key][2]) + '\t' + str(ecs[key][3]) + '\t' + str(ecs[key][4]) + '\t' + str(ecs[key][5]) + '\t' + str(ecs[key][6]) + '\t' + str(ecs[key][7]) + '\t' + str(
        ecs[key][8]) + '\t' + str(ecs[key][9]) + '\t' + str(ecs[key][10]) + '\t' + str(ecs[key][11]) + '\t' + str(ecs[key][12]) + '\t' + str(ecs[key][13]) + '\t' + str(ecs[key][14]) + '\t' + str(ecs[key][15]) + '\t' + str(ecs[key][16]) + '\t' + str(ecs[key][17]) + '\n')

ec_expression_data_all.close()
ec_expression_data.close()
os.remove('ec.expression_data_all')

# do Pearson's correlation for data sets (linked/unlinked genes and linked/unlinked ecs)

gene_expression_data = open('gene.expression_data', 'r')
ec_expression_data = open('ec.expression_data', 'r')
gene_linked_pairs = open('gene.linked', 'r')
gene_unlinked_pairs = open('gene.unlinked', 'r')
ec_linked_pairs = open('ec.linked', 'r')
ec_unlinked_pairs = open('ec.unlinked', 'r')
gene_linked_correlation_data = open('gene.linked.correlation_data', 'w')
gene_unlinked_correlation_data = open('gene.unlinked.correlation_data', 'w')
ec_linked_correlation_data = open('ec.linked.correlation_data', 'w')
ec_unlinked_correlation_data = open('ec.unlinked.correlation_data', 'w')

gene_data = {}
for line in gene_expression_data:
    cols = line.strip().split('\t')
    if '#' in line:
        continue
    else:
        gene_data[cols[0]] = np.array([float(e) for e in cols[1:]])

for line in gene_linked_pairs:
    cols = line.strip().split('\t')
    if cols[0] in gene_data and cols[1] in gene_data:
        r = sp.pearsonr(gene_data[cols[0]], gene_data[cols[1]])
        gene_linked_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

for line in gene_unlinked_pairs:
    cols = line.strip().split('\t')
    if cols[0] in gene_data and cols[1] in gene_data:
        r = sp.pearsonr(gene_data[cols[0]], gene_data[cols[1]])
        gene_unlinked_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

ec_data = {}
for line in ec_expression_data:
    cols = line.strip().split('\t')
    if '#' in line:
        continue
    else:
        ec_data[cols[0]] = np.array([float(e) for e in cols[1:]])

for line in ec_linked_pairs:
    cols = line.strip().split('\t')
    if cols[0] in ec_data and cols[1] in ec_data:
        r = sp.pearsonr(ec_data[cols[0]], ec_data[cols[1]])
        ec_linked_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

for line in ec_unlinked_pairs:
    cols = line.strip().split('\t')
    if cols[0] in ec_data and cols[1] in ec_data:
        r = sp.pearsonr(ec_data[cols[0]], ec_data[cols[1]])
        ec_unlinked_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

gene_expression_data.close()
ec_expression_data.close()
gene_linked_pairs.close()
gene_unlinked_pairs.close()
ec_linked_pairs.close()
ec_unlinked_pairs.close()
gene_linked_correlation_data.close()
gene_unlinked_correlation_data.close()
ec_linked_correlation_data.close()
ec_unlinked_correlation_data.close()

# plot correlation data in skyline format

a = np.loadtxt('ec.linked.correlation_data', usecols=[2])
b = np.loadtxt('ec.unlinked.correlation_data', usecols=[2])

bins = linspace(-1, 1, 22)
bins2 = np.arange(-1, 1.1, 0.1)

hista, bin_edges = np.histogram(a, bins=bins)
histb, bin_edges = np.histogram(b, bins=bins)

percentsa = []
for x in hista:
    percentsa.append((float(x) / sum(hista)) * 100)
line1 = plt.plot(bins2, percentsa, linestyle='-',
                 label='FLOs (145 pairs) r=0.54', color='black', drawstyle='steps')

percentsb = []
for x in histb:
    percentsb.append((float(x) / sum(histb)) * 100)
line2 = plt.plot(bins2, percentsb, linestyle=':',
                 label='Background (1,586 pairs) r=0.32', color='black', drawstyle='steps')

plt.ylabel('Frequency (%)')
plt.xlabel('Correlation')
plt.xticks([-1, -0.5, 0, 0.5, 1], ['-1', '-0.5', '0', '0.5', '1'])
plt.yticks([0, 5, 10, 15, 20, 25], ['0', '5', '10', '15', '20', '25'])
plt.xlim(-1.1, 1.1)
plt.ylim(-1, 26)
plt.legend(loc=2)
plt.savefig('ec.svg', format='svg')

a = np.loadtxt('gene.linked.correlation_data', usecols=[2])
b = np.loadtxt('gene.unlinked.correlation_data', usecols=[2])

bins = linspace(-1, 1, 22)
bins2 = np.arange(-1, 1.1, 0.1)

hista, bin_edges = np.histogram(a, bins=bins)
histb, bin_edges = np.histogram(b, bins=bins)

percentsa = []
for x in hista:
    percentsa.append((float(x) / sum(hista)) * 100)
line1 = plt.plot(bins2, percentsa, linestyle='-',
                 label='FLOs (1,813 pairs) r=0.32', color='black', drawstyle='steps')

percentsb = []
for x in histb:
    percentsb.append((float(x) / sum(histb)) * 100)
line2 = plt.plot(bins2, percentsb, linestyle=':',
                 label='Background (10,295 pairs) r=0.22', color='black', drawstyle='steps')

plt.ylabel('Frequency (%)')
plt.xlabel('Correlation')
plt.xticks([-1, -0.5, 0, 0.5, 1], ['-1', '-0.5', '0', '0.5', '1'])
plt.yticks([0, 5, 10, 15, 20, 25], ['0', '5', '10', '15', '20', '25'])
plt.xlim(-1.1, 1.1)
plt.ylim(-1, 26)
plt.legend(loc=2)
plt.savefig('gene.svg', format='svg')

# print mann whitney statistic

a = np.loadtxt('gene.linked.correlation_data', usecols=[2])
b = np.loadtxt('gene.unlinked.correlation_data', usecols=[2])
c = np.loadtxt('ec.linked.correlation_data', usecols=[2])
d = np.loadtxt('ec.unlinked.correlation_data', usecols=[2])

print 'gene:', sp.mannwhitneyu(a, b, use_continuity=True)
print 'ec:', sp.mannwhitneyu(c, d, use_continuity=True)

# ---PART 2-----------------------------------------------------
# ---toxic versus nontoxic pairs--------------------------------

# separate toxic and nontoxic ecs

ec_pairs = open('ec.pairs', 'r')
ec_compound = open('ec.compound', 'r')
compound_names = open('compound.names', 'r')
compound_toxicity = open('compound.toxicity', 'r')
ec_toxic = open('ec.toxic', 'w')
ec_nontoxic = open('ec.nontoxic', 'w')

ec_pairs = open('ec.pairs', 'r')
ec_compound = open('ec.compound', 'r')
compound_names = open('compound.names', 'r')
compound_toxicity = open('compound.toxicity', 'r')
ec_toxic = open('ec.toxic', 'w')
ec_nontoxic = open('ec.nontoxic', 'w')

ec_compound_dict = {}
for line in ec_compound:
    cols = line.strip().split('\t')
    ec_compound_dict[cols[0]] = cols[1]

compound_names_dict = {}
for line in compound_names:
    cols = line.strip().split('\t')
    compound_names_dict.setdefault(cols[0], []).append(cols[1])

compound_toxicity_dict = {}
for line in compound_toxicity:
    cols = line.strip().split('\t')
    compound_toxicity_dict[cols[0]] = cols[1]

all_toxicity_dict = {}
for i in compound_names_dict:
    for j in compound_toxicity_dict:
        if str(j) in str(compound_names_dict[i]):
            all_toxicity_dict.setdefault(i, []).append(
                compound_toxicity_dict[j])

max_toxicity_dict = {}
for i in all_toxicity_dict:
    if '2' in all_toxicity_dict[i]:
        max_toxicity_dict[i] = '2'
    elif '1' in all_toxicity_dict[i]:
        max_toxicity_dict[i] = '1'
    else:
        max_toxicity_dict[i] = '?'

for line in ec_pairs:
    cols = line.strip().split('\t')
    if cols[0] in ec_compound_dict and cols[1] in ec_compound_dict:
        if ec_compound_dict[cols[0]] in max_toxicity_dict and ec_compound_dict[cols[1]] in max_toxicity_dict:
            if max_toxicity_dict[ec_compound_dict[cols[0]]] == '1':
                ec_toxic.write(line)
            elif max_toxicity_dict[ec_compound_dict[cols[0]]] == '2':
                ec_toxic.write(line)
            else:
                ec_nontoxic.write(line)
        else:
            ec_nontoxic.write(line)
    else:
        ec_nontoxic.write(line)

ec_pairs.close()
ec_compound.close()
compound_names.close()
compound_toxicity.close()
ec_toxic.close()
ec_nontoxic.close()

# make toxic and nontoxic genes

entrez_ec = open('entrez.ec', 'r')
ec_toxic = open('ec.toxic', 'r')
ec_nontoxic = open('ec.nontoxic', 'r')
gene_toxic = open('gene.toxic', 'w')
gene_nontoxic = open('gene.nontoxic', 'w')

entrez_ec_dict = {}
for line in entrez_ec:
    cols = line.strip().split('\t')
    entrez_ec_dict.setdefault(cols[1], []).append(cols[0])
for line in ec_toxic:
    cols = line.strip().split('\t')
    if cols[0] in entrez_ec_dict and cols[1] in entrez_ec_dict:
        genes1 = entrez_ec_dict[cols[0]]
        genes2 = entrez_ec_dict[cols[1]]
        for i in genes1:
            for j in genes2:
                gene_toxic.write('%s\t%s\n' % (min(i, j), max(i, j)))
for line in ec_nontoxic:
    cols = line.strip().split('\t')
    if cols[0] in entrez_ec_dict and cols[1] in entrez_ec_dict:
        genes1 = entrez_ec_dict[cols[0]]
        genes2 = entrez_ec_dict[cols[1]]
        for i in genes1:
            for j in genes2:
                gene_nontoxic.write('%s\t%s\n' % (min(i, j), max(i, j)))

entrez_ec.close()
ec_toxic.close()
ec_nontoxic.close()
gene_toxic.close()
gene_nontoxic.close()

# do Pearson's correlation for data sets (toxic/nontoxic genes and
# toxic/nontoxic ecs)

gene_expression_data = open('gene.expression_data', 'r')
ec_expression_data = open('ec.expression_data', 'r')
gene_toxic_pairs = open('gene.toxic', 'r')
gene_nontoxic_pairs = open('gene.nontoxic', 'r')
gene_linked_pairs = open('gene.linked', 'r')
gene_unlinked_pairs = open('gene.unlinked', 'r')
ec_toxic_pairs = open('ec.toxic', 'r')
ec_nontoxic_pairs = open('ec.nontoxic', 'r')
ec_linked_pairs = open('ec.linked', 'r')
ec_unlinked_pairs = open('ec.unlinked', 'r')
gene_toxic_correlation_data = open('gene.toxic.correlation_data', 'w')
gene_nontoxic_correlation_data = open('gene.nontoxic.correlation_data', 'w')
gene_linked_toxic_correlation_data = open(
    'gene.linked.toxic.correlation_data', 'w')
gene_linked_nontoxic_correlation_data = open(
    'gene.linked.nontoxic.correlation_data', 'w')
gene_unlinked_toxic_correlation_data = open(
    'gene.unlinked.toxic.correlation_data', 'w')
gene_unlinked_nontoxic_correlation_data = open(
    'gene.unlinked.nontoxic.correlation_data', 'w')
ec_toxic_correlation_data = open('ec.toxic.correlation_data', 'w')
ec_nontoxic_correlation_data = open('ec.nontoxic.correlation_data', 'w')
ec_linked_toxic_correlation_data = open(
    'ec.linked.toxic.correlation_data', 'w')
ec_linked_nontoxic_correlation_data = open(
    'ec.linked.nontoxic.correlation_data', 'w')
ec_unlinked_toxic_correlation_data = open(
    'ec.unlinked.toxic.correlation_data', 'w')
ec_unlinked_nontoxic_correlation_data = open(
    'ec.unlinked.nontoxic.correlation_data', 'w')

gene_data = {}
for line in gene_expression_data:
    cols = line.strip().split('\t')
    if '#' in line:
        continue
    else:
        gene_data[cols[0]] = np.array([float(e) for e in cols[1:]])

gene_linked_pairs_dict = {}
for line in gene_linked_pairs:
    gene_linked_pairs_dict[line] = 1

gene_unlinked_pairs_dict = {}
for line in gene_unlinked_pairs:
    gene_unlinked_pairs_dict[line] = 1

for line in gene_toxic_pairs:
    cols = line.strip().split('\t')
    if cols[0] in gene_data and cols[1] in gene_data:
        r = sp.pearsonr(gene_data[cols[0]], gene_data[cols[1]])
        gene_toxic_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')
        if line in gene_linked_pairs_dict:
            gene_linked_toxic_correlation_data.write(
                str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')
        else:
            gene_unlinked_toxic_correlation_data.write(
                str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

for line in gene_nontoxic_pairs:
    cols = line.strip().split('\t')
    if cols[0] in gene_data and cols[1] in gene_data:
        r = sp.pearsonr(gene_data[cols[0]], gene_data[cols[1]])
        gene_nontoxic_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')
        if line in gene_linked_pairs_dict:
            gene_linked_nontoxic_correlation_data.write(
                str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')
        else:
            gene_unlinked_nontoxic_correlation_data.write(
                str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

ec_data = {}
for line in ec_expression_data:
    cols = line.strip().split('\t')
    if '#' in line:
        continue
    else:
        ec_data[cols[0]] = np.array([float(e) for e in cols[1:]])

ec_linked_pairs_dict = {}
for line in ec_linked_pairs:
    ec_linked_pairs_dict[line] = 1

ec_unlinked_pairs_dict = {}
for line in ec_unlinked_pairs:
    ec_unlinked_pairs_dict[line] = 1

for line in ec_toxic_pairs:
    cols = line.strip().split('\t')
    if cols[0] in ec_data and cols[1] in ec_data:
        r = sp.pearsonr(ec_data[cols[0]], ec_data[cols[1]])
        ec_toxic_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')
        if line in ec_linked_pairs_dict:
            ec_linked_toxic_correlation_data.write(
                str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')
        else:
            ec_unlinked_toxic_correlation_data.write(
                str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

for line in ec_nontoxic_pairs:
    cols = line.strip().split('\t')
    if cols[0] in ec_data and cols[1] in ec_data:
        r = sp.pearsonr(ec_data[cols[0]], ec_data[cols[1]])
        ec_nontoxic_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')
        if line in ec_linked_pairs_dict:
            ec_linked_nontoxic_correlation_data.write(
                str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')
        else:
            ec_unlinked_nontoxic_correlation_data.write(
                str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

gene_expression_data.close()
ec_expression_data.close()
gene_toxic_pairs.close()
gene_nontoxic_pairs.close()
gene_linked_pairs.close()
gene_unlinked_pairs.close()
ec_toxic_pairs.close()
ec_nontoxic_pairs.close()
ec_linked_pairs.close()
ec_unlinked_pairs.close()
gene_toxic_correlation_data.close()
gene_nontoxic_correlation_data.close()
gene_linked_toxic_correlation_data.close()
gene_unlinked_toxic_correlation_data.close()
gene_linked_nontoxic_correlation_data.close()
gene_unlinked_nontoxic_correlation_data.close()
ec_toxic_correlation_data.close()
ec_nontoxic_correlation_data.close()
ec_linked_toxic_correlation_data.close()
ec_unlinked_toxic_correlation_data.close()
ec_linked_nontoxic_correlation_data.close()
ec_unlinked_nontoxic_correlation_data.close()

# plot correlation data in skyline format

a = np.loadtxt('gene.linked.nontoxic.correlation_data', usecols=[2])
b = np.loadtxt('gene.unlinked.nontoxic.correlation_data', usecols=[2])

bins = linspace(-1, 1, 22)
bins2 = np.arange(-1, 1.1, 0.1)

hista, bin_edges = np.histogram(a, bins=bins)
histb, bin_edges = np.histogram(b, bins=bins)

percentsa = []
for x in hista:
    percentsa.append((float(x) / sum(hista)) * 100)
line1 = plt.plot(bins2, percentsa, linestyle='-',
                 label='linked unknown toxicity genes (1,179 pairs)', color='black', drawstyle='steps-mid')

percentsb = []
for x in histb:
    percentsb.append((float(x) / sum(histb)) * 100)
line2 = plt.plot(bins2, percentsb, linestyle=':',
                 label='unlinked unknown toxicity genes (8,244 pairs)', color='black', drawstyle='steps-mid')

plt.ylabel('Frequency (%)')
plt.xlabel('Correlation')
plt.xticks([-1, -0.5, 0, 0.5, 1], ['-1', '-0.5', '0', '0.5', '1'])
plt.yticks([0, 5, 10, 15, 20, 25, 30], [
           '0', '5', '10', '15', '20', '25', '30'])
plt.xlim(-1.1, 1.1)
plt.ylim(-1, 31)
plt.legend(loc=2)
plt.savefig('gene.nontoxic.svg', format='svg')

a = np.loadtxt('gene.linked.toxic.correlation_data', usecols=[2])
b = np.loadtxt('gene.unlinked.toxic.correlation_data', usecols=[2])

bins = linspace(-1, 1, 22)
bins2 = np.arange(-1, 1.1, 0.1)

hista, bin_edges = np.histogram(a, bins=bins)
histb, bin_edges = np.histogram(b, bins=bins)

percentsa = []
for x in hista:
    percentsa.append((float(x) / sum(hista)) * 100)
line1 = plt.plot(bins2, percentsa, linestyle='-',
                 label='linked toxic genes (634 pairs)', color='black', drawstyle='steps-mid')

percentsb = []
for x in histb:
    percentsb.append((float(x) / sum(histb)) * 100)
line2 = plt.plot(bins2, percentsb, linestyle=':',
                 label='unlinked toxic genes (2,051 pairs)', color='black', drawstyle='steps-mid')

plt.ylabel('Frequency (%)')
plt.xlabel('Correlation')
plt.xticks([-1, -0.5, 0, 0.5, 1], ['-1', '-0.5', '0', '0.5', '1'])
plt.yticks([0, 5, 10, 15, 20, 25, 30], [
           '0', '5', '10', '15', '20', '25', '30'])
plt.xlim(-1.1, 1.1)
plt.ylim(-1, 31)
plt.legend(loc=2)
plt.savefig('gene.toxic.svg', format='svg')

a = np.loadtxt('ec.linked.nontoxic.correlation_data', usecols=[2])
b = np.loadtxt('ec.unlinked.nontoxic.correlation_data', usecols=[2])

bins = linspace(-1, 1, 22)
bins2 = np.arange(-1, 1.1, 0.1)

hista, bin_edges = np.histogram(a, bins=bins)
histb, bin_edges = np.histogram(b, bins=bins)

percentsa = []
for x in hista:
    percentsa.append((float(x) / sum(hista)) * 100)
line1 = plt.plot(bins2, percentsa, linestyle='-',
                 label='linked unknown toxicity ECs (123 pairs)', color='black', drawstyle='steps-mid')

percentsb = []
for x in histb:
    percentsb.append((float(x) / sum(histb)) * 100)
line2 = plt.plot(bins2, percentsb, linestyle=':',
                 label='unlinked unknown toxicity ECs (1,358 pairs)', color='black', drawstyle='steps-mid')

plt.ylabel('Frequency (%)')
plt.xlabel('Correlation')
plt.xticks([-1, -0.5, 0, 0.5, 1], ['-1', '-0.5', '0', '0.5', '1'])
plt.yticks([0, 5, 10, 15, 20, 25, 30], [
           '0', '5', '10', '15', '20', '25', '30'])
plt.xlim(-1.1, 1.1)
plt.ylim(-1, 31)
plt.legend(loc=2)
plt.savefig('ec.nontoxic.svg', format='svg')

a = np.loadtxt('ec.linked.toxic.correlation_data', usecols=[2])
b = np.loadtxt('ec.unlinked.toxic.correlation_data', usecols=[2])

bins = linspace(-1, 1, 22)
bins2 = np.arange(-1, 1.1, 0.1)

hista, bin_edges = np.histogram(a, bins=bins)
histb, bin_edges = np.histogram(b, bins=bins)

percentsa = []
for x in hista:
    percentsa.append((float(x) / sum(hista)) * 100)
line1 = plt.plot(bins2, percentsa, linestyle='-',
                 label='linked toxic ECs (22 pairs)', color='black', drawstyle='steps-mid')

percentsb = []
for x in histb:
    percentsb.append((float(x) / sum(histb)) * 100)
line2 = plt.plot(bins2, percentsb, linestyle=':',
                 label='unlinked toxic ECs (228 pairs)', color='black', drawstyle='steps-mid')

plt.ylabel('Frequency (%)')
plt.xlabel('Correlation')
plt.xticks([-1, -0.5, 0, 0.5, 1], ['-1', '-0.5', '0', '0.5', '1'])
plt.yticks([0, 5, 10, 15, 20, 25, 30], [
           '0', '5', '10', '15', '20', '25', '30'])
plt.xlim(-1.1, 1.1)
plt.ylim(-1, 31)
plt.legend(loc=2)
plt.savefig('ec.toxic.svg', format='svg')


# print mann whitney statistic

a = np.loadtxt('gene.linked.toxic.correlation_data', usecols=[2])
b = np.loadtxt('gene.linked.nontoxic.correlation_data', usecols=[2])
c = np.loadtxt('ec.linked.toxic.correlation_data', usecols=[2])
d = np.loadtxt('ec.linked.nontoxic.correlation_data', usecols=[2])

print 'gene:', sp.mannwhitneyu(a, b, use_continuity=True)
print 'ec:', sp.mannwhitneyu(c, d, use_continuity=True)

# plot where genes in pairs are most highly expressed

b = 0
c = 0
h = 0
k = 0
l = 0
t = 0
bb = 0
cc = 0
hh = 0
kk = 0
ll = 0
tt = 0

tissues = {1: 'brain', 2: 'brain', 3: 'brain', 4: 'brain', 5: 'brain', 6: 'brain', 7: 'cerebellum', 8: 'cerebellum', 9: 'heart',
           10: 'heart', 11: 'heart', 12: 'kidney', 13: 'kidney', 14: 'kidney', 15: 'liver', 16: 'liver', 17: 'testes', 18: 'testes'}
max_tissue = {}
gene_tissue_keys = open('gene.tissues', 'r')
gene_correlations = open('gene.correlations', 'r')

for line in gene_tissue_keys:
    cols = line.strip().split('\t')
    max_tissue[cols[0]] = cols[1]

for line in gene_correlations:
    cols = line.strip().split('\t')
    if max_tissue[cols[0]] == '1' or max_tissue[cols[0]] == '2' or max_tissue[cols[0]] == '3' or max_tissue[cols[0]] == '4' or max_tissue[cols[0]] == '5' or max_tissue[cols[0]] == '6':
        bb += 1
    elif max_tissue[cols[0]] == '7' or max_tissue[cols[0]] == '8':
        cc += 1
    elif max_tissue[cols[0]] == '9' or max_tissue[cols[0]] == '10' or max_tissue[cols[0]] == '11':
        hh += 1
    elif max_tissue[cols[0]] == '12' or max_tissue[cols[0]] == '13' or max_tissue[cols[0]] == '14':
        kk += 1
    elif max_tissue[cols[0]] == '15' or max_tissue[cols[0]] == '16':
        ll += 1
    elif max_tissue[cols[0]] == '17' or max_tissue[cols[0]] == '18':
        tt += 1
    if float(cols[1]) >= 0.7:
        if max_tissue[cols[0]] == '1' or max_tissue[cols[0]] == '2' or max_tissue[cols[0]] == '3' or max_tissue[cols[0]] == '4' or max_tissue[cols[0]] == '5' or max_tissue[cols[0]] == '6':
            b += 1
        elif max_tissue[cols[0]] == '7' or max_tissue[cols[0]] == '8':
            c += 1
        elif max_tissue[cols[0]] == '9' or max_tissue[cols[0]] == '10' or max_tissue[cols[0]] == '11':
            h += 1
        elif max_tissue[cols[0]] == '12' or max_tissue[cols[0]] == '13' or max_tissue[cols[0]] == '14':
            k += 1
        elif max_tissue[cols[0]] == '15' or max_tissue[cols[0]] == '16':
            l += 1
        else:
            t += 1

total = b + c + h + k + l + t
total2 = bb + cc + hh + kk + ll + tt

fig = [float(h) / float(total) * 100, float(b) / float(total) * 100, float(k) / float(total) * 100,
       float(c) / float(total) * 100, float(t) / float(total) * 100, float(l) / float(total) * 100]
fig2 = [float(hh) / float(total2) * 100, float(bb) / float(total2) * 100, float(kk) / float(total2) *
        100, float(cc) / float(total2) * 100, float(tt) / float(total2) * 100, float(ll) / float(total2) * 100]

width = 0.4
N = 6
ind = np.arange(N)

rects1 = plt.bar(ind, fig2, color='black', width=width, edgecolor='none')
rects2 = plt.bar(ind + width, fig, color='gray', width=width, edgecolor='none')
plt.legend((rects1[0], rects2[0]), ('All', 'Correlation >= 0.7'), loc=2)
plt.ylabel('Frequency (%)')
plt.xlabel('Tissue')
plt.xlim(-0.2, 6)
plt.ylim(0, 70)
plt.xticks(ind + width, ['heart', 'brain', 'kidney',
                         'cerebellum', 'testes', 'liver'])
plt.savefig('tissue.expression.svg', format='svg')
gene_tissue_keys.close()
gene_correlations.close()

# find gene pairs corresponding to highly correlated ec pairs >= 0.7

ec_correlated = open('ec.correlated', 'r')
entrez_ec = open('entrez.ec', 'r')
gene_correlated = open('gene.correlated', 'w')
gene_linked_pairs = open('gene.linked.correlation_data', 'r')
gene_unlinked_pairs = open('gene.unlinked.correlation_data', 'r')

entrez_ec_dict = {}
for line in entrez_ec:
    cols = line.strip().split('\t')
    entrez_ec_dict.setdefault(cols[1], []).append(cols[0])
linked = {}
for line in gene_linked_pairs:
    cols = line.strip().split('\t')
    linked[cols[0] + '\t' + cols[1]] = cols[2]
unlinked = {}
for line in gene_unlinked_pairs:
    cols = line.strip().split('\t')
    unlinked[cols[0] + '\t' + cols[1]] = cols[2]
for line in ec_correlated:
    line.strip()
    cols = line.strip().split('\t')
    if cols[0] in entrez_ec_dict and cols[1] in entrez_ec_dict:
        mylist0 = entrez_ec_dict[cols[0]]
        mylist1 = entrez_ec_dict[cols[1]]
        for i in mylist0:
            for j in mylist1:
                string = min(i, j) + '\t' + max(i, j)
                if string in linked:
                    gene_correlated.write('%s\t%s\n' %
                                          (string, linked[string]))
                elif string in unlinked:
                    gene_correlated.write('%s\t%s\n' %
                                          (string, unlinked[string]))


ec_correlated.close()
entrez_ec.close()
gene_correlated.close()
gene_linked_pairs.close()
gene_unlinked_pairs.close()

# plot gene pairs corresponding to highly correlated ec pairs

a = np.loadtxt('gene.correlated', usecols=[2])
b = np.loadtxt('ec.correlated', usecols=[2])

bins = linspace(-1, 1, 22)
bins2 = np.arange(-1, 1.1, 0.1)

hista, bin_edges = np.histogram(a, bins=bins)
histb, bin_edges = np.histogram(b, bins=bins)

percentsa = []
for x in hista:
    percentsa.append((float(x) / sum(hista)) * 100)
line1 = plt.plot(bins2, percentsa, linestyle='-',
                 label='genes (1,114 pairs)', color='black', drawstyle='steps-mid')

percentsb = []
for x in histb:
    percentsb.append((float(x) / sum(histb)) * 100)
line2 = plt.plot(bins2, percentsb, linestyle=':',
                 label='ecs (67 pairs)', color='black', drawstyle='steps-mid')

plt.ylabel('Frequency (%)')
plt.xlabel('Correlation')
plt.xticks([-1, -0.5, 0, 0.5, 1], ['-1', '-0.5', '0', '0.5', '1'])
plt.xlim(-1.1, 1.1)
plt.ylim(-1, 41)
plt.legend(loc=2)
plt.savefig('correlated.svg', format='svg')

# ---PART 3-----------------------------------------------------
# interspecies analysis using normalized data-------------------

# sum gene data to create EC data

entrez_ec = open('entrez.ec', 'r')
gene_expression_data = open('gene.expression_data', 'r')
ec_expression_data_all = open('ec.expression_data_all', 'w')

entrez_ec_dict = {}
for line in entrez_ec:
    cols = line.strip().split('\t')
    entrez_ec_dict.setdefault(cols[0], []).append(cols[1])
for line in gene_expression_data:
    cols = line.strip().split('\t')
    if cols[0] in entrez_ec_dict:
        ec_expression_data_all.write(line.replace(
            str(cols[0]), str(entrez_ec_dict[cols[0]])[2:-2], 1))

entrez_ec.close()
gene_expression_data.close()
ec_expression_data_all.close()

ec_expression_data_all = open('ec.expression_data_all', 'r')
ec_expression_data = open('ec.expression_data', 'w')

ecs = {}
for line in ec_expression_data_all:
    cols = line.strip().split('\t')
    if cols[0] in ecs:
        ecs[cols[0]] += np.array([float(e) for e in cols[1:]])
    else:
        ecs[cols[0]] = np.array([float(e) for e in cols[1:]])
for key in ecs:
    ec_expression_data.write(key + '\t' + str(ecs[key][0]) + '\t' + str(ecs[key][1]) + '\t' + str(ecs[key][2]) + '\t' + str(ecs[key][3]) + '\t' + str(
        ecs[key][4]) + '\t' + str(ecs[key][5]) + '\t' + str(ecs[key][6]) + '\t' + str(ecs[key][7]) + '\t' + str(ecs[key][8]) + '\n')

ec_expression_data_all.close()
ec_expression_data.close()
os.remove('ec.expression_data_all')

# do Pearson's correlation for data sets (linked/unlinked genes and
# linked/unlinked ecs)

gene_expression_data = open('gene.expression_data', 'r')
ec_expression_data = open('ec.expression_data', 'r')
gene_linked_pairs = open('gene.linked', 'r')
gene_unlinked_pairs = open('gene.unlinked', 'r')
ec_linked_pairs = open('ec.linked', 'r')
ec_unlinked_pairs = open('ec.unlinked', 'r')
gene_linked_correlation_data = open('gene.linked.correlation_data', 'w')
gene_unlinked_correlation_data = open('gene.unlinked.correlation_data', 'w')
ec_linked_correlation_data = open('ec.linked.correlation_data', 'w')
ec_unlinked_correlation_data = open('ec.unlinked.correlation_data', 'w')

gene_data = {}
for line in gene_expression_data:
    cols = line.strip().split('\t')
    if '#' in line:
        continue
    else:
        gene_data[cols[0]] = np.array([float(e) for e in cols[1:]])

for line in gene_linked_pairs:
    cols = line.strip().split('\t')
    if cols[0] in gene_data and cols[1] in gene_data:
        r = sp.pearsonr(gene_data[cols[0]], gene_data[cols[1]])
        gene_linked_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

for line in gene_unlinked_pairs:
    cols = line.strip().split('\t')
    if cols[0] in gene_data and cols[1] in gene_data:
        r = sp.pearsonr(gene_data[cols[0]], gene_data[cols[1]])
        gene_unlinked_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

ec_data = {}
for line in ec_expression_data:
    cols = line.strip().split('\t')
    if '#' in line:
        continue
    else:
        ec_data[cols[0]] = np.array([float(e) for e in cols[1:]])

for line in ec_linked_pairs:
    cols = line.strip().split('\t')
    if cols[0] in ec_data and cols[1] in ec_data:
        r = sp.pearsonr(ec_data[cols[0]], ec_data[cols[1]])
        ec_linked_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

for line in ec_unlinked_pairs:
    cols = line.strip().split('\t')
    if cols[0] in ec_data and cols[1] in ec_data:
        r = sp.pearsonr(ec_data[cols[0]], ec_data[cols[1]])
        ec_unlinked_correlation_data.write(
            str(cols[0]) + '\t' + str(cols[1]) + '\t' + str(r[0]) + '\n')

gene_expression_data.close()
ec_expression_data.close()
gene_linked_pairs.close()
gene_unlinked_pairs.close()
ec_linked_pairs.close()
ec_unlinked_pairs.close()
gene_linked_correlation_data.close()
gene_unlinked_correlation_data.close()
ec_linked_correlation_data.close()
ec_unlinked_correlation_data.close()

# plot correlation data in skyline format

a = np.loadtxt('ec.linked.correlation_data', usecols=[2])
b = np.loadtxt('ec.unlinked.correlation_data', usecols=[2])

bins = linspace(-1, 1, 22)
bins2 = np.arange(-1, 1.1, 0.1)

hista, bin_edges = np.histogram(a, bins=bins)
histb, bin_edges = np.histogram(b, bins=bins)

percentsa = []
for x in hista:
    percentsa.append((float(x) / sum(hista)) * 100)
line1 = plt.plot(bins2, percentsa, linestyle='-',
                 label='ptr linked ECs (107 pairs)', color='black', drawstyle='steps-mid')

percentsb = []
for x in histb:
    percentsb.append((float(x) / sum(histb)) * 100)
line2 = plt.plot(bins2, percentsb, linestyle=':',
                 label='ptr unlinked ECs (1,123 pairs)', color='black', drawstyle='steps-mid')

plt.ylabel('Frequency (%)')
plt.xlabel('Correlation')
plt.xticks([-1, -0.5, 0, 0.5, 1], ['-1', '-0.5', '0', '0.5', '1'])
plt.yticks([0, 5, 10, 15, 20, 25], ['0', '5', '10', '15', '20', '25'])
plt.xlim(-1.1, 1.1)
plt.ylim(-1, 26)
plt.legend(loc=2)
plt.savefig('ec.svg', format='svg')

# plot all primate data

a = np.loadtxt('hsa.gene.linked', usecols=[2])
b = np.loadtxt('hsa.gene.unlinked', usecols=[2])
c = np.loadtxt('ptr.gene.linked', usecols=[2])
d = np.loadtxt('ptr.gene.unlinked', usecols=[2])
e = np.loadtxt('ggo.gene.linked', usecols=[2])
f = np.loadtxt('ggo.gene.unlinked', usecols=[2])
g = np.loadtxt('ppy.gene.linked', usecols=[2])
h = np.loadtxt('ppy.gene.unlinked', usecols=[2])
i = np.loadtxt('mml.gene.linked', usecols=[2])
j = np.loadtxt('mml.gene.unlinked', usecols=[2])

bins = linspace(-1, 1, 21)
bins2 = np.arange(-0.95, 1.05, 0.1)

hista, bin_edges = np.histogram(a, bins=bins)
histb, bin_edges = np.histogram(b, bins=bins)
histc, bin_edges = np.histogram(c, bins=bins)
histd, bin_edges = np.histogram(d, bins=bins)
histe, bin_edges = np.histogram(e, bins=bins)
histf, bin_edges = np.histogram(f, bins=bins)
histg, bin_edges = np.histogram(g, bins=bins)
histh, bin_edges = np.histogram(h, bins=bins)
histi, bin_edges = np.histogram(i, bins=bins)
histj, bin_edges = np.histogram(j, bins=bins)

percentsa = []
for x in hista:
    percentsa.append((float(x) / sum(hista)) * 100)
line1 = plt.plot(bins2, percentsa, linestyle='-', label='hsa',
                 color='red', drawstyle='steps-mid')

percentsb = []
for x in histb:
    percentsb.append((float(x) / sum(histb)) * 100)
line2 = plt.plot(bins2, percentsb, linestyle=':',
                 color='red', drawstyle='steps-mid')

percentsc = []
for x in histc:
    percentsc.append((float(x) / sum(histc)) * 100)
line3 = plt.plot(bins2, percentsc, linestyle='-', label='ptr',
                 color='green', drawstyle='steps-mid')

percentsd = []
for x in histd:
    percentsd.append((float(x) / sum(histd)) * 100)
line4 = plt.plot(bins2, percentsd, linestyle=':',
                 color='green', drawstyle='steps-mid')

percentse = []
for x in histe:
    percentse.append((float(x) / sum(histe)) * 100)
line5 = plt.plot(bins2, percentse, linestyle='-', label='ggo',
                 color='blue', drawstyle='steps-mid')

percentsf = []
for x in histf:
    percentsf.append((float(x) / sum(histf)) * 100)
line6 = plt.plot(bins2, percentsf, linestyle=':',
                 color='blue', drawstyle='steps-mid')

percentsg = []
for x in histg:
    percentsg.append((float(x) / sum(histg)) * 100)
line7 = plt.plot(bins2, percentsg, linestyle='-', label='ppy',
                 color='orange', drawstyle='steps-mid')

percentsh = []
for x in histh:
    percentsh.append((float(x) / sum(histh)) * 100)
line8 = plt.plot(bins2, percentsh, linestyle=':',
                 color='orange', drawstyle='steps-mid')

percentsi = []
for x in histi:
    percentsi.append((float(x) / sum(histi)) * 100)
line9 = plt.plot(bins2, percentsi, linestyle='-', label='mml',
                 color='purple', drawstyle='steps-mid')

percentsj = []
for x in histj:
    percentsj.append((float(x) / sum(histj)) * 100)
line10 = plt.plot(bins2, percentsj, linestyle=':',
                  color='purple', drawstyle='steps-mid')

plt.ylabel('Frequency (%)')
plt.xlabel('Correlation')
plt.xticks([-1, -0.5, 0, 0.5, 1], ['-1', '-0.5', '0', '0.5', '1'])
plt.yticks([0, 5, 10, 15, 20], ['0', '5', '10', '15', '20'])
plt.xlim(-1.06, 1.06)
plt.ylim(-1, 21)
plt.legend(loc=2)
plt.savefig('species.svg', format='svg')

# flux analysis

linked = open('ec.linked', 'r')
unlinked = open('ec.unlinked', 'r')
data = open('flux.tab', 'r')
data2 = open('Waern.correlation_data', 'r')
keys = open('sce.ec', 'r')
outfile = open('ec.linked.data', 'w')
outfile2 = open('ec.unlinked.data', 'w')
outfile3 = open('ec.linked.flux', 'w')
outfile4 = open('ec.unlinked.flux', 'w')

dict = {}
linkedflux = {}
unlinkedflux = {}
flux = []
dictkeys = {}
dictkeys2 = {}

for line in data:
    cols = line.strip().split('\t')
    dict[cols[0]] = abs(float(cols[1]))
n = 0
for line in linked:
    cols = line.strip().split('\t')
    ec1 = cols[0]
    ec2 = cols[1]
    if ec1 in dict and ec2 in dict:
        flux = min(dict[ec1], dict[ec2])
        pair = str(ec1) + ', ' + str(ec2)
        n += float(flux)
        linkedflux[pair] = flux
        outfile3.write(pair + '\t' + str(flux) + '\n')
m = 0
for line in unlinked:
    cols = line.strip().split('\t')
    ec1 = cols[0]
    ec2 = cols[1]
    if ec1 in dict and ec2 in dict:
        flux = min(dict[ec1], dict[ec2])
        pair = str(ec1) + ', ' + str(ec2)
        m += float(flux)
        unlinkedflux[pair] = flux
        outfile4.write(pair + '\t' + str(flux) + '\n')


print 'linked flux = ' + str(n / 94)
print 'unlinked flux = ' + str(m / 826)

for line in keys:
    cols = line.strip().split('\t')
    dictkeys[cols[1]] = cols[0]

n = 0
genes = data2.readline().split('\t')
for line in data2:
    n += 1
    vals = line.strip('\n').split('\t')
    gene = vals[0]
    for i in range(n + 1, len(vals)):
        if gene in dictkeys and genes[i] in dictkeys:
            pair = str(dictkeys[gene]) + ', ' + str(dictkeys[genes[i]])
            if pair in linkedflux:
                outfile.write(dictkeys[gene] + '\t' + dictkeys[genes[i]] +
                              '\t' + vals[i] + '\t' + str(linkedflux[pair]) + '\n')
            elif pair in unlinkedflux:
                outfile2.write(dictkeys[gene] + '\t' + dictkeys[genes[i]] +
                               '\t' + vals[i] + '\t' + str(unlinkedflux[pair]) + '\n')


linked.close()
unlinked.close()
data.close()
data2.close()
keys.close()
outfile.close()
outfile2.close()
outfile3.close()
outfile4.close()

# network component analysis

keys = open('ec.network.component', 'r')
linked = open('ec.linked.correlation_data', 'r')
unlinked = open('ec.unlinked.correlation_data', 'r')
outfile = open('ec.linked.network', 'w')
outfile2 = open('ec.unlinked.network', 'w')

dict = {}
for line in keys:
    cols = line.strip().split('\t')
    dict[cols[0]] = cols[1]

for line in linked:
    line.strip()
    cols = line.strip().split('\t')
    if cols[0] in dict and cols[1] in dict:
        if dict[cols[0]] == 'IN' and dict[cols[1]] == 'IN':
            outfile.write(line.strip() + '\t' + 'IN' + '\n')
        if dict[cols[0]] == 'OUT' and dict[cols[1]] == 'OUT':
            outfile.write(line.strip() + '\t' + 'OUT' + '\n')
        if dict[cols[0]] == 'GSC' and dict[cols[1]] == 'GSC':
            outfile.write(line.strip() + '\t' + 'GSC' + '\n')
        else:
            outfile.write(line.strip() + '\t' + 'MIXED' + '\n')

for line in unlinked:
    line.strip()
    cols = line.strip().split('\t')
    if cols[0] in dict and cols[1] in dict:
        if dict[cols[0]] == 'IN' and dict[cols[1]] == 'IN':
            outfile2.write(line.strip() + '\t' + 'IN' + '\n')
        if dict[cols[0]] == 'OUT' and dict[cols[1]] == 'OUT':
            outfile2.write(line.strip() + '\t' + 'OUT' + '\n')
        if dict[cols[0]] == 'GSC' and dict[cols[1]] == 'GSC':
            outfile2.write(line.strip() + '\t' + 'GSC' + '\n')
        else:
            outfile2.write(line.strip() + '\t' + 'MIXED' + '\n')

keys.close()
linked.close()
unlinked.close()
outfile.close()
outfile2.close()
