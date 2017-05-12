# -*- coding: utf-8 -*-
"""
Created on Wed May  3 16:51:44 2017

@author: Keith Tilley
"""

import re
from Bio import SeqIO

def get_PredictBias_tag(first_gene, feat_list):
    for feat in feat_list:
        if "old_locus_tag" in feat.qualifiers:
            if first_gene in feat.qualifiers["old_locus_tag"]:
                if len(feat.qualifiers["old_locus_tag"]) == 1:
                    return ("old_locus_tag",0)
                else:
                    for (i,gene) in enumerate(feat.qualifiers["old_locus_tag"]):
                        if first_gene == gene:
                            return ("old_locus_tag",i)
        if "locus_tag" in feat.qualifiers:
            if first_gene in feat.qualifiers["locus_tag"]:
                if len(feat.qualifiers["locus_tag"]) == 1:
                    return ("locus_tag",0)
                else:
                    for (i,gene) in enumerate(feat.qualifiers["locus_tag"]):
                        if first_gene == gene:
                            return ("locus_tag",i)
    
def parse_PredictBias(genbankfile, PBresultfile, outputfile):
    record = next(SeqIO.parse(genbankfile,"genbank"))
    feat_list = record.features[1:]
    feat_list = [feat for feat in feat_list if feat.type == "gene"]    
    output = ""
    o = open(outputfile, "w")
    result = open(PBresultfile,"r").read()
    if len(result):
        # gi_genes will be a list of lists
        # Each entry a list of length two, containing the start and end gene names.
        gi_genes = re.findall("([A-Za-z0-9]+_?[A-Za-z0-9]+)\r?\t([A-Za-z0-9]+_?[A-Za-z0-9]+)\r?\t[YN]\t[YN]\t[YN]\t.*",result)
        if len(gi_genes) == 0:
            print(PBresultfile + ": no regex matches")
        gi_genes = sorted(gi_genes)
        # First find out if the results use old_locus_tag or locus_tag and the index
        get_tag = None
        i = 0
        while get_tag == None:
            get_tag = get_PredictBias_tag(gi_genes[i][0], feat_list)
            i += 1
        (tag, tag_i) = get_tag
        # create a dictionary of each gene with its start and stop to simplify indexing
        gbk_genes = dict()
        for feat in feat_list:
            if tag in feat.qualifiers and len(feat.qualifiers[tag]) > tag_i:
                gbk_genes[feat.qualifiers[tag][tag_i]] = (int(feat.location.start+1), int(feat.location.end))
        gis = []    
        for gene in gi_genes:
            # some gene results from PB are not present in the submitted genbank file
            if gene[0] in gbk_genes and gene[1] in gbk_genes:
                # some gene results are placed backwards by PB
                if gbk_genes[gene[0]][0] < gbk_genes[gene[1]][0]:
                    gis.append((gbk_genes[gene[0]][0],gbk_genes[gene[1]][1]))
                else:
                    gis.append((gbk_genes[gene[1]][0],gbk_genes[gene[0]][1]))
        # merge overlaps
        i = 0
        while i < len(gis)-1:
            if gis[i][1] > gis[i+1][0]:
                if gis[i+1][1] > gis[i][1]:
                    gis[i] = (gis[i][0],gis[i+1][1])
                del gis[i+1]
            else:
                i += 1
        for (i,gi) in enumerate(gis):
            output += ("PredictBias_"+str(i+1)+"\t"+str(gi[0])+"\t"+str(gi[1])+"\n")
    o.write(output)
    o.close()