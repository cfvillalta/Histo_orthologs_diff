#!/usr/bin/python
"""
Given log2 Y/H expression ratios and domain annotations for 4 Hc
strains, find domains for which there are both Y _and_ H enriched
genes, and write lists of (gene, enriched state) for each domain for
downstream composition analysis.
"""

import sys
from CdtFile import CdtFile, CdtRow
from GenomeFactory import GenomeFactory
from optparse import OptionParser

#Copy code in from CompositionHeatmap.py to get amino acid freqs.

# N.B.: Don't want silent out-of-range bugs, so not doing ascii math
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
aa2offset = dict((aa,i) for (i,aa) in enumerate(amino_acids))
for (key,val) in aa2offset.items():
    aa2offset[key.lower()] = val

def composition_vector(s):
    v = [0.]*20
    for i in s:
        try:
            v[aa2offset[i]] += 1.
        except KeyError:
            pass
    return v

ivywrel_aa = "IVYWREL"
ivywrel2offset = dict((aa,i) for (i,aa) in enumerate(ivywrel_aa))
for (key,val) in ivywrel2offset.items():
    ivywrel2offset[key.lower()] = val

    
def ivywrel_comp_vector(s):
    v = [0.]*7
    for i in s:
        try:
            v[ivywrel2offset[i]] += 1.
        except KeyError:
            pass
    return v


if(__name__=="__main__"):
    # Input and find common pfam domains between Yeast and Mycelia pfam
    # domain lists for gets up regulated log 1.5 fold or more.

    m_pfam_in = open("/home/cfvillalta/ThermalAdaptation/M_1.5_domains.txt")
    y_pfam_in = open("/home/cfvillalta/ThermalAdaptation/Y_1.5_domains.txt")

    m_pfam = set(i.strip() for i in m_pfam_in.readlines())
    y_pfam = set(i.strip() for i in y_pfam_in.readlines())

    shared_domains = dict((i,[]) for i in m_pfam.intersection(y_pfam))

    # Search through yeast and mycelia sig up reg .cdt files and get
    # gene names for genes that have the common domains found in shared_domains
    cdt_file = "/home/cfvillalta/ThermalAdaptation/s13cYsort.cdt"
    cdt_in = CdtFile.fromCdt(open(cdt_file))
    up_m_genes_in = open("/home/cfvillalta/ThermalAdaptation/M_1.5_predictions.HcG217B.txt")
    up_y_genes_in = open("/home/cfvillalta/ThermalAdaptation/Y_1.5_predictions.HcG217B.txt")

    cdt_gene_pfam = {} 
    pfam_indices = [n for (n,i) in enumerate(cdt_in.extranames)
                    if(i.endswith("Pfam name"))]
    for gene in cdt_in:
        pfam_name_split = set()
        for i in pfam_indices:
            pfam_name_split.update(gene.extra[i].split("|"))
        for name in gene.Name().split("|"):
            #some genes have two names separated by |
            cdt_gene_pfam[name]=pfam_name_split
            #gene.extra[2] prints out the transcirpt name.
            #print cdt_gene_pfam
    up_m_genes = dict((gene.strip(), cdt_gene_pfam[gene.strip()])
                      for gene in up_m_genes_in)

    up_y_genes = dict((gene.strip(), cdt_gene_pfam[gene.strip()])
                      for gene in up_y_genes_in)

    # Determine if the genes with common domains are present in the list
    # of genes sig up in yeast or mycelia.

    shared_domain_counts = {}
    shared_domain_genes = {}
    
    for domain in sorted(shared_domains):
        d=0
        for (gene, domains) in up_m_genes.items():
            if domain in domains:
                d = d+1
                #print '%s\t%s\t%s\tMycelia' %(domain,gene,domain_name)
                shared_domains[domain].append([gene,'mycelia'])
        for (gene, domains) in up_y_genes.items():
            if domain in domains:
                d = d+1
                shared_domains[domain].append([gene,'yeast'])
                #print '%s\t%s\t%s\tYeast' %(domain,gene,domain_name)
        shared_domain_counts[domain] = d
    
#        print "%s\t%s" %(shared_domains[domain],shared_domains[domain][-1])

        if shared_domain_counts[domain] >= 2:
            y = 0
            m = 0
            for (gene, state) in shared_domains[domain]:
                if state == 'yeast':
                    y = y +1
                else:
                    assert(state == 'mycelia')
                    m = m +1

            if y>=1 and m >=1:
                shared_domain_genes[domain]=shared_domains[domain]
    #print shared_domain_genes
    f=GenomeFactory()
    genes = []
    uid=0
    #write domains to txt file.
    domain_file = open("yeast_mycelia_pfam_up.txt", "w")
    domain_gene_file = open("yeast_mycelia_pfam_up_domain_seq.txt", "w")
    for domain in shared_domain_genes:
        domain_file.write("%s\n" %(domain))
        domain_fa = open('%s_G217B_up.fa' %(domain),'w')
        for g in shared_domain_genes[domain]: 
            gene = f.getGene(g[0])
            seq = gene.ProteinSequence()
            #print ">%s\n%s" %(g[0],seq)
            domain_fa.write(">%s\n%s\n" %(g[0],seq))
            uid = uid+1
            ivywrel = sum(ivywrel_comp_vector(seq))
            #print seq
            domain_gene_file.write("%s\t%s\t%s\n" %(gene,domain,g[1]))
            if g[1] == 'yeast':
                genes.append([gene,seq,g[1],domain,uid, 4,ivywrel])
            else:
                assert(g[1] == 'mycelia')
                genes.append([gene,seq,g[1],domain,uid, -4,ivywrel])
    domain_file.close()
    domain_gene_file.close()



'''
    #print seq
            #print composition_vector(seq)

            #ratios = [-4, len(seq)] + composition_vector(seq) + [ivywrel]

            #print ratios
    fieldnames = ["morph","length"]+list(amino_acids)+ ["IVYWREL"]
#    print len(fieldnames)
    cdt = CdtFile(
        probes = [CdtRow(gid = gene[0].Name(),uniqid = gene[4],name = gene[0].Name(),extra=[gene[2],gene[3]],ratios = [float(gene[5]),len(gene[1])]+composition_vector(gene[1])+[gene[6]]) for gene in genes if(len(gene[1]) > 0)],
        fieldnames = fieldnames,
        eweights = [1.]*23,
        extranames = ["morphology", "domain"])
     
     #Normalize counts to length
    cdt = CdtFile.fromPrototype(cdt,probes = [CdtRow.fromPrototype(i,ratios = i.ratios+[j/float(i.ratios[1]) for j in i[2:]])for i in cdt],eweights = cdt.eweights + [1.]*21,fieldnames = cdt.fieldnames + [j.lower() for j in cdt.fieldnames[2:]])
    
    # Write unsorted data
    filename = cdt_file.split(".cdt")
    cdt.writeCdt(open("%s_aa_ratio.cdt" %(filename[0]),"w"))
    
    #Cluster on normalized counts
    tree = cdt.cluster(dist = "e", method = "m", cols = range(22,42))
    cdt.writeCdtGtr("%s.composition.norm_em" %(filename[0]), tree)

    #Cluster on normalized counts
    tree = cdt.cluster(dist = "e", method = "m", cols = [43])
    cdt.writeCdtGtr("%s.composition.ivywrel" %(filename[0]), tree)
'''
