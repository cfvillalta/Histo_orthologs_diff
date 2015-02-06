#!/usr/bin/python

"""
Part I
Input and find common pfam domains between Yeast and Mycelia pfam domain lists for gets up regulated log 1.5 fold or more.

Part II
Search through yeast and mycelia sig up reg .cdt files and get gene names for genes that have the common domains found in Part I

Part III
Determine if the genes found in Part II are present in the list of genes sig up in yeast or mycelia.
"""

#Part I

if(__name__=="__main__"):
    m_pfam_in = open("/home/cfvillalta/ThermalAdaptation/M_1.5_domains.txt")
    y_pfam_in = open("/home/cfvillalta/ThermalAdaptation/Y_1.5_domains.txt")

    m_pfam = m_pfam_in.readlines()
    y_pfam = y_pfam_in.readlines()

    shared_domains = {}
    
    for domain in m_pfam:
        if any(domain in y for y in y_pfam):    
            domain=domain.strip()
            shared_domains[domain]=[]
        else:
            pass
    #PART II
    
    cdt_in = open("/home/cfvillalta/ThermalAdaptation/s13cMsort.cdt")
    up_m_genes_in = open("/home/cfvillalta/ThermalAdaptation/M_1.5_predictions.HcG217B.txt")
    up_y_genes_in = open("/home/cfvillalta/ThermalAdaptation/Y_1.5_predictions.HcG217B.txt")

    cdt_gene_pfam = {} 
    up_m_genes = {}
    up_y_genes = {}

    for gene in cdt_in.readlines():
        gene=gene.strip()
        gene_split = gene.split("\t")
        pfam_name_split = gene_split[14].split("|")
        cdt_gene_pfam[gene_split[2]]=pfam_name_split
    
    for gene in up_m_genes_in.readlines():
        gene = gene.strip()
        if gene in cdt_gene_pfam:
            up_m_genes[gene]=cdt_gene_pfam[gene]

    for gene in up_y_genes_in.readlines():
        gene = gene.strip()
        if gene in cdt_gene_pfam:
            up_y_genes[gene]=cdt_gene_pfam[gene]

    
    
    
