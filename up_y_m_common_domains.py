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

    #PART III
    for domain in shared_domains:
        d=0
        for gene in up_m_genes:
            for domain_name in up_m_genes[gene]:
                if domain == domain_name:
                    d = d+1
                    #print '%s\t%s\t%s\tMycelia' %(domain,gene,domain_name)
                    shared_domains[domain].append([gene,'mycelia'])
        for gene in up_y_genes:
            for domain_name in up_y_genes[gene]:
                if domain == domain_name:
                    d = d+1
                    shared_domains[domain].append([gene,'yeast'])
                    #print '%s\t%s\t%s\tYeast' %(domain,gene,domain_name)
        shared_domains[domain].append(d)
        d=0
    
#        print "%s\t%s" %(shared_domains[domain],shared_domains[domain][-1])

        final_domain_list = {}
        if shared_domains[domain][-1] >= 2:
            list = range(0,shared_domains[domain][-1])
            y = 0
            m = 0
            for gene in list:
                if 'yeast' in shared_domains[domain][gene]:
                    y = y +1
                if 'mycelia' in shared_domains[domain][gene]:
                    m = m +1

            if y>=1 and m >=1:
                final_domain_list[domain]=shared_domains[domain]
               # print "%s\t%s" %(domain,shared_domains[domain])
                    
            y = 0
            m = 0
            for domain in final_domain_list:
                del final_domain_list[domain][-1]
                domain_out = open('%s_gene_up_list.txt' %(domain), 'w')
                for gene in final_domain_list[domain]:
                    domain_out.write('%s\n' %("\t".join(gene)))
                domain_out.close()
