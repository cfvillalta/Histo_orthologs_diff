#!/usr/bin/python
"""
Given log2 Y/H expression ratios and domain annotations for 4 Hc
strains, find domains for which there are both Y _and_ H enriched
genes, and write lists of (gene, enriched state) for each domain for
downstream composition analysis.
"""

if(__name__=="__main__"):
    # Input and find common pfam domains between Yeast and Mycelia pfam
    # domain lists for gets up regulated log 1.5 fold or more.

    m_pfam_in = open("/home/mvoorhie/data/Chris/ThermalAdaptation/M_1.5_domains.txt")
    y_pfam_in = open("/home/mvoorhie/data/Chris/ThermalAdaptation/Y_1.5_domains.txt")

    m_pfam = set(i.strip() for i in m_pfam_in.readlines())
    y_pfam = set(i.strip() for i in y_pfam_in.readlines())

    shared_domains = dict((i,[]) for i in m_pfam.intersection(y_pfam))

    # Search through yeast and mycelia sig up reg .cdt files and get
    # gene names for genes that have the common domains found in shared_domains

    from CdtFile import CdtFile
    cdt_in = CdtFile.fromCdt(open("/home/mvoorhie/data/Chris/ThermalAdaptation/s13cMsort.cdt"))
    up_m_genes_in = open("/home/mvoorhie/data/Chris/ThermalAdaptation/M_1.5_predictions.HcG217B.txt")
    up_y_genes_in = open("/home/mvoorhie/data/Chris/ThermalAdaptation/Y_1.5_predictions.HcG217B.txt")

    cdt_gene_pfam = {} 
    pfam_indices = [n for (n,i) in enumerate(cdt_in.extranames)
                    if(i.endswith("Pfam name"))]
    for gene in cdt_in:
        pfam_name_split = set()
        for i in pfam_indices:
            pfam_name_split.update(gene.extra[i].split("|"))
        for name in gene.Name().split("|"):
            cdt_gene_pfam[name]=pfam_name_split

    up_m_genes = dict((gene.strip(), cdt_gene_pfam[gene.strip()])
                      for gene in up_m_genes_in)

    up_y_genes = dict((gene.strip(), cdt_gene_pfam[gene.strip()])
                      for gene in up_y_genes_in)

    # Determine if the genes with common domains are present in the list
    # of genes sig up in yeast or mycelia.

    for domain in shared_domains:
        d=0
        for (gene, domains) in up_m_genes.items():
            for domain_name in domains:
                if domain == domain_name:
                    d = d+1
                    #print '%s\t%s\t%s\tMycelia' %(domain,gene,domain_name)
                    shared_domains[domain].append([gene,'mycelia'])
        for (gene, domains) in up_y_genes.items():
            for domain_name in domains:
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
                #del final_domain_list[domain][-1]
                domain_out = open('%s_gene_up_list.txt' %(domain), 'w')
                for gene in final_domain_list[domain][:-1]:
                    domain_out.write('%s\n' %("\t".join(gene)))
                domain_out.close()
