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

#Part II
    m_up_cdt = open("/home/cfvillalta/ThermalAdaptation/s13cMsort.cdt")
    y_up_cdt = open("/home/cfvillalta/ThermalAdaptation/s13cYsort.cdt")


#    x=0
#    for domain in shared_domains:
#        x=x+1
#        print domain

#    print x
#    print m_pfam
#    print y_pfam
