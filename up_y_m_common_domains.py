#!/usr/bin/python

"""
Part I
Input and find commpn pfam domains between Yeast and Mycelia pfam domain lists for gets up regulated log 1.5 fold or more.

Part II
Search through yeast and myclia sig up reg .cdt files and get gene names for genes that have the common domains found in Part I

Part III
Determine if the genes found in Part II are present in the list of genes sig up in yeast or mycelia.
"""

#Part I

if(__name__=="__main__"):
    m_pfam_in = open("/home/cfvillalta/ThermalAdaptation/M_1.5_domains.txt")
    y_pfam_in = open("/home/cfvillalta/ThermalAdaptation/Y_1.5_domains.txt")

    m_pfam = m_pfam_in.readlines()
    y_pfam = y_pfam_in.readlines()

    print m_pfam
    print y_pfam
