#!/usr/bin/env python
from hmmalign2norm import aa_ratio_data

'''
input domains of interest and hmm alignment files. Output in tab delimited format: seqname domain amino acid ratios  amino acid norm rations mean median number of positions yeast/mycelia/none
'''

if(__name__=="__main__"):
    domain_file = "/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/Lucien_eval0.01/yeast_mycelia_pfam_up.txt"
    #domain_file = "/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/pfam_test.txt"
    domain_list = open(domain_file)
    domains = [i.strip() for i in domain_list.readlines()]
    
    amino_acids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    aa_ch = ["Amino_Acid_Counts_%s" %aa for aa in amino_acids]
    aa_rh = ["Amino_Acid_Ratios_%s" %aa for aa in amino_acids]
    n_aa_rh = ["Normalized_Amino_Acid_Ratios_%s" %aa for aa in amino_acids]
    file_out = open("yeast_mycelia_pfam_up_data.txt", "w")
    header = "#Gene\tDomain\t%s\t%s\t%s\n" %("\t".join(aa_ch),"\t".join(aa_rh),"\t".join(n_aa_rh))
    print header
    file_out.write(header)
    for x in range(len(domains)):
        print domains[x]
        data = aa_ratio_data(domains[x])
        print data
        for y in range(len(data)):
            print "%s\t%s"%(x,y)
            #print len(data[0])
        
            col_2 = "\t".join([str(i) for i in data[y][2]])
            col_3 = "\t".join([str(i) for i in data[y][3]])
            col_4 = "\t".join([str(i) for i in data[y][4]])
            data_out = "%s\t%s\t%s\t%s\t%s\n" %(data[y][0],data[y][1],col_2,col_3,col_4)
            file_out.write(data_out)
            print data_out
        
