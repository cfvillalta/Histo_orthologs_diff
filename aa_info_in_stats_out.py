#!/usr/bin/env python


"""
input data for genes of interest and their domain info and amino acid composition info.

Run a wilcoxon text on norm amino acid ratios for each domain.

"""
import re
import numpy as np
import matplotlib as mpl
mpl.use('agg')
from scipy.stats import mannwhitneyu
import random
import matplotlib.pyplot as plt

if(__name__=="__main__"):
    data_file = "/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/Lucien_eval0.01/yeast_mycelia_pfam_up_data.txt"
    data_in = open(data_file)


    genes_in = open("/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/yeast_mycelia_pfam_up_domain_seq.txt")
    genes = [i.strip() for i in genes_in.readlines()]
    genes_s = [i.split("\t") for i in genes]
    #print genes_s        

    data_in = [i.split("\t") for i in [i.strip() for i in data_in.readlines()]]
    #print data_in
    header = data_in[0]
    n=0
    for x in header:        
        #print "%s\t%s" % (n,x)
        n=n+1

    data_in.remove(data_in[0])
    domains = set(data_in[x][1] for x in range(len(data_in)))
    # print domains
    
    aa_range = range(42,62)
    goi = [g[0] for g in genes_s]
    #print goi
    h = re.compile("|".join(goi))
    
           
    aa_rnm = [[] for _ in aa_range]
    aa_rny = [[] for _ in aa_range]
       

    for d in domains:
        #print d
        aa_rn = [[] for _ in aa_range]
        median = []
        for x in range(len(data_in)):
            if data_in[x][1] == d:
                aa_num = 0
                for aa in aa_range:
                    aa_rn[aa_num].append(float(data_in[x][aa]))
                    aa_num = aa_num+1
                find_h = h.search(data_in[x][0])
                aa_num = 0
                if find_h:
                    #print data_in[x][0]
                    for g in genes_s:
                        if g[0] == data_in[x][0] and g[1] == d:
                            if g[2] == 'mycelia':
                                for aa in aa_range:
                                    aa_rnm[aa_num].append(float(data_in[x][aa]))
                                    aa_num = aa_num+1
                            if g[2] == 'yeast':
                                for aa in aa_range:
                                    aa_rny[aa_num].append(float(data_in[x][aa]))
                                    aa_num = aa_num+1
                            #print "%s\t%s" %(g[1],g[0])
                            
        [median.append(np.median(np.array(aa))) for aa in aa_rn]

        #print median
    amino_acids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    f1,ax1 = plt.subplots(1,1,sharex=True,sharey=True)
    f1.set_size_inches(10,10)
    plt.ylim(-0.2,1)
    plt.xlim(0,21)
    plt.setp(ax1,xticks=range(len(amino_acids)))
    ax1.set_xticklabels(amino_acids,visible=True)
    for aa in range(len(aa_range)):
        
        mw_test =  mannwhitneyu(aa_rnm[aa],aa_rny[aa])
        print header[aa_range[aa]]
        print "median mycelia  %s" %(np.median(np.array(aa_rnm[aa])))
        print "median yeast  %s" %(np.median(np.array(aa_rny[aa])))
        print mw_test

        a=aa-.15
        b=aa+.15

        x_axis_m = [random.uniform(a,b) for p in range(0, len(aa_rnm[aa]))]
        x_axis_y = [random.uniform(a,b) for p in range(0, len(aa_rny[aa]))]

        ax1.plot(x_axis_m,aa_rnm[aa],'bo') 
        ax1.plot(x_axis_y,aa_rny[aa],'ro') 


    f1.savefig('aa_ratios_all_domains_dp.pdf')
