#!/usr/bin/env python
'''
Will input list of domains, run lucien with histo transcriptome data for each domain, from each hmm domain alignment will get amino acid ratios and normalize by mean. Plot unnormalized and normalized data boxplots.

-L force runs lucien.
'''
import sys
import re
from hmmalign2norm import aa_mean_normilization
from hmmalign2norm import aa_mean_gene
import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
import os
from ClustalTools import MultipleAlignment
import random

if(__name__=="__main__"):
    domain_in = "/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/yeast_mycelia_pfam_up.txt"
    domain_in = "/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/pfam_test.txt"

    domain_list = open(domain_in)

    domains = [i.strip() for i in domain_list.readlines()]
    domains_aa_norm = []
    domains_len=len(domains)
    num=0
    for domain in domains:
        num=num+1
        #if this option present force run Lucien.py
        force_lucien = re.compile('-L')
        run_lucien = force_lucien.search(str(sys.argv))
        if run_lucien: 
            print 'running Lucien.py on %s, %s out of %s' %(domain,num,domains_len)
            Lucien = Popen(['Lucien.py','-a', '-D', domain, '-G', '-g', '/home/cfvillalta/AssembledTranscriptomes/ForLucien/'],stdout=PIPE)
            Lucien.communicate()
            print 'Lucien run done.'

        if os.path.isfile("%s.hmmalign" %(domain)):
            print 'hmmalign file present. Lucien will not run on %s, %s out of %s' %(domain,num,domains_len)
            if os.stat("%s.hmmalign" %(domain)).st_size == 0:
                print "Hmmalign file for %s is empty or misformatted" %(domain)
                domains_aa_norm.append('None')
            else:
                hmmalign =open("%s.hmmalign" %(domain))
                aa_mean_norm = aa_mean_normilization(hmmalign)
                domains_aa_norm.append(aa_mean_norm)            
        else:
            print 'running Lucien.py on %s, %s out of %s' %(domain,num,domains_len)
            Lucien = Popen(['Lucien.py','-a', '-D', domain, '-G', '-g', '/home/cfvillalta/AssembledTranscriptomes/ForLucien/'],stdout=PIPE)
            Lucien.communicate()
            print 'Lucien run done.'

    #domain seqs graphing portion.

    bp_range = range(10)
    bp_cord = []
    for x in bp_range:
        for y in bp_range:
            bp_cord.append([x,y])
    

    seq_counts=[]
    for x in range(domains_len):
        if domains_aa_norm[x] != 'None':
            stockholm = open("%s.hmmalign" %(domains[x])) 
            alignment = MultipleAlignment.fromStockholm(stockholm)
#            print len(alignment.seqnames)
            seq_counts.append(len(alignment.seqnames))
        else:
            seq_counts.append(0)
#    print seq_counts
    f, ax = plt.subplots(10,10,sharex=True,sharey=True)
    f.set_size_inches(40,40)
    plt.tight_layout()
    amino_acids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    for x in range(domains_len):
        if domains_aa_norm[x] != 'None':
            ax[bp_cord[x][0],bp_cord[x][1]].boxplot(domains_aa_norm[x][2])
 #           print domains[x]
            ax[bp_cord[x][0],bp_cord[x][1]].set_title('%s - %s seqs' %(domains[x],seq_counts[x]))
            ax[bp_cord[x][0],bp_cord[x][1]].set_xticklabels(amino_acids)
        else:
#            print domains[x]
            ax[bp_cord[x][0],bp_cord[x][1]].set_title('%s' %(domains[x]))
            ax[bp_cord[x][0],bp_cord[x][1]].set_xticklabels(amino_acids,visible=True)
    f.savefig('mean_normalized.pdf')

    #up in yeast and mycelia graphing portion.
    
    genes_in = open("/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/yeast_mycelia_pfam_up_domain_seq.txt")
    genes = [i.strip() for i in genes_in.readlines()]
    genes_s = [i.split("\t") for i in genes]

    #list of lists name of gene and amino acid frequencies.

    f1, ax1 = plt.subplots(10,10,sharex=True,sharey=True)
    f1.set_size_inches(50,50)
    plt.setp(ax1,xticks=range(len(amino_acids)))
    plt.ylim(-0.2,1)
    
    #ax1=plt.gca()
    #f1.tight_layout()


    for x in range(domains_len):
        goi = []
        #print x
        print domains[x]
        if domains_aa_norm[x] != 'None':
            #print hmmalign
            medians = []
            for aa in range(len(amino_acids)):
                #print amino_acids[aa]
                #print domains_aa_norm[x][2][aa]
                median = np.median(np.array(domains_aa_norm[x][2][aa]))
                # print median
                medians.append(median)

            domains_aa_norm[x].append(medians)

            for gene in range(len(genes_s)):

                #print "%s\t%s" %(domains[x],gene[1])
                if domains[x] == genes_s[gene][1]:
                    #print "should be present %s" %( genes_s[gene][0])
                    #print domains_aa_norm[x][0]
                    for n in range(len(domains_aa_norm[x][0])):
                        
                        #print domains_aa_norm[x][0][n]
                        if genes_s[gene][0] == domains_aa_norm[x][0][n]:
                            genes_s[gene].append(domains[x])
                            goi.append(genes_s[gene]) 
        if not goi:            
            print 'None'
        else:
            gene_aa_ratios = aa_mean_gene(goi)
            for aa in range(len(amino_acids)):
                a=aa-.15
                b=aa+.15
                x_axis = [random.uniform(a,b) for p in range(0, len(gene_aa_ratios))]
                ax1[bp_cord[x][0],bp_cord[x][1]].plot(aa,(domains_aa_norm[x][4][aa]),'g_',ms=10,mew=5)
                for g in range(len(gene_aa_ratios)):
                    if gene_aa_ratios[g][2] == 'mycelia':
                        #print 'mycelia'
                        ax1[bp_cord[x][0],bp_cord[x][1]].plot(x_axis[g],(gene_aa_ratios[g][1][aa]-domains_aa_norm[x][3][aa]),'bo')
                    elif gene_aa_ratios[g][2] == 'yeast':  
                        #print 'yeast'
                        ax1[bp_cord[x][0],bp_cord[x][1]].plot(x_axis[g],(gene_aa_ratios[g][1][aa]-domains_aa_norm[x][3][aa]),'ro')
        ax1[bp_cord[x][0],bp_cord[x][1]].set_title('%s' %(domains[x]))
        ax1[bp_cord[x][0],bp_cord[x][1]].set_xticklabels(amino_acids,visible=True)
    f1.savefig('mean_normalized_goi.pdf')


