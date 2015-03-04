#!/usr/bin/python
'''
Will input list of domains, run lucien with histo transcriptome data for each domain, from each hmm domain alignment will get amino acid ratios and normalize by mean. Plot unnormalized and normalized data boxplots.
'''

from hmmalign2norm import aa_mean_normilization
import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
import os


if(__name__=="__main__"):
    domain_in = "/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/yeast_mycelia_pfam_up.txt"
    domain_list = open(domain_in)

    domains = [i.strip() for i in domain_list.readlines()]
    domains_aa_norm = []
    domains_len=len(domains)
    num=0
    for domain in domains:
        num=num+1
        print 'running Lucien.py on %s, %s out of %s' %(domain,num,domains_len)
        Lucien = Popen(['Lucien.py', '-D', domain, '-G', '-g', '/home/cfvillalta/AssembledTranscriptomes/ForLucien/'],stdout=PIPE)
        
        Lucien.communicate()
        print 'Lucien run done.'

        hmmalign =open("%s.hmmalign" %(domain))
        aa_mean_norm = aa_mean_normilization(hmmalign)
        domains_aa_norm.append(aa_mean_norm)
    #dir =  os.path.dirname(os.path.abspath(domain_in))

#    print domains_len
    bp_range = range(10)
    bp_cord = []
    for x in bp_range:
        for y in bp_range:
            bp_cord.append([x,y])
    #print bp_cord

    f, ax = plt.subplots(10,10,sharex=True,sharey=True)
    f.set_size_inches(40,40)
    plt.tight_layout()
    amino_acids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    for x in range(domains_len):
        ax[bp_cord[x][0],bp_cord[x][1]].boxplot(domains_aa_norm[x][2])
        #print domains[x]
        ax[bp_cord[x][0],bp_cord[x][1]].set_title('%s' %(domains[x]))
        ax[bp_cord[x][0],bp_cord[x][1]].set_xticklabels(amino_acids)
    #ax.plot(x, y)
    #ax.set_title('Simple plot')
    #fig=plt.subplots(1,3)
    #plt.show()
    #print 'done'
    #print plt.get_backend()
    f.savefig('mean_normalized.pdf')
#print domains

'''
hmmalign = open("/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/Lucien_Zn_clu/Zn_clus.hmmalign")

x =  aa_mean_normilization(hmmalign)

print x 
'''