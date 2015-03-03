#!/usr/bin/python
'''
Will input list of domains, run lucien with histo transcriptome data for each domain, from each hmm domain alignment will get amino acid ratios and normalize by mean. Plot unnormalized and normalized data boxplots.
'''

from hmmalign2norm import aa_mean_normilization
#import numpy as np
#import matplotlib as mpl
#mpl.use('agg')
#import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
import os


if(__name__=="__main__"):
    domain_in = "/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/pfam_test.txt"
    domain_list = open(domain_in)

    domains = [i.strip() for i in domain_list.readlines()]
    for domain in domains:
        print 'running Lucien.py on %s' %(domain)
#        Lucien = Popen(['Lucien.py', '-D', domain, '-G', '-g', '/home/cfvillalta/AssembledTranscriptomes/ForLucien/'],stdout=PIPE)
        
        Lucien.communicate()
        print 'Luncien run done.'

        hmmalign =open("%s.hmmalign" %(domain))
        aa_mean_norm = aa_mean_normilization(hmmalign)
    
    #dir =  os.path.dirname(os.path.abspath(domain_in))
    
    
#print domains

'''
hmmalign = open("/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/Lucien_Zn_clu/Zn_clus.hmmalign")

x =  aa_mean_normilization(hmmalign)

print x 
'''
