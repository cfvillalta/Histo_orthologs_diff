#!/usr/bin/python
'''
Will input list of domains, run lucien with histo transcriptome data for each domain, from each hmm domain alignment will get amino acid ratios and normalize by mean. Plot unnormalized and normalized data boxplots.
'''
from hmmalign2norm import aa_mean_normilization

hmmalign = open("/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/Lucien_Zn_clu/Zn_clus.hmmalign")

x =  aa_mean_normilization(hmmalign)

print x 
