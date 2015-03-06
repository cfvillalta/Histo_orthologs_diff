#!/usr/bin/python 

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from ClustalTools import MultipleAlignment

# N.B.: Don't want silent out-of-range bugs, so not doing ascii math
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
aa2offset = dict((aa,i) for (i,aa) in enumerate(amino_acids))
for (key,val) in aa2offset.items():
    aa2offset[key.lower()] = val
    
    def composition_vector(s):
        v = [0.]*20
        for i in s:
            try:
                v[aa2offset[i]] += 1.
            except KeyError:
                pass
        return v
#module input is hmmaling stockholm file, and outputs a list of sublists: name, aa_ratios, aa_norm_ratios.                         
def aa_mean_normilization(hmmalign):
#if(__name__=="__main__"):
#    hmmalign = open("/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/Lucien_Zn_clu/Zn_clus.hmmalign"
#    hmmalign = open("/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/Lucien_Zn_clu/Zn_clus_test.stockholm")
    alignment = MultipleAlignment.fromStockholm(hmmalign)
#    print alignment.seqnames
    ref= alignment.colAnnotations
#    print ref.keys()
    x_indices = [n for (n,i) in enumerate(ref['RF']) if (i=="x")]
    #print x_indices
    names = []
    aa_counts =[]
    aa_ratios = []
    for (seq,name) in zip(alignment, alignment.seqnames):
        #print "%s\n%s" %(name, "".join(seq))
        aa_c = composition_vector(seq)
        aa_r = [float(aa/sum(aa_c)) for aa in aa_c]
        #print aa_ratios
        #print aa_counts
        #print sum(aa_counts)
        names.append(name)
        aa_counts.append(aa_c)
        aa_ratios.append(aa_r)

    aa_cm = np.array(aa_counts)
    aa_rm = np.array(aa_ratios)
    aa_mean = []
    aa_range=range(0,len(amino_acids))
    #for x in aa_range:
    [aa_mean.append(np.mean(aa_rm[:,x])) for x in aa_range]
       #need to mean normalize
    bp_aa_ratios_norm = [[] for _ in aa_range]
    bp_aa_ratios = [[] for _ in aa_range]
    #print bp_aa_ratios
    for x in aa_range:
        for aa in aa_ratios:  
            #print amino_acids[x]
            #print aa[x]
            #print aa_mean[x]
            aa_norm=aa[x]-aa_mean[x]
            #print aa_norm
           # print x
            bp_aa_ratios_norm[x].append(aa_norm)
            bp_aa_ratios[x].append(aa[x])
    return(names, bp_aa_ratios, bp_aa_ratios_norm)
'''    
    #print bp_aa_ratios
    fig=plt.figure()
    ax=fig.add_subplot(211)
    bp=ax.boxplot(bp_aa_ratios_norm)
    ax.set_xticklabels(amino_acids)
    #fig.savefig('test.pdf')

    ax=fig.add_subplot(212)
    bp=ax.boxplot(bp_aa_ratios)
    ax.set_xticklabels(amino_acids)
    fig.savefig('test.pdf')
'''    
