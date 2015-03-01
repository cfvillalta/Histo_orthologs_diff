#!/usr/bin/python 

import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

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
                                                


from ClustalTools import MultipleAlignment
if(__name__=="__main__"):
    hmmalign = open("/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/Lucien_Zn_clu/Zn_clus.hmmalign")
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
    [aa_mean.append(np.mean(aa_cm[:,x])) for x in aa_range]
       #need to mean normalize
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
            bp_aa_ratios[x].append(aa_norm)
    
    print bp_aa_ratios
    fig=plt.figure()
    ax=fig.add_subplot(111)
    bp=ax.boxplot(bp_aa_ratios)
    fig.savefig('test.pdf')
    


'''
#print aa_data
A = []
for aa in aa_data:
    #print aa
    A.append(aa[2][0])

A_np = np.array(A)
mean = np.mean(A_np)
median = np.median(A_np)
print "mean = %s" %(mean)
print "median = %s" %(median)
x = np.random.normal(0,1,1000)
#print x
fig = plt.figure()
ax = fig.add_subplot(111)
numBins = 100
histo=ax.hist(A,numBins,color='green',alpha=0.8)
fig.savefig('fig1.pdf')

fig2 =plt.figure()
ax2=fig2.add_subplot(111)
bp=ax2.boxplot(A)
fig2.savefig('fig2.pdf')
'''
