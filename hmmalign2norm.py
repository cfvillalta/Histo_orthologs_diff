#!/usr/bin/python 

import numpy as np
import matplotlib as mpl
mpl.use('agg')
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
    aa_mean = [[] for _ in aa_range]
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
            aa_mean[x].append(aa_mean[x])
    return([names, bp_aa_ratios, bp_aa_ratios_norm,aa_mean])

def aa_mean_gene(goi):
    goi_out = []
    for g in goi:
        hmmalign=open("%s.hmmalign" %(g[1]))
        alignment = MultipleAlignment.fromStockholm(hmmalign)
        #print 'stockholm uploaded'
        #print hmmalign
        ref= alignment.colAnnotations
        #    print ref.keys()
        x_indices = [n for (n,i) in enumerate(ref['RF']) if (i=="x")]
    #print x_indices
        names = ''
        aa_counts =[]
        aa_ratios = []
        for (seq,name) in zip(alignment, alignment.seqnames):
            #print name
            #print len(genes)
            if name == g[0]:
                #print name            
                #print "%s\n%s" %(name, "".join(seq))
                aa_c = composition_vector(seq)
                aa_r = [float(aa/sum(aa_c)) for aa in aa_c]
                names=name
                aa_counts.append(aa_c)
                aa_ratios.append(aa_r)

        aa_cm = np.array(aa_counts)
        aa_rm = np.array(aa_ratios)
        aa_mean = []
        aa_range=range(0,len(amino_acids))
    #for  in aa_range:
        
        #need to mean normalize
        
        bp_aa_ratios = [[] for _ in aa_range]
        #print bp_aa_ratios
        for x in aa_range:
            for aa in aa_ratios:  
                #print amino_acids[x]
                #print aa[x]
                #print aa_mean[x]
                
                #print aa_norm
                #print x
                
                bp_aa_ratios[x].append(aa[x])
        goi_out.append([names, bp_aa_ratios,g[2]])
    
    return(goi_out)

#hmmalign =open("SMC_N.hmmalign")
#print aa_mean_gene('HISTO_DM.Contig933.Fgenesh_Aspergillus.157.final_new', hmmalign)
#print aa_mean_gene('HISTO_DA.Contig93-snap.315.final_new', hmmalign)
