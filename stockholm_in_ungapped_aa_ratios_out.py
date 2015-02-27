#!/usr/bin/python 

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
    hmmalign = open("/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/Lucien_Zn_clu/Zn_clus_test.stockholm")
    alignment = MultipleAlignment.fromStockholm(hmmalign)
#    print alignment.seqnames
    ref= alignment.colAnnotations
#    print ref.keys()
    x_indices = [n for (n,i) in enumerate(ref['RF']) if (i=="x")]
    #print x_indices

    for (seq,name) in zip(alignment, alignment.seqnames):
        print "%s\n%s" %(name, "".join(seq))
        aa_counts = composition_vector(seq)
        aa_ratios = [float(aa/sum(aa_counts)) for aa in aa_counts]
        print aa_ratios
        print aa_counts
        print sum(aa_counts)
        
        
        
        
