#!/usr/bin/python 

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
        #print x.seqnames
        
        
        
        
