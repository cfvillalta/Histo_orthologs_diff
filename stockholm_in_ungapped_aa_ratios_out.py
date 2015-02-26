#!/usr/bin/python 

from ClustalTools import MultipleAlignment
if(__name__=="__main__"):
    hmmalign = open("/home/cfvillalta/ThermalAdaptation/thermal_domain_ratio_20150226/Lucien_Zn_clu/Zn_clus.hmmalign")
    alignment = MultipleAlignment.fromStockholm(hmmalign)

    print alignment
