#!/usr/bin/env python
# Time-stamp: <CompositionHeatmap.py 2013-01-10 17:16:48 Mark Voorhies>
import sys
from CdtFile import CdtFile, CdtRow
from GenomeFactory import GenomeFactory

from optparse import OptionParser

# N.B.: Don't want silent out-of-range bugs, so not doing ascii math
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
aa2offset = dict((aa,i) for (i,aa) in enumerate(
    amino_acids))
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

'''def seqiter(genome):
    for (name, gene) in genome.Genes():
        yield (gene, gene.ProteinSequence())
'''
if(__name__ == "__main__"):
    if(len(sys.argv) != 2):
        sys.stderr.write("""Usage: %s id_in type_in > id_in.fasta
        where:
        id_in is the systemic gene name
        type_in is the annotation to be passed to the output FASTA file"""%(sys.argv[0]))
        sys.exit(0)
    else:    
        f = GenomeFactory()
        #genome = f.getGenome("HcG217B")
        cdt_in = CdtFile.fromCdt(open(sys.argv[1]))
        genes = []
        
        for g in cdt_in:
            print g.Name()
            gene=f.getGene(g.name)      
            seq =gene.ProteinSequence()
            genes.append([gene,seq])
        for gene in genes:
            print len(gene[1])
       
            #Raw lengths and aa counts
        cdt = CdtFile(probes = [CdtRow(
            gid = gene[0].Name(),
            uniqid = gene[0].Name(),
            name = gene[0].Name(),
            ratios = [len(gene[1])]+composition_vector(gene[1]))
                    for gene in genes
                        if(len(gene[1]) > 0)],
                      fieldnames = ["length"]+list(amino_acids),
                      eweights = [1.]*21)

            #Normalize counts to length
        cdt = CdtFile.fromPrototype(
            cdt,
            probes = [CdtRow.fromPrototype(i,ratios = i.ratios+[
                j/float(i.ratios[0]) for j in i[1:]])
                      for i in cdt],
            eweights = cdt.eweights + [1.]*20,
            fieldnames = cdt.fieldnames + [j.lower() for j in cdt.fieldnames[1:]])

    # Write unsorted data
        filename = sys.argv[1].split(".cdt")
        cdt.writeCdt(open("%s_aa_ratio.cdt" %(filename[0]),"w"))

    # Cluster on normalized counts
    tree = cdt.cluster(dist = "e", method = "m", cols = range(21,41))
    cdt.writeCdtGtr("HcG217B.composition.norm_em", tree)
                                    
