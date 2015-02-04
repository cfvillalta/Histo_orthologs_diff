#!/usr/bin/python 
"""
The script can have a gene id input and will output the protein sequence extracted from histobase.
"""

import sys
from GenomeFactory import GenomeFactory

if(__name__ == "__main__"):
    if(len(sys.argv) != 3):
        sys.stderr.write("""Usage: %s id_in type_in > id_in.fasta
        where:
        id_in is the systemic gene name
        type_in is the annotation to be passed to the output FASTA file"""%(sys.argv[0]))
        sys.exit(0)
    else:    
        id_in = sys.argv[1]
        type_in =  sys.argv[2]
        f = GenomeFactory()
        gene = f.getGene(id_in)
        prot_seq = gene.ProteinSequence()
        #BUGFIX: fixed variable name
        print prot_seq.FormatFasta(name=gene.Name(), annotation = type_in)
