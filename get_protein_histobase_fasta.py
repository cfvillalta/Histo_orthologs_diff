#!/usr/bin/python 
#The script can have a gene id input and will output the protein sequence extracted from histobase.

import sys

id_in = sys.argv[1]
type_in =  sys.argv[2]

from GenomeFactory import GenomeFactory

f = GenomeFactory()

gene = f.getGene(id_in)

prot_seq = gene.ProteinSequence()

print seq.FormatFasta(name=gene.Name(), annotation = type_in)
