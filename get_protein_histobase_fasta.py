#!/usr/bin/python
#  chmod "a+x" get_protein_histobase_fasta.py to take advantage of the shebang.
#
#  Top level documentation as string.  Doesn't matter for script,
#  but for an importable module, this is used as the top-level
#  docstring.
"""The script can have a gene id input and will output the protein
sequence extracted from histobase.
"""

# Imports at top of file, external dependencies followed by local dependencies.
#  [This is the convention in the Python world.  I often reverse this order
#   in HistoBase, but I am wrong =) ]
import sys

from GenomeFactory import GenomeFactory

# Distinguish "main" function from local functions and variables.
# (Aids readability and makes it easier to turn scripts into library modules)

if(__name__ == "__main__"):
    # Script with wrong arguments should print documentation.
    # Where possible, script with no arguments should also print documentation.
    if(len(sys.argv) != 3):
        sys.stderr.write("""Usage: %s id_in type_in > id_in.fasta
   where:
      id_in is the systematic gene name and
      type_in is the annotation to be passed to the output FASTA file
""" % sys.argv[0])
        # exit with error code
        sys.exit(0)
        
    id_in = sys.argv[1]
    # I would make type_in optional and provide a sane default (e.g.,
    # no annotation), updating the documentation above.
    type_in =  sys.argv[2]

    f = GenomeFactory()
    gene = f.getGene(id_in)
    prot_seq = gene.ProteinSequence()

    # BUGFIX: fixed variable name
    print prot_seq.FormatFasta(name=gene.Name(), annotation = type_in)
