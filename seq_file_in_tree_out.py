#!/usr/bin/python

'''
1. Open text files in folders with list of genes that share a domain and where at least one of each is sigificantly upregulated in yeast and mycelia. Get prote sequences from gene list.

2. Run HMMsearch using HMM of domain type (e.g. if file is for "Catalase") run HMM search with that domain HMM from pfam. Run a second HMMsearch using the HMM made above G217B sequences. Build HMM from the results of the two HMM searches.

3. Run Lucien using the new HMM.
'''

import sys
import os
from GenomeFactory import GenomeFactory
from subprocess import Popen, PIPE

#1.
if(__name__=="__main__"):
    if(len(sys.argv) != 2):
        sys.stderr.write("""Usage: %s directory_in
        where:
        directory_in is the directory with .txt files named after domain with gene lists
        """)
        sys.exit(0)
    else:
        dir_in = sys.argv[1]
        for filename in os.listdir(dir_in):
            file_name_split = filename.split("_")
            domain = file_name_split[0]
            #print domain
            list_open = open("%s/%s" %(dir_in,filename), "rU")
            gene_list = list_open.readlines()
            #print gene_list
            prot_fasta = open('%s%s_G217B_up.fa' %(dir_in,domain), 'w')
            for gene in gene_list:
                gene=gene.strip()
                #print gene
                gene_split = gene.split("\t")
                #f=GenomeFactory(assembly_caching="eager")
                f=GenomeFactory()
                gene=f.getGene(gene_split[0])
                prot_seq = gene.ProteinSequence()
                prot_fasta.write(prot_seq.FormatFasta(name=gene.Name()))
            prot_fasta.close()
            clustalo = Popen(['time', 'clustalo', '-i', '%s%s_G217B_up.fa' %(dir_in,domain), '-o', '%s%s_G217B_up_clustalo.fa' %(dir_in,domain), '--force', '--threads=4'])
            clustalo.communicate()
            hmmbuild = Popen(['hmmbuild', '--cpu', '4', '%s%s_G217B_up.hmm' %(dir_in,domain), '%s\
            %s_G217B_up_clustalo.fa' %(dir_in,domain)])
            hmmbuild.communicate()
        

