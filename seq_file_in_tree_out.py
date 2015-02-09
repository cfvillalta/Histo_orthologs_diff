#!/usr/bin/python

'''
1. Open text files in folders with list of genes that share a domain and where at least one of each is sigificantly upregulated in yeast and mycelia.

2. Run HMMsearch using HMM of domain type (e.g. if file is for "Catalase") run HMM search with that domain HMM from pfam. Run a second HMMsearch using the HMM made above G217B sequences. Build HMM from the results of the two HMM searches.

3. Run Lucien using the new HMM.
'''

if(__name__=="__main__"):
    if(len(sys.argv) != 2):
        sys.stderr.write("""Usage: %s directory_in
        where:
        directory_in is the directory with .txt files named after domain with gene lists
        """)
        sys.exit(0)
    else:
        for filename in os.listdir(os.getcwd()):
            print filename
    

