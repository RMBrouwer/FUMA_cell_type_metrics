# This script runs cellex in order to prepare the input for cellect
# Author: Tanya Phung (t.n.phung@vu.nl) - adapted to match the Brainscapes_VTA_human datafile by Rachel Brouwer 19/12/23

import pandas as pd
import argparse
import os



def get_top_perc(ctmatrix,perc, magma_outdir, method,output):
        """ ctmatrix is a pandas dataframe containing a parameter of interest; rows (celltypes) and columns (genes)
            the top "perc" percentage of genes for each celltype will be written to the outputfile 
        """ 
        # Initialize outfile for magma
        magma_outfile = open(os.path.join(magma_outdir, output), "w")
        for ct in ctmatrix.index:
            sorted = ctmatrix.loc[ct].sort_values(ascending=False)
            # remove the nonzero values
            sorted=sorted[sorted>0]
            numexpressed = len(sorted)
            maxincl=sorted[round(perc*numexpressed)-1]
            ct_fmt0 = ct.replace(' ',"_")
            ct_fmt = ct_fmt0.replace('/',"_")
            print(ct)
            print(ct_fmt)
            # write outfile for magma
            magma_out = [ct_fmt]
            for i in ctmatrix.columns[ctmatrix.loc[ct] >= maxincl]:
                magma_out.append(i)
            print("\t".join(magma_out), file=magma_outfile)
        magma_outfile.close()

def main(args):

    spec_cepo_orig = pd.read_table(args.spec_linear, sep="\t")
    spec_cepo_orig.set_index('GENE', inplace=True)
    print(spec_cepo_orig.head)
    spec_cepo = spec_cepo_orig.transpose()
    print(spec_cepo.head)
 
    get_top_perc(spec_cepo, 0.1, args.outdir_magma,"cepo",args.output)
    

def parse_args():
    parser = argparse.ArgumentParser(description='Generate top-10% from continuous data.')
    parser.add_argument('--spec_linear', required=True, help='Path to the specificity file.')
    parser.add_argument('--outdir_magma', required=True, help='Path to the output directory where magma genelists are stored')
    parser.add_argument('--output', required=True, help='Name of the output file')
    return parser.parse_args()


main(parse_args())


