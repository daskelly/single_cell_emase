#!/usr/bin/env python

"""
Convert a *.total.gene.counts matrix from single cell EMASE
into a single matrix of total counts.
Our output format is MTX (for sparse matrices). This follows
a common convention in scRNA-Seq.
By default we round EMASE fractional counts to the nearest int.
"""
import sys, argparse
import numpy as np
from scipy import sparse
from scipy.io import mmwrite
import mygene

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("countsfile", metavar="file.total.gene.counts")
    parser.add_argument("outdir", help="dir to write barcodes, features, matrix.mtx")
    parser.add_argument("--round_ndigits", metavar="N", default=0, type=int)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    
    # Read through the file once and get all barcodes, plus pointers
    # to where in the file the reads start
    f = open(args.countsfile, 'rb')   # rb since using tell and a for loop
    positions = dict()
    this_pos = f.tell()
    for line in f:
        if line.startswith(b'#'):
            line = line.decode('utf-8')
            if line.startswith('#target_id'):
                assert line.rstrip() =="#target_id\tA\tB\tC\tD\tE\tF\tG\tH\ttotal"
            else:
                assert line.startswith('#sample_id')
                fields = line.rstrip().split(' ')
                assert len(fields) == 2, "fields was {}".format(fields)
                barcode = fields[1]
                positions[barcode] = this_pos
        this_pos = f.tell()
    f.close()
    
    if args.verbose:
        print("Reformatting file into a matrix of {} cells".format(len(positions)))
    
    # Now read through and get the counts for each cell
    counts_list = dict()
    for i, barcode in enumerate(positions):
        cell_data = read_sample(args.countsfile, positions[barcode], args.round_ndigits)
        if i == 0:
            expected_genes = cell_data[0]
        else:
            # double check the genes are in same order for every cell
            assert cell_data[0] == expected_genes
            counts_list[barcode] = cell_data[1]
    
    # Get the gene symbol for each gene ID
    mg = mygene.MyGeneInfo()
    symb = mg.querymany(expected_genes, scopes='ensembl.gene', returnall=True)
    symbols = dict()
    for i in range(len(symb['out'])):
        gene_id = symb['out'][i]['query']
        try:
            gene_symbol = symb['out'][i]['symbol']
        except KeyError:
            assert gene_id in symb['missing']
            gene_symbol = gene_id
        symbols[gene_id] = gene_symbol
        
    # Print out a matrix of all the counts
    # First make a matrix, then write it as a sparse MTX to args.outdir
    barcodes = sorted(counts_list.keys())
    bc_out = open('{}/barcodes.tsv'.format(args.outdir), 'w')
    for b in barcodes:
        bc_out.write('{}\n'.format(b))
    bc_out.close()
    feat_out = open('{}/features.tsv'.format(args.outdir), 'w')
    for gene in expected_genes:
        feat_out.write('{}\t{}\n'.format(gene, symbols[gene]))
    feat_out.close()
    tmat = np.transpose(np.array([counts_list[b] for b in barcodes]))   # genes in rows, cells in columns
    tmat_sparse = sparse.csc_matrix(tmat)
    mmwrite("{}/matrix.mtx".format(args.outdir), tmat_sparse)

def read_sample(f, pos, round_digits):
    genes, counts = list(), list()
    if round_digits == 0:
        def to_num(x):
            return int(round(float(x)))
    else:
        def to_num(x):
            return round(float(x), digits=round_digits)
    fh = open(f)
    fh.seek(pos)
    headerline = fh.readline()
    line = fh.readline()
    while line:
        if line.startswith('#'):
            break   # done reading for this cell
        else:
            fields = line.rstrip('\n').split('\t')
            genes.append(fields[0])
            counts.append(to_num(fields[9]))
        line = fh.readline()
    fh.close()
    return (genes, np.array(counts))

if __name__ == "__main__":
    sys.exit(main())
