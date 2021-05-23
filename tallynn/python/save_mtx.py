import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
import Levenshtein
import pandas as pd
import scipy.sparse as sparse
import scipy.io as io
import os


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("save_mtx.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--data", default=None, type=str,
                    help='counts data from umi_tools')
parser.add_argument("--dir", default=None, type=str,
                    help='dir for output mtx')
args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #
def save_mtx(data, destination, cell_names=None, gene_names=None):
    """Save a mtx file - taken from https://www.programcreek.com/python/?code=KrishnaswamyLab%2Fscprep%2Fscprep-master%2Fscprep%2Fio%2Fmtx.py
    
    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data, saved to destination/matrix.mtx
    destination : str
        Directory in which to save the data
    cell_names : list-like, shape=[n_samples], optional (default: None)
        Cell names associated with rows, saved to destination/cell_names.tsv.
        If `data` is a pandas DataFrame and `cell_names` is None,
        these are autopopulated from `data.index`.
    gene_names : list-like, shape=[n_features], optional (default: None)
        Cell names associated with rows, saved to destination/gene_names.tsv.
        If `data` is a pandas DataFrame and `gene_names` is None,
        these are autopopulated from `data.columns`.

    Examples
    --------
    >>> import scprep
    >>> scprep.io.save_mtx(data, destination="my_data")
    >>> reload = scprep.io.load_mtx("my_data/matrix.mtx",
    ...                             cell_names="my_data/cell_names.tsv",
    ...                             gene_names="my_data/gene_names.tsv")
    """
    if isinstance(data, pd.DataFrame):
        if cell_names is None:
            cell_names = data.index
        if gene_names is None:
            gene_names = data.columns
    data = sparse.coo_matrix(data)
    # handle ~/ and relative paths
    #print(cell_names)
    destination = os.path.expanduser(destination)
    if not os.path.isdir(destination):
        os.mkdir(destination)
    if cell_names is not None:
        with open(os.path.join(destination, "genes.barcodes.txt"), "w") as handle:
            for name in cell_names:
                handle.write("{}\n".format(name))
    if gene_names is not None:
        with open(os.path.join(destination, "genes.genes.txt"), "w") as handle:
            for name in gene_names:
                handle.write("{}\n".format(name))
    io.mmwrite(os.path.join(destination, "genes.mtx"), data) 


infile = pd.read_table(args.data, sep="\t", header=0)
infile = infile[infile['count'] > 0]

infile = infile.pivot(index='cell', columns='gene', values='count')
infile.fillna(0, inplace=True)

save_mtx(infile, args.dir)
