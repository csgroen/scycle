import anndata
import pickle as pkl
import zipfile
import os
import shutil
import re

def write(adata, fname, verbose = True):
    """ Write AnnData file with scycle analysis

    Parameters
    -------------
    adata: AnnData
        The Anndata object to be saved
    fname: String
        Path of the file to be saved. 'zip' extension will
        automatically be added.
    verbose: bool
        If True, the function will print messages.
    """
    adata_save = adata.copy()
    # Write analysis results
    fname_uns = f'{fname}_uns.pkl'
    _write(adata.uns, fname_uns)
    adata_save.uns = {}

    # Write AnnData
    adata_save.write(f'{fname}.h5ad')

    # Zip results
    if verbose: print('Saving AnnData and analyses...')
    scycle_files = [f'{fname}_uns.pkl', f'{fname}.h5ad']
    with zipfile.ZipFile(f'{fname}.zip', 'w') as zipF:
        for file in scycle_files:
            base_name = re.sub('^.*\/', '', file)
            zipF.write(file, compress_type=zipfile.ZIP_DEFLATED, arcname=base_name)
    os.remove(f'{fname}_uns.pkl'); os.remove(f'{fname}.h5ad')
    if verbose: print("Done")

def read(fname, verbose = True):
    """Read scycle analysis zip file

    Paramaters
    ---------------
    fname: string
        A string pointing to the path where the zipped file exported
        using cc.tl.write is stored.
    verbose: bool
        If True, the function will print messages
    """
    #-- Get file name variables
    fdir = fname.replace('.zip', '')
    base_name = re.sub('^.*\/', '', fdir)

    #-- Unzip
    with zipfile.ZipFile(fname, 'r') as zip_ref:
        zip_ref.extractall(fdir)

    #-- Read
    adata = anndata.read_h5ad(f'{fdir}/{base_name}.h5ad')
    uns = _read(f'{fdir}/{base_name}_uns.pkl')
    adata.uns = uns

    #-- Delete decompressed dir and contents
    shutil.rmtree(fdir)

    #-- Return
    return adata

def _write(data, fname):
    with open(fname, 'wb') as pickle_file:
        pkl.dump(data, pickle_file)

    return True

def _read(fname):
    with open(fname, 'rb') as pickle_file:
        obj = pkl.load(pickle_file)
    return obj
