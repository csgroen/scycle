#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import requests
import os
import anndata
import re

#-- Data locations
# sc200_CCLE ref = SCP542
datasets = ['CHLA9, sc200_CCLE']
url_chla9 = 'https://xfer.curie.fr/get/nil/zXRLaZtAPi8/CHLA9.loom'
url_sc200 = 'http://xfer.curie.fr/get/nil/8yt0mkg9OEF/sc200CL_pp.h5ad'

#-- Main function
def get_data(dataset):
    """Download scycle example dataset

    Parameters
    ----------
    dataset: str
        dataset is a string with the name of the dataset to be downloaded.
        Must be one of: 'CHLA9' or 'sc200_CCLE'
    """
    #-- Get cache location
    cache_dir = os.path.dirname(os.path.realpath(__file__))

    #-- Check if cached, download otherwise
    #------ CHLA9
    if dataset == 'CHLA9':
        fname = cache_dir + '/chla9.loom'
        if 'chla9.loom' not in os.listdir(cache_dir):
            print('-- Downloading CHLA9 data from Xfer...')
            _download_scdata(url_chla9, fname)
            print('-- Download concluded.')
    elif dataset == 'sc200_CCLE':
        fname = cache_dir + '/sc200_ccle.h5ad'
        if 'sc200_ccle.h5ad' not in os.listdir(cache_dir):
            print('-- Downloading sc200_ccle data from Xfer...')
            _download_scdata(url_sc200, fname)
            print('-- Download concluded.')
    else:
        print("Dataset not in list of supported datasets. Must be one of:" + datasets)
        return None

    #-- Load from cache
    print('-- Loading data from cache...')
    if len(re.findall('loom$', fname)) > 0:
        scdata = anndata.read_loom(fname)
        scdata.var_names_make_unique()
        print('Done.')
    elif len(re.findall('h5ad$', fname)) > 0:
        scdata = anndata.read_h5ad(fname)




    return scdata

#-- Download data
def _download_scdata(url, fname):
    r = requests.get(url, stream = True)
    f = open(fname, 'wb')
    for chunk in r.iter_content(chunk_size=1024):
        f.write(chunk)
