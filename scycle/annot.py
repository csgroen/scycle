import pkg_resources
import pickle

DATA_PATH = pkg_resources.resource_filename('scycle', 'data/')

def cellcycle_signatures():
    with open(DATA_PATH + 'cellcycle_signatures.pkl', 'rb') as f:
        cc_sigs = pickle.load(f)
    return cc_sigs
