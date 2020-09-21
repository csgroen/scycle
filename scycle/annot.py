import pickle
def cellcycle_signatures():
    with open('data/cellcycle_signatures.pkl', 'rb') as f:
        cc_sigs = pickle.load(f)
    return cc_sigs
