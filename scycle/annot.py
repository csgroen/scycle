import pkg_resources
import pickle

DATA_PATH = pkg_resources.resource_filename("scycle", "data/")


def cellcycle_signatures():
    with open(DATA_PATH + "cellcycle_signatures.pkl", "rb") as f:
        cc_sigs = pickle.load(f)
    return cc_sigs


def pathway_annotation(path_name):
    # -- get source
    if path_name == "REACTOME":
        fpath = "reactome_v7_2.pkl"
    elif path_name == "BIOCARTA":
        fpath = "biocarta_v7_2.pkl"
    elif path_name == "KEGG":
        fpath = "kegg_v7_2.pkl"
    # -- return
    with open(DATA_PATH + fpath, "rb") as f:
        path_annot = pickle.load(f)
    return path_annot


# =============================================================================
# def _read_gmt(fname, cname):
#     with open(DATA_PATH + fname) as file:
#         paths = {}
#         for line in file:
#           path = line
#           path_info = path.split('\t')
#           path_info[len(path_info)-1] = path_info[len(path_info)-1].replace('\n', '')
#           path_name = path_info[0].replace(cname+'_', '')
#           paths[path_name] = path_info[2:len(path_info)]
#     return paths
# =============================================================================
