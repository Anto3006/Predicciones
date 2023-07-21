import pandas
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
from descriptoresRDKit import calcularDescriptoresRDKit
from descriptoresJazzy import calcularDescriptoresJazzy


pandas2ri.activate()


def calcularDescriptoresCDK(smiles):
    rcdk = rpackages.importr("rcdk")
    base = rpackages.importr("base")
    descCategories = rcdk.get_desc_categories()
    getDescNames = rcdk.get_desc_names

    descNames = base.unique(base.unlist(base.sapply(descCategories,getDescNames)))

    mols = rcdk.parse_smiles(ro.StrVector(smiles))
    descriptorsCDK = rcdk.eval_desc(mols,descNames)

    descriptorsCDK = ro.conversion.rpy2py(descriptorsCDK)
    descriptorsCDK.index = [i for i in range(len(descriptorsCDK.index))]

    return descriptorsCDK




def calcularDescriptoresObabel(smiles):
    rcpi = rpackages.importr('Rcpi')
    descriptorsObabel = pandas.DataFrame()
    for smile in smiles:
        x = rcpi.extractDrugDescOB(smile, type = 'smile')
        y = ro.conversion.rpy2py(x)
        descriptorsObabel = pandas.concat([descriptorsObabel,y],ignore_index=True)
    descriptorsObabel.drop(columns=["cansmi","cansmiNS","formula","title","InChI"],inplace=True)
    return descriptorsObabel

def calcularDescriptores(smiles):
    descriptorsCDK = calcularDescriptoresCDK(smiles)
    descriptorsObabel = calcularDescriptoresObabel(smiles)
    descriptoresRDKit = calcularDescriptoresRDKit(smiles)
    descriptoresJazzy = calcularDescriptoresJazzy(smiles)
    descriptors = pandas.concat([descriptoresJazzy,descriptorsCDK,descriptoresRDKit,descriptorsObabel],axis=1)
    return descriptors


