from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
import pandas as pd


def calculateDescriptorsMolRDKit(molecule):
    descriptors = Descriptors.CalcMolDescriptors(molecule)
    descriptors["NumAmideBonds"] = rdMolDescriptors.CalcNumAmideBonds(molecule)
    descriptors["NumSpiroAtoms"] = rdMolDescriptors.CalcNumSpiroAtoms(molecule)
    descriptors["NumBridgeheadAtoms"] = rdMolDescriptors.CalcNumBridgeheadAtoms(molecule)
    PEOE_VSA = rdMolDescriptors.PEOE_VSA_(molecule)
    for i in range(1,15):
        descriptors["PEOE_VSA_"+str(i)] = PEOE_VSA[i-1]

    SMR_VSA = rdMolDescriptors.SMR_VSA_(molecule)
    for i in range(1,11):
        descriptors["SMR_VSA_"+str(i)] = SMR_VSA[i-1]

    SlogP_VSA = rdMolDescriptors.SlogP_VSA_(molecule)
    for i in range(1,11):
        descriptors["SlogP_VSA_"+str(i)] = SlogP_VSA[i-1]
    
    
    MQNs = rdMolDescriptors.MQNs_(molecule)
    for i in range(1,43):
        descriptors["MQNs_"+str(i)] = MQNs[i-1]
    return descriptors
    


def calculateDescriptorsRDKit(smiles):
    descriptorsRDKit = {}
    for smile in smiles:
        try:
            molecule = Chem.MolFromSmiles(smile)
            descriptors = calculateDescriptorsMolRDKit(molecule)
            for descriptor in descriptors:
                if (descriptor+"_rdkit") in descriptorsRDKit:
                    descriptorsRDKit[descriptor+"_rdkit"].append(descriptors[descriptor])
                else:
                    descriptorsRDKit[descriptor+"_rdkit"] = [descriptors[descriptor]]
        except:
            for descriptor in descriptorsRDKit:
                descriptorsRDKit[descriptor+"_rdkit"].append("")
    return pd.DataFrame(descriptorsRDKit)