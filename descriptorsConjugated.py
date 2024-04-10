from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
import subprocess
import os
import pandas as pd


def smiles2sdf(smiles, output_file):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    writer = Chem.SDWriter(output_file)
    writer.write(mol)
    writer.close()



def calculateConjugatedDescriptors(smiles_list):
    for i, smiles in enumerate(smiles_list):
        output_file = f'molecule_{i+1}.sdf'
        smiles2sdf(smiles, output_file)
    subprocess.call("Rscript conjugaR.R", shell=True)
    files = os.listdir()
    sdf_files = [file for file in files if file.split(".")[-1] == "sdf"]
    csv_files = [file for file in files if file.split(".")[-1] == "csv"]
    csv_files.remove("cjsystems.csv")
    descriptors_file_path = 'conjugated_descriptors.csv'
    descriptors = pd.read_csv(descriptors_file_path)
    descriptors.drop(columns=[descriptors.columns[0]],axis=1,inplace=True)
    for file in sdf_files:
        os.remove(file)
    for i in range(len(smiles_list)):
        file = f"molecule_{i+1}_conjugate_counts.csv"
        os.remove(file)
    os.remove(descriptors_file_path)
    return descriptors
