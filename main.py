from calcularDescriptores import calcularDescriptores
import pickle
import pandas as pd
import os
import sys
import warnings

def cargarModelo(nombreModelo):
    modelo = pickle.load(open("Modelos//"+nombreModelo, 'rb'))
    return modelo

def predecir(nombreModelo,descriptores):
    modelo = cargarModelo(nombreModelo)
    descriptoresModelo = modelo.feature_names_in_
    descriptores = descriptores[descriptoresModelo]
    descriptores.columns = descriptores.columns.astype(str)
    descriptores = descriptores.loc[:, ~descriptores.columns.duplicated()]
    predicciones = modelo.predict(descriptores)
    return predicciones

def realizarPredicciones(smiles):
    nombreModelos = os.listdir("Modelos")
    predicciones = {}
    descriptores = calcularDescriptores(smiles)
    for nombreModelo in nombreModelos:
        p = predecir(nombreModelo,descriptores)
        predicciones[nombreModelo] = p
    return pd.DataFrame(predicciones,index=smiles)


if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    modo = sys.argv[1]
    if modo == "-s":
        smiles = [sys.argv[i] for i in range(2,len(sys.argv))]
        print(realizarPredicciones(smiles))
    elif modo == "-f":
        nombreArchivo = sys.argv[2]
        smiles = pd.read_csv(nombreArchivo)["smiles"]
        realizarPredicciones(smiles).to_csv("predicciones_"+nombreArchivo)