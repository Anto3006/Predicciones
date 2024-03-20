from calcularDescriptores import calcularDescriptores
import pickle
import pandas as pd
import os
import sys
import warnings

def cargarModelo(nombre_modelo):
    modelo = pickle.load(open("Modelos//"+nombre_modelo, 'rb'))
    return modelo

def predecir(tipo,descriptores):
    nombre_modelos = os.listdir(f"Modelos//{tipo}")
    prediccion = 0
    for nombre_modelo in nombre_modelos:
        modelo = cargarModelo(nombre_modelo)
        descriptores_modelo = modelo.feature_names_in_
        descriptores = descriptores[descriptores_modelo]
        descriptores.columns = descriptores.columns.astype(str)
        prediccion_modelo = modelo.predict(descriptores)
        prediccion += prediccion_modelo
    return prediccion/len(nombre_modelos)

def realizarPredicciones(smiles):
    descriptores = calcularDescriptores(smiles)
    prediccion_logKoa = predecir("logKoa",descriptores)
    prediccion_logP = predecir("logP",descriptores)
    predicciones = pd.DataFrame(data=[smiles,prediccion_logKoa,prediccion_logP],columns=["smiles","logKoa","logP"])
    return predicciones


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