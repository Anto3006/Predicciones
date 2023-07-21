from calcularDescriptores import calcularDescriptores
import pickle
import pandas as pd
import os

def cargarModelo(nombreModelo):
    modelo = pickle.load(open("Modelos//"+nombreModelo, 'rb'))
    return modelo

def predecir(nombreModelo,descriptores):
    modelo = cargarModelo(nombreModelo)
    descriptoresModelo = modelo.feature_names_in_
    descriptores = descriptores[descriptoresModelo]
    descriptores.columns = descriptores.columns.astype(str)
    descriptores = descriptores.loc[:, ~descriptores.columns.duplicated()]
    for i in range(len(descriptoresModelo)):
        print(descriptoresModelo[i],descriptores.columns[i])
    predicciones = modelo.predict(descriptores)
    return predicciones

def realizarPredicciones(smiles):
    nombreModelos = os.listdir("Modelos")
    print(nombreModelos)
    predicciones = {}
    descriptores = calcularDescriptores(smiles)
    for nombreModelo in nombreModelos:
        p = predecir(nombreModelo,descriptores)
        predicciones[nombreModelo] = p
    return pd.DataFrame(predicciones,index=smiles)

archivoSmiles = input("ingrese direccion del archivo csv con los datos\n")
smiles = pd.read_csv(archivoSmiles)["smiles"]
realizarPredicciones(smiles).to_csv("predicciones_"+archivoSmiles)