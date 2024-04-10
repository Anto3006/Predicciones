from calculateDescriptors import calculateDescriptors
import pickle
import pandas as pd
import os
import sys
import warnings

def mean(list_values):
  total = sum(list_values)
  if len(list_values) > 0:
    return total/len(list_values)
  else:
    return np.nan

def load_model(model_name):
    model = pickle.load(open("Models//"+model_name, 'rb'))
    return model

def predict(predict_type,descriptors):
    model_names = os.listdir(f"{predict_type}")
    predictions = {smiles:[] for smiles in descriptors["smiles"]}
    col_na = []
    for model_name in model_names:
      if not model_name.startswith('.'):
        model = load_model(predict_type,model_name)
        na_smiles = []
        non_na_smiles = []
        descriptors_model = list(model.feature_names_in_)
        descriptors_prediction = descriptors[["smiles"] + descriptors_model]
        if descriptors_prediction.isna().any().any():
          na_rows = descriptors_prediction[descriptors_prediction.isna().any(axis=1)]
          na_smiles = list(na_rows["smiles"])
          print(f"At least one value is missing for model {model_name} for molecules {na_smiles}")
          descriptors_prediction.dropna(axis=0,inplace=True)
        non_na_smiles = descriptors_prediction["smiles"]
        
        descriptors_prediction.drop(columns=["smiles"],inplace=True)
        descriptors_prediction.columns = descriptors_prediction.columns.astype(str)
        prediccion_modelo = model.predict(descriptors_prediction)
        for smiles,prediction in zip(non_na_smiles,prediccion_modelo):
          predictions[smiles].append(prediction)
    return [mean(predictions[smiles]) for smiles in descriptors["smiles"]]

def make_predictions(smiles):
    descriptors = calculateDescriptors(smiles)
    prediccion_logKoa = predict("logKoa",descriptors)
    prediccion_logP = predict("logP",descriptors)
    predictions = pd.DataFrame(data=[smiles,prediccion_logKoa,prediccion_logP],columns=["smiles","logKoa","logP"])
    return predictions


if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    mode = sys.argv[1]
    if mode == "-s":
        smiles = [sys.argv[i] for i in range(2,len(sys.argv))]
        print(make_predictions(smiles))
    elif mode == "-f":
        file_name = sys.argv[2]
        smiles = pd.read_csv(file_name)["smiles"]
        make_predictions(smiles).to_csv("predictions_"+file_name)