#acofl-x-flga-cpst

import pandas as pd
import os

# Lista dos arquivos CSV a serem unidos
arquivos = [
    "comparacao_ACO_GA_Cycles.csv",
    "comparacao_ACO_GA_Paths.csv",
    "comparacao_ACO_GA_Stars.csv",
    "comparacao_ACO_GA_Trees.csv"
]

# Lista para armazenar os DataFrames de cada arquivo
dataframes = []

# Ler cada arquivo CSV e adicionar ao DataFrame
for arquivo in arquivos:
    if os.path.exists(arquivo):  # Verifica se o arquivo existe
        df = pd.read_csv(arquivo)
        dataframes.append(df)
    else:
        print(f"Arquivo {arquivo} não encontrado. Verifique o nome e o caminho do arquivo.")

# Concatenar todos os DataFrames em um único DataFrame
df_final = pd.concat(dataframes, ignore_index=True)

# Salvar o DataFrame final em um novo arquivo CSV
df_final.to_csv("comparacao_ACO_GA_completo.csv", index=False)

print("Arquivos unidos e salvos em 'comparacao_ACO_GA_completo.csv'.")
