import pandas as pd
from tabulate import tabulate

# Caminho do arquivo CSV
csv_file = "/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/Tabelas/ACOFLxFLGA-RG/comparacao_ACO_GA.csv"

# Lendo o arquivo CSV
df = pd.read_csv(csv_file)

# Criando a tabela formatada
tabela_formatada = tabulate(df, headers='keys', tablefmt='grid')

# Salvando a tabela em um arquivo de texto
output_file = "tabela.txt"
with open(output_file, "w") as f:
    f.write(tabela_formatada)

print(f"Tabela salva em {output_file}")
