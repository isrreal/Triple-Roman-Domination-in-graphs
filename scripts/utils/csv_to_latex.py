import pandas as pd

# Carrega o arquivo CSV
df = pd.read_csv("/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/Tabelas/Tabelas Restantes/Tabela MN/comparacao_ACO_GA_MN.csv")

# Converte para LaTeX
latex_table = df.to_latex(index=False, 
	column_format="|l|r|r|r|r|r|r|r|r|r|r|"
)

# Salva em um arquivo .tex
with open('tabela.tex', 'w') as f:
    f.write(latex_table)

