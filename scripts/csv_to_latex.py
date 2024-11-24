import pandas as pd

# Carrega o arquivo CSV
df = pd.read_csv('saida_unificada.csv')

# Converte para LaTeX
latex_table = df.to_latex(index=False, 
	column_format="|c|c|c|c|c|c|c|c|c|"
)

# Salva em um arquivo .tex
with open('tabela.tex', 'w') as f:
    f.write(latex_table)

