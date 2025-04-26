import os
import shutil
import networkx as nx

# Caminho da pasta principal onde estão as subpastas com os arquivos de grafos
pasta_principal = "."

# Pasta onde serão salvos os grafos filtrados, mantendo a estrutura original
pasta_destino = "./atmost333/"

# Percorre todas as subpastas dentro da pasta principal
for subpasta in os.listdir(pasta_principal):
    caminho_subpasta = os.path.join(pasta_principal, subpasta)
    
    if os.path.isdir(caminho_subpasta):  # Verifica se é uma pasta
        # Criar a subpasta correspondente na pasta de destino
        pasta_destino_subpasta = os.path.join(pasta_destino, subpasta)
        os.makedirs(pasta_destino_subpasta, exist_ok=True)

        # Percorre os arquivos dentro da subpasta
        for arquivo in os.listdir(caminho_subpasta):
            if arquivo.endswith(".txt"):  # Apenas arquivos de grafos
                caminho_arquivo = os.path.join(caminho_subpasta, arquivo)

                try:
                    # Ler o grafo usando networkx
                    G = nx.read_edgelist(caminho_arquivo, nodetype=int)

                    # Verificar a ordem do grafo
                    if len(G.nodes()) <= 333:
                        shutil.copy(caminho_arquivo, os.path.join(pasta_destino_subpasta, arquivo))
                        print(f"Copiado para {pasta_destino_subpasta}: {arquivo}")

                except Exception as e:
                    print(f"Erro ao processar {arquivo}: {e}")

print("Processamento concluído.")

