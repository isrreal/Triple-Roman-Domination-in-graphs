import os

def calculate_graph_orders(directory):
    """
    Lê grafos de arquivos .txt em um diretório, calcula a ordem de cada grafo e exibe o resultado ordenado.

    :param directory: Caminho para o diretório contendo os arquivos .txt.
    """
    results = []
    for file_name in os.listdir(directory):
        if file_name.endswith('.txt'):
            file_path = os.path.join(directory, file_name)
            with open(file_path, 'r') as file:
                vertices = set()
                for line in file:
                    u, v = map(int, line.split())  # Lê as arestas
                    vertices.add(u)
                    vertices.add(v)
                order = len(vertices)  # Calcula o número de vértices
                results.append((file_name, order))
    
    # Ordena os resultados pelo número de vértices (ordem do grafo)
    results.sort(key=lambda x: x[1])

    # Exibe os resultados ordenados
    for file_name, order in results:
        print(f"{file_name} ----> {order}")

# Diretório contendo os arquivos extraídos
directory_path = "/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/input_files/Harwell-Boeing"  # Substitua pelo caminho correto

# Calcular e exibir as ordens dos grafos
calculate_graph_orders(directory_path)

