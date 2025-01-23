from networkx.generators.random_graphs import erdos_renyi_graph

order = (25, 50, 75, 100, 125, 150, 175, 200, 225, 250)
edge_probabilities = (0.2, 0.5, 0.8)

graphs = []
graph_counter = 1  # Contador para nomear os arquivos

for n in order:
    for p in edge_probabilities:
        # Gera o grafo
        graph = erdos_renyi_graph(n, p)
        graphs.append(graph)
        
        # Salva as arestas em um arquivo
        filename = f"random_graph{graph_counter}-order{n}-edge_probability-{p}.txt"
        with open(filename, "w") as file:
            for u, v in graph.edges():
                file.write(f"{u} {v}\n")
        
        graph_counter += 1











