import pyomo.environ as pyo
import networkx as nx
from pathlib import Path
import time
import csv

# Caminho da pasta
input_dir = Path("input_files")
output_dir = Path("output_files")

# Criar a pasta de saída se não existir
output_dir.mkdir(exist_ok=True)

# Tempo máximo de execução em segundos
TIME_LIMIT = 900  # 900s = 15 minutos 

# Percorre todas as subpastas e arquivos dentro de 'input_files'
for subfolder in input_dir.iterdir():
    if subfolder.is_dir(): # verifica se é pasta

        # Criar diretório da subpasta dentro de 'output_files'
        subfolder_output_dir = output_dir / subfolder.name
        subfolder_output_dir.mkdir(exist_ok=True, parents=True)

        output_file = subfolder_output_dir / "resultados.csv"

        # Criar arquivo CSV e escrever o cabeçalho (se ainda não existir)
        if not output_file.exists():
            with open(output_file, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["filename", "vertex", "edge", "density", "TRD_number", "time(seconds)", "status"])


        for arquivo in subfolder.iterdir():
            if arquivo.is_file():
                # Ler o grafo do arquivo de lista de arestas
                G = nx.read_edgelist(arquivo, nodetype=int)
                
                # Remover arestas de loop
                G.remove_edges_from([(u, v) for u, v in G.edges() if u == v])
                
                edge_count = len(G.edges)
                vertex_count = len(G.nodes)
                density = nx.density(G)
                
                # Criar modelo Pyomo
                model = pyo.ConcreteModel()
                
                # Definir variáveis binárias
                nodes = list(G.nodes())
                model.q = pyo.Var(nodes, within=pyo.Binary)
                model.r = pyo.Var(nodes, within=pyo.Binary)
                model.s = pyo.Var(nodes, within=pyo.Binary)
                
                # Restrições
                def domination_rule(model, i):
                    sum_labels = 2 * model.q[i] + 3 * model.r[i] + 4 * model.s[i]
                    for neighbor in G.neighbors(i):
                        sum_labels += model.q[neighbor] + 2 * model.r[neighbor] + 3 * model.s[neighbor]
                    return sum_labels >= 3
                
                model.domination_constraints = pyo.Constraint(nodes, rule=domination_rule)
                
                def single_assignment_rule(model, i):
                    return model.q[i] + model.r[i] + model.s[i] <= 1
                
                model.single_assignment_constraints = pyo.Constraint(nodes, rule=single_assignment_rule)
                
                # Função objetivo
                model.objective = pyo.Objective(
                    expr=sum(2 * model.q[i] + 3 * model.r[i] + 4 * model.s[i] for i in nodes),
                    sense=pyo.minimize
                )
                
                try:
                    # Resolver modelo com limite de tempo
                    if(len(G.nodes()) <= 333):
                        solver = pyo.SolverFactory("cplex")
                        solver.options['timelimit'] = TIME_LIMIT
                    else:
                        solver = pyo.SolverFactory("cbc")
                        solver.options["seconds"] = TIME_LIMIT

                    
                    start_time = time.time()
                    results = solver.solve(model, tee=False)
                    end_time = time.time()
                    elapsed_time = (end_time - start_time)  # em segundos
                    
                    # Verificar status
                    status = results.solver.termination_condition
                    if status == pyo.TerminationCondition.optimal:
                        best_value = pyo.value(model.objective)
                        result_line = [arquivo.name, vertex_count, edge_count, f"{density:.5f}", f"{best_value:.2f}", f"{elapsed_time:.2f}", "Optimal"]
                    elif status in {pyo.TerminationCondition.feasible, pyo.TerminationCondition.maxTimeLimit}:
                        best_value = pyo.value(model.objective)
                        result_line = [arquivo.name, vertex_count, edge_count, f"{density:.5f}", f"{best_value:.2f}", f"{elapsed_time:.2f}", "Best Found"]
                    else:
                        result_line = [arquivo.name, vertex_count, edge_count, f"{density:.5f}", "No feasible solution", f"{elapsed_time:.2f}", "Infeasible"]
                except Exception as e:
                    result_line = [arquivo.name, "Error", "Error", "Error", "Error", "Error", f"erro:: {e}"]
                
                # Escrever no arquivo CSV
                with open(output_file, "a", newline="", encoding="utf-8") as f:
                    writer = csv.writer(f)
                    writer.writerow(result_line)