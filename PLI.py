import pyomo.environ as pyo
import networkx as nx
from pathlib import Path
import time

# Caminho da pasta
pasta = Path("input_files")

# Tempo máximo de execução em segundos
TIME_LIMIT = 900  # Defina o limite desejado. 900s = 15 minutos 

print("filename,vertex,edge,density,TRD_number,time(miliseconds)")

# Listar todos os arquivos na pasta
for arquivo in pasta.iterdir():
    if arquivo.is_file():
        # Ler o grafo do arquivo de lista de arestas
        G = nx.read_edgelist(arquivo, nodetype=int)
        
        # Remover arestas de loop
        G.remove_edges_from([(u, v) for u, v in G.edges() if u == v])
        
        edge_count = len(G.edges())
        vertex_count = len(G.nodes())
        density = (2.0 * edge_count) / (vertex_count * (vertex_count - 1))
        
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
        
        # Resolver modelo com limite de tempo
        solver = pyo.SolverFactory("cplex")  # usando o solver comercial CPLEX
        start_time = time.time()
        solver.options['timelimit'] = TIME_LIMIT
        results = solver.solve(model, tee=False)
        end_time = time.time()
        elapsed_time = (end_time - start_time) * 1000  # Convertendo para milissegundos
        
        # Verificar status
        if results.solver.termination_condition == pyo.TerminationCondition.optimal:
            optimal_value = pyo.value(model.objective)
            print(f"{arquivo},{vertex_count},{edge_count},{density:.5f},{optimal_value},{elapsed_time:.2f}")
        else:
            print(f"{arquivo},{vertex_count},{edge_count},{density:.5f},No optimal solution found,{elapsed_time:.2f}")
