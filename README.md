# Meta-heuristic Algorithms for Triple Roman Domination Problem

In this project, we implement two algorithms based on meta-heuristics for solving the Triple Roman Domination Problem, an NP-Complete Combinatorial Optimization Problem, implementation of the content of my [monagraphy](https://drive.google.com/file/d/1vaPgLHMxWFpXXzfMNvM1Z-oMuRDs-Nt8/view). These meta-heuristics are based on Genetic Algorithms and Ant Colony Optimization.

# What is Triple Roman Domination Problem?

1. If $h(v) = 0$, then $v$ has a neighbor $u$ such that $h(u) = 4$; or $v$ has two neighbors $u_1$ and $u_2$ such that $h(u_1) = 2$ and $h(u_2) = 3$; or $v$ has three neighbors $u_1$, $u_2$, and $u_3$ such that $h(u_1) = h(u_2) = h(u_3) = 2$.

2. If $h(v) = 1$, then $v$ has at least one neighbor $u$ such that $h(u) \geq 3$; or two neighbors $u_1$ and $u_2$ such that $h(u_1) = h(u_2) = 2$.

3. If $h(v) = 2$, then $v$ has a neighbor $u$ such that $h(u) \geq 2$.

The weight of the triple Roman domination function is denoted by $h(V(G)) = \sum_{v \in V(G)} h(v)$. A TRDF is *optimal* if it has the smallest possible weight. The *triple Roman domination number* of a graph $G$ is the smallest weight of a TRDF of $G$ and is denoted by $\gamma_{3R}(G)$. The *Triple Roman Domination Problem* (TRDP) consists of determining $\gamma_{3R}(G)$ for an arbitrary graph $G$.

The decision version of the TRDP consists of, given a graph $G$ and an integer $\ell$, deciding if $G$ has a TRDF with weight less than or equal to $\ell$. This decision problem is NP-Complete.

# What is a Optimization Problem?

According to *Papadimitriou e Steiglitz (1982)*, an *optimization problem (OP)* consists of three core components: a set of *instances* (inputs), a set of *feasible solutions* $Sol(I)$ for each instance $I$, and an *objective function* $val(S)$, which assigns a numerical cost to each feasible solution $S \in Sol(I)$. If the feasible solution set $Sol(I)$ for an instance $I$ is empty, $I$ is deemed *infeasible*; otherwise, it is *feasible*. The goal of an OP is to find the best solution among all feasible solutions. These problems are generally classified as *minimization* or *maximization* problems, focusing on the solutions with the minimum or maximum cost, respectively. A feasible solution with optimal cost is called an *optimal solution*, and its value for any instance $I$ is denoted as $OPT(I)$, where $OPT(I) = val(S^\*)$ for an optimal solution $S^*$. Note that $OPT(I)$ is defined only if $I$ is feasible.

Optimization problems are often computationally challenging. Complexity theory classifies these problems as follows:

- **P**: Problems solvable in polynomial time.
- **NP**: Decision problems (YES/NO answers) whose solutions are verifiable in polynomial time.
- **NP-hard**: Problems to which every problem in NP can be reduced in polynomial time, indicating high computational difficulty.
- **NP-complete**: Problems that are both in NP and NP-hard, representing many real-world problems without known polynomial-time solutions.

# What is a Meta-Heuristic?
A meta-heuristic is a general exploration method, often stochastic, that applies similarly to various different combinatorial optimization problems. It typically begins with any feasible solution and iteratively improves it over several iterations until a stopping criterion is met. Examples of meta-heuristics include Memetic Genetic Algorithms and Ant Colony Optimization.


## Genetic Algorithm

The stochastic metaheuristic proposed by Holland (1992a), known as the genetic algorithm (HOLLAND, 1992b), aims to provide an optimal or near-optimal solution, inspired by Charles Darwin's evolutionary theory (DARWIN, 2017). The central idea is that, in a population, only the individuals most adapted to the environment will be able to pass on their genes across generations through successive reproduction processes. This process is subject to gene mutation rates and elitism, which refers to the selection of individuals with genetic traits better suited to the environment for reproduction. This results from environmental changes, ensuring the survival of individuals best adapted to them.

The GA proposed in my monography represents solutions, or chromosomes, as vectors of size $|V|$, where each position represents a specific vertex and its label, according to the constraints of the Triple Roman Domination Function (TRDF). The initial population is generated using three heuristics (the fourth heuristic is a mix of all others heuristics). Each heuristic utilizes an auxiliary graph $G$, a copy of the original graph. During each iteration, vertices are removed from $G$, with the stopping condition being $V(G) = \emptyset$.

### Heuristics for Initial Population

1. **First Heuristic**: A greedy approach is applied. In each iteration, a vertex $v$ is randomly selected from the auxiliary graph and assigned the label 2 (3 if it has degree 0). The adjacent vertices are labeled 0. If remains isolated vertices in the auxiliary graph, it is labeled 3; after this, feasibility check routine is called for fix eventuals errors in TRDF, ensuring that all vertices adhere its constraints.

2. **Second Heuristic**:  A greedy approach is applied. In each iteration, a vertex $v$ is randomly selected from the auxiliary graph and assigned the label 4 (3 if it has degree 0). The adjacent vertices are labeled 0. If remains isolated vertices in the auxiliary graph, it is labeled 3; after this, toggle vertices routine is called to try minimize the labels assigned in TRDF, trying improve the *triple roman domination number*.

3. **Third Heuristic**: Vertices are sorted in descending order by degree. The highest-degree vertexs are selected and labeled 4 (3 if it has degree 0), while its neighbors in $N(v)$ are labeled 0. If any remaining vertices in the auxiliary graph have degree 0, they are labeled 3.

4. **Fourth heuristic**: Create the population using all the above heuristics

### One Point Mutation Operator

The `mutation` function introduces random variations in a random gene of the chromosome, allowing for greater exploration of the solution space and preventing premature convergence. This function operates as follows:

1. **Random Number Generation**: A random seed generator (`std::random_device`) and a pseudo-random number generator (`std::mt19937`) are initialized, along with two distributions: `gap`, which generates a probability between 0 and 1, and `randomLabel`, which selects an integer between 0 and 4 to represent a possible gene value.

2. **Mutation Process**: For each gene in the chromosome:
   - A random probability is generated using getTandomInt method.
   - If this probability is less than or equal to the mutation rate (`mutationRate`), the gene is replaced with a randomly selected value from `randomLabel`, representing a new gene label between 0 and 4.
   
3. **Feasibility Check**: After each mutation, the `feasibilityCheck` function is called to ensure that the chromosome maintains valid properties, as required by the problem's constraints.

This mutation operator, therefore, enhances genetic diversity within the population, helping the algorithm to explore new solutions by occasionally altering gene values.

### Elitism Selection Operator

The `elitism` function implements an elitism strategy, ensuring that the best solution in the current population is preserved across generations. This approach prevents the loss of the most optimal solutions and is outlined as follows:

1. **Identifying the Best Solution**: The function identifies the best chromosome in the population using `getBestChromosome`, which helps maintain high-quality solutions.

2. **Elitism Rate and Iterations**: The number of elite individuals to retain is determined by `elitismRate`. Specifically, this rate is multiplied by the population size to calculate the number of copies of the best solution to retain in the next generation.

3. **Population Replacement**: The selected elite chromosomes are copied into a temporary vector, `temp`, and then swapped with the original population, ensuring that the best solutions are propagated forward.

This elitism operator helps maintain strong solutions across generations, balancing exploration and exploitation by retaining top performers, thus enhancing convergence towards an optimal solution.

### One Point Cross Over Operator

The crossover operator combines two solutions, $S_1$ and $S_2$, from the current population of $n$ chromosomes. This is achieved through the following steps:

1. **Selection of Crossover Point**: Random indice, $R_1$, is chosen within the chromosome length.
2. **Gene Exchange**: The labels (genes) after $R_1$ index in $S_1$ and $S_2$ are swapped, producing two new offspring solutions.
3. **Feasibility Check**: The resulting offspring solutions are verified to ensure they satisfy the constraints of the Triple Roman Domination Function (TRDF). Solutions that violate constraints are discarded or adjusted.

In this proposed algorithm, elitism and mutation rate do not influence the crossover operator directly.


### Two Point Cross Over Operator

The crossover operator combines two solutions, $S_1$ and $S_2$, from the current population of 1000 chromosomes. This is achieved through the following steps:

1. **Selection of Crossover Points**: Random indices, $R_1$ and $R_2$, are chosen within the chromosome length.
2. **Gene Exchange**: The labels (genes) between indices $R_1$ and $R_2$ in $S_1$ and $S_2$ are swapped, producing two new offspring solutions.
3. **Feasibility Check**: The resulting offspring solutions are verified to ensure they satisfy the constraints of the Triple Roman Domination Function (TRDF). Solutions that violate constraints are discarded or adjusted.

In this proposed algorithm, elitism and mutation rate do not influence the crossover operator directly.

## Ant Colony Optimization (ACO)

The Ant Colony Optimization algorithm simulates the behavior of ants searching for food, where ants deposit pheromones along their paths. The best path is the one with the highest pheromone concentration, guiding future ants. In this algorithm, each vertex in the graph is associated with an initial pheromone value of 0.5, and each solution is represented as a vector of size $|V|$. The weight of the solution is the sum of the labels in the solution vector, with the goal of finding a solution with the smallest possible weight. The algorithm runs until a maximum number of iterations is reached or there are no improvements in the solution for a set number of iterations.

The type of Ant Colony Optimization proposed is based on two variants: **MAX-MIN Ant System** (MMAS) and **Hyper-Cube Framework** (HCF).

### ACO Sub-routines

1. **ConstructSolution**: A random vertex is selected and labeled as 4, while its neighbors are labeled as 0. These vertices are then removed from the auxiliary graph. If there are any remaining vertices with degree 0, they are labeled as 2. After this, the algorithm checks if there are any vertices labeled 2 that have only neighbors labeled as 0; if so, a random neighbor is assigned the label 2 to ensure the constraints of the TRDF.

2. **ExtendSolution**: A vertex currently labeled 0, 2, or 3 is selected and upgraded to a label of 4 in the current solution.

3. **ReduceSolution**: Redundant vertices are removed from the current solution. A vertex $v$ is considered redundant if all vertices in its closed neighborhood $N[v]$ are already dominated by other vertices in the solution.

4. **RVNS (Reduced Variable Neighborhood Search)**: This function enhances the best result obtained from the ACO algorithm by exploring nearby solution spaces. RVNS increases the likelihood of finding a solution closer to the optimal one by modifying the solution in the following steps:
    - The **DestroySolution** sub-routine randomly selects a vertex labeled 0, 2, or 3 and also unlabeled it, setting its value to -1.
    - **ConstructSolution**, **ExtendSolution**, and **ReduceSolution** are then reapplied to form a new solution, which is compared with the previous one. If the new solution has a lower weight, it replaces the current solution.

## How to Use

1. Clone the repository
   ```bash
   git clone https://github.com/isrreal/Double-Roman-Domination-In-Graphs-meta-heuristics.git
   
2. Run the Code

To execute the application, use the following syntax:

```bash
# Syntax:
# ./app <graph_file> <graph_name> <heuristic 1-4>

# Parameters:
#   population_size               - Defines the size of the initial population. A larger size increases solution diversity but may increase computation time.
#   generations                   - Number of generations (iterations) the Genetic Algorithm will run to evolve solutions.
#   chromosome_creation_heuristic - Heuristic used for initializing chromosomes:
#                                     1 - heuristic 1
#                                     2 - heuristic 2
#                                     3 - heuristic 3
#                                     4 - heuristic 4

# Example:
./app graph.txt graph_name heuristic
```

