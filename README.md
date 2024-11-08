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


## Memetic Genetic Algorithm

A **Memetic Genetic Algorithm (MGA)** is an evolutionary algorithm that combines the principles of genetic algorithms (GAs) with local search techniques to refine solutions. In MGAs, after the crossover and mutation processes, each new chromosome (solution) undergoes a local improvement phase (or "meme"). This additional step helps to explore the search space more effectively and improve convergence towards optimal solutions, making MGAs particularly useful for complex optimization problems. The local search, or meme, can vary based on the problem context and may include heuristics, gradient-based optimization, or other methods to improve solution quality. MGAs are widely used in combinatorial optimization, including problems in NP-hard domains.

The MGA proposed by the authors represents solutions, or chromosomes, as vectors of size $|V|$, where each position represents a specific vertex and its label, according to the constraints of the Triple Roman Domination Function (TRDF). The initial population is generated using three heuristics, producing 1000 chromosomes. Each heuristic utilizes an auxiliary graph $G$, a copy of the original graph. During each iteration, vertices are removed from $G$, with the stopping condition being $V(G) = \emptyset$.

### Heuristics for Initial Population

1. **First Heuristic**: A greedy approach is applied. In each iteration, a vertex $v$ is randomly selected from the auxiliary graph and assigned the label 4. The corresponding position in the solution vector is labeled 4, while adjacent vertices are labeled 0. If only one vertex remains in the auxiliary graph, it is labeled 3, ensuring that all vertices adhere to the TRDF constraints.

2. **Second Heuristic**: Similar to the first, but after labeling vertices and removing them from the auxiliary graph, any isolated vertices (degree 0) are labeled 2. These isolated vertices' positions in the solution vector are updated, and the vertices are then removed from $G$. If a vertex labeled 2 has only neighbors labeled 0, one of these neighbors is randomly assigned a label of 2, maintaining connectivity under TRDF requirements.

3. **Third Heuristic**: Vertices are sorted in descending order by degree. The highest-degree vertex is selected and labeled 4, while its neighbors in $N(v)$ are labeled 0. If any remaining vertices in the auxiliary graph have degree 0, they are labeled 2. The auxiliary graph is then checked for any vertices labeled 2 that have only neighbors labeled 0; if such vertices are found, a random neighbor is assigned the label 2.

### Memetic Local Search

After the crossover operation, each offspring undergoes a local search process to improve solution quality. The local search involves:

- **ExtendSolution**: Any vertex currently labeled 0, 2, or 3 may be upgraded to a label of 4 to strengthen domination within the solution.
- **ReduceSolution**: Redundant vertices are identified and removed from the solution, preserving TRDF constraints. A vertex $v$ is redundant if all vertices in its closed neighborhood $N[v]$ are dominated by other vertices in the solution.

By using these local refinement techniques, the algorithm effectively balances exploration with focused improvements, increasing the chances of finding optimal or near-optimal solutions for the TRDF problem.

This approach leverages both the global exploration of genetic algorithms and the precision of local search, making it highly effective for Triple Roman Domination problems in conex graph structures.

### Mutation Operator

The `mutation` function introduces random variations in the genes of a chromosome, allowing for greater exploration of the solution space and preventing premature convergence. This function operates as follows:

1. **Random Number Generation**: A random seed generator (`std::random_device`) and a pseudo-random number generator (`std::mt19937`) are initialized, along with two distributions: `gap`, which generates a probability between 0 and 1, and `randomLabel`, which selects an integer between 0 and 4 to represent a possible gene value.

2. **Mutation Process**: For each gene in the chromosome:
   - A random probability is generated using `gap`.
   - If this probability is less than or equal to the mutation rate (`mutationRate`), the gene is replaced with a randomly selected value from `randomLabel`, representing a new gene label between 0 and 4.
   
3. **Feasibility Check**: After each mutation, the `feasibilityCheck` function is called to ensure that the chromosome maintains valid properties, as required by the problem's constraints.

This mutation operator, therefore, enhances genetic diversity within the population, helping the algorithm to explore new solutions by occasionally altering gene values.


### Elitism Selection Operator

The `elitism` function implements an elitism strategy, ensuring that the best solution in the current population is preserved across generations. This approach prevents the loss of the most optimal solutions and is outlined as follows:

1. **Identifying the Best Solution**: The function identifies the best chromosome in the population using `getBestChromosome`, which helps maintain high-quality solutions.

2. **Elitism Rate and Iterations**: The number of elite individuals to retain is determined by `elitismRate`. Specifically, this rate is multiplied by the population size to calculate the number of copies of the best solution to retain in the next generation.

3. **Population Replacement**: The selected elite chromosomes are copied into a temporary vector, `temp`, and then swapped with the original population, ensuring that the best solutions are propagated forward.

This elitism operator helps maintain strong solutions across generations, balancing exploration and exploitation by retaining top performers, thus enhancing convergence towards an optimal solution.


### Crossover Operator

The crossover operator combines two solutions, $S_1$ and $S_2$, from the current population of 1000 chromosomes. This is achieved through the following steps:

1. **Selection of Crossover Points**: Random indices, $R_1$ and $R_2$, are chosen within the chromosome length.
2. **Gene Exchange**: The labels (genes) between indices $R_1$ and $R_2$ in $S_1$ and $S_2$ are swapped, producing two new offspring solutions.
3. **Feasibility Check**: The resulting offspring solutions are verified to ensure they satisfy the constraints of the Double Roman Domination Function (DRDF). Solutions that violate constraints are discarded or adjusted.

In this proposed algorithm, elitism and mutation rate do not influence the crossover operator directly.

### RVNS (Random Variable Neighborhood Search)

The RVNS method further refines the best solution found by the genetic algorithm, seeking solutions in nearby solution spaces. This neighborhood exploration enhances the chances of finding an improved solution, closer to the optimal. The RVNS process involves the following steps:

1. **DestroySolution Sub-Routine**: A randomly chosen vertex labeled 0, 2, or 3 is selected and "destroyed" by unlabeled it (setting its label to -1). This introduces a controlled disruption in the current solution.

2. **Solution Reconstruction**: After "destroying" part of the solution, the following sub-routines are applied to form a new solution:
   - **ConstructSolution**: Partially reconstructs the solution, assigning feasible labels to vertices based on DRDF rules.
   - **ExtendSolution**: Expands the labeling where possible to improve coverage or reduce weights.
   - **ReduceSolution**: Optimizes the labeling to reduce the overall weight while satisfying DRDF constraints.

3. **Solution Comparison**: The newly formed solution is compared to the original. If the new solution has a lower weight, it replaces the current solution as the best candidate.

This iterative improvement through RVNS allows for fine-tuning of the solution, making it a valuable addition to the genetic algorithm's exploration-exploitation balance.


## Ant Colony Optimization (ACO)

The Ant Colony Optimization algorithm simulates the behavior of ants searching for food, where ants deposit pheromones along their paths. The best path is the one with the highest pheromone concentration, guiding future ants. In this algorithm, each vertex in the graph is associated with an initial pheromone value of 0.5, and each solution is represented as a vector of size $|V|$. The weight of the solution is the sum of the labels in the solution vector, with the goal of finding a solution with the smallest possible weight. The algorithm runs until a maximum number of iterations is reached or there are no improvements in the solution for a set number of iterations.

The type of Ant Colony Optimization proposed is based on two variants: **MAX-MIN Ant System** (MMAS) and **Hyper-Cube Framework** (HCF).

### ACO Sub-routines

1. **ConstructSolution**: A random vertex is selected and labeled as 4, while its neighbors are labeled as 0. These vertices are then removed from the auxiliary graph. If there are any remaining vertices with degree 0, they are labeled as 2. After this, the algorithm checks if there are any vertices labeled 2 that have only neighbors labeled as 0; if so, a random neighbor is assigned the label 2 to ensure the constraints of the TRDF.

2. **ExtendSolution**: A vertex currently labeled 0, 2, or 3 is selected and upgraded to a label of 4 in the current solution.

3. **ReduceSolution**: Redundant vertices are removed from the current solution. A vertex $v$ is considered redundant if all vertices in its closed neighborhood $N[v]$ are already dominated by other vertices in the solution.

4. **RVNS (Random Variable Neighborhood Search)**: This function enhances the best result obtained from the ACO algorithm by exploring nearby solution spaces. RVNS increases the likelihood of finding a solution closer to the optimal one by modifying the solution in the following steps:
    - The **DestroySolution** sub-routine randomly selects a vertex labeled 0, 2, or 3 and also unlabeled it, setting its value to -1.
    - **ConstructSolution**, **ExtendSolution**, and **ReduceSolution** are then reapplied to form a new solution, which is compared with the previous one. If the new solution has a lower weight, it replaces the current solution.


## How to Use

1. Clone the repository
   ```bash
   git clone https://github.com/isrreal/Double-Roman-Domination-In-Graphs-meta-heuristics.git
### 2. Edit the File "graph.txt":

The `graph.txt` file is essential for defining the structure of the graph. It should be formatted in a way that clearly outlines the relationships between the vertices.

1. **File Structure**:
   **Edge List**:
      - Subsequent lines must list pairs of integers, each representing an edge between two vertices. Each pair indicates a direct connection between the specified vertices.

**Example of `graph.txt`**:

By default, the Roman Empire graph is represented as follows:

![Graph Representation](https://github.com/user-attachments/assets/d7ec6492-6090-4e5a-8173-01f1ace57b5a)

#### Vertex Mapping:
- **Britannia** - 0
- **Gaul** - 1
- **Iberia** - 2
- **Rome** - 3
- **North Africa** - 4
- **Constantinople** - 5
- **Asia Minor** - 6
- **Egypt** - 7

The content below illustrates the Roman Empire represented as a graph with **8 vertices**:

```plaintext
0 1 0 2
1 2 1 3
2 3 2 4 
3 4 3 5 3 7
4 7 5 6
5 7 6 7
```
   
5. Run the Code

To execute the application, use the following syntax:

```bash
# Syntax:
# ./app <population_size> <generations> <chromosome_creation_heuristic> <mutation_rate> <elitism_rate> <number_of_ants> <iterations>

# Parameters:
#   population_size               - Defines the size of the initial population. A larger size increases solution diversity but may increase computation time.
#   generations                   - Number of generations (iterations) the Genetic Algorithm will run to evolve solutions.
#   chromosome_creation_heuristic - Heuristic used for initializing chromosomes:
#                                     1 - Basic randomized approach
#                                     2 - Weighted by vertex degree
#                                     3 - Weighted by neighborhood constraints
#   mutation_rate                 - Rate at which mutation occurs, altering labels on selected vertices to introduce diversity and escape local optima.
#   elitism_rate                  - Percentage of top-performing solutions retained in each generation, ensuring the best solutions are carried over.
#   number_of_ants                - Number of ants used in the Ant Colony Optimization (ACO) process, influencing solution exploration and reinforcement.
#   iterations                    - Number of iterations for the ACO phase, allowing for finer convergence in conjunction with the Genetic Algorithm.

# Example:
./app 100 200 1 0.05 0.1 50 100
```

