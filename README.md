# Metaheuristics and Integer Linear Programming for the Triple Roman Domination Problem in Graphs

This repository contains exact and metaheuristic approaches for solving the **Triple Roman Domination Problem (TRDP)**, an NP-complete combinatorial optimization problem defined on graphs.

The project was developed as part of my undergraduate thesis in Computer Science at the **Federal University of Ceará — UFC**.

## Undergraduate Thesis

**Original title**

> Meta-heurísticas e Programação Linear Inteira para o Problema da Dominação Romana Tripla em Grafos

**English title**

> Metaheuristics and Integer Linear Programming for the Triple Roman Domination Problem in Graphs

- **Author:** Israel Souza Ferreira
- **Advisor:** Prof. Dr. Atilio Gomes Luiz
- **Institution:** Federal University of Ceará — UFC
- **Campus:** Quixadá
- **Approval date:** February 28, 2025

📄 **[Read the full undergraduate thesis](https://drive.google.com/file/d/1IxBw_odR99kklE6b3ANmi-9bubvGC4mp/view?usp=sharing)**

---

## Overview

The project investigates three approaches for the Triple Roman Domination Problem:

1. an **Integer Linear Programming formulation** for exact optimization;
2. the **Ferreira–Luiz Genetic Algorithm**, referred to as **FLGA**;
3. the **Ferreira–Luiz Ant Colony Optimization algorithm**, referred to as **ACO-FL**.

According to the literature review conducted in the thesis, this work presents:

- the first Genetic Algorithm specifically designed for the TRDP;
- the first Ant Colony Optimization algorithm specifically designed for the TRDP;
- the first correct Integer Linear Programming formulation for the problem in the literature;
- a counterexample showing that a previously published ILP formulation does not correctly model the TRDP.

The metaheuristics were evaluated against exact or best-known solutions obtained through Integer Linear Programming.

---

## Main Contributions

- Identified an error in a previously published Integer Linear Programming formulation for the TRDP.
- Constructed a counterexample using the cycle $C_{10}$, showing that the previous model could return an invalid Triple Roman Dominating Function.
- Developed a new Integer Linear Programming formulation with:
  - $3\lvert V\rvert$ binary decision variables;
  - $2\lvert V\rvert$ constraints.
- Implemented the exact model using Pyomo and NetworkX.
- Solved the model with CPLEX and CBC.
- Designed and implemented the FLGA metaheuristic.
- Developed four initial-population strategies for FLGA.
- Implemented tournament selection, elitism, mutation, and one-point and two-point crossover.
- Implemented specialized feasibility-repair and label-reduction procedures.
- Designed and implemented the ACO-FL metaheuristic.
- Combined MAX-MIN Ant System with the Hyper-Cube Framework.
- Integrated Reduced Variable Neighborhood Search into ACO-FL.
- Automatically calibrated the metaheuristic parameters using `irace`.
- Evaluated the approaches on 362 graph instances.
- Generated 30 Erdős–Rényi random graphs for the experimental study.
- Produced graph instances, reference values, raw outputs, and experimental tables that can support future research.

---

## Triple Roman Domination Problem

Let

$$
G=(V,E)
$$

be a graph, and let

$$
h:V(G)\rightarrow\{0,1,2,3,4\}
$$

be a vertex-labeling function.

A vertex $v\in V(G)$ is called **active** when

$$
h(v)>0.
$$

The active neighborhood of $v$, denoted by $AN(v)$, is the set of active neighbors of $v$.

A function $h$ is a **Triple Roman Dominating Function**, or TRDF, when, for every vertex $v\in V(G)$,

$$
\sum_{u\in N[v]}h(u)\geq 3+\lvert AN(v)\rvert.
$$

Equivalently, the following conditions must hold.

### Vertices labeled 0

If $h(v)=0$, at least one of these conditions must be satisfied:

- $v$ has a neighbor labeled 4;
- $v$ has one neighbor labeled 2 and another labeled 3;
- $v$ has at least three neighbors labeled 2.

### Vertices labeled 1

If $h(v)=1$, at least one of these conditions must be satisfied:

- $v$ has a neighbor labeled at least 3;
- $v$ has at least two neighbors labeled 2.

### Vertices labeled 2

If $h(v)=2$, then $v$ must have a neighbor labeled at least 2.

The weight of a TRDF is

$$
w(h)=\sum_{v\in V(G)}h(v).
$$

The **Triple Roman domination number** of $G$, denoted by $\gamma_{3R}(G)$, is the minimum possible weight of a TRDF:

$$
\gamma_{3R}(G)
=
\min_h
\sum_{v\in V(G)}h(v).
$$

The Triple Roman Domination Problem consists of determining $\gamma_{3R}(G)$ for an arbitrary graph $G$.

Its decision version is NP-complete, including when restricted to some specific graph classes.

---

# Integer Linear Programming

## Analysis of a Previous Formulation

A previously published Integer Linear Programming formulation for the TRDP used six binary variables for each graph vertex.

During this research, the formulation was found to be incorrect.

A counterexample based on the cycle $C_{10}$ showed that the model could produce a labeling satisfying all of its linear constraints while violating the definition of a Triple Roman Dominating Function.

The issue was associated with the constraint intended to represent the requirements for vertices labeled 0.

Because the same incorrect constraint was retained in refinements of the original formulation, those models were also unable to represent the TRDP correctly.

## Proposed Formulation

The corrected formulation uses the result that every connected nontrivial graph admits an optimal TRDF without vertices labeled 1.

Three binary variables are therefore associated with each vertex $v$:

$$
q_v=
\begin{cases}
1, & \text{if } h(v)=2,\\
0, & \text{otherwise},
\end{cases}
$$

$$
r_v=
\begin{cases}
1, & \text{if } h(v)=3,\\
0, & \text{otherwise},
\end{cases}
$$

and

$$
s_v=
\begin{cases}
1, & \text{if } h(v)=4,\\
0, & \text{otherwise}.
\end{cases}
$$

The objective is to minimize the total weight:

$$
\min
\left(
\sum_{v\in V(G)}2q_v
+
\sum_{v\in V(G)}3r_v
+
\sum_{v\in V(G)}4s_v
\right).
$$

The main domination constraint is derived directly from the definition of a Triple Roman Dominating Function:

$$
(2q_v+3r_v+4s_v)
+
\sum_{u\in N(v)}
(2q_u+3r_u+4s_u)
\geq
3+
\sum_{u\in N(v)}
(q_u+r_u+s_u),
\qquad
\forall v\in V(G).
$$

Each vertex receives at most one positive label:

$$
q_v+r_v+s_v\leq 1,
\qquad
\forall v\in V(G).
$$

The decision variables are binary:

$$
q_v,r_v,s_v\in\{0,1\},
\qquad
\forall v\in V(G).
$$

The resulting formulation contains:

- $3\lvert V\rvert$ binary variables;
- $2\lvert V\rvert$ constraints.

## Implementation

The exact model was implemented using:

- Python;
- Pyomo;
- NetworkX;
- IBM CPLEX;
- CBC.

CPLEX and CBC were executed with a maximum time limit of 900 seconds per instance.

When the solver proved optimality, the returned value was used as the exact Triple Roman domination number. Otherwise, the best solution found within the execution limit was retained as a reference value.

---

# Ferreira–Luiz Genetic Algorithm

The **Ferreira–Luiz Genetic Algorithm**, or **FLGA**, represents a candidate solution as a chromosome of size $\lvert V\rvert$.

Each gene corresponds to one graph vertex and contains a label from

$$
\{0,2,3,4\}.
$$

The fitness of a chromosome is its total labeling weight:

$$
\operatorname{fitness}(h)
=
\sum_{v\in V(G)}h(v).
$$

The objective is to minimize the fitness while preserving all TRDF constraints.

## Main Components

FLGA includes:

- specialized initial-population heuristics;
- tournament selection;
- elitism;
- one-point crossover;
- two-point crossover;
- mutation;
- feasibility checking;
- infeasible-solution repair;
- label-reduction procedures;
- early stopping based on iterations without improvement.

---

## Initial-Population Heuristics

Four strategies were implemented for generating the initial population.

### Heuristic 1 — Randomized Label-2 Construction

An auxiliary graph is created as a copy of the original graph.

At each step:

1. a vertex is randomly selected;
2. the selected vertex receives label 2, or label 3 when isolated;
3. its neighbors receive label 0;
4. the selected vertex and its neighbors are removed from the auxiliary graph;
5. remaining isolated vertices receive label 3;
6. a feasibility-repair procedure is applied.

### Heuristic 2 — Randomized Label-4 Construction

At each step:

1. a vertex is randomly selected;
2. the selected vertex receives label 4, or label 3 when isolated;
3. its neighbors receive label 0;
4. the selected vertex and its neighbors are removed;
5. remaining isolated vertices receive label 3;
6. the label-reduction procedure attempts to decrease the total solution weight.

### Heuristic 3 — Degree-Based Greedy Construction

Vertices are ordered by decreasing degree.

At each step:

1. a highest-degree valid vertex is selected;
2. the selected vertex receives label 4, or label 3 when isolated;
3. its neighbors receive label 0;
4. the selected vertex and its neighbors are removed;
5. the label-reduction procedure is applied.

### Heuristic 4 — Mixed Population

The fourth strategy combines chromosomes produced by the three previous heuristics.

Approximately one-third of the initial population is generated by each construction method.

This mixed strategy was the strongest overall FLGA initialization method in the experiments.

---

## Tournament Selection

FLGA uses tournament selection.

A fixed number of chromosomes is randomly sampled from the population, and the chromosome with the smallest fitness is selected.

This strategy favors high-quality solutions while retaining population diversity.

---

## Elitism

The elitism procedure preserves the best chromosomes from the current generation.

The number of preserved solutions is

$$
k=
\left\lceil
\text{population size}\cdot\text{elitism rate}
\right\rceil.
$$

The population is sorted by fitness, and the $k$ best chromosomes are transferred directly to the next generation.

---

## One-Point Crossover

Given two parent chromosomes, a random crossover index is selected.

The genes after that index are exchanged, producing two offspring.

A feasibility procedure is then applied to repair assignments that violate the TRDF constraints.

---

## Two-Point Crossover

Two random indices are selected.

The chromosome segments between these indices are exchanged.

The resulting offspring are checked and repaired before being inserted into the population.

---

## Mutation

Mutation is applied according to the configured mutation probability.

When mutation occurs:

1. a random gene is selected;
2. its value is replaced by a random label from $\{0,2,3,4\}$;
3. the feasibility-repair procedure is applied.

---

## Feasibility Repair

Crossover, mutation, and some population heuristics can generate infeasible chromosomes.

The feasibility procedure inspects the labels and neighborhoods of the graph vertices and repairs invalid assignments.

Examples include:

- changing an unsupported vertex labeled 0 to label 2 or 3;
- changing a vertex labeled 2 to label 3 when it has no neighboring vertex labeled at least 2;
- re-evaluating the chromosome fitness after repairs.

---

## Label Reduction

The label-reduction procedure attempts to decrease labels without compromising feasibility.

For each positively labeled vertex, the method tests progressively smaller values.

A reduction is retained only when:

- the modified vertex remains feasible;
- its active neighborhood remains feasible;
- the complete labeling continues satisfying the TRDF constraints.

---

# Ferreira–Luiz Ant Colony Optimization

The **Ferreira–Luiz Ant Colony Optimization algorithm**, or **ACO-FL**, is based on:

- Ant Colony Optimization;
- MAX-MIN Ant System;
- Hyper-Cube Framework;
- Reduced Variable Neighborhood Search.

A candidate solution is represented as a vector of size $\lvert V\rvert$, where each position stores the label assigned to one graph vertex.

Each vertex $v$ has an associated pheromone value:

$$
\tau_v\in[0,1].
$$

The initial pheromone value is

$$
\tau_v=0.5.
$$

Higher pheromone values increase the probability of including a vertex as an important component of a constructed solution.

---

## Main ACO-FL Stages

Each ACO-FL iteration applies the following stages:

1. initialize or read the current pheromone vector;
2. construct candidate solutions;
3. extend the candidate solutions;
4. reduce redundant or excessive labels;
5. apply RVNS;
6. update the best solution;
7. update pheromone values;
8. compute the convergence factor;
9. reinitialize pheromones when convergence becomes excessive.

---

## Vertex Selection

Vertex selection combines graph structure and pheromone information.

The evaluation function is

$$
f(v)=\deg(v)\tau_v.
$$

Depending on the configured probability, the algorithm either:

- selects the vertex maximizing $f(v)$; or
- applies roulette-wheel selection using values proportional to $f(v)$.

---

## Construct Solution

The construction procedure operates over an auxiliary graph.

At each iteration:

1. a vertex is selected;
2. the selected vertex receives label 4;
3. its neighbors receive label 0;
4. the selected vertex and its neighbors are removed from the auxiliary graph;
5. remaining isolated vertices receive label 3.

The process continues until the auxiliary graph becomes empty.

---

## Extend Solution

The extension procedure selects vertices labeled below 4 and upgrades some of them to label 4.

The number of selected vertices is determined by the configured extension rate.

This operation introduces additional domination capacity before the reduction phase.

---

## Reduce Solution

Vertices are processed in decreasing degree order.

For every vertex labeled 2, 3, or 4, the procedure attempts to reduce its label while preserving feasibility.

This stage decreases the total solution weight after construction and extension.

---

## Pheromone Update

The pheromone update is influenced by:

- the best solution found in the current iteration;
- the best global solution found so far;
- the current convergence factor.

During early exploration, greater importance is assigned to the current iteration solution.

As convergence increases, the global best solution receives greater influence.

Pheromones are restricted to

$$
\tau_{\min}=0.001
$$

and

$$
\tau_{\max}=0.999.
$$

When the convergence factor becomes greater than 0.99, pheromone values are reinitialized to encourage further exploration.

---

# Reduced Variable Neighborhood Search

ACO-FL integrates a **Reduced Variable Neighborhood Search**, or RVNS, to improve constructed solutions.

The RVNS procedure applies four main operations:

1. destroy part of the current solution;
2. reconstruct the destroyed assignments;
3. extend the reconstructed solution;
4. reduce labels while preserving feasibility.

## Destroy Solution

Some vertices labeled 0, 2, or 3 are temporarily assigned the value `-1`.

The destruction intensity changes according to the current neighborhood:

$$
\mu
=
d_{\min}
+
(k_{\mathrm{current}}-1)
\left(
\frac{d_{\max}-d_{\min}}
{k_{\max}-1}
\right).
$$

Smaller neighborhoods produce limited perturbations and favor intensification.

Larger neighborhoods produce stronger perturbations and increase diversification.

## Acceptance

When the reconstructed solution has a smaller weight, it replaces the current solution and the neighborhood index is reset.

Otherwise, the neighborhood index is increased.

The process stops when:

- the maximum number of RVNS iterations is reached; or
- the maximum number of iterations without improvement is reached.

---

# Experimental Evaluation

## Execution Environment

The experiments were conducted using:

- **Processor:** Intel Core i5-8265U at 1.60 GHz;
- **Memory:** 8 GB DDR4 at 2400 MHz;
- **Operating system:** Ubuntu 22.04.5 LTS;
- **Compiler:** G++ 11.4.0;
- **Language standard:** C++17.

The C++ implementations were compiled using optimization flags including:

```bash
-std=c++17 -Wall -Wextra -Ofast -finline-functions -march=native
```

The Integer Linear Programming formulation was implemented in Python with Pyomo and NetworkX.

---

## Graph Collections

The complete experimental study used **362 graph instances**:

| Collection | Instances |
|:---|---:|
| BAI sparse matrices | 50 |
| Harwell–Boeing collection | 186 |
| Miscellaneous Networks | 56 |
| Classical graph families | 40 |
| Erdős–Rényi random graphs | 30 |
| **Total** | **362** |

The classical graph families included:

- cycles;
- paths;
- stars;
- trees.

The Erdős–Rényi graphs had:

- orders between 25 and 250 vertices;
- connection probabilities of 0.2, 0.5, and 0.8.

---

## Evaluation Criteria

The algorithms were compared according to:

- execution time;
- best solution weight;
- fitness;
- relative optimality gap;
- proximity to exact or best-known ILP solutions;
- behavior across different graph densities;
- behavior across different graph families.

---

## Parameter Tuning

The FLGA and ACO-FL parameters were automatically calibrated using `irace`.

The tuning process used a representative subset of **115 graph instances**, containing small, medium, and large graphs from the experimental collections.

### FLGA Parameters

The calibrated FLGA configuration included:

| Parameter | Tuned value |
|:---|---:|
| Population divisor | 2 |
| Generations | 532 |
| Elitism rate | 0.3288 |
| Mutation rate | 0.4886 |
| Crossover rate | 0.4862 |
| Tournament size | 2 |
| Maximum iterations without improvement | 86 |

The population size was defined according to the graph order divided by the tuned population divisor.

### ACO-FL Parameters

The calibrated ACO-FL configuration with RVNS included:

| Parameter | Tuned value |
|:---|---:|
| Number of ants | 1 |
| Iterations | 5 |
| Evaporation rate | 0.2871 |
| Minimum destruction rate | 0.243 |
| Maximum destruction rate | 0.7 |
| Maximum RVNS neighborhoods | 5 |
| Maximum RVNS iterations | 50 |
| Maximum RVNS iterations without improvement | 20 |
| Construct-solution selection rate | 0.1 |
| Extend-solution selection rate | 0.9 |
| Extend-solution addition rate | 0.05 |

---

# Key Findings

## Integer Linear Programming

The research identified that a previous ILP formulation did not correctly represent the TRDP.

A new formulation was proposed using fewer decision variables and constraints:

- $3\lvert V\rvert$ binary variables instead of $6\lvert V\rvert$;
- $2\lvert V\rvert$ constraints instead of $4\lvert V\rvert$.

The corrected model was used to obtain exact or best-known reference solutions for the experimental instances.

---

## Genetic Algorithm

Four initial-population strategies were evaluated.

The fourth strategy, which combines the three individual construction heuristics, produced the best overall FLGA performance.

It generally achieved:

- lower fitness values;
- competitive execution times;
- better population diversity;
- greater robustness across different random graph densities.

The configuration using the fourth heuristic was selected as the final FLGA algorithm.

---

## Ant Colony Optimization

ACO was evaluated with and without RVNS.

The version with RVNS generally produced better solutions, particularly on several smaller and less dense instances.

This demonstrated that the local-search stage was an important component of the final ACO-FL method.

---

## ACO-FL vs FLGA

ACO-FL obtained the strongest overall performance in most graph collections.

It was generally more effective at finding lower-weight solutions, especially in:

- larger graph instances;
- sparse graph collections;
- cycles and paths;
- several structurally complex instances.

FLGA remained competitive in specific scenarios, including:

- small and dense graphs;
- star graphs;
- instances favorable to population diversity;
- some sparse or structured graphs;
- scenarios where execution time was the primary concern.

In the classical graph-family experiments:

- ACO-FL frequently found optimal or near-optimal solutions for cycles and paths;
- both algorithms found the optimal value for star graphs;
- FLGA was generally faster on star graphs;
- ACO-FL usually produced better-quality solutions for trees.

Across the complete experimental study, only six graph instances produced a relative gap greater than 50%.

The remaining instances stayed below that threshold, indicating that both proposed metaheuristics were competitive on the evaluated datasets.

---

## Approach Summary

| Approach | Type | Main purpose | Main outcome |
|:---|:---|:---|:---|
| **Integer Linear Programming** | Exact optimization | Compute exact or best-known reference solutions | Correct formulation with $3\lvert V\rvert$ variables and $2\lvert V\rvert$ constraints |
| **FLGA** | Population-based metaheuristic | Explore the solution space through evolutionary operators | Competitive on specific small, dense, and structured graphs |
| **ACO-FL** | Constructive hybrid metaheuristic | Combine pheromone learning with RVNS | Best overall solution quality across most evaluated collections |

> These conclusions apply to the graph instances, execution limits, parameter configurations, and computational environment used in this study.

---

# Repository Structure

```text
.
├── PLI/           # Integer Linear Programming implementation and experiments
├── Tabelas/       # Generated experimental tables
├── inc/           # C++ header files
├── input_files/   # Graph instances
├── output_files/  # Raw algorithm outputs
├── results/       # Consolidated experimental results
├── scripts/       # Experiment and result-processing scripts
├── src/           # C++ source code
├── makefile
└── README.md
```

---

# Building the Project

## Requirements

### Metaheuristics

- GNU/Linux or another compatible Unix-like environment;
- G++ with C++17 support;
- GNU Make;
- Bash;
- Python 3 for experiment-processing scripts.

### Integer Linear Programming

- Python 3;
- Pyomo;
- NetworkX;
- CBC or IBM CPLEX.

### Parameter Tuning

- `irace`;
- R.

---

## Clone the Repository

```bash
git clone https://github.com/isrreal/Triple-Roman-Domination-in-graphs.git
cd Triple-Roman-Domination-in-graphs
```

---

## Compile

```bash
make
```

---

# Running FLGA

The basic execution format is:

```bash
./app <graph_file> <graph_name> <initialization_heuristic>
```

Example:

```bash
./app input_files/graph.txt graph_name 4
```

The initialization heuristic is selected using:

| Value | Strategy |
|---:|:---|
| `1` | Randomized construction primarily using label 2 |
| `2` | Randomized construction primarily using label 4 |
| `3` | Degree-based greedy construction |
| `4` | Mixed population using all construction strategies |

Additional parameters such as population size, generations, mutation rate, crossover rate, and elitism rate are defined by the implementation or the experimental scripts.

---

# Running ACO-FL and the ILP Model

The scripts used to execute ACO-FL, RVNS experiments, parameter tuning, and Integer Linear Programming experiments are organized in:

```text
scripts/
PLI/
```

Because some commands depend on the selected graph collection, solver, and experiment configuration, consult the scripts in these directories for the exact commands used in the thesis experiments.

---

# Reproducibility

The repository contains:

- C++ implementations of FLGA and ACO-FL;
- the Python Integer Linear Programming implementation;
- graph instances used in the experiments;
- randomly generated graph instances;
- exact or best-known ILP reference values;
- raw metaheuristic outputs;
- experiment-execution scripts;
- parameter-tuning configurations;
- result-processing scripts;
- consolidated experimental tables.

Because FLGA and ACO-FL are stochastic, multiple independent executions should be performed for each graph instance and parameter configuration.

For reproducible comparisons, record:

- random seed;
- compiler version;
- compiler flags;
- solver version;
- solver time limit;
- algorithm parameters;
- graph source and preprocessing steps;
- number of independent executions.

---

# Limitations

- FLGA and ACO-FL do not guarantee optimal solutions.
- Their performance depends on parameter configuration and graph structure.
- Solver results may not be proven optimal when the time limit is reached.
- Some experimental scripts depend on the original GNU/Linux environment.
- The experimental evaluation was conducted on a single hardware configuration.
- The current implementation prioritizes research experimentation over a library-style API.
- The conclusions should not be generalized to every graph class without additional experiments.

---

# Future Work

## Genetic Algorithm

Potential extensions include:

- alternative selection strategies;
- additional crossover operators;
- adaptive mutation rates;
- alternative elitism policies;
- new initial-population heuristics;
- hybridization with local-search algorithms;
- graph-class-specific parameter configurations;
- automated parameter tuning over larger training collections.

## Ant Colony Optimization

Potential extensions include:

- alternative pheromone-update policies;
- adaptive pheromone limits;
- different candidate-selection mechanisms;
- additional construction heuristics;
- different RVNS neighborhood structures;
- alternative destruction and reconstruction procedures;
- broader evaluation of the interaction between ACO and local search.

## Integer Linear Programming

Potential extensions include:

- valid inequalities;
- symmetry-breaking constraints;
- graph preprocessing;
- graph-reduction rules;
- alternative variable formulations;
- decomposition approaches;
- comparisons between additional optimization solvers;
- expansion of the collection of instances with proven optimal values.

## Performance Engineering

Implementation improvements may include:

- profiling computational bottlenecks;
- more efficient graph data structures;
- reducing dynamic memory allocations;
- parallel execution of independent experiments;
- parallel chromosome evaluation;
- parallel ant-solution construction;
- additional compiler and architecture optimizations;
- improved experiment orchestration.

---

# Thesis Citation

```bibtex
@misc{ferreira2025triple,
  author       = {Israel Souza Ferreira},
  title        = {Meta-Heurísticas e Programação Linear Inteira para o Problema da Dominação Romana Tripla em Grafos},
  year         = {2025},
  howpublished = {Undergraduate thesis, Federal University of Ceará},
  address      = {Quixadá, Brazil},
  note         = {Advisor: Atilio Gomes Luiz}
}
```

---

# Technologies and Methods

- C++;
- Python;
- Bash;
- GNU Make;
- Linux;
- Pyomo;
- NetworkX;
- IBM CPLEX;
- CBC;
- `irace`;
- Integer Linear Programming;
- Genetic Algorithms;
- Tournament Selection;
- Ant Colony Optimization;
- MAX-MIN Ant System;
- Hyper-Cube Framework;
- Reduced Variable Neighborhood Search;
- Graph Theory;
- Combinatorial Optimization;
- Automated Parameter Tuning.

---

# Author

**Israel Souza Ferreira**

- GitHub: [@isrreal](https://github.com/isrreal)
- LinkedIn: [Israel Souza](https://www.linkedin.com/in/israel-souza-84b8102b0)
- Email: [israel.souza.dev@gmail.com](mailto:israel.souza.dev@gmail.com)
