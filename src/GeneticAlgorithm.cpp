#include "GeneticAlgorithm.hpp"
  
/**
 * @brief Creates a population of chromosomes with a specific number of genes.
 * 
 * If a heuristic function is provided, the chromosomes are initialized using this heuristic.
 * Otherwise, chromosomes are initialized with random genes.
 * 
 * @param heuristic A pointer to a function that generates chromosomes based on a graph.
 * @param graph Object to the graph used to initialize the chromosomes.
 */

void GeneticAlgorithm::createPopulation(
    std::vector<std::function<Chromosome(const Graph&)>> generateChromosomeHeuristics, 
    const Graph& graph, 
    size_t heuristic) {
    
	if (!generateChromosomeHeuristics.empty()) {
        if (heuristic == 4) { 
            size_t portion_size { static_cast<size_t>(population_size / 3) };
			size_t remainder { population_size % 3 }; 
			size_t index {0};

			for (size_t i {0}; i < 3; ++i) {
				for (size_t j {0}; j < portion_size; ++j) {
					this->population[index++] = generateChromosomeHeuristics[i](graph);
				}
			}

			for (size_t i {0}; i < remainder; ++i) {
				this->population[index++] = generateChromosomeHeuristics[i](graph);
			}
			
			std::random_device rd;
			std::mt19937 random_engine(rd());
			std::shuffle(population.begin(), population.end(), random_engine);		
        } 
        
        else if (heuristic > 0 && heuristic < 4) {
            Chromosome func { generateChromosomeHeuristics[heuristic - 1](graph) } ;
            for (size_t i {0}; i < population_size; ++i) {
                this->population[i] = func;          
            }
        }
        
     	else {
		    throw population;
    	}
    } 
}


/**
 * @brief Generates a new population by crossing over Chromosomes from the current population.
 *     
 * 
 * @return std::vector<Chromosome>& A new population of Chromosomes.
 */

std::vector<Chromosome>& GeneticAlgorithm::createNewPopulation() {
    std::vector<Chromosome> temp {population};
   	
   	this->elitism(this->elitism_rate);
    
    temp.reserve(population_size);
    
    Chromosome selected1;
    Chromosome selected2;
    Chromosome offspring;
    
    while (population.size() < population_size) {          
        selected1 = tournamentSelection(temp); 
        selected2 = tournamentSelection(temp);      	
       	
       	if (getRandomFloat(0.0, 1.0) <= crossover_rate) {  		
        	offspring = this->twoPointCrossOver(selected1, selected2);
        } 
        
        else {      
        	offspring = this->onePointCrossOver(selected1, selected2);     
        }
        
	    mutation(offspring);	    
		
        population.emplace_back(offspring);       
    }
	
    population.swap(temp);
    
    return population;
}

std::vector<Chromosome>& GeneticAlgorithm::elitism(float elitism_rate) {
    size_t iterations { static_cast<size_t>(std::ceil(this->population_size * elitism_rate)) };

    Chromosome best_solution = [&]() {
		Chromosome best_solution { population[getRandomInt(0, population_size - 1)] };
		Chromosome solution;
		for (const auto& chromosome: population) {
			solution = chromosome;		
			if (chromosome.fitness > best_solution.fitness) {
				best_solution = chromosome;
			}
		}
		return best_solution;
	}();
    
    population.clear();
  	population.reserve(iterations);
  	
    for (size_t i {0}; i < iterations; ++i) {
    	this->population.emplace_back(best_solution);
	}
	
    return this->population;
}
   
Chromosome& GeneticAlgorithm::mutation(Chromosome& chromosome) {	
	if (getRandomFloat(0.0, 1.0) < this->mutation_rate) {
		std::vector<size_t> labels {0, 2, 3, 4};
		size_t randomIndex { getRandomInt(0, genes_size - 1) };
		short random_label { static_cast<short>(getRandomInt(0, labels.size() - 1)) };
			
		chromosome.genes[randomIndex] = labels[random_label];
		feasibilityCheck(this->graph, chromosome);
	}

	return chromosome;   
}

Chromosome& GeneticAlgorithm::onePointCrossOver(const Chromosome& chromosome1, const Chromosome& chromosome2) {
  
   size_t index { getRandomInt(0, genes_size - 1) };
   
   Chromosome solution1 { chromosome1 };
   Chromosome solution2 { chromosome2 };
    
   for(size_t i {index + 1}; i <= genes_size - 1; ++i) {
		std::swap(solution1.genes[i], solution2.genes[i]);
   }
   
   feasibilityCheck(this->graph, solution1);
   feasibilityCheck(this->graph, solution2);
   
   fitness(solution1);
   fitness(solution2);
    
   return chooseBestSolution(solution1, solution2);
}

Chromosome& GeneticAlgorithm::twoPointCrossOver(const Chromosome& chromosome1, const Chromosome& chromosome2) {
   size_t range1 { getRandomInt(0, genes_size - 1) };
   size_t range2 { getRandomInt(0, genes_size - 1) };
   
   Chromosome solution1 { chromosome1 };
   Chromosome solution2 { chromosome2 };

   if (range1 > range2) {
        std::swap(range1, range2);
   }
   
   for (size_t i {range1} ; i <= range2; ++i) {
		std::swap(solution1.genes[i], solution2.genes[i]);
   }
        
   feasibilityCheck(this->graph, solution1);
   feasibilityCheck(this->graph, solution2);
   
   fitness(solution1);
   fitness(solution2);
    
   return chooseBestSolution(solution1, solution2);
}

Chromosome& GeneticAlgorithm::tournamentSelection(const std::vector<Chromosome>& population) {
    constexpr float parameter { 0.75f }; 
    
    const Chromosome& c1 { population[getRandomInt(0, population.size() - 1)] };
    const Chromosome& c2 { population[getRandomInt(0, population.size() - 1)] };

    float probability { getRandomFloat(0.0, 1.0) };

    if (probability < parameter) {
        return chooseBestSolution(const_cast<Chromosome&>(c1), const_cast<Chromosome&>(c2));
    } 
    else {
        return chooseWorstSolution(const_cast<Chromosome&>(c1), const_cast<Chromosome&>(c2));
    }
}


Chromosome& GeneticAlgorithm::chooseBestSolution(Chromosome& chromosome1, Chromosome& chromosome2) {
    return (chromosome1.fitness < chromosome2.fitness ? chromosome1 : chromosome2);
}

Chromosome& GeneticAlgorithm::chooseWorstSolution(Chromosome& chromosome1, Chromosome& chromosome2) {
    return (chromosome1.fitness > chromosome2.fitness ? chromosome1 : chromosome2);
}
		
// public methods 

std::vector<Chromosome> GeneticAlgorithm::getPopulation() { return population; }

size_t GeneticAlgorithm::getPopulationSize() { return population_size; }

size_t GeneticAlgorithm::getGenesSize() { return genes_size; }

size_t GeneticAlgorithm::getGenerations() { return generations; }

double GeneticAlgorithm::getMutationRate() { return mutation_rate; }

double GeneticAlgorithm::getElitismRate() { return elitism_rate; }

Graph GeneticAlgorithm::getGraph() { return graph; }

std::vector<int> GeneticAlgorithm::getBestSolution() { return best_solution; }

void GeneticAlgorithm::run(size_t generations, std::vector<std::function<Chromosome(const Graph&)>> heuristics, size_t chosen_heuristic) { 
   this->createPopulation(heuristics, graph, chosen_heuristic);
   
   Chromosome current_best_solution { tournamentSelection(this->population) };   
   
   Chromosome best_solution { current_best_solution };
   
   size_t iteration {0};
   size_t current_no_improvement_iteration {0};
   
   while ((iteration < generations) && (current_no_improvement_iteration < max_no_improvement_iterations)) {
   
		this->population.swap(createNewPopulation());
       	
        current_best_solution = tournamentSelection(population);                                       
		
        ++iteration;
        ++current_no_improvement_iteration;

        if (best_solution.fitness > current_best_solution.fitness) {
            best_solution = current_best_solution;       
            current_no_improvement_iteration = 1;
        }
   }

    this->best_solution = best_solution.genes;
}
