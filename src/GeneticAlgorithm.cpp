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
	
		Chromosome func;
		
        if (heuristic == 4) { 
            size_t portion_size { static_cast<size_t>(population_size / 3) };
			size_t remainder { population_size % 3 }; 
			size_t index {0};
			
			for (size_t i {0}; i < 3; ++i) {		
		        func = generateChromosomeHeuristics[i](graph);
		        
				for (size_t j {0}; j < portion_size; ++j) {
					this->population[index++] = func;
				}
			}
			
			for (size_t i {0}; i < remainder; ++i) {
				this->population[index++] = generateChromosomeHeuristics[getRandomInt(0, 2)](graph);
			}
			
			std::random_device rd;
			std::mt19937 random_engine(rd());
			std::shuffle(population.begin(), population.end(), random_engine);	

        } 
        
        else if (heuristic > 0 && heuristic < 4) {
        
            func = generateChromosomeHeuristics[heuristic - 1](graph);
            
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
	std::vector<Chromosome> old_population = population;
   	this->elitism(population, elitism_rate);
    
    Chromosome selected1;
    Chromosome selected2;
    Chromosome offspring;
    
    while (population.size() < population_size) {          
        selected1 = tournamentSelection(old_population, tournament_population_size); 
        
        selected2 = tournamentSelection(old_population, tournament_population_size);      

       	if (getRandomFloat(0.0, 1.0) <= crossover_rate) {  		
        	offspring = this->twoPointCrossOver(selected1, selected2);
		} 
        
		else {      
        	offspring = this->onePointCrossOver(selected1, selected2);     
    	}

	    mutation(offspring);	    
    	
        population.emplace_back(offspring);       
    }
    
    return population;
}

void GeneticAlgorithm::elitism(std::vector<Chromosome>& population, float elitism_rate) {
    size_t iterations { static_cast<size_t>(std::ceil(population.size() * elitism_rate)) };
	
	std::sort(population.begin(), population.end(), 
            [](const Chromosome& a, const Chromosome& b){
            	return a.fitness < b.fitness;
	});
	
	population.resize(iterations);
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

Chromosome GeneticAlgorithm::onePointCrossOver(const Chromosome& chromosome1, const Chromosome& chromosome2) {
  
   size_t index { getRandomInt(0, genes_size - 1) };
   
   Chromosome solution1 { chromosome1 };
   Chromosome solution2 { chromosome2 };
    
   for(size_t i {index + 1}; i <= genes_size - 1; ++i) {
		std::swap(solution1.genes[i], solution2.genes[i]);
   }
   
   feasibilityCheck(this->graph, solution1);
   feasibilityCheck(this->graph, solution2);
    
   return chooseBestSolution(solution1, solution2);
}

Chromosome GeneticAlgorithm::twoPointCrossOver(const Chromosome& chromosome1, const Chromosome& chromosome2) {
	size_t range1 { getRandomInt(0, genes_size - 1) };
	size_t range2 { getRandomInt(0, genes_size - 1) };

	Chromosome solution1 { chromosome1 };
	Chromosome solution2 { chromosome2 };

	if (range1 > range2) {
		std::swap(range1, range2);
	}

	for (size_t i {range1}; i <= range2; ++i) {
		std::swap(solution1.genes[i], solution2.genes[i]);
	}
	
	feasibilityCheck(this->graph, solution1);
	feasibilityCheck(this->graph, solution2);

	return chooseBestSolution(solution1, solution2);
}

const Chromosome& GeneticAlgorithm::tournamentSelection(const std::vector<Chromosome>& population, size_t individuals_size) {  
    std::vector<size_t> tournament_individuals_indices; 
    tournament_individuals_indices.reserve(individuals_size);
    size_t random_index {0};
    
    for (size_t i {0}; i < individuals_size; ++i) {
        random_index = getRandomInt(0, population.size() - 1);
        tournament_individuals_indices.push_back(random_index);
    }
	
    size_t best_index = tournament_individuals_indices[0];
	
    for (auto index : tournament_individuals_indices) {
        if (population[index].fitness < population[best_index].fitness) {
            best_index = index;
        }
    }
		
    return population[best_index];
}

Chromosome& GeneticAlgorithm::chooseBestSolution(Chromosome& chromosome1, Chromosome& chromosome2) {
    return (chromosome1.fitness < chromosome2.fitness ? chromosome1 : chromosome2);
}

Chromosome GeneticAlgorithm::findBestSolution(const std::vector<Chromosome>& population) { 
	Chromosome best_solution { population[0] };

	for (const auto& chromosome: population) {
		if (best_solution.fitness > chromosome.fitness) {
			best_solution = chromosome;
		}
	}

	return best_solution; 
}
		
// public methods 

size_t GeneticAlgorithm::getGenerations() { return generations; }

std::vector<int> GeneticAlgorithm::getBestSolution() { return best_solution; }

void GeneticAlgorithm::run(size_t generations, std::vector<std::function<Chromosome(const Graph&)>> heuristics, size_t chosen_heuristic) { 
   	this->createPopulation(heuristics, graph, chosen_heuristic);
   
	Chromosome current_best_solution { findBestSolution(population) };  
   
   	Chromosome best_solution { current_best_solution };
   
   	size_t iteration {0};
   	size_t current_no_improvement_iteration {0};
   
  	while ((iteration < generations) && (current_no_improvement_iteration < max_no_improvement_iterations)) {
   
		this->population.swap(createNewPopulation());
		
        current_best_solution = findBestSolution(population);  
		
        ++iteration;
        ++current_no_improvement_iteration;

        if (best_solution.fitness > current_best_solution.fitness) {
            best_solution = current_best_solution;       
            current_no_improvement_iteration = 1;
        }
        
   	}

    this->best_solution.swap(best_solution.genes);
}
