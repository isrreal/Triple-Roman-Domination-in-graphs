#include "GeneticAlgorithm.hpp"
 
std::vector<Chromosome> GeneticAlgorithm::getPopulation() { return this->population; }

size_t GeneticAlgorithm::getPopulationSize() { return this->populationSize; }

size_t GeneticAlgorithm::getGenesSize() { return this->genesSize; }

size_t GeneticAlgorithm::getGenerations() { return this->generations; }

double GeneticAlgorithm::getMutationRate() { return this->mutationRate; }

double GeneticAlgorithm::getElitismRate() { return this->elitismRate; }

Graph GeneticAlgorithm::getGraph() { return this->graph; }

std::vector<int> GeneticAlgorithm::getBestSolution() { return this->bestSolution; } 

Chromosome GeneticAlgorithm::getBestChromosome(std::vector<Chromosome> population) {
	std::random_device randomNumber;
	std::mt19937 seed(randomNumber());
	std::uniform_int_distribution<int> gap(0, populationSize - 1);
	
	Chromosome bestSolution = GeneticAlgorithm::fitness(population[gap(seed)], nullptr);
	Chromosome solution;
	
	for (size_t i = 0; i < populationSize; ++i) {
		solution = fitness(population[i], nullptr);
		if (solution.fitnessValue > bestSolution.fitnessValue)
			bestSolution = solution;
	}

	return bestSolution;
}


/**
 * @brief Creates a population of chromosomes with a specific number of genes.
 * 
 * If a heuristic function is provided, the chromosomes are initialized using this heuristic.
 * Otherwise, chromosomes are initialized with random genes.
 * 
 * @param heuristic A pointer to a function that generates chromosomes based on a graph.
 * @param graph A pointer to the graph used to initialize the chromosomes.
 */

void GeneticAlgorithm::createPopulation(Chromosome(*generateChromosomeHeuristic)(Graph), Graph graph) {
    if (generateChromosomeHeuristic) {  
       Chromosome func = (*generateChromosomeHeuristic)(graph);  
       for (size_t i = 0; i < populationSize; ++i) {
            this->population[i] = Chromosome(func);
            this->population[i].indexRemove = i;
       }
   } 
        
   else {
       for (size_t i = 0; i < populationSize; ++i) {
           this->population[i] = Chromosome(genesSize); 
           this->population[i].indexRemove = i;
       }
   }
}

/**
 * @brief Calculates the fitness score for a given chromosome.
 * 
 * If a fitness heuristic is provided, it is used to calculate the fitness score.
 * Otherwise, returns nullptr.
 * 
 * @param chromosome Reference to the chromosome to evaluate.
 * @param fitnessHeuristic A pointer to a function that evaluates the fitness of the chromosome.
 * @return Chromosome The chromosome with its fitness score calculated, or nullptr if no heuristic is provided.
 */

Chromosome GeneticAlgorithm::fitness(Chromosome& chromosome, Chromosome(*fitnessHeuristic)(Chromosome&) = nullptr) { 
    if (fitnessHeuristic)
        return (*fitnessHeuristic)(chromosome);

    for (size_t i = 0; i < chromosome.genesSize; ++i)
        chromosome.fitnessValue += chromosome.genes[i];
    return chromosome;
}

/**
 * @brief Selects the chromosome from the population using tournament selection.
 * 
 * This method randomly selects two chromosomes from the population, evaluates their fitness, 
 * and choses a random number r between 0 and 1. If r < k (where k is a parameters, for example 0.75),
 * fitter of the two parameters is selected to be a parent; otherwise, the one with the higher fitness value is selected.
 * 
 * @param population The vector of chromosomes in the current population.
 * @return Chromosome The chromosome with the highest or lowest fitness value.
 */

Chromosome GeneticAlgorithm::tournamentSelection(std::vector<Chromosome> population) { 
    constexpr float parameter = 0.75f; 
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_int_distribution<int> gap(0, population.size() - 1);
    std::uniform_real_distribution<float> probability(0, 1); 
    
    
    Chromosome c1 = GeneticAlgorithm::fitness(population[gap(seed)], nullptr);
    Chromosome c2 = GeneticAlgorithm::fitness(population[gap(seed)], nullptr);
   
    if (probability(seed) < parameter) 
       return GeneticAlgorithm::chooseBestSolution(c1, c2);
    else 
       return GeneticAlgorithm::chooseWorstSolution(c1, c2); 
}

/**
 * @brief Selects a chromosome from the population using roulette wheel selection.
 * 
 * This method randomly selects and returns a chromosome from the population.
 * 
 * @param population The vector of chromosomes in the current population.
 * @return Chromosome The randomly selected chromosome.
 */

Chromosome GeneticAlgorithm::rouletteWheelSelection(std::vector<Chromosome> population) {
    size_t totalFitness = 0;

    for (size_t i = 0; i < population.size(); ++i) {
        fitness(population[i], nullptr);
        totalFitness += population[i].fitnessValue;
    }

    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_int_distribution<> distribution(0, totalFitness - 1);
    size_t randomValue = distribution(seed);

    size_t cumulativeFitness = 0;
    for (auto& chromosome: population) {
        cumulativeFitness += chromosome.fitnessValue;
        if (cumulativeFitness >= randomValue) 
            return chromosome;
    }
    
    return population.back();
}



/**
 * @brief Selects a chromosome from the population using a specific selection heuristic and removes it.
 * 
 * Uses the provided selection heuristic to select the best chromosome from the population and removes it
 * to prevent duplications in future selections.
 * 
 * @param selectionHeuristic A pointer to a function that implements the selection heuristic.
 * @return Chromosome The selected chromosome, or defualt Chromosome object if no valid selection is made.
 */

Chromosome GeneticAlgorithm::selectionMethod(Chromosome(*selectionHeuristic)(std::vector<Chromosome>)) {
    if (!selectionHeuristic) 
        return Chromosome(); 
  
    Chromosome selected = (*selectionHeuristic)(this->population);
    
    if (selected.indexRemove >= 0 && selected.indexRemove < this->population.size()) 
    	this->population.erase(this->population.begin() + selected.indexRemove);

    return selected; 
}


/**
 * @brief Chooses the best Chromosome between two options based on their fitness values.
 * 
 * @param chromosome1 constant reference to the first Chromosome.
 * @param chromosome2 constant reference  to the second Chromosome.
 * @return Chromosome The Chromosome with the higher fitness value.
 */

Chromosome GeneticAlgorithm::chooseBestSolution(const Chromosome& chromosome1, const Chromosome& chromosome2) {
	    return (chromosome1.fitnessValue > chromosome2.fitnessValue ? chromosome1 : chromosome2);
}

/**
 * @brief Chooses the worst Chromosome between two options based on their fitness values.
 * 
 * @param chromosome1 constant reference to the first Chromosome.
 * @param chromosome2 constant reference to the second Chromosome.
 * @return Chromosome The Chromosome with the lower fitness value.
 */


Chromosome GeneticAlgorithm::chooseWorstSolution(const Chromosome& chromosome1, const Chromosome& chromosome2) {
	    return (chromosome1.fitnessValue < chromosome2.fitnessValue ? chromosome1 : chromosome2);
}

/**
 * @brief Performs a crossover between two chromosomes to create an offspring.
 * 
 * Combines genes from two parent chromosomes to create a new chromosome with genes 
 * inherited from both parents.
 * 
 * @param chromosome1 Reference to the first parent chromosome.
 * @param chromosome2 Reference to the second parent chromosome.
 * @param crossOverHeuristic Optional pointer to a function that performs the crossover.
 * @return Chromosome A new chromosome offspring.
 */


Chromosome GeneticAlgorithm::crossOver(Chromosome& chromosome1, Chromosome& chromosome2,
 	Chromosome(*crossOverHeuristic)(Chromosome&, Chromosome&)) {
 	
   if (crossOverHeuristic)
        return (*crossOverHeuristic)(chromosome1, chromosome2);
     
   std::random_device randomNumber;
   std::mt19937 seed(randomNumber()); 
   std::uniform_int_distribution<> gap(0, genesSize - 1);
   size_t range1 = gap(seed);
   size_t range2 = gap(seed);
   
   std::vector<int> x, y;

   if (range1 > range2) 
        std::swap(range1, range2);
    
   for (size_t i = range1; i <= range2; ++i) { 
        x.push_back(chromosome1.genes[i]);
        y.push_back(chromosome2.genes[i]);
   }

   for (size_t i = range1, j = 0; i <= range2; ++i, ++j) {  
       chromosome1.genes[i] = y[j];
       chromosome2.genes[i] = x[j]; 
   }
 
   Chromosome solution1 = chromosome1;
   Chromosome solution2 = chromosome2;

   feasibilityCheck(solution1);
   feasibilityCheck(solution2);
    
   return chooseBestSolution(solution1, solution2);
}


Chromosome GeneticAlgorithm::mutation(Chromosome& chromosome) {
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber()); 
    std::uniform_real_distribution<> gap(0.0, 1.0);
    
	double probability = 0.0;
	
    for (size_t i = 0; i < chromosome.genes.size(); ++i) {
        probability = gap(seed);
        if (probability <= this->mutationRate) {
        	if (chromosome.genes[i] == 0)
        		chromosome.genes[i] = 2;
        		
        	else if (chromosome.genes[i] == 3)
        		chromosome.genes[i] = 2;
        	
        	feasibilityCheck(chromosome);
        }

    }

    return  chromosome;
}


std::vector<Chromosome>& GeneticAlgorithm::elitism(float elitismRate) {
    std::vector<Chromosome> temp;

    Chromosome bestSolution = this->getBestChromosome(this->population);
    
    size_t iterations = static_cast<size_t>(this->populationSize * elitismRate);
	
    for (size_t i = 0; i < iterations; ++i)
        temp.push_back(bestSolution);
		
	population.swap(temp);
	
    return population; 
}


/**
 * @brief Checks the feasibility of a chromosome.
 * 
 * Adjusts genes based on the adjacency list of the graph, ensuring that constraints of Double Roman Domination are met.
 * 
 * @param chromosome A reference to the chromosome to be checked.
 * @return Chromosome The adjusted chromosome.
 */
 
 
Chromosome GeneticAlgorithm::feasibilityCheck(Chromosome& chromosome) {	
    bool hasNeighborWith4 = false; 
    bool hasTwoNeighbors3or2[] = {false, false}; 
    bool hasThreeNeighbors2[] = {false, false, false};
    bool hasNeighborAtLeast2 = false;

    for (size_t i = 0; i < genesSize; ++i) {
        if (chromosome.genes[i] == 0) {
            for (auto& neighbor : this->graph.getAdjacencyList(i)) {
                if (chromosome.genes[neighbor] == 4) {
                    hasNeighborWith4 = true;
                    break; 
                }

                else if (chromosome.genes[neighbor] == 2) {
                    if (chromosome.genes[neighbor] == 3) {
                        hasTwoNeighbors3or2[0], hasTwoNeighbors3or2[1] = true;
                        break;
                    }

                    else if (chromosome.genes[neighbor] == 2) {
                        if (chromosome.genes[neighbor] == 2) {
                            hasThreeNeighbors2[0], hasThreeNeighbors2[0], hasThreeNeighbors2[0] = true;
                            break;
                        } 
                    } 
                }
                
                if (!hasNeighborWith4 && !hasTwoNeighbors3or2[0] && !hasTwoNeighbors3or2[1]
                        && !hasThreeNeighbors2[0] && !hasThreeNeighbors2[1] && !hasThreeNeighbors2[2])
                    chromosome.genes[i] = 4;
            }
        }

        else if (chromosome.genes[i] == 2) {
            for (auto& neighbor : this->graph.getAdjacencyList(i)) {
                if (chromosome.genes[neighbor] == 4) {
                    hasNeighborAtLeast2 = true;
                    break; 
                }
            }

            if (!hasNeighborAtLeast2) 
                chromosome.genes[i] = 3;
        }
    }

    return chromosome;	
}

/**
 * @brief Generates a new population by crossing over Chromosomes from the current population.
 * 
 *  Produces a new genetically superior population.
 * 
 * @return std::vector<Chromosome> A new population of Chromosomes.
 */

std::vector<Chromosome>& GeneticAlgorithm::createNewPopulation() {
    this->population = this->elitism(this->elitismRate);
    
    std::vector<Chromosome> temp;
      
    while (temp.size() < populationSize) {
        Chromosome selected1 = this->selectionMethod(tournamentSelection);
        Chromosome selected2 = this->selectionMethod(rouletteWheelSelection);
        Chromosome offspring = this->crossOver(selected1, selected2, nullptr);

        temp.push_back(offspring);
    }
    
    population.swap(temp);
    
    return population;
}

/**
 * @brief Runs the genetic algorithm for a specified number of generations.
 * 
 * @param generations Number of generations to evolve.
 * @param heuristic Pointer to a function that generates initial Chromosomes from a graph.
 * @param graph Object to the graph used to generate initial solutions.
 * @return Chromosome The best solution found after all generations.
 */

void GeneticAlgorithm::run(size_t generations, Chromosome(*heuristic)(Graph)) { 

   this->createPopulation(heuristic, graph);

   Chromosome currentBestSolution = this->tournamentSelection(this->population);                                         
   Chromosome bestSolution = currentBestSolution;

   for (size_t i = 0; i < generations; ++i) {        
   
		this->population.swap(this->createNewPopulation());
       	
        currentBestSolution = this->tournamentSelection(this->population);                                       

        if (bestSolution.fitnessValue > currentBestSolution.fitnessValue)
            bestSolution = currentBestSolution; 
   }

    this->bestSolution = bestSolution.genes;
}
