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
	Chromosome bestSolution = GeneticAlgorithm::fitness(population[getRandomInt(0, populationSize - 1)]);
	Chromosome solution;
	
	for (size_t i = 0; i < populationSize; ++i) {
		solution = fitness(population[i]);
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
 * @param graph Object to the graph used to initialize the chromosomes.
 */

void GeneticAlgorithm::createPopulation(
    std::vector<std::function<Chromosome(const Graph&)>> generateChromosomeHeuristics, 
    const Graph& graph, 
    size_t heuristic) {

    if (!generateChromosomeHeuristics.empty()) {
        if (heuristic == 4) {
        
            for (size_t i = 0; i <= populationSize - 5; ++i) {
                this->population[i] = generateChromosomeHeuristics[0](graph);
            }

            this->population[populationSize - 4] = generateChromosomeHeuristics[1](graph);
            this->population[populationSize - 3] = generateChromosomeHeuristics[1](graph);
            this->population[populationSize - 2] = generateChromosomeHeuristics[2](graph);
            this->population[populationSize - 1] = generateChromosomeHeuristics[2](graph);
				
			std::random_device randomNumber; 
    		std::mt19937 seed(randomNumber());
    		
			std::shuffle(population.begin(), population.end(), randomNumber);
        } 
        
        else if (heuristic == 1) {
            Chromosome func = generateChromosomeHeuristics[0](graph);
            for (size_t i = 0; i < populationSize; ++i) {
                this->population[i] = func;
            }
        }
        
         else if (heuristic == 2) {
            Chromosome func = generateChromosomeHeuristics[1](graph);
            for (size_t i = 0; i < populationSize; ++i) {
                this->population[i] = func;
            }
        }
        
        else if (heuristic == 3) {
            Chromosome func = generateChromosomeHeuristics[2](graph);
            for (size_t i = 0; i < populationSize; ++i) {
                this->population[i] = func;
            }
        }
        
    } 
    
    else {
    
        for (size_t i = 0; i < populationSize; ++i) {
            this->population[i] = Chromosome(genesSize);
        }
    }
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


Chromosome GeneticAlgorithm::twoPointCrossOver(Chromosome& chromosome1, Chromosome& chromosome2) {
   size_t range1 = GeneticAlgorithm::getRandomInt(0, genesSize - 1);
   size_t range2 = GeneticAlgorithm::getRandomInt(0, genesSize - 1);
   
    Chromosome solution1 = chromosome1;
    Chromosome solution2 = chromosome2;

   if (range1 > range2) 
        std::swap(range1, range2);
    
   for (size_t i = range1; i <= range2; ++i) 
		std::swap(solution1.genes[i], solution2.genes[i]);
        
   feasibilityCheck(solution1);
   feasibilityCheck(solution2);
   
   fitness(solution1);
   fitness(solution2);
    
   return chooseBestSolution(solution1, solution2);
}

Chromosome GeneticAlgorithm::onePointCrossOver(Chromosome& chromosome1, Chromosome& chromosome2) {
  
   size_t index = GeneticAlgorithm::getRandomInt(0, genesSize - 1);
   
   Chromosome solution1 = chromosome1;
   Chromosome solution2 = chromosome2;
    
   for(size_t i = index + 1; i <= genesSize - 1; ++i) 
		std::swap(solution1.genes[i], solution2.genes[i]);
        
   feasibilityCheck(solution1);
   feasibilityCheck(solution2);
   
   fitness(solution1);
   fitness(solution2);
    
   return chooseBestSolution(solution1, solution2);
}

std::vector<Chromosome>& GeneticAlgorithm::elitism(float elitismRate) {
    size_t iterations = static_cast<size_t>(std::ceil(this->populationSize * elitismRate));

    Chromosome bestSolution = this->getBestChromosome(this->population);
    
	this->population.clear();
    this->population.reserve(iterations);
    
    for (size_t i = 0; i < iterations; ++i) 
        this->population.emplace_back(bestSolution);
	
    return this->population;
}


/**
 * @brief Checks the feasibility of a chromosome.
 * 
 * Verify if the constraints of Triple Roman Domination are met.
 * 
 * @param chromosome A reference to the chromosome to be checked.
 * @return True if the constraints of Triple Roman Domination are met, otherwise, False 
 */
 
bool GeneticAlgorithm::feasible(Chromosome& chromosome) {
    bool isValid = false;
    bool hasNeighborAtLeast2 = false;
    
    for (size_t i = 0; i < genesSize; ++i) {

        isValid = false;

        if (chromosome.genes[i] == 0) {                                 
            size_t countNeighbors2 = 0;
            size_t countNeighbors3 = 0;
            
            for (auto& neighbor : this->graph.getAdjacencyList(i)) {          

                if ((countNeighbors2 == 1 && chromosome.genes[neighbor] >= 3) ||
                    (countNeighbors2 == 2 && chromosome.genes[neighbor] >= 2)) {
                    isValid = true;
                    break;
                }

                if (countNeighbors3 == 1 && chromosome.genes[neighbor] >= 2) {
                     isValid = true;
                     break;
                }

                if (chromosome.genes[neighbor] == 4) {
                    isValid = true;
                    break;
                }
                
                if (chromosome.genes[neighbor] == 3) 
                    ++countNeighbors3;
                    
                if (chromosome.genes[neighbor] == 2) 
                    ++countNeighbors2;
            }
            
            if (!isValid)
                return false;
        } 

        else if (chromosome.genes[i] == 2) {
            hasNeighborAtLeast2 = false;
            for (auto& neighbor : this->graph.getAdjacencyList(i)) {
                if (chromosome.genes[neighbor] >= 2) {
                    hasNeighborAtLeast2 = true;
                    break; 
                }
            }

            if (!hasNeighborAtLeast2) 
                return false;           
        }
    }
    
    return true;    
}

/**
 * @brief Checks the feasibility of a chromosome and adjusts if it isn't feasible.
 * 
 * Adjusts genes based on the adjacency list of the graph, ensuring that constraints of Triple Roman Domination are met.
 * 
 * @param chromosome A reference to the chromosome to be checked.
 * @return Chromosome The adjusted chromosome.
 */
 
 
Chromosome GeneticAlgorithm::feasibilityCheck(Chromosome& chromosome) {  
    bool isValid = false;
    bool hasNeighborAtLeast2 = false;
                                                           
    for (size_t i = 0; i < genesSize; ++i) {
    	isValid = false;
    	 
        if (chromosome.genes[i] == 0) {                                 
            size_t countNeighbors2 = 0;
            size_t countNeighbors3 = 0;
            
            for (auto& neighbor : this->graph.getAdjacencyList(i)) {          
                                                                                 
                if ((countNeighbors2 == 1 && chromosome.genes[neighbor] >= 3) ||
                    (countNeighbors2 == 2 && chromosome.genes[neighbor] >= 2)) {
                    isValid = true;
                    break;
                }
                                                                                 
                if (countNeighbors3 == 1 && chromosome.genes[neighbor] >= 2) {
                     isValid = true;
                     break;
                }
                                                                                 
                if (chromosome.genes[neighbor] == 4) {
                    isValid = true;
                    break;
                }
                
                if (chromosome.genes[neighbor] == 3) 
                    ++countNeighbors3;
                    
                if (chromosome.genes[neighbor] == 2) 
                    ++countNeighbors2;
            }
            
            if (!isValid) {
            	if (countNeighbors3 == 0) {
            		if (countNeighbors2 == 0) 
            			chromosome.genes[i] = 3;
            			
            		else if (countNeighbors2 > 0)
            			chromosome.genes[i] = 2;
            	} 
            	
            	else if (countNeighbors3 == 1)
            		chromosome.genes[i] = 2;
            }
        }

        else if (chromosome.genes[i] == 2) {
            hasNeighborAtLeast2 = false;
            for (auto& neighbor : this->graph.getAdjacencyList(i)) {
                if (chromosome.genes[neighbor] >= 2) {
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
 * 
 * @return std::vector<Chromosome>& A new population of Chromosomes.
 */

std::vector<Chromosome>& GeneticAlgorithm::createNewPopulation() {
    std::vector<Chromosome> temp = population;
   	
   	this->elitism(this->elitismRate);
    
    temp.reserve(populationSize);
    
    Chromosome selected1;
    Chromosome selected2;
    Chromosome offspring;
    
    while (population.size() < populationSize) {          
        selected1 = tournamentSelection(temp); 
        selected2 = tournamentSelection(temp);      	
       	
       	if (getRandomFloat(0.0, 1.0) <= crossOverRate)   		
        	offspring = this->twoPointCrossOver(selected1, selected2); 
        else      
        	offspring = this->onePointCrossOver(selected1, selected2);     
        
	    mutation(offspring);	    
		
        population.emplace_back(offspring);       
    }
	
    population.swap(temp);
    
    return population;
}

void GeneticAlgorithm::run(size_t generations, std::vector<std::function<Chromosome(const Graph&)>> heuristics, size_t chosenHeuristic) { 
   this->createPopulation(heuristics, graph, chosenHeuristic);
	
   Chromosome currentBestSolution = this->tournamentSelection(this->population);                                         
   Chromosome bestSolution = currentBestSolution;
   size_t iteration = 0;
   size_t currentNoImprovementIteration = 0;
   
   while ((iteration < generations) && (currentNoImprovementIteration < maxNoImprovementIterations)) {
   
		this->population.swap(this->createNewPopulation());
       	
        currentBestSolution =  this->tournamentSelection(this->population);                                       
		
        ++iteration;
        ++currentNoImprovementIteration;

        if (bestSolution.fitnessValue > currentBestSolution.fitnessValue) {
            bestSolution = currentBestSolution;       
            currentNoImprovementIteration = 1;
        }
   }

    this->bestSolution = bestSolution.genes;
}
