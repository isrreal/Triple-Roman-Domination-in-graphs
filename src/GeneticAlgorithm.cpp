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

size_t GeneticAlgorithm::chooseVertex(Graph& temp) {
    float selectionVertexRateConstructSolution = 0.7f;    
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_real_distribution<float> probabilityGap(0.0, 1.0);
    std::uniform_int_distribution<int> gap(0, this->graph.getOrder() - 1);
    
    float number = probabilityGap(seed); 
    
    size_t vertex = gap(seed);
    float value = 0.0;
    
    while (!temp.vertexExists(vertex)) 
        vertex = gap(seed);

    if (selectionVertexRateConstructSolution < number) 
        vertex = rouletteWheelSelection(temp);
        
    return vertex;
}

// returns a index of vertex randomly selected
size_t GeneticAlgorithm::chooseVertex(std::vector<int> twoOrZeroOrThreeLabeledVertices) {
    constexpr float selectionVertexRateExtendSolution = 0.9f;

    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_real_distribution<float> probabilityGap(0.0, 1.0);
    std::uniform_int_distribution<int> randomIndex(0, twoOrZeroOrThreeLabeledVertices.size() - 1);                                                               

    float randomNumberGenerated = probabilityGap(seed); 

    size_t choosenIndex = randomIndex(seed);

    if (randomNumberGenerated > selectionVertexRateExtendSolution)
    	choosenIndex = rouletteWheelSelection(twoOrZeroOrThreeLabeledVertices);
    	
    return choosenIndex;
}


Chromosome GeneticAlgorithm::destroySolution(Chromosome& chromosome) {
    float destructionRate = minDestructionRate + ((currentRVNSnumber - 1) *
                ((maxDestructionRate - minDestructionRate)) 
                / (maxRVNSfunctions - 1));
    Graph temp = this->graph;
    size_t itr = genesSize * destructionRate; 
    size_t vertex = 0;
    
    while (itr != 0 && (temp.getOrder() > 0)) {
       vertex = chooseVertex(temp);
       if ((chromosome.genes[vertex] == 0) || 
       	  (chromosome.genes[vertex] == 2)||
       	  (chromosome.genes[vertex] == 3))
            	chromosome.genes[vertex] = -1;
       else
            ++itr;
       temp.deleteVertex(vertex);
       --itr;
    }

    return chromosome;
}

Chromosome GeneticAlgorithm::extendSolution(Chromosome& chromosome) {
    constexpr float addVerticesRate = 0.05f;
    size_t itr = 0;
    std::vector<int> twoOrZeroOrThreeLabeledVertices;
    
    for (size_t i = 0; i < genesSize; ++i) {
        if ((chromosome.genes[i] == 0) || (chromosome.genes[i] == 2) || (chromosome.genes[i] == 3))
            twoOrZeroOrThreeLabeledVertices.push_back(i);
    }

    itr = addVerticesRate * twoOrZeroOrThreeLabeledVertices.size();
    
    size_t vertex = 0;

    while (itr != 0 && !twoOrZeroOrThreeLabeledVertices.empty()) {          
        vertex = chooseVertex(twoOrZeroOrThreeLabeledVertices);         
        chromosome.genes[vertex] = 4; 
        twoOrZeroOrThreeLabeledVertices.erase(twoOrZeroOrThreeLabeledVertices.begin() + vertex);                                                                
        --itr;
    }

    return chromosome;
}



Chromosome GeneticAlgorithm::reduceSolution(Chromosome& chromosome) {
    Graph temp = this->graph;
    std::vector<int> sortedVertices;
    int initLabel = -1;

    for (size_t i = 0; i < temp.getOrder(); ++i)
        sortedVertices.push_back(i);

    std::sort(sortedVertices.begin(), sortedVertices.end(),
        [&](size_t a, size_t b) {
            return temp.getVertexDegree(a) < temp.getVertexDegree(b);                                                                             
        });
    
    size_t choosenVertex = 0;
    
    while ((temp.getOrder() > 0) && (choosenVertex < sortedVertices.size())) {
        if (choosenVertex >= sortedVertices.size()) break;
        
        while (choosenVertex < sortedVertices.size() && 
                (!temp.vertexExists(sortedVertices[choosenVertex]))) {
            ++choosenVertex;
        }
       
        if (choosenVertex >= sortedVertices.size()) break;

        if (chromosome.genes[sortedVertices[choosenVertex]] == 4 || chromosome.genes[sortedVertices[choosenVertex]] == 3) {         
    		initLabel = chromosome.genes[sortedVertices[choosenVertex]];
    		chromosome.genes[sortedVertices[choosenVertex]] = 0;
    		
    		if (!feasible(chromosome)) {
        		chromosome.genes[sortedVertices[choosenVertex]] = 2;
        		
        		if (!feasible(chromosome)) {
        			chromosome.genes[sortedVertices[choosenVertex]] = 3;
            		
            		if (!feasible(chromosome)) 
                		chromosome.genes[sortedVertices[choosenVertex]] = initLabel;        
        		}
    		}
		}
            	
        temp.deleteAdjacencyList(sortedVertices[choosenVertex++]);
    }

    return chromosome;
}

Chromosome GeneticAlgorithm::RVNS(Chromosome& chromosome, Chromosome(*heuristic)(Graph)) {
	size_t currentNoImprovementIteration = 0;
	size_t currentRVNSnumber = 1;
    Chromosome temp = chromosome;
    
    while ((currentNoImprovementIteration < maxRVNSnoImprovementIterations) && (maxRVNSiterations > 0)) {
        destroySolution(temp);
        
        temp = (*heuristic)(graph);
        
        extendSolution(temp);
        reduceSolution(temp);

        if (temp.fitnessValue < chromosome.fitnessValue) {
            chromosome = temp;
            currentRVNSnumber = 1;
            currentNoImprovementIteration = 0;
        }

        else {
            ++currentRVNSnumber;
            ++currentNoImprovementIteration;

            if (currentRVNSnumber > maxRVNSiterations)
                currentRVNSnumber = 1;
        }

        --maxRVNSiterations;
    }	

    return chromosome;
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
    std::uniform_int_distribution<> randomLabel(0, 4);

	double probability = 0.0;
	
    for (size_t i = 0; i < chromosome.genes.size(); ++i) {
        probability = gap(seed);
        if (probability <= this->mutationRate) {
            chromosome.genes[i] = randomLabel(seed); 
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
 * Verify if the constraints of Triple Roman Domination are met.
 * 
 * @param chromosome A reference to the chromosome to be checked.
 * @return True if the constraints of Triple Roman Domination are met, otherwise, False 
 */
 
bool GeneticAlgorithm::feasible(Chromosome& chromosome) {
    bool isValid = false;

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
            bool hasNeighborAtLeast2 = false;
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
                                                                                 
    for (size_t i = 0; i < genesSize; ++i) {
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
                    isValid == true;
                    break;
                }
                
                if (chromosome.genes[neighbor] == 3) 
                    ++countNeighbors3;
                    
                if (chromosome.genes[neighbor] == 2) 
                    ++countNeighbors2;
            }
            
            if (!isValid) {
            	if (countNeighbors3 == 0 && countNeighbors2 == 0)
                	chromosome.genes[i] = 3; 
                else if (countNeighbors3 == 1 || countNeighbors2 == 2) 
                    chromosome.genes[i] = 2;
            }
        }

        else if (chromosome.genes[i] == 2) {
            bool hasNeighborAtLeast2 = false;
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
 * @return std::vector<Chromosome> A new population of Chromosomes.
 */

std::vector<Chromosome>& GeneticAlgorithm::createNewPopulation1() {
    this->population = this->elitism(this->elitismRate);
    
    std::vector<Chromosome> temp;
      
    while (temp.size() < populationSize) {
        Chromosome selected1 = this->selectionMethod(tournamentSelection);
        Chromosome selected2 = this->selectionMethod(rouletteWheelSelection);
        Chromosome offspring = this->crossOver(selected1, selected2, nullptr);
        
		mutation(offspring);
		
        temp.push_back(offspring);       
    }
    
    population.swap(temp);
    
    return population;
}

/**
 * @brief Generates a new population by crossing over Chromosomes from the current population.
 *     
 * 
 * @return std::vector<Chromosome> A new population of Chromosomes.
 */

std::vector<Chromosome>& GeneticAlgorithm::createNewPopulation2() {
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

void GeneticAlgorithm::run1(size_t generations, Chromosome(*heuristic)(Graph)) { 

   this->createPopulation(heuristic, graph);
	
    
   Chromosome currentBestSolution = this->tournamentSelection(this->population);                                         
   Chromosome bestSolution = currentBestSolution;

   for (size_t i = 0; i < generations; ++i) {        
   
		this->population.swap(this->createNewPopulation1());
       	
        currentBestSolution = this->tournamentSelection(this->population);                                       
		
        if (bestSolution.fitnessValue > currentBestSolution.fitnessValue)
            bestSolution = currentBestSolution; 
        
         RVNS(bestSolution, heuristic);
   }

    this->bestSolution = bestSolution.genes;
}


void GeneticAlgorithm::run2(size_t generations, Chromosome(*heuristic)(Graph)) { 

   this->createPopulation(heuristic, graph);
	
    
   Chromosome currentBestSolution = this->tournamentSelection(this->population);                                         
   Chromosome bestSolution = currentBestSolution;

   for (size_t i = 0; i < generations; ++i) {        
   
		this->population.swap(this->createNewPopulation2());
       	
        currentBestSolution = this->tournamentSelection(this->population);                                       
		
        if (bestSolution.fitnessValue > currentBestSolution.fitnessValue)
            bestSolution = currentBestSolution;       
   }

    this->bestSolution = bestSolution.genes;
}

// selects a vertex of graph by roulette wheel selection method embiased by vertex degree and total fitness

size_t GeneticAlgorithm::rouletteWheelSelection(Graph& temp) {
    float totalFitness = 0.0f;
    std::vector<std::pair<size_t, float>> probabilities;

    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_real_distribution<float> gap(0.0, 1.0);

    for (size_t i = 0; i < graph.getOrder(); ++i) 
        if (temp.vertexExists(i)) 
            totalFitness += population[i].fitnessValue;

    for (size_t i = 0; i < graph.getOrder(); ++i)
        if (temp.vertexExists(i))
            probabilities.push_back({i, (temp.getVertexDegree(i) / totalFitness)});
    
    float randomValue = gap(seed);

    float cumulativeSum = 0.0f;
    for (const auto& prob : probabilities) {
        cumulativeSum += prob.second;
        if (randomValue <= cumulativeSum) 
            return prob.first;
    }

    return probabilities.back().first;
}

// selects a vertex of graph that has 0, 2 or 3 label, by roulette wheel selection method embiased by vertex degree and total fitness
size_t GeneticAlgorithm::rouletteWheelSelection(std::vector<int> twoOrZeroOrThreeLabeledVertices) {
    float totalFitness = 0.0f;
    std::vector<std::pair<size_t, float>> probabilities;
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_real_distribution<float> gap(0.0, 1.0);

    for (size_t i = 0; i < twoOrZeroOrThreeLabeledVertices.size(); ++i) 
        totalFitness += population[i].fitnessValue;
    
    for (size_t i = 0; i < twoOrZeroOrThreeLabeledVertices.size(); ++i)
        probabilities.push_back({i, (graph.getVertexDegree(twoOrZeroOrThreeLabeledVertices[i]) / totalFitness)});

    float randomValue = gap(seed);
                                             
    float cumulativeSum = 0.0f;
    for (const auto& prob : probabilities) {
        cumulativeSum += prob.second;
        if (randomValue <= cumulativeSum) 
            return prob.first;
    }
                                             
    return probabilities.back().first;
}
