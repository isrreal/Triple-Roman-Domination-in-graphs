#ifndef GENETIC_ALGORITHM_HPP
#define GENETIC_ALGORITHM_HPP

#include <iostream>
#include <random>
#include <vector>
#include "Chromosome.hpp"

class GeneticAlgorithm {
	private:
		size_t populationSize;
		size_t genesSize;
		std::vector<Chromosome> population;
	    size_t generations;
        Graph graph;
        std::vector<int> bestSolution;
        
        float mutationRate;
        float elitismRate;
        float crossOverRate;
		size_t maxNoImprovementIterations;

		inline void createPopulation(std::vector<std::function<Chromosome(const Graph&)>> generateChromosomeHeuristics,
		 	const Graph& graph, size_t heuristic);
	
        inline std::vector<Chromosome>& createNewPopulation();    
        
		inline Chromosome onePointCrossOver(Chromosome& chromosome1, Chromosome& cromossomo2); 
                	
    	inline Chromosome twoPointCrossOver(Chromosome& chromosome1, Chromosome& cromossomo2);
                	
        inline Chromosome& mutation(Chromosome& chromosome) {
        	std::vector<int> labels = {0, 2, 3, 4};
        	
			if (getRandomFloat(0.0, 1.0) <= this->mutationRate) {
				size_t randomIndex = getRandomInt(0, genesSize - 1);
				short randomLabel = getRandomInt(0, labels.size() - 1);
				std::cout << labels[randomLabel] << std::endl;
				chromosome.genes[randomIndex] = labels[randomLabel];
    			feasibilityCheck(chromosome);
    		}
    		
			return chromosome;   
        }
        
        inline std::vector<Chromosome>& elitism(float elitismRate);
        
        bool feasible(Chromosome& chromosome);
                
		Chromosome feasibilityCheck(Chromosome& chromosome);

        Chromosome getBestChromosome(std::vector<Chromosome> population);
        
        inline static Chromosome fitness(Chromosome& chromosome) {
        	for (size_t i = 0; i < chromosome.genesSize; ++i) {
        		chromosome.fitnessValue += chromosome.genes[i];
        	}
        		
   		 	return chromosome; 
        }
        
		inline static Chromosome tournamentSelection(const std::vector<Chromosome>& population) {
			constexpr float parameter = 0.75f; 
				
			Chromosome c1 = population[getRandomInt(0, population.size() - 1)]; 
			Chromosome c2 = population[getRandomInt(0, population.size() - 1)];
		   	
		   	float probability = GeneticAlgorithm::getRandomFloat(0.0, 1.0);
			if (probability < parameter) 
			   return GeneticAlgorithm::chooseBestSolution(c1, c2);
			else 
			   return GeneticAlgorithm::chooseWorstSolution(c1, c2); 
		}
		
		static Chromosome chooseBestSolution(const Chromosome& chromosome1, const Chromosome& chromosome2);
        static Chromosome chooseWorstSolution(const Chromosome& chromosome1, const Chromosome& chromosome2);
        
        inline static int getRandomInt(int start, int end) {
        	std::random_device randomNumber; 
    		std::mt19937 seed(randomNumber()); 
    		std::uniform_int_distribution<> gap(start, end); 
    
    		return gap(seed);  
        }
        
        inline static float getRandomFloat(float start, float end) {
        	std::random_device randomNumber; 
    		std::mt19937 seed(randomNumber()); 
    		std::uniform_real_distribution<> gap(start, end); 
    
    		return gap(seed);
    	}

	public:
		GeneticAlgorithm(Graph& graph, size_t populationSize, size_t genesSize, size_t generations,
			float mutationRate, float elitismRate, float crossOverRate):
			  populationSize(populationSize), genesSize(genesSize), 
			  population(populationSize), generations(generations), 
			  graph(graph), bestSolution(), 
			  mutationRate(mutationRate), elitismRate(elitismRate),
			  crossOverRate(crossOverRate),
              maxNoImprovementIterations(100) {}               

		~GeneticAlgorithm() {}
		
        Graph getGraph();
		std::vector<Chromosome> getPopulation();     
		size_t getPopulationSize();    
		size_t getGenesSize();
		size_t getGenerations();   
		double getMutationRate();
		double getElitismRate();
        std::vector<int> getBestSolution();		      
        void run(size_t generations, std::vector<std::function<Chromosome(const Graph&)>>, size_t chosenHeuristic);
};	

#endif
