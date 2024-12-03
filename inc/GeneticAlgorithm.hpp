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
        
     	inline std::vector<Chromosome>& elitism(float elitismRate);	
        
		inline Chromosome& onePointCrossOver(const Chromosome& chromosome1, const Chromosome& cromossomo2); 
                	
    	inline Chromosome& twoPointCrossOver(const Chromosome& chromosome1, const Chromosome& cromossomo2);
                	
        inline Chromosome& mutation(Chromosome& chromosome);
        
        bool feasible(const Chromosome& chromosome);
        
		inline static Chromosome& tournamentSelection(const std::vector<Chromosome>& population);
		
		static Chromosome& chooseBestSolution(Chromosome& chromosome1, Chromosome& chromosome2);
        static Chromosome& chooseWorstSolution(Chromosome& chromosome1, Chromosome& chromosome2);

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
		
        static Chromosome& fitness(Chromosome& chromosome);
        
        static size_t getRandomInt(size_t start, size_t end);
        
        static float getRandomFloat(float start, float end);
        
		static Chromosome& feasibilityCheck(const Graph& graph, Chromosome& chromosome);
		
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
