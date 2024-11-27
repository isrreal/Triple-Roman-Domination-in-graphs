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

		void createPopulation(Chromosome(*heuristic)(Graph), Graph graph);
	
        inline std::vector<Chromosome>& createNewPopulation();    
		inline Chromosome onePointCrossOver(Chromosome& chromosome1, Chromosome& cromossomo2,
                	Chromosome(*crossOverHeuristic)(Chromosome&, Chromosome&)); 
                	
    	inline Chromosome twoPointCrossOver(Chromosome& chromosome1, Chromosome& cromossomo2,
        			Chromosome(*crossOverHeuristic)(Chromosome&, Chromosome&));
                	
        Chromosome& mutation(Chromosome& chromosome);
        
        inline std::vector<Chromosome>& elitism(float elitismRate);
        
        bool feasible(Chromosome& chromosome);
                
		Chromosome feasibilityCheck(Chromosome& chromosome);

        Chromosome getBestChromosome(std::vector<Chromosome> population);
        
        static Chromosome fitness(Chromosome& chromosome, Chromosome(*fitnessHeuristic)(Chromosome&));
        
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
		
		
		inline static Chromosome rouletteWheelSelection(std::vector<Chromosome> population); 
		static Chromosome chooseBestSolution(const Chromosome& chromosome1, const Chromosome& chromosome2);
        static Chromosome chooseWorstSolution(const Chromosome& chromosome1, const Chromosome& chromosome2);
        
        static int getRandomInt(int start, int end);
        static float getRandomFloat(float start, float end);
	
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
        void run(size_t generations, Chromosome(*heuristic)(Graph));
};	

#endif
