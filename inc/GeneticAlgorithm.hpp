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

		void createPopulation(Chromosome(*heuristic)(Graph), Graph graph);
		
		Chromosome chooseBestSolution(const Chromosome& chromosome1, const Chromosome& chromosome2);
		
		Chromosome crossOver(Chromosome& chromosome1, Chromosome& cromossomo2,
                	Chromosome(*crossOverHeuristic)(Chromosome&, Chromosome&)); 
                
		Chromosome feasibilityCheck(Chromosome& chromosome);
		
		std::vector<Chromosome>& createNewPopulation();
		
        Chromosome selectionMethod(Chromosome(*selectionHeuristic)(std::vector<Chromosome>)); 
        	
        static Chromosome fitness(Chromosome& chromosome, Chromosome(*fitnessHeuristic)(Chromosome&));
		static Chromosome tournamentSelection(std::vector<Chromosome> population);
		static Chromosome rouletteWheelSelection(std::vector<Chromosome> population); 
		
	public:
		GeneticAlgorithm(Graph& graph, size_t populationSize, size_t genesSize, size_t generations):
        				populationSize(populationSize), genesSize(genesSize),
                        generations(generations), population(populationSize),
                        graph(graph) {}                 

		~GeneticAlgorithm();
        Graph getGraph();
		std::vector<Chromosome> getPopulation();     
		size_t getPopulationSize();    
		size_t getGenesSize();
		size_t getGenerations();   
        std::vector<int> getBestSolution();		      

		void run(size_t generations, Chromosome(*heuristic)(Graph));
};	

#endif
