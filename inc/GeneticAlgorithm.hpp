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
        
        double mutationRate;
        double elitismRate;

		void createPopulation(Chromosome(*heuristic)(Graph), Graph graph);
		
		Chromosome chooseBestSolution(const Chromosome& chromosome1, const Chromosome& chromosome2);
	   
        Chromosome chooseWorstSolution(const Chromosome& chromosome1, const Chromosome& chromosome2);

		Chromosome crossOver(Chromosome& chromosome1, Chromosome& cromossomo2,
                	Chromosome(*crossOverHeuristic)(Chromosome&, Chromosome&)); 
                	
        Chromosome mutation(Chromosome& chromosome);
        
        std::vector<Chromosome>& elitism();
                
		Chromosome feasibilityCheck(Chromosome& chromosome);
		
		std::vector<Chromosome>& createNewPopulation();
		
        Chromosome selectionMethod(Chromosome(*selectionHeuristic)(std::vector<Chromosome>)); 
        	
        static Chromosome fitness(Chromosome& chromosome, Chromosome(*fitnessHeuristic)(Chromosome&));
		static Chromosome tournamentSelection(std::vector<Chromosome> population);
		static Chromosome rouletteWheelSelection(std::vector<Chromosome> population); 
		
	public:
		GeneticAlgorithm(Graph& graph, size_t populationSize, size_t genesSize, size_t generations,
			double mutationRate, double elitismRate):
        				populationSize(populationSize), genesSize(genesSize),
                        generations(generations), population(populationSize),
                        mutationRate(mutationRate), elitismRate(elitismRate),
                        graph(graph) {}                 

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
