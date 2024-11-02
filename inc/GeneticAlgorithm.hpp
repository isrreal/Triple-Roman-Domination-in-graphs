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
        
        float minDestructionRate;
        float maxDestructionRate;
        size_t maxRVNSnoImprovementIterations;
        size_t maxRVNSiterations;
        size_t currentRVNSnumber;
        size_t maxRVNSfunctions;
        
        float mutationRate;
        float elitismRate;

		void createPopulation(Chromosome(*heuristic)(Graph), Graph graph);
		
		Chromosome crossOver(Chromosome& chromosome1, Chromosome& cromossomo2,
                	Chromosome(*crossOverHeuristic)(Chromosome&, Chromosome&)); 
                	
        Chromosome mutation(Chromosome& chromosome);
        
        std::vector<Chromosome>& elitism(float elitismRate);
        
        bool feasible(Chromosome& chromosome);
                
		Chromosome feasibilityCheck(Chromosome& chromosome);
		
		std::vector<Chromosome>& createNewPopulation();
		
        Chromosome selectionMethod(Chromosome(*selectionHeuristic)(std::vector<Chromosome>)); 
        
        Chromosome getBestChromosome(std::vector<Chromosome> population);
        
        Chromosome RVNS(Chromosome& chromosome, Chromosome(*heuristic)(Graph));
        
        Chromosome destroySolution(Chromosome& chromosome);
        
        Chromosome extendSolution(Chromosome& chromosome);
        
        Chromosome reduceSolution(Chromosome& chromosome);
        
        size_t chooseVertex(Graph& graph);
        size_t chooseVertex(std::vector<int> twoOrZeroOrThreeLabeledVertices);
        
        size_t rouletteWheelSelection(Graph& graph);
        size_t rouletteWheelSelection(std::vector<int> twoOrZeroOrThreeLabeledVertices);
        
        static Chromosome fitness(Chromosome& chromosome, Chromosome(*fitnessHeuristic)(Chromosome&));
		static Chromosome tournamentSelection(std::vector<Chromosome> population);
		static Chromosome rouletteWheelSelection(std::vector<Chromosome> population); 
		static Chromosome chooseBestSolution(const Chromosome& chromosome1, const Chromosome& chromosome2);
        static Chromosome chooseWorstSolution(const Chromosome& chromosome1, const Chromosome& chromosome2);
	
	public:
		GeneticAlgorithm(Graph& graph, size_t populationSize, size_t genesSize, size_t generations,
			double mutationRate, double elitismRate):
        				populationSize(populationSize), genesSize(genesSize),
                        generations(generations), population(populationSize),
                        mutationRate(mutationRate), elitismRate(elitismRate),
                        graph(graph), maxRVNSiterations(150), 
                        maxRVNSnoImprovementIterations(10),
                        currentRVNSnumber(1), maxRVNSfunctions(5),
                        minDestructionRate(0.2), maxDestructionRate(0.5) {}                 

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
