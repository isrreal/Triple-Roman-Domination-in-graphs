#ifndef TRIPLE_ROMAN_DOMINATION_HPP
#define TRIPLE_ROMAN_DOMINATION_HPP

#include "Graph.hpp"
#include "GeneticAlgorithm.hpp"
#include "Chromosome.hpp"
#include "AntColonyOptimization.hpp"
#include <vector>
#include <random>

class TripleRomanDomination {
	private:
    	Graph graph;
    	GeneticAlgorithm* geneticAlgorithm;
    	AntColonyOptimization* ACO;
		std::vector<int> solutionACO;
        std::vector<int> solutionGeneticAlgorithm;
        size_t geneticAlgorithmBestFitness;		
        size_t ACOBestFitness;        
	public:
		TripleRomanDomination(Graph& graph, size_t populationSize, size_t genesSize, size_t generations, 
			short int heuristic, double mutationRate, double elitismRate,
			size_t numberOfAnts, size_t iterations) 
    			: graph(graph), geneticAlgorithmBestFitness(0), ACOBestFitness(0),   
                geneticAlgorithm(new GeneticAlgorithm(graph, populationSize, genesSize, generations, mutationRate, elitismRate)),
    		    ACO(new AntColonyOptimization(graph, iterations, numberOfAnts)) {}

        ~TripleRomanDomination();
        Graph& getGraph();
        std::vector<int> getSolutionACO();
        std::vector<int> getSolutionGeneticAlgorithm();
        bool feasible(std::vector<int> solution);
        
        size_t getGeneticAlgorithmBestFitness();
        size_t getACOBestFitness();

        void runGeneticAlgorithm(short int heuristic);
        void runACO();

        static Chromosome heuristic1(Graph graph);
        static Chromosome heuristic2(Graph graph);
        static Chromosome heuristic3(Graph graph);   
};
#endif
