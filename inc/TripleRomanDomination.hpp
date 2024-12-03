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
    	static GeneticAlgorithm* geneticAlgorithm;
    	static AntColonyOptimization* ACO;
		std::vector<int> solutionACO;
        std::vector<int> solutionGeneticAlgorithm;
        size_t geneticAlgorithmBestFitness;		
        size_t ACOBestFitness;   
              
	public:
		TripleRomanDomination(Graph& graph, size_t populationSize, size_t genesSize, size_t generations, 
		float mutationRate, float elitismRate, float crossOverRate, 
		size_t numberOfAnts, size_t iterations);

        ~TripleRomanDomination();
        
        Graph& getGraph();
        std::vector<int> getSolutionACO();
        std::vector<int> getSolutionGeneticAlgorithm();
        
        size_t getGeneticAlgorithmBestFitness();
        size_t getACOBestFitness();

        void runGeneticAlgorithm(short int heuristic);
        void runACO();
        
        static bool feasible(const Graph& graph, const std::vector<int>& solution);
		
        static Chromosome heuristic1(const Graph& graph);
        static Chromosome heuristic2(const Graph& graph);
        static Chromosome heuristic3(const Graph& graph);   
};
#endif
