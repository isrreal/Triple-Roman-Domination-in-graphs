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
        static void setNeighbor2(const Graph& graph, Chromosome& solution);     
        static void toggleLabels(const Graph& graph, Chromosome& solution);
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
        
        
        size_t getGeneticAlgorithmBestFitness();
        size_t getACOBestFitness();

        void runGeneticAlgorithm(short int heuristic, bool hasRVNS);
        void runACO();
		
		static bool feasible(const Graph& graph, std::vector<int> solution);
		static std::vector<int> feasibilityCheck(const Graph& graph, std::vector<int> solution);
		
        static Chromosome heuristic1(Graph graph);
        static Chromosome heuristic2(Graph graph);
        static Chromosome heuristic3(Graph graph);   
        
        static Chromosome heuristic1RVNS(Graph graph, Chromosome& chromosome);
        static Chromosome heuristic2RVNS(Graph graph, Chromosome& chromosome);
        static Chromosome heuristic3RVNS(Graph graph, Chromosome& chromosome);   
};
#endif
