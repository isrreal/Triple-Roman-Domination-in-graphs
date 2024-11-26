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
        static void toggleLabels(const Graph& graph, Chromosome& solution);
        static int getRandomInt(int start, int end);
        static float getRandomFloat(float start, float end);
        static void fitness(Chromosome& chromosome);
        
	public:
		TripleRomanDomination(Graph& graph, size_t populationSize, size_t genesSize, size_t generations, 
		float mutationRate, float elitismRate, float crossOverRate, 
		size_t numberOfAnts, size_t iterations):
		  graph(graph), 
		  geneticAlgorithm(new GeneticAlgorithm(graph, populationSize, genesSize, generations, mutationRate, elitismRate, crossOverRate)),
		  ACO(new AntColonyOptimization(graph, iterations, numberOfAnts)),
		  solutionACO(), 
		  solutionGeneticAlgorithm(),
		  geneticAlgorithmBestFitness(0), 
		  ACOBestFitness(0) {}

        ~TripleRomanDomination();
        Graph& getGraph();
        std::vector<int> getSolutionACO();
        std::vector<int> getSolutionGeneticAlgorithm();
        
        size_t getGeneticAlgorithmBestFitness();
        size_t getACOBestFitness();

        void runGeneticAlgorithm(short int heuristic);
        void runACO();
		
		static bool feasible(const Graph& graph, std::vector<int> solution);
		static std::vector<int> feasibilityCheck(const Graph& graph, std::vector<int>& solution);
		
        inline static Chromosome heuristic1(Graph graph);
        static Chromosome heuristic2(Graph graph);
        static Chromosome heuristic3(Graph graph);   
};
#endif
