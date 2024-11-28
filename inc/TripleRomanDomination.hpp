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
        
        inline static void toggleLabels(const Graph& graph, Chromosome& solution) {
        	size_t initLabel = 0;
	
			for (size_t i = 0; i < graph.getOrder(); ++i) {
				if (solution.genes[i] == 4 || solution.genes[i] == 3) { 
					initLabel = solution.genes[i];
					solution.genes[i] = 0;
						
					if (!feasible(graph, solution.genes)) {
						solution.genes[i] = 2;
						
						if (!feasible(graph, solution.genes)) {
							solution.genes[i] = 3;
							
							if (!feasible(graph, solution.genes)) 
					    		solution.genes[i] = initLabel;        
						}
					}
				}
			}
        }
        
        inline static int getRandomInt(int start, int end) {
        	std::random_device randomNumber; 
    		std::mt19937 seed(randomNumber()); 
    		std::uniform_int_distribution<> gap(start, end); 
    
    		return gap(seed);
        }
        
        inline static void fitness(Chromosome& chromosome) {
        	for (size_t i = 0; i < chromosome.genesSize; ++i)
        		chromosome.fitnessValue += chromosome.genes[i];
        }
        
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
		
        static Chromosome heuristic1(const Graph& graph);
        static Chromosome heuristic2(const Graph& graph);
        static Chromosome heuristic3(const Graph& graph);   
};
#endif
