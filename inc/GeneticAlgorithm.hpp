#ifndef GENETIC_ALGORITHM_HPP
#define GENETIC_ALGORITHM_HPP

#include <iostream>
#include <random>
#include <vector>
#include "Chromosome.hpp"
#include "util_functions.hpp"

class GeneticAlgorithm {
	private:
		size_t population_size;
		size_t genes_size;
		std::vector<Chromosome> population;
	    size_t generations;
        Graph graph;
        std::vector<int> best_solution;
        
        float mutation_rate;
        float elitism_rate;
        float crossover_rate;
		float selection_chromosome_rate;
		size_t max_no_improvement_iterations;


		inline void createPopulation(std::vector<std::function<Chromosome(const Graph&)>> generateChromosomeHeuristics,
		 	const Graph& graph, size_t heuristic);
	
        inline std::vector<Chromosome>& createNewPopulation();    
        
     	inline std::vector<Chromosome>& elitism(float elitism_rate);	
        
		inline Chromosome& onePointCrossOver(const Chromosome& chromosome1, const Chromosome& cromossomo2); 
                	
    	inline Chromosome& twoPointCrossOver(const Chromosome& chromosome1, const Chromosome& cromossomo2);
                	
        inline Chromosome& mutation(Chromosome& chromosome);
        
		inline Chromosome& tournamentSelection(const std::vector<Chromosome>& population);
		
		inline static Chromosome& chooseBestSolution(Chromosome& chromosome1, Chromosome& chromosome2);
        inline static Chromosome& chooseWorstSolution(Chromosome& chromosome1, Chromosome& chromosome2);

	public:
		GeneticAlgorithm(Graph& graph, size_t population_size, size_t genes_size, size_t generations,
			float mutation_rate, float elitism_rate, float crossover_rate,
			float selection_chromosome_rate, size_t max_no_improvement_iterations):
			  population_size(population_size), genes_size(genes_size), 
			  population(population_size), generations(generations), 
			  graph(graph), best_solution(), 
			  mutation_rate(mutation_rate), elitism_rate(elitism_rate),
			  crossover_rate(crossover_rate),
			  selection_chromosome_rate(selection_chromosome_rate),
              max_no_improvement_iterations(max_no_improvement_iterations) {}               

		~GeneticAlgorithm() {}
		
        Graph getGraph();
		std::vector<Chromosome> getPopulation();     
		size_t getPopulationSize();    
		size_t getGenesSize();
		size_t getGenerations();   
		double getMutationRate();
		double getElitismRate();
        std::vector<int> getBestSolution();		      
        void run(size_t generations, std::vector<std::function<Chromosome(const Graph&)>>, size_t chosen_heuristic);
};	

#endif
