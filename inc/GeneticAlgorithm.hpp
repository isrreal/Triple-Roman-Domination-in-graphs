#ifndef GENETIC_ALGORITHM_HPP
#define GENETIC_ALGORITHM_HPP

#include <iostream>
#include <random>
#include <vector>
#include <functional>
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
        size_t tournament_population_size;
		size_t max_no_improvement_iterations;


		inline void createPopulation(std::vector<std::function<Chromosome(const Graph&)>> generateChromosomeHeuristics,
		 	const Graph& graph, size_t heuristic);
	
        inline std::vector<Chromosome>& createNewPopulation();    
        
     	inline void elitism(std::vector<Chromosome>& population, float elitism_rate);	
        
		inline Chromosome onePointCrossOver(const Chromosome& chromosome1, const Chromosome& cromossomo2); 
                	
    	inline Chromosome twoPointCrossOver(const Chromosome& chromosome1, const Chromosome& cromossomo2);
                	
        inline Chromosome& mutation(Chromosome& chromosome);
        
		inline const Chromosome& tournamentSelection(const std::vector<Chromosome>& population, size_t individuals_size);
		
		inline static Chromosome& chooseBestSolution(Chromosome& chromosome1, Chromosome& chromosome2);

	public:
		GeneticAlgorithm(Graph& graph, size_t population_size, size_t genes_size, size_t generations,
			float mutation_rate, float elitism_rate, float crossover_rate,
			size_t tournament_population_size, size_t max_no_improvement_iterations):
			  population_size(population_size), genes_size(genes_size), 
			  population(population_size), generations(generations), 
			  graph(graph), best_solution(), 
			  mutation_rate(mutation_rate), elitism_rate(elitism_rate),
			  crossover_rate(crossover_rate),
			  tournament_population_size(tournament_population_size),
              max_no_improvement_iterations(max_no_improvement_iterations) {}               

		~GeneticAlgorithm() {}
		
		size_t getGenerations();   
        std::vector<int> getBestSolution();		      
        void run(size_t generations, std::vector<std::function<Chromosome(const Graph&)>>, size_t chosen_heuristic);
};	

#endif
