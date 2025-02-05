#include "GeneticAlgorithm.hpp"  
#include "TripleRomanDomination.hpp"
#include "Graph.hpp"             
#include "AntColonyOptimization.hpp"
#include "util_functions.hpp"
#include <thread>
#include <chrono>

void printGeneticAlgorithmLog(short heuristic) {
	std::cout << "graph_name,graph_order,graph_size,graph_min_degree,graph_max_degree,fitness_heuristic_" << heuristic;
    std::cout << ",lower_bound,upper_bound,graph_density,elapsed_time(seconds)\n";
}

void printAntColonyOptimizationLog(size_t number_of_ants, size_t iterations) {
	std::cout << "graph_name,graph_order,graph_size,graph_min_degree,graph_max_degree,fitness_" << number_of_ants << '_' << iterations;
    std::cout << ",lower_bound,upper_bound,graph_density,elapsed_time(seconds)\n";
}

void computeAntColonyOptimization(TripleRomanDomination& trd, int upper_bound, int lower_bound, double graph_density) {   
    std::chrono::duration<double> elapsed_time;
   	
	auto start = std::chrono::high_resolution_clock::now();
	
	trd.runACO();
	
	auto end = std::chrono::high_resolution_clock::now();
	elapsed_time = end - start;
  	    
	std::cout << trd.getACOBestFitness() << ',';
	std::cout << lower_bound << ',';
    std::cout << upper_bound << ',';  	
    std::cout << graph_density << ',';
	std::cout << elapsed_time.count() << '\n';
}

void computeGeneticAlgorithm(TripleRomanDomination& trd, short heuristic, int upper_bound, int lower_bound, double graph_density) {
    std::chrono::duration<double> elapsed_time;
  
	auto start = std::chrono::high_resolution_clock::now();
	
	trd.runGeneticAlgorithm(heuristic);
	
	auto end = std::chrono::high_resolution_clock::now();
	
	elapsed_time = end - start;
  	    
	std::cout << trd.getGeneticAlgorithmBestFitness() << ',';
	std::cout << lower_bound << ',';
    std::cout << upper_bound << ',';  	
    std::cout << graph_density << ',';
	std::cout << elapsed_time.count() << '\n';
}

auto main(int argc, char** argv) -> int {
    if (argc > 4) {
        Graph graph(argv[1]);
        
        if (graph.getOrder() == 0) {
            return -1;
        }
                     
    	constexpr size_t trial {10};
    	
    	// Genetic Algorithm parameters
    	
        size_t population_size { static_cast<size_t>(graph.getOrder() / 4) };
        constexpr size_t generations {623};
        short heuristic = std::stoi(argv[3]);
        constexpr float mutation_rate {0.2248};
        constexpr float elitism_rate {0.3095};
        constexpr float cross_over_rate {0.4995};
        size_t tournament_population_size { static_cast<size_t>(graph.getOrder() / 4) };
        constexpr size_t max_no_improvement_iterations {33};
        
        // ACO parameters
        
        constexpr size_t number_of_ants {6};
        constexpr size_t iterations {1};
        constexpr float evaporation_rate {0.4922};
        constexpr float min_destruction_rate {0.214};
        constexpr float max_destruction_rate {0.9};
        constexpr size_t max_rvns_functions {10};
        constexpr size_t max_rvns_iterations {150};
        constexpr size_t max_rvns_no_improvement_iterations {95};
        constexpr float selection_vertex_rate_extend_solution {0.3};
        constexpr float selection_vertex_rate_construct_solution {0.5};
        constexpr float add_vertices_rate_extend_solution {0.05};
        int upper_bound {0};
    	int lower_bound {0};
    	
		upper_bound = computeRightUpperBound(graph, upper_bound);

		lower_bound = computeRightLowerBound(graph, lower_bound);

        size_t Delta { graph.getMaxDegree() };
       	size_t delta { graph.getMinDegree() };
       	
   	   	double graph_density { static_cast<double>(2 * graph.getSize()) / (graph.getOrder() * (graph.getOrder() - 1)) };
       	
       // graph, populationSize, genesSize, generations, heuristic, mutation rate, elitism rate, numberOfAnts, iterations
        TripleRomanDomination trd(graph, population_size, graph.getOrder(), generations,
            	mutation_rate, elitism_rate, cross_over_rate, tournament_population_size,
            	max_no_improvement_iterations,
            	
            	number_of_ants, iterations, evaporation_rate,
			  	min_destruction_rate, max_destruction_rate, 
			  	max_rvns_functions, max_rvns_iterations, max_rvns_no_improvement_iterations,
			  	selection_vertex_rate_extend_solution, selection_vertex_rate_construct_solution,
			  	add_vertices_rate_extend_solution);
		
        if (std::stoi(argv[4]) == 1) {
        	printAntColonyOptimizationLog(number_of_ants, iterations);
        }
        else {
        	printGeneticAlgorithmLog(heuristic);
        }
        
       	for (size_t i {0}; i < trial; ++i) {         	
        	std::cout << argv[2] << ',';
	        std::cout << graph.getOrder() << ',';
	        std::cout << graph.getSize() << ',';
	        std::cout << delta << ',';
	        std::cout << Delta << ',';
	        
		 	if (std::stoi(argv[4]) == 1) {
                computeAntColonyOptimization(trd, upper_bound, lower_bound, graph_density);
            } 
            else {
                computeGeneticAlgorithm(trd, heuristic, upper_bound, lower_bound, graph_density);
            }
 		}
 		
    	return EXIT_SUCCESS;
	}
	
	return -1;
	
}
