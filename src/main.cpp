#include "GeneticAlgorithm.hpp"  
#include "TripleRomanDomination.hpp"
#include "Graph.hpp"             
#include "AntColonyOptimization.hpp"
#include "util_functions.hpp"
#include <thread>
#include <chrono>

int computeRightLowerBound(const Graph& graph, int lowerBound) {
    lowerBound = -1; 
    if (graph.getMaxDegree() >= 3 && graph.getOrder() >= 2) {
        lowerBound = std::ceil(static_cast<size_t>(4.0 * graph.getOrder() / (graph.getMaxDegree() + 1.0)));
	}
	
    return lowerBound;
}

int computeRightUpperBound(Graph& graph, int upperBound) {
    upperBound = -1;
    auto components { graph.connectedComponents() };
    
    if(components.size() == 1 && graph.getMinDegree() >= 2) {
        upperBound = std::floor(3.0 * graph.getOrder() / 2.0);
    }
    else if (components.size() == 1 && graph.getOrder() >= 3) {
        upperBound = std::floor(7.0 * graph.getOrder() / 4.0);
    }
    else if (components.size() > 1 && graph.getOrder() >= 3) {
        upperBound = 0;
        for (auto& par: components) {
            if(par.first == 1) {
                upperBound += 3;
            }
            else if(par.first == 2) {
                upperBound += 4;
            }            
            else if(par.first >= 3 && par.second >= 2) {
                upperBound += std::floor(static_cast<size_t>(3.0 * par.first / 2.0));
            }           
            else if(par.first >= 3) {
                upperBound += std::floor(static_cast<size_t>(7.0 * par.first / 4.0));
            }
        }
    }

    return upperBound;
}

auto main(int argc, char** argv) -> int {
    if (argc > 3) {
        Graph graph(argv[1], false);
        
        if (graph.getOrder() == 0) {
            return -1;
        }
                     
    	constexpr size_t trial {10};
    	
    	// Genetic Algorithm parameters
    	
        size_t population_size = static_cast<size_t>(graph.getOrder() / 3);
        constexpr size_t generations {144};
        short heuristic = std::stoi(argv[3]);
        constexpr float mutation_rate {0.5013};
        constexpr float elitism_rate {0.2534};
        constexpr float cross_over_rate {0.9719};
        constexpr float selection_chromosome_rate {0.6035};
        constexpr size_t max_no_improvement_iterations {75};
        
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
        int upperBound {0};
    	int lowerBound {0};
        
        size_t Delta { graph.getMaxDegree() };
       	size_t delta { graph.getMinDegree() };
       	
       // graph, populationSize, genesSize, generations, heuristic, mutation rate, elitism rate, numberOfAnts, iterations
        TripleRomanDomination trd(graph, population_size, graph.getOrder(), generations,
            	mutation_rate, elitism_rate, cross_over_rate, selection_chromosome_rate,
            	max_no_improvement_iterations,
            	
            	number_of_ants, iterations, evaporation_rate,
			  	min_destruction_rate,  max_destruction_rate, 
			  	max_rvns_functions, max_rvns_iterations, max_rvns_no_improvement_iterations,
			  	selection_vertex_rate_extend_solution, selection_vertex_rate_construct_solution,
			  	add_vertices_rate_extend_solution);
                 
        std::cout << "graph_name,graph_order,graph_size,graph_min_degree,graph_max_degree,GA_fitness_heuristic" << heuristic;
	    std::cout << ",ACO_fitness_" << number_of_ants << "_" << iterations; 
        std::cout << ",lower_bound,upper_bound,elapsed_time_GA(seconds)" << ",elapsed_time_ACO(seconds),is_3RDF_GA,is_3RDF_ACO\n";

       for (size_t i {0}; i < trial; ++i) {         	
	        std::cout << argv[2] << ",";
	        std::cout << graph.getOrder() << ",";
	        std::cout << graph.getSize() << ",";
	        std::cout << delta << ",";
	        std::cout << Delta << ",";
	        
	        std::chrono::duration<double> elapsedGA;
	    	std::chrono::duration<double> elapsedACO;
	   	
	    	std::thread gaThread([&]() {
	    		auto startGA = std::chrono::high_resolution_clock::now();
	    		trd.runGeneticAlgorithm(heuristic);
	    		auto endGA = std::chrono::high_resolution_clock::now();
	    		elapsedGA = endGA - startGA;
	    	});		
	    	
    		std::thread acoThread([&]() {
	    		auto startACO = std::chrono::high_resolution_clock::now();
	    		trd.runACO();
	    		auto endACO = std::chrono::high_resolution_clock::now();
	    		elapsedACO = endACO - startACO;
	    	});
	    	
	    	gaThread.join();
  	    	acoThread.join();
  	    	
        	upperBound = computeRightUpperBound(graph, upperBound);

        	lowerBound = computeRightLowerBound(graph, lowerBound);

	    	std::cout << trd.getGeneticAlgorithmBestFitness() << ",";
	     	std::cout << trd.getACOBestFitness() << ",";
	    	std::cout << lowerBound << ",";
            std::cout << upperBound << ","; 
	    	
	    	std::cout << elapsedGA.count() << ",";
	    	std::cout << elapsedACO.count() << ",";	    	
	    	
	    	std::cout << feasible(graph, trd.getSolutionGeneticAlgorithm()) << ",";
	    	std::cout << feasible(graph, trd.getSolutionACO()) << '\n';	
 		}
	}
	
    return EXIT_SUCCESS;
}
