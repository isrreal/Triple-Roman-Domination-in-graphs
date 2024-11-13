#include "GeneticAlgorithm.hpp"  
#include "TripleRomanDomination.hpp"
#include "Graph.hpp"             
#include "AntColonyOptimization.hpp"
#include <thread>
#include <chrono>

int main(int argc, char** argv) {
    if (argc > 1) {
        Graph graph(argv[1], false);
        
        if (graph.getOrder() == 0)
            return -1;
        
    	constexpr size_t trial = 10;
        constexpr size_t populationSize = 1000;
        size_t generations = graph.getOrder() / 2;
        constexpr short int heuristic = 1;
        constexpr bool hasRVNS  = false;
        constexpr float mutationRate = 0.05;
        constexpr float elitismRate = 0.15;
        constexpr size_t numberOfAnts = 20;
        constexpr size_t iterations = 10;
        
        size_t Delta = graph.getMaxDegree();
       	size_t delta = graph.getMinDegree();
        
        // graph, populationSize, genesSize, generations, heuristic, mutation rate, elitism rate, numberOfAnts, iterations
        TripleRomanDomination* trd = new TripleRomanDomination(graph, populationSize, graph.getOrder(), generations, heuristic,
                mutationRate, elitismRate, numberOfAnts, iterations); 
        std::cout << "graph_name,graph_order,graph_size,graph_min_degree,graph_max_degree,GA_fitness_heuristic" << heuristic << ",";
        std::cout << "ACO_fitness_" << numberOfAnts << "_" << iterations << ",lower_bound,upper_bound,elapsed_time_GA(seconds),elapsed_time_ACO(seconds)" << std::endl;
             
        for (size_t i = 0; i < trial; ++i) {         	
		    std::cout << argv[1] << ",";
		    std::cout << graph.getOrder() << ",";
		    std::cout << graph.getSize() << ",";
		    std::cout << delta << ",";
		    std::cout << Delta << ",";
		    
		    std::chrono::duration<double> elapsedGA;
			std::chrono::duration<double> elapsedACO;
		
			std::thread gaThread([&]() {
				auto startGA = std::chrono::high_resolution_clock::now();
				trd->runGeneticAlgorithm(heuristic, hasRVNS);
				auto endGA = std::chrono::high_resolution_clock::now();
				elapsedGA = endGA - startGA;
			});			
		
			std::thread acoThread([&]() {
				auto startACO = std::chrono::high_resolution_clock::now();
				trd->runACO();
				auto endACO = std::chrono::high_resolution_clock::now();
				elapsedACO = endACO - startACO;
			});
			
			gaThread.join();
			acoThread.join();
	
			std::cout << trd->getGeneticAlgorithmBestFitness() << ",";
			std::cout << trd->getACOBestFitness() << ",";
			std::cout << std::ceil(static_cast<double>(4 * graph.getOrder()) / (Delta + 1)) << ",";
			std::cout << static_cast<size_t>((3 * graph.getOrder()) / 2) << ",";

			std::cout << elapsedGA.count() << ",";
			std::cout << elapsedACO.count() << std::endl;	
    	}
    	
		delete trd;
		return EXIT_SUCCESS;
    }
}
