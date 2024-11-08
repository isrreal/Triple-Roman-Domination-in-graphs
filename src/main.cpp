#include "GeneticAlgorithm.hpp"  
#include "TripleRomanDomination.hpp"
#include "Graph.hpp"             
#include "AntColonyOptimization.hpp"


int main(int argc, char** argv) {
    if (argc > 7) {
        Graph graph("graph.txt", false);
        
        if (graph.getOrder() == 0)
        	return -1;
        	
        // graph, populationSize, genesSize, generations, heuristic, mutation rate, elitism rate, numberOfAnts, iterations, 
        TripleRomanDomination* drd = new TripleRomanDomination(graph, std::stoi(argv[1]), graph.getOrder(), std::stoi(argv[2]), std::stoi(argv[3]),
                std::stof(argv[4]), std::stof(argv[5]), std::stoi(argv[6]), std::stoi(argv[7])); 
        std::cout << "Triple Roman Domination Number computed by Genetic Algorithm: " << drd->getGamma3rGeneticAlgorithm() << std::endl;
        std::cout << "Triple Roman Domination Number computed by ACO: " << drd->getGamma3rACO() << std::endl;
        
        std::cout << "\nGenetic Algorithm solution: " << std::endl;
        for (const auto& it: drd->getSolutionGeneticAlgorithm())
           std::cout << it << " ";
        std::cout << std::endl;

        std::cout << "\nTriple Roman Domination Function: \n" << std::endl;
        std::cout << "ACO solution: " << std::endl;
        for (const auto& it: drd->getSolutionACO())
           std::cout << it << " ";
        std::cout << std::endl;

        delete drd;
    }

	return 0;
}
