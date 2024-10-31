#ifndef ANT_COLONY_OPTIMIZATION_HPP
#define ANT_COLONY_OPTIMIZATION_HPP

#include "Graph.hpp"
#include <vector>

class AntColonyOptimization {
    private:
        Graph graph;
        std::vector<float> graphPheromone;
        std::vector<int> solution;
        size_t numberOfAnts;
        size_t iterations;

        std::vector<int> currentBestSolution;
        std::vector<int> bestSolution;
         
        float convergenceFactor;
        float evaporationRate;

        float minDestructionRate;
        float maxDestructionRate;
        size_t currentRVNSnumber;
        size_t maxRVNSfunctions;
        size_t maxRVNSiterations;
        size_t maxRVNSnoImprovementIterations;


        void initializePheromones(std::vector<float>& graphPheromone);
        std::vector<int> constructSolution(std::vector<int> solution);
        std::vector<int> extendSolution(std::vector<int> solution);
        std::vector<int> reduceSolution(std::vector<int> solution);
        std::vector<int> RVNS(std::vector<int> solution);
                                                                       
        std::vector<int> destroySolution(std::vector<int> solution);

        size_t chooseVertex(Graph& temp);
        size_t chooseVertex(std::vector<int> twoOrZeroLabeledVertices);

        bool feasible(std::vector<int> solution);

        size_t summation(std::vector<int> solution);
        
        float getMaxPheromoneValue(std::vector<float> graphPheromone);
        
        float getMinPheromoneValue(std::vector<float> graphPheromone);
        
        bool delta(std::vector<int> solution, size_t vertex);

        void updatePheromones(std::vector<int>& currentBestSolution,
                std::vector<int>& bestSolution,
                std::vector<float>& graphPheromone);

        float computeConvergence(std::vector<float> graphPheromone);
        
        size_t rouletteWheelSelection(Graph& temp);
        size_t rouletteWheelSelection(std::vector<int> twoOrZeroLabeledVertices);

    public:

        AntColonyOptimization(Graph& graph, size_t iterations, size_t numberOfAnts):
             graph(graph), solution(graph.getOrder(), -1),
             graphPheromone(graph.getOrder(), 0.0),
             numberOfAnts(numberOfAnts), iterations(iterations),
             convergenceFactor(0), evaporationRate(0.2),
             minDestructionRate(0.2), maxDestructionRate(0.5),
             currentRVNSnumber(1), maxRVNSfunctions(5), maxRVNSiterations(150),
             maxRVNSnoImprovementIterations(10), 
             currentBestSolution(graph.getOrder(), 3), bestSolution(graph.getOrder(), 3) {}

        ~AntColonyOptimization() {} 
        std::vector<int> getBestSolution();


        void run();
};

#endif
