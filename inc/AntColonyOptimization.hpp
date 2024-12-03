#ifndef ANT_COLONY_OPTIMIZATION_HPP
#define ANT_COLONY_OPTIMIZATION_HPP

#include "Graph.hpp"
#include <vector>

class AntColonyOptimization {
    private:
        Graph graph;
        std::vector<float> graphPheromone;
        size_t numberOfAnts;
        size_t iterations;
    
        std::vector<int> bestSolution;
         
        float convergenceFactor;
        float evaporationRate;

        float minDestructionRate;
        float maxDestructionRate;
        size_t currentRVNSnumber;
        size_t maxRVNSfunctions;
        size_t maxRVNSiterations;
        size_t maxRVNSnoImprovementIterations;
        float addVerticesRate;
        float selectionVertexRateExtendSolution;
			
        void initializePheromones(std::vector<float>& graphPheromone);
        void constructSolution(std::vector<int>& solution);
        
        void extendSolution(std::vector<int>& solution);
        void reduceSolution(std::vector<int>& solution);
        void RVNS(std::vector<int>& solution);
        
        size_t chooseVertex(const Graph& temp);
        size_t chooseVertex(const std::vector<int>& twoOrZeroLabeledVertices);
        
        size_t rouletteWheelSelection(const Graph& temp);
        size_t rouletteWheelSelection(const std::vector<int>& twoOrZeroLabeledVertices);
                                                                                      
        void destroySolution(std::vector<int>& solution);
        
        void updatePheromones(std::vector<int>& currentBestSolution, std::vector<int>& bestSolution);

        float computeConvergence(const std::vector<float>& graphPheromone);

        static bool feasible(const Graph& graph, const std::vector<int>& solution);

        size_t summation(const std::vector<int>& solution);
        
     	bool delta(const std::vector<int>& solution, size_t vertex);
        
        float getMinPheromoneValue(const std::vector<float>& graphPheromone);
        
        float getMaxPheromoneValue(const std::vector<float>& graphPheromone);   
        
    public:

        AntColonyOptimization(Graph& graph, size_t iterations, size_t numberOfAnts):
             graph(graph), graphPheromone(graph.getOrder(), 0.0),
             numberOfAnts(numberOfAnts), iterations(iterations),
             convergenceFactor(0), evaporationRate(0.2),
             minDestructionRate(0.2), maxDestructionRate(0.5),
             currentRVNSnumber(1), maxRVNSfunctions(5), maxRVNSiterations(150),
             maxRVNSnoImprovementIterations(50),
             addVerticesRate(0.05),
             selectionVertexRateExtendSolution(0.7) {}

        ~AntColonyOptimization() {} 
        
        static void toggleLabels(const Graph& graph, std::vector<int>& solution);
        std::vector<int> getBestSolution();
       
        void run();
};

#endif
