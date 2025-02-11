#ifndef ANT_COLONY_OPTIMIZATION_HPP
#define ANT_COLONY_OPTIMIZATION_HPP

#include "Graph.hpp"
#include "util_functions.hpp"
#include <vector>

class AntColonyOptimization {
    private:
        Graph graph;
        std::vector<float> graph_pheromones;
        size_t number_of_ants;
        size_t iterations;
    
        std::vector<int> best_solution;
         
        float convergence_factor;
        float evaporation_rate;

        float min_destruction_rate;
        float max_destruction_rate;
        size_t current_rvns_number;
        size_t max_rvns_functions;
        size_t max_rvns_iterations;
        size_t max_rvns_no_improvement_iterations;
     	float selection_vertex_rate_construct_solution;
        float selection_vertex_rate_extend_solution;
        float add_vertices_rate_extend_solution;
        			
        inline void initializePheromones(std::vector<float>& graph_pheromones);
        void constructSolution(std::vector<int>& solution);
        
        inline void extendSolution(std::vector<int>& solution);
        inline void reduceSolution(std::vector<int>& solution);
        inline void RVNS(std::vector<int>& solution);
        
        size_t chooseVertex(const Graph& temp);
        size_t chooseVertex(const std::vector<int>& two_or_zero_labeled_vertices);
        
        size_t rouletteWheelSelection(const Graph& temp);
        size_t rouletteWheelSelection(const std::vector<int>& two_or_zero_labeled_vertices);
                                                                                      
        inline void destroySolution(std::vector<int>& solution);
        
        inline void updatePheromones(std::vector<int>& current_best_solution, std::vector<int>& best_solution);

        inline float computeConvergence(const std::vector<float>& graph_pheromones);

        inline size_t summation(const std::vector<int>& solution);
        
     	inline bool delta(const std::vector<int>& solution, size_t vertex);
        
        inline float getMinPheromoneValue(const std::vector<float>& graph_pheromones);
        
        inline float getMaxPheromoneValue(const std::vector<float>& graph_pheromones);   
        
    public:

        AntColonyOptimization(Graph& graph, size_t number_of_ants, size_t iterations, float evaporation_rate,
			  	float min_destruction_rate,  float max_destruction_rate, 
			  	size_t max_rvns_functions, size_t max_rvns_iterations, size_t max_rvns_no_improvement_iterations,
			  	float selection_vertex_rate_extend_solution, float selection_vertex_rate_construct_solution,
			  	float add_vertices_rate_extend_solution):
			  	
             graph(graph), graph_pheromones(graph.getOrder(), 0.0),
             number_of_ants(number_of_ants), iterations(iterations),
             convergence_factor(0), evaporation_rate(evaporation_rate),
             min_destruction_rate(min_destruction_rate), max_destruction_rate(max_destruction_rate),
             current_rvns_number(1), max_rvns_functions(max_rvns_functions), max_rvns_iterations(max_rvns_iterations),
             max_rvns_no_improvement_iterations(max_rvns_no_improvement_iterations),
             selection_vertex_rate_construct_solution(selection_vertex_rate_construct_solution),  
             selection_vertex_rate_extend_solution(selection_vertex_rate_extend_solution),
             add_vertices_rate_extend_solution(add_vertices_rate_extend_solution) {}

        ~AntColonyOptimization() {} 
        
        std::vector<int> getBestSolution();
       
        void run(bool with_RVNS);
};

#endif
