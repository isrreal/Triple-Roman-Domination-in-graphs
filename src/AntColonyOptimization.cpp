#include "AntColonyOptimization.hpp"

void AntColonyOptimization::run() {
    size_t temp = iterations;
    while (temp > 0) {
        for (size_t i = 0; i < numberOfAnts; ++i) {
            solution = constructSolution(solution);
            solution = extendSolution(solution);
            solution = reduceSolution(solution);
            solution = RVNS(solution);
            if (summation(solution) < summation(currentBestSolution))
                currentBestSolution = solution;
        }    
                                                           
        if (summation(currentBestSolution) < summation(bestSolution))
              bestSolution = currentBestSolution;
        
        updatePheromones(currentBestSolution, bestSolution, graphPheromone);
                                                               
        convergenceFactor = computeConvergence(graphPheromone);
                                 
        if (convergenceFactor > 0.99)
            initializePheromones(graphPheromone);

        --temp;
    }
}

void AntColonyOptimization::initializePheromones(std::vector<float>& graphPheromone) {
    for (size_t i = 0; i < graphPheromone.size(); ++i)
        graphPheromone[i] = 0.5;
}



/**
 * @brief Constructs a feasible solution for the Triple Roman Domination problem.
 * 
 * @details This function selects a random vertex and labels it with a value of 4, setting the labels 
 * of all its adjacent vertices to 0. These vertices are then removed from the temporary graph.
 * If any remaining vertices in the temporary graph have a degree of 0, they are labeled with 3 and 
 * also removed. This process continues until there are no vertices left in the temporary graph.
 *
 * @return A vector of integers representing the computed solution for the graph, where each element 
 * indicates the label assigned to each vertex according to the Triple Roman Domination rules.
 */

std::vector<int> AntColonyOptimization::constructSolution(std::vector<int> solution) {
    Graph temp(this->graph);
    size_t vertex = 0;
    
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_int_distribution<int> gap(0, graph.getOrder() - 1);
	
	size_t graphOrder = graph.getOrder();
    while (temp.getOrder() > 0) {
        vertex = chooseVertex(temp);  
        solution[vertex] = 4;

        for (const auto& it: temp.getAdjacencyList(vertex)) 
            if (solution[it] == -1) 
                solution[it] = 0;

        temp.deleteAdjacencyList(vertex);
        
        for (size_t i = 0; i < graphOrder; ++i) {
            if (temp.vertexExists(i) && temp.getVertexDegree(i) == 0) {
                    solution[i] = 2;
                    temp.deleteVertex(i);
            }
        }
    }
  
  	  
    for (size_t i = 0; i < graphOrder; ++i) {
        if (solution[i] == 2) {
            bool onlyLabel0 = true;
            for (const auto& neighbor: graph.getAdjacencyList(i)) {	
                if (solution[neighbor] >= 2) {
                    onlyLabel0 = false;
                    break;
                }
            }

            if (onlyLabel0) {
                auto neighbors = graph.getAdjacencyList(i);
                size_t randomNeighbor = gap(seed) % neighbors.size();
                auto it = std::next(neighbors.begin(), randomNeighbor);
                solution[*it] = 2;
            }
        }
    }
    
    return solution;
}

/**
 * @brief Increases the labels of selected vertices in the input solution.
 * 
 * @details This function iteratively selects random vertices from the input graph that are labeled 
 * with 0, 2, or 3, and increases their label to 4. The number of vertices to be updated is determined 
 * by the constant `addVerticesRate` and the count of vertices with labels 0, 2, or 3. This process 
 * continues until the target number of vertices is updated or there are no more eligible vertices.
 * 
 * @return A vector of integers representing the updated solution, where selected labels have been 
 * increased to 4 based on the specified rate.
 */

std::vector<int> AntColonyOptimization::extendSolution(std::vector<int> solution) {
    constexpr float addVerticesRate = 0.05f;
    size_t itr = 0;
    std::vector<int> twoOrZeroOrThreeLabeledVertices;
    
    for (size_t i = 0; i < solution.size(); ++i) {
        if ((solution[i] == 0) || (solution[i] == 2) || (solution[i] == 3))
            twoOrZeroOrThreeLabeledVertices.push_back(i);
    }

    itr = addVerticesRate * twoOrZeroOrThreeLabeledVertices.size();
    
    size_t vertex = 0;

    while (itr != 0 && !twoOrZeroOrThreeLabeledVertices.empty()) {          
        vertex = chooseVertex(twoOrZeroOrThreeLabeledVertices);         
        solution[vertex] = 4; 
        twoOrZeroOrThreeLabeledVertices.erase(twoOrZeroOrThreeLabeledVertices.begin() + vertex);                                                                
        --itr;
    }

    return solution;
}

/**	@brief Tries to reduce the labels of vertices 
*	@details Sorts the vertices of graph on descendent order based on its respective degree,
*			and selects its first vertex (that has higher degree value). After this, 
			tries to reduce the values of vertices in solution that have labels 4 or 3, to 2,
			if the solution are feasible, then the label is updated to 0, otherwise,
			the original label is recovered.
*
* @return A solution that can have its label values reduced
*/

std::vector<int> AntColonyOptimization::reduceSolution(std::vector<int> solution) {
    Graph temp = this->graph;
    std::vector<int> sortedVertices;
    int initLabel = -1;

    for (size_t i = 0; i < temp.getOrder(); ++i)
        sortedVertices.push_back(i);

    std::sort(sortedVertices.begin(), sortedVertices.end(),
        [&](size_t a, size_t b) {
            return temp.getVertexDegree(a) < temp.getVertexDegree(b);                                                                             
        });
    
    size_t choosenVertex = 0;
    
    while ((temp.getOrder() > 0) && (choosenVertex < sortedVertices.size())) {
        if (choosenVertex >= sortedVertices.size()) break;
        
        while (choosenVertex < sortedVertices.size() && 
                (!temp.vertexExists(sortedVertices[choosenVertex]))) {
            ++choosenVertex;
        }
       
        if (choosenVertex >= sortedVertices.size()) break;

        if (solution[sortedVertices[choosenVertex]] == 4 || solution[sortedVertices[choosenVertex]] == 3) {         
    		initLabel = solution[sortedVertices[choosenVertex]];
    		solution[sortedVertices[choosenVertex]] = 0;
    		
    		if (!feasible(solution)) {
        		solution[sortedVertices[choosenVertex]] = 2;
        		
        		if (!feasible(solution)) {
        			solution[sortedVertices[choosenVertex]] = 3;
            		
            		if (!feasible(solution)) 
                		solution[sortedVertices[choosenVertex]] = initLabel;        
        		}
    		}
		}
         	
        temp.deleteAdjacencyList(sortedVertices[choosenVertex++]);
    }
    return solution;
}

/**
 * @brief Attempts to reduce the labels of vertices in the solution.
 * 
 * @details This function sorts the vertices of the graph in descending order based on their degree 
 * and selects the vertex with the highest degree. For each selected vertex with a label of 4 or 3 
 * in the solution, it attempts to reduce the label to 2 and checks if the solution remains feasible.
 * If the solution remains feasible, the label is updated to 0; otherwise, the original label is restored.
 * This process continues until all vertices have been processed.
 * 
 * @return A vector of integers representing the solution with potentially reduced labels, while 
 * maintaining feasibility according to the problem's constraints.
 */

std::vector<int> AntColonyOptimization::destroySolution(std::vector<int> solution) {
    float destructionRate = minDestructionRate + ((currentRVNSnumber - 1) *
                ((maxDestructionRate - minDestructionRate)) 
                / (maxRVNSfunctions - 1));
    Graph temp = this->graph;
    size_t itr = solution.size() * destructionRate; 
    size_t vertex = 0;
    while (itr != 0 && (temp.getOrder() > 0)) {
       vertex = chooseVertex(temp);
       if ((solution[vertex] == 0) || (solution[vertex] == 2) || (solution[vertex] == 3))
            solution[vertex] = -1;
       else
            ++itr;
       temp.deleteVertex(vertex);
       --itr;
    }

    return solution;
} 

/**
 * @brief Attempts to improve the final solution computed by the algorithm.
 * 
 * @details This function takes a previously computed solution and applies the following subroutines:
 * - @subroutine destroySolution
 * - @subroutine constructSolution
 * - @subroutine extendSolution
 * - @subroutine reduceSolution
 * 
 * After each iteration, the function compares the weight of the newly generated solution with the original one.
 * If the new solution has a lower weight, it updates the best solution to this new value.
 * The process repeats until the maximum number of iterations without improvement is reached.
 * 
 * @param solution A vector of integers representing the previously computed by the algorithm.
 * 
 * @return A vector of integers representing the solution with the least weight found by the algorithm.
 */

std::vector<int> AntColonyOptimization::RVNS(std::vector<int> solution) {
    size_t currentNoImprovementIteration = 0;
    std::vector<int> temp = solution;
    currentRVNSnumber = 1;

    while ((currentNoImprovementIteration < maxRVNSnoImprovementIterations) && (maxRVNSiterations > 0)) {
        destroySolution(solution);
        constructSolution(solution);
        extendSolution(solution);
        reduceSolution(solution);

        if (summation(temp) < summation(solution)) {
            solution = temp;
            currentRVNSnumber = 1;
            currentNoImprovementIteration = 0;
        }

        else {
            ++currentRVNSnumber;
            ++currentNoImprovementIteration;

            if (currentRVNSnumber > maxRVNSiterations)
                currentRVNSnumber = 1;
        }

        --maxRVNSiterations;
    }


    return solution;
}

size_t AntColonyOptimization::chooseVertex(Graph& temp) {
    float selectionVertexRateConstructSolution = 0.7f;    
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_real_distribution<float> probabilityGap(0.0, 1.0);
    std::uniform_int_distribution<int> gap(0, this->graph.getOrder() - 1);
    
    float number = probabilityGap(seed); 
    
    size_t vertex = gap(seed);
    float value = 0.0;
    
    while (!temp.vertexExists(vertex)) 
        vertex = gap(seed);

    if (selectionVertexRateConstructSolution < number) {
        for (size_t i = 0; i < this->graph.getOrder(); ++i) {
            if (temp.vertexExists(i)) { 
                if (value < (temp.getVertexDegree(i) * graphPheromone[i])) {
                    value = temp.getVertexDegree(i) * graphPheromone[i];
                	vertex = i;
																																																																											                                  																																																																											          																																																																					                                                                                               
                }
            }
        }
    }

    else
        vertex = rouletteWheelSelection(temp);

    return vertex;
}

/**
 * @brief Selects the index of a vertex that maximizes the objective function or selects one at random.
 * 
 * @details This function either chooses a vertex index that maximizes the product of its degree 
 * and pheromone level (objective function) or randomly selects an index from the 
 * `twoOrZeroOrThreeLabeledVertices` vector. The choice between these strategies is determined by the 
 * `selectionVertexRateExtendSolution`, which acts as a threshold probability. If a randomly generated 
 * number exceeds this rate, the function will maximize the objective function; otherwise, a random 
 * index is selected.
 * 
 * @param twoOrZeroOrThreeLabeledVertices A vector containing indices of vertices labeled as 0, 2, or 3.
 * @return The index within `twoOrZeroOrThreeLabeledVertices` that either maximizes the objective function
 *         or is randomly selected.
 */

size_t AntColonyOptimization::chooseVertex(std::vector<int> twoOrZeroOrThreeLabeledVertices) {
    constexpr float selectionVertexRateExtendSolution = 0.9f;

    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_real_distribution<float> probabilityGap(0.0, 1.0);
    std::uniform_int_distribution<int> randomIndex(0, twoOrZeroOrThreeLabeledVertices.size() - 1);                                                               

    float randomNumberGenerated = probabilityGap(seed); 

    size_t choosenIndex = randomIndex(seed);

    if (randomNumberGenerated > selectionVertexRateExtendSolution) {
        float maxObjectiveValue = -1.0f;
        float objectiveValue = 0.0;
        size_t vertex = 0;

        for (size_t i = 0; i < twoOrZeroOrThreeLabeledVertices.size(); ++i) {
            vertex = twoOrZeroOrThreeLabeledVertices[i];

            objectiveValue = this->graph.getVertexDegree(vertex) * graphPheromone[vertex];

            if (objectiveValue > maxObjectiveValue) {
                maxObjectiveValue = objectiveValue;
                choosenIndex = i;
            }
        }
    }

    else 
        choosenIndex = rouletteWheelSelection(twoOrZeroOrThreeLabeledVertices);

    return choosenIndex;
}

/**
 * @brief Evaluates if the given solution is feasible according to the constraints of the Triple Roman Domination Function.
 * 
 * @details This function checks if each vertex in the solution meets the required conditions for 
 * Triple Roman Domination. For each vertex labeled as 0, it verifies if there is at least one adjacent 
 * vertex labeled as 4, two adjacent vertices labeled as 3 or 2, or three adjacent vertices labeled as 2.
 * For vertices labeled as 2, it ensures that they have at least one neighboring vertex labeled as 4.
 * 
 * @param solution A vector representing the labeled solution of the graph.
 * @return True if the input solution meets the constraints of the Triple Roman Domination Function; False otherwise.
 */

bool AntColonyOptimization::feasible(std::vector<int> solution) {
    bool isValid = false;

    for (size_t i = 0; i < solution.size(); ++i) {

        isValid = false;

        if (solution[i] == 0) {                                 
            size_t countNeighbors2 = 0;
            size_t countNeighbors3 = 0;
            
            for (auto& neighbor : this->graph.getAdjacencyList(i)) {          

                if ((countNeighbors2 == 1 && solution[neighbor] >= 3) ||
                    (countNeighbors2 == 2 && solution[neighbor] >= 2)) {
                    isValid = true;
                    break;
                }

                if (countNeighbors3 == 1 && solution[neighbor] >= 2) {
                     isValid = true;
                     break;
                }

                if (solution[neighbor] == 4) {
                    isValid = true;
                    break;
                }
                
                if (solution[neighbor] == 3) 
                    ++countNeighbors3;
                    
                if (solution[neighbor] == 2) 
                    ++countNeighbors2;
            }
            
            if (!isValid)
                return false;
        } 

        else if (solution[i] == 2) {
            bool hasNeighborAtLeast2 = false;
            for (auto& neighbor : this->graph.getAdjacencyList(i)) {
                if (solution[neighbor] >= 2) {
                    hasNeighborAtLeast2 = true;
                    break; 
                }
            }

            if (!hasNeighborAtLeast2) 
                return false;           
        }
    }
    
    return true;    
}

size_t AntColonyOptimization::summation(std::vector<int> solution) {
    size_t temp = 0;
    for (size_t i = 0; i < solution.size(); ++i)
        temp += solution[i];
    return temp;
}


/* @brief checks if the vertex was choosen
 * The critery of selection vertex is that maximize the function degree(v) * pheromone(v)
 * if 4, so the vertex was choosen
 * 
 *@return the vertex was choosen or not
 *
 */ 
 
bool AntColonyOptimization::delta(std::vector<int> solution, size_t vertex) {
    return solution[vertex] == 4 ? true : false; 
}

float AntColonyOptimization::getMaxPheromoneValue(std::vector<float> graphPheromone) {
    float value = graphPheromone[0];

    for (const auto& it: graphPheromone)
        if (value < it)
            value = it;
    return value;
}

float AntColonyOptimization::getMinPheromoneValue(std::vector<float> graphPheromone) {
    float value = graphPheromone[0];                                     
    for (const auto& it: graphPheromone) 
        if (value > it)
            value = it;
    return value;                   
}


void AntColonyOptimization::updatePheromones(std::vector<int>& currentBestSolution,                                                               
                                             std::vector<int>& bestSolution,        
                                             std::vector<float>& pheromoneValues) {                                                               
    size_t weightCurrentBestSolution = summation(currentBestSolution);
    size_t weightBestSolution = summation(bestSolution); 
    float equation = 0.0;        

    for (size_t i = 0; i < graphPheromone.size(); ++i) {  
        equation = (((weightCurrentBestSolution * delta(currentBestSolution, i) + 
                      weightBestSolution * delta(bestSolution, i))) 
                    / (weightCurrentBestSolution + weightBestSolution)) - graphPheromone[i];                                                      
         
        graphPheromone[i] += evaporationRate * equation;
    }
}


float AntColonyOptimization::computeConvergence(std::vector<float> graphPheromone) {
    float maxPheromone = getMaxPheromoneValue(graphPheromone); 
    float minPheromone = getMinPheromoneValue(graphPheromone); 
    float temp = 0.0;

    size_t numberOfAnts = graphPheromone.size();

    for (size_t i = 0; i < numberOfAnts; ++i) {
        temp += std::max(maxPheromone - graphPheromone[i], graphPheromone[i] - minPheromone);
    }

    this->convergenceFactor = 2 * ((temp / (numberOfAnts * (maxPheromone + minPheromone)))) - 1;

    return this->convergenceFactor;
}

size_t AntColonyOptimization::rouletteWheelSelection(Graph& temp) {
    float totalFitness = 0.0f;
    std::vector<std::pair<size_t, float>> probabilities;

    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_real_distribution<float> gap(0.0, 1.0);

    for (size_t i = 0; i < graph.getOrder(); ++i) 
        if (temp.vertexExists(i)) 
            totalFitness += temp.getVertexDegree(i) * graphPheromone[i];

    for (size_t i = 0; i < graph.getOrder(); ++i)
        if (temp.vertexExists(i))
            probabilities.push_back({i, ((temp.getVertexDegree(i) * graphPheromone[i]) / totalFitness)});
    
    float randomValue = gap(seed);

    float cumulativeSum = 0.0f;
    for (const auto& prob : probabilities) {
        cumulativeSum += prob.second;
        
        if (randomValue <= cumulativeSum) 
            return prob.first;
    }

    return probabilities.back().first;
}



size_t AntColonyOptimization::rouletteWheelSelection(std::vector<int> twoOrZeroOrThreeLabeledVertices) {
    float totalFitness = 0.0f;
    std::vector<std::pair<size_t, float>> probabilities;
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_real_distribution<float> gap(0.0, 1.0);

    for (size_t i = 0; i < twoOrZeroOrThreeLabeledVertices.size(); ++i) 
        totalFitness += graph.getVertexDegree(twoOrZeroOrThreeLabeledVertices[i]) * graphPheromone[twoOrZeroOrThreeLabeledVertices[i]];
    
    for (size_t i = 0; i < twoOrZeroOrThreeLabeledVertices.size(); ++i)
        probabilities.push_back({i, ((graph.getVertexDegree(twoOrZeroOrThreeLabeledVertices[i]) * graphPheromone[twoOrZeroOrThreeLabeledVertices[i]]) / totalFitness)});

    float randomValue = gap(seed);
                                             
    float cumulativeSum = 0.0f;
    for (const auto& prob : probabilities) {
        cumulativeSum += prob.second;
        if (randomValue <= cumulativeSum) 
            return prob.first;
    }
                                             
    return probabilities.back().first;
}


std::vector<int> AntColonyOptimization::getBestSolution() { return this->bestSolution; }
