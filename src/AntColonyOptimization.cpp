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

std::vector<int> AntColonyOptimization::constructSolution(std::vector<int> solution) {
    Graph temp(this->graph);
    size_t vertex = 0;

    while (temp.getOrder() > 0) {
        vertex = chooseVertex(temp);  
        solution[vertex] = 4;

        for (const auto& it: temp.getAdjacencyList(vertex)) 
            if (solution[it] == -1) 
                solution[it] = 0;

        temp.deleteAdjacencyList(vertex);
    }

    return solution;
}

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

        if (solution[choosenVertex] == 3 || solution[choosenVertex] == 2) {
            initLabel = solution[sortedVertices[choosenVertex]];
            solution[sortedVertices[choosenVertex]] = 0;
        
            if (!feasible(solution)) {
                solution[sortedVertices[choosenVertex]] = 2;
            
                if (!feasible(solution)) 
                    solution[sortedVertices[choosenVertex]] = initLabel;
            }
        }

        temp.deleteAdjacencyList(sortedVertices[choosenVertex++]);
    }
    
    return solution;
}



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
                	vertex = 																																																																											 i;                                                                                                                   
                }
            }
        }
    }

    else
        vertex = rouletteWheelSelection(temp);

    return vertex;
}

/**
 * @brief Chooses the index of a vertex that maximizes the objective function or selects one randomly.
 * 
 * This function randomly selects an index from the `twoOrZeroLabeledVertices` vector.
 * Depending on the selection rate (`selectionVertexRateExtendSolution`), it may opt to choose
 * the vertex that maximizes the product of its degree and pheromone level.
 * 
 * @param twoOrZeroLabeledVertices Vector containing indices of vertices labeled as 0 or 2.
 * @return The index of the `twoOrZeroLabeledVertices` vector that maximizes the objective function, 
 *         or one selected randomly.
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
 
bool AntColonyOptimization::feasible(std::vector<int> solution) {	
    bool hasNeighborWith4 = false; 
    bool hasTwoNeighbors3or2[] = {false, false}; 
    bool hasThreeNeighbors2[] = {false, false, false};
    bool hasNeighborAtLeast2 = false;

    for (size_t i = 0; i < solution.size(); ++i) {
        if (solution[i] == 0) {
            for (auto& neighbor : this->graph.getAdjacencyList(i)) {
                if (solution[neighbor] == 4) {
                    hasNeighborWith4 = true;
                    break; 
                }

                else if (solution[neighbor] == 2) {
                    if (solution[neighbor] == 3) {
                        hasTwoNeighbors3or2[0], hasTwoNeighbors3or2[1] = true;
                        break;
                    }

                    else if (solution[neighbor] == 2) {
                        if (solution[neighbor] == 2) {
                            hasThreeNeighbors2[0], hasThreeNeighbors2[0], hasThreeNeighbors2[0] = true;
                            break;
                        } 
                    } 
                }
                
                if (!hasNeighborWith4 && !hasTwoNeighbors3or2[0] && !hasTwoNeighbors3or2[1]
                        && !hasThreeNeighbors2[0] && !hasThreeNeighbors2[1] && !hasThreeNeighbors2[2])
                    return false;
            }
        }

        else if (solution[i] == 2) {
            for (auto& neighbor : this->graph.getAdjacencyList(i)) {
                if (solution[neighbor] == 4) {
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
 * if 3, so the vertex was choosen
 * 
 *@return the vertex was choosen or not
 *
 */ 
bool AntColonyOptimization::delta(std::vector<int> solution, size_t vertex) {
    return solution[vertex] == 3 ? true : false; 
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
