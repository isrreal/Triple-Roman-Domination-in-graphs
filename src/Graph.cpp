#include "Graph.hpp"

Graph::Graph(const std::string& filename): 
	order(0), size(0) {
	
    std::ifstream file(filename);
    
    if (!file) {
        throw std::runtime_error("Error opening the file!");
    }
    
    size_t source, destination {0};
    std::string line {};
    
    while (std::getline(file, line)) {
        std::stringstream ssEdges(line);

        if (ssEdges >> source >> destination) {
            addVertex(source);
            addVertex(destination);
            addEdge(source, destination);
        }
        
    }
    
    file.close();
}

Graph::Graph(const Graph& graph): order(graph.order), size(graph.size), adjList(graph.adjList) {}

void Graph::addVertex(size_t source) {
	if (!vertexExists(source)) {
        adjList[source] = {};
        ++order;
    }
}

void Graph::addEdge(size_t source, size_t destination) {
	if (source == destination) { return; }
	
    this->adjList[source].push_back(destination);
    this->adjList[destination].push_back(source);            
    this->size += 1;
}

void Graph::DFSVisit(size_t u, std::vector<bool>& discovered, size_t& numberOfVertices, size_t& minDegree) {
    discovered[u] = true;

    if (getVertexDegree(u) < static_cast<size_t>(minDegree)) { 
        minDegree = getVertexDegree(u);
    }

    for (const size_t& v : getAdjacencyList(u)) {
        if (!discovered[v]) {
            ++numberOfVertices;
            if (getVertexDegree(v) < static_cast<size_t>(minDegree)) {
                minDegree = getVertexDegree(v);
            }
            
            DFSVisit(v, discovered, numberOfVertices, minDegree);
        }
    }
}

size_t Graph::computeMaxVertexDegree() const {
	if (getAdjacencyList().empty()) { return 0; }
	
	size_t maxDegree { getVertexDegree(getAdjacencyList().begin()->first) };
    size_t temp {0};
    
    for (const auto& [i, _] : getAdjacencyList()) {
        temp = getVertexDegree(i);
        if (maxDegree < temp) {
            maxDegree = temp;
        }
    }
            
    return maxDegree;
}


size_t Graph::computeMinVertexDegree() const {
	if (getAdjacencyList().empty()) { return 0; }
	
 	size_t minDegree { getVertexDegree(getAdjacencyList().begin()->first) };
    size_t temp {0};
    
    for (const auto& [i, _] : getAdjacencyList()) {
        temp = getVertexDegree(i);
        if (minDegree > temp) {
            minDegree = temp;
        }
    }
            
    return minDegree;
}

size_t Graph::getOrder() const { return this->order; }

size_t Graph::getSize() const { return this->size; }

size_t Graph::getVertexDegree(size_t vertex) const {
    return !vertexExists(vertex) ? 0 : const_cast<const std::unordered_map<size_t, std::vector<size_t>>&>(adjList).at(vertex).size();
}

size_t Graph::getMinDegree() const { return computeMinVertexDegree(); }

size_t Graph::getMaxDegree() const { return computeMaxVertexDegree(); }

const std::unordered_map<size_t, std::vector<size_t>>& Graph::getAdjacencyList() const { return this->adjList; }

const std::vector<size_t>& Graph::getAdjacencyList(size_t vertex) const { return this->adjList.at(vertex); }




bool Graph::edgeExists(size_t u, size_t v) const {
    return std::find(adjList.at(u).begin(), adjList.at(u).end(), v) != adjList.at(u).end();
}

bool Graph::vertexExists(size_t vertex) const { return adjList.find(vertex) != adjList.end(); }

std::vector<std::pair<int, int>> Graph::connectedComponents() {
    std::vector<bool> discovered(getOrder(), false);
    std::vector<std::pair<int, int>> components;

    for (const auto& [u, _] : adjList) {
        size_t vertex = u;
        if (!discovered[vertex]) {
            size_t numberOfVertices = 1;
            size_t minDegree = std::numeric_limits<int>::max();
            DFSVisit(vertex, discovered, numberOfVertices, minDegree);
            components.push_back({numberOfVertices, minDegree});
        }
    }
    return components;
}

void Graph::deleteVertex(size_t vertex) {
    for (auto& [v, neighbors] : this->adjList) {
        neighbors.erase(std::remove(neighbors.begin(), neighbors.end(), vertex), neighbors.end());
    }

    this->size -= this->adjList[vertex].size();
    this->adjList.erase(vertex);
    --this->order;
}

void Graph::deleteAdjacencyList(size_t vertex) {
    if (adjList.find(vertex) == adjList.end()) { return; }

    std::queue<size_t> toDelete;
    toDelete.push(vertex);
    
    int currentVertex {-1};

    for (const auto& it: this->adjList[vertex]) {
    	toDelete.push(it);
    }

    while (!toDelete.empty()) {
        currentVertex = toDelete.front();
        toDelete.pop();

        for (auto& [v, neighbors] : this->adjList) {
        	neighbors.erase(std::remove(neighbors.begin(), neighbors.end(), vertex), neighbors.end());
    	}
       	
       	deleteVertex(currentVertex);
    }
}

std::ostream& operator<< (std::ostream& os, const Graph& graph) {
    for (const auto& [u, v] : graph.adjList) {
        size_t vertex = u;  

        os << vertex << " ----> ";
        for (const auto& neighbor : v) { 
            os << neighbor << " ";  
    	}
    	        
        os << '\n';
    	
 	}
 	   
    return os;
}
