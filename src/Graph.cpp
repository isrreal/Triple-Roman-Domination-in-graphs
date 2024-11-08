#include "Graph.hpp"

Graph::Graph(size_t order, bool isDirected, float probabilityOfEdge) {
    this->isDirected = isDirected;
    this->order = order;
    this->size = 0;

    size_t connectedVertex = 0;
    float probability = 0.0;

    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_int_distribution<int> gap(0, order - 1);
    std::uniform_real_distribution<float> probabilityGap(0.0, 1.0);

    for (size_t i = 0; i < order; ++i)
        adjList[i] = {};

    for (size_t i = 0; i < order; ++i) {
        connectedVertex = gap(seed);

        while (i == connectedVertex)
            connectedVertex = gap(seed);

        if (!edgeExists(i, connectedVertex))
            addEdge(i, connectedVertex);

        for (size_t j = i + 1; j < order; ++j) {
            if (!edgeExists(i, j)) {
                probability = probabilityGap(seed);
                if (probabilityOfEdge >= probability)
                    addEdge(i, j);
           }
        }
    }
}

Graph::Graph(const std::string& filename, bool isDirected) : isDirected(isDirected), size(0) {
    std::ifstream file(filename, std::fstream::in);

    if (!file) {
        std::cerr << "Error opening the file!" << std::endl;
        throw std::runtime_error("File not found");
    }

    std::string line = "";
    size_t source, destination = 0;

    while (std::getline(file, line)) {
        std::stringstream ssEdges(line);

        while (ssEdges >> source >> destination) {
        	if (source != destination) {
        		addVertex(source);
	    		addVertex(destination);
	    		
		    	if ((!edgeExists(source, destination))) 		    		
		        	addEdge(source, destination);		       	
		    }
    	}
    }

    file.close();
}

Graph::Graph(const Graph& graph) {
    this->adjList = graph.adjList;
    this->order = graph.order;
    this->size = graph.size;
    this->isDirected = graph.isDirected;
}

Graph::Graph() {}

void Graph::addVertex(size_t source) {
	if (!vertexExists(source)) {
        adjList[source] = {};
        ++order;
    }
}

void Graph::addEdge(size_t source, size_t destination) {
    if (this->isDirected == false) {
        this->adjList[source].push_back(destination);
        this->adjList[destination].push_back(source);
        this->size += 2;      
    } 
    
    else {
        this->adjList[source].push_back(destination);
        this->size += 1;
    }
}

bool Graph::edgeExists(size_t u, size_t v) const {
    return std::find(adjList.at(u).begin(), adjList.at(u).end(), v) != adjList.at(u).end();
}

size_t Graph::getVertexDegree(size_t vertex) const {
    return const_cast<const std::unordered_map<size_t, std::list<size_t>>&>(adjList).at(vertex).size();
}

size_t Graph::getSize() const { return this->size; }

size_t Graph::getOrder() const { return this->order; }

std::unordered_map<size_t, std::list<size_t>> Graph::getAdjacencyList() const { return this->adjList; }

const std::list<size_t>& Graph::getAdjacencyList(size_t vertex) const { return this->adjList.at(vertex); }

bool Graph::vertexExists(size_t vertex) const { return adjList.find(vertex) != adjList.end(); }

/*
void Graph::breadthFirstSearch() {
    std::vector<bool> visited(this->order, false);
    std::queue<int> fila;
    visited[vertices[0]->identificator] = true;
    fila.push(vertices[0]->identificator);

    while (!fila.empty()) {
        int temp = fila.front();
        fila.pop();
        std::cout << temp << std::endl;

        for (const auto& it : edges[temp]) {
            if (!visited[it]) {
                visited[it] = true;
                fila.push(it);
            }
        }
    }
}

void Graph::depthFirstSearch() {
    std::vector<bool> visited(this->order, false);
    std::stack<int> stack;
    stack.push(vertices[0]->identificator);

    while (!stack.empty()) {
        int temp = stack.top();
        stack.pop();

        if (!visited[temp]) {
            visited[temp] = true;
            std::cout << temp << std::endl;

            for (const auto& neighbor : edges[temp]) {
                if (!visited[neighbor])
                    stack.push(neighbor);
            }
        }
    }
}

*/

void Graph::deleteAdjacencyList(size_t vertex) {
    if (adjList.find(vertex) == adjList.end())
    	return;

    std::queue<size_t> toDelete;
    toDelete.push(vertex);
    int currentVertex = -1;

    for (const auto& it: this->adjList[vertex])
    	toDelete.push(it);

    while (!toDelete.empty()) {
        currentVertex = toDelete.front();
        toDelete.pop();

        for (const auto& it: this->adjList[currentVertex])
        	this->adjList[it].remove(currentVertex);
       	deleteVertex(currentVertex);
    }
}

void Graph::deleteVertex(size_t vertex) {
    this->size -= this->adjList[vertex].size();
    this->adjList.erase(vertex);
    --this->order;
}

std::ostream& operator<< (std::ostream& os, const Graph& graph) {
    for (const auto& pair : graph.adjList) {
        size_t vertex = pair.first;  

        os << vertex << " ----> ";
        for (const auto& neighbor : pair.second) 
            os << neighbor << " ";  
        os << std::endl;
    }
    return os;
}
