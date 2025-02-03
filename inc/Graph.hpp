#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <random>
#include <algorithm>
#include <list>
#include <unordered_map>
#include <queue>
#include <sstream>
#include <fstream>

class Graph {
private:
    size_t order;
    size_t size;
    
    std::unordered_map<size_t, std::vector<size_t>> adjList;
    void addVertex(size_t source);
    void addEdge(size_t source, size_t destination);
    void DFSVisit(size_t u, std::vector<bool>& discovered, size_t& numberOfVertices, size_t& minDegree);
    size_t computeMaxVertexDegree() const;
    size_t computeMinVertexDegree() const;
    
public:	
    Graph(const std::string& filename);
    Graph(size_t order);	
    Graph(const Graph& graph);

    Graph() = default;
    ~Graph() = default;
 	
    size_t getOrder() const;
    size_t getSize() const;
    size_t getVertexDegree(size_t vertex) const;
    size_t getMinDegree() const;
    size_t getMaxDegree() const;
    const std::unordered_map<size_t, std::vector<size_t>>& getAdjacencyList() const;
    const std::vector<size_t>& getAdjacencyList(size_t vertex) const;
    
    bool edgeExists(size_t u, size_t v) const;
    
    bool vertexExists(size_t vertex) const;	

    std::vector<std::pair<int, int>> connectedComponents();
    
    void deleteVertex(size_t vertex);
    
    void deleteAdjacencyList(size_t vertex);
    
    friend std::ostream& operator<< (std::ostream& os, const Graph& graph);
};

#endif
