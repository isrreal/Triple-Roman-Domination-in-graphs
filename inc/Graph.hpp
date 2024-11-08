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
    bool isDirected;  
    std::unordered_map<size_t, std::list<size_t>> adjList;
    void addVertex(size_t source);
    void addEdge(size_t source, size_t destination);
    Graph readGraph(const std::string& filename);
public:	

    Graph(size_t order, bool isDirected, float probabilityOfEdge);	
    Graph(const std::string& filename, bool isDirected);
    Graph(size_t order, bool isDirected);	
    Graph(const Graph& graph);
    Graph();
    ~Graph() = default;
 
    size_t getSize() const;
    size_t getOrder() const;
    size_t getVertexDegree(size_t vertex) const;
    
    std::unordered_map<size_t, std::list<size_t>> getAdjacencyList() const;
    
    const std::list<size_t>& getAdjacencyList(size_t vertex) const;
    
    bool edgeExists(size_t u, size_t v) const;
    
    bool vertexExists(size_t vertex) const;	

    void breadthFirstSearch();  
    
    void setVertexLabel(size_t vertex, int label);
    
    void setAdjacenciesLabel(size_t vertex, int label);
    
    void depthFirstSearch();
    
    void deleteAdjacencyList(size_t vertex);
    
    void deleteVertex(size_t vertex);
    
    friend std::ostream& operator<< (std::ostream& os, const Graph& graph);
};

#endif
