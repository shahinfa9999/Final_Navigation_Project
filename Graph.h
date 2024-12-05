#ifndef GRAPH_H
#define GRAPH_H

#include <unordered_map>
#include <vector>
#include <string>
#include <queue>
#include <limits>

struct Node {
    double lat;
    double lon;
    // Other node attributes...
};

struct Edge {
    long long target;
    double weight;
};

class Graph {
public:
    void loadFromOSMFile(const std::string& filename);
    std::vector<long long> findSecondShortestPath(long long startNode, long long endNode);
    long long findNearestNode(double lat, double lon);
    void addEdge(long long source, long long target, double weight);
    std::vector<Edge>& getNeighbors(long long node);
    const std::vector<Edge>& getNeighbors(long long node) const;
    const std::unordered_map<long long, Node>& getNodes() const;
    const std::unordered_map<long long, std::vector<Edge>>& getAdjacencyList() const;

    // Calculate the total distance of a given path in the graph
    double calculatePathDistance(const std::vector<long long>& path) const;

    // Calculate the distance between two OSM nodes
    double calculateDistance(double lat1, double lon1, double lat2, double lon2) const;

    // Find the shortest path using Dijkstra's algorithm
    std::vector<long long> findShortestPath(long long startNode, long long endNode) const;

    // Find the second shortest path using Yin's algorithm
    std::vector<long long> findSecondShortestPath(long long startNode, long long endNode) const;

private:
    std::unordered_map<long long, std::vector<Edge>> adjacencyList;
    std::unordered_map<long long, Node> nodes;
};

#endif // GRAPH_H