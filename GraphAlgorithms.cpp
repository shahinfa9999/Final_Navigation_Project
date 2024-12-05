#include "GraphAlgorithms.h"
#include <queue>
#include <unordered_map>
#include <limits>
#include <algorithm>

/*
std::vector<int> secondShortestPath(const Graph& graph, int start, int end) {
    std::vector<int> shortestPath = dijkstra(graph, start, end);
    double secondShortestDist = std::numeric_limits<double>::infinity();
    std::vector<int> secondShortestPath;

    for (size_t i = 0; i < shortestPath.size() - 1; ++i) {
        int u = shortestPath[i];
        int v = shortestPath[i + 1];

        Graph modifiedGraph = graph;
        auto& edges = modifiedGraph.getNeighbors(u);
        edges.erase(std::remove_if(edges.begin(), edges.end(), [&](Edge& e) { return e.target == v; }), edges.end());
        //edges.erase(std::remove_if(edges.begin(), edges.end(),
            //[&](Edge& e) { return e.target == v; }), edges.end());

        std::vector<int> newPath = dijkstra(modifiedGraph, start, end);
        double newDist = calculatePathDistance(modifiedGraph, newPath);

        if (newDist < secondShortestDist) {
            secondShortestDist = newDist;
            secondShortestPath = newPath;
        }
    }
    return secondShortestPath;
}
*/

// Yin's Algorithm
std::vector<int> secondShortestPath(const Graph& graph, double start, double end) {
    std::vector<int> shortestPath = dijkstra(graph, start, end);
    std::vector<int> secondShortestPath;
    double secondShortestDist = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < shortestPath.size() - 1; ++i) {
        int u = shortestPath[i];
        int v = shortestPath[i + 1];

        Graph modifiedGraph = graph;
        std::vector<Edge>& edges = modifiedGraph.getNeighbors(u);
        edges.erase(std::remove_if(edges.begin(), edges.end(), [&](Edge& e) { return e.target == v; }), edges.end());

        std::vector<int> newPath = dijkstra(modifiedGraph, start, end);
        double newDist = calculatePathDistance(modifiedGraph, newPath);

        if (newDist < secondShortestDist) {
            secondShortestDist = newDist;
            secondShortestPath = newPath;
        }
    }
    return secondShortestPath;
}

double calculatePathDistance(const Graph& graph, const std::vector<int>& path) {
    double totalDistance = 0.0;

    for (size_t i = 0; i < path.size() - 1; ++i) {
        int currentNode = path[i];
        int nextNode = path[i + 1];
        const std::vector<Edge>& neighbors = graph.getNeighbors(currentNode);

        bool edgeFound = false;
        for (const Edge& edge : neighbors) {
            if (edge.target == nextNode) {
                totalDistance += edge.weight;
                edgeFound = true;
                break;
            }
        }

        if (!edgeFound) {
            // If there is no edge between consecutive nodes, return an error value or handle it appropriately
            return -1.0; // Indicating that the path is invalid
        }
    }

    return totalDistance;
}

// Dijkstra's Algorithm
std::vector<int> dijkstra(const Graph& graph, int start, int end) {
    std::unordered_map<int, double> distances;
    std::unordered_map<int, int> previous;
    //auto compare = [&](int a, int b) { return distances[a] > distances[b]; };
    auto compare = [&](int a, int b) -> bool { return distances[a] > distances[b]; };
    std::priority_queue<int, std::vector<int>, decltype(compare)> minHeap(compare);

    for (const auto& pair : graph.getNeighbors(start)) {
        distances[pair.target] = std::numeric_limits<double>::infinity();
    }
    distances[start] = 0;
    minHeap.push(start);

    while (!minHeap.empty()) {
        int current = minHeap.top();
        minHeap.pop();

        if (current == end) break;

        for (const auto& edge : graph.getNeighbors(current)) {
            double newDist = distances[current] + edge.weight;
            if (newDist < distances[edge.target]) {
                distances[edge.target] = newDist;
                previous[edge.target] = current;
                minHeap.push(edge.target);
            }
        }
    }

    std::vector<int> path;
    for (int at = end; at != start; at = previous[at]) {
        path.push_back(at);
    }
    path.push_back(start);
    std::reverse(path.begin(), path.end());
    return path;
}

// Yin's Algorithm


