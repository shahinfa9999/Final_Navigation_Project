#include "Graph.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

void outputNodes(const Graph& graph, const std::string& filename) {
    std::ofstream file(filename);
    for (const auto& node : graph.getNodes()) {
        file << node.first << " " << node.second.lat << " " << node.second.lon << std::endl;
    }
}

void outputWays(const Graph& graph, const std::string& filename) {
    std::ofstream file(filename);
    for (const auto& pair : graph.getAdjacencyList()) {
        for (const auto& edge : pair.second) {
            file << pair.first << " " << edge.target << std::endl;
        }
    }
}

int main() {
    Graph graph;

    // Load graph data from the OSM file
    try {
        std::cout << "loading file" << std::endl;
        graph.loadFromOSMFile("meeker-roads.osm");
        std::cout << "loaded files" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error loading graph: " << e.what() << std::endl;
        return 1;
    }

    // Output nodes and ways to files
    outputNodes(graph, "nodes.txt");
    outputWays(graph, "ways.txt");

    // Input start and destination nodes
    long long startNode, endNode;
    std::cout << "Enter start node ID: ";
    std::cin >> startNode;
    std::cout << "Enter end node ID: ";
    std::cin >> endNode;

    if (graph.getNodes().find(startNode) == graph.getNodes().end() || graph.getNodes().find(endNode) == graph.getNodes().end()) {
        std::cerr << "Error: Invalid node ID encountered!" << std::endl;
        return 1;
    }

    // Find the shortest path
    std::vector<long long> shortestPath = graph.findShortestPath(startNode, endNode);
    std::ofstream shortestPathFile("shortest_path.txt");
    for (long long node : shortestPath) {
        shortestPathFile << node << " ";
    }
    shortestPathFile << std::endl;

    // Find the second shortest path
    std::vector<long long> secondShortestPath = graph.findSecondShortestPath(startNode, endNode);
    std::ofstream secondShortestPathFile("second_shortest_path.txt");
    for (long long node : secondShortestPath) {
        secondShortestPathFile << node << " ";
    }
    secondShortestPathFile << std::endl;

    return 0;
}