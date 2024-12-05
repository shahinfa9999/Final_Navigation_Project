#ifndef GRAPH_ALGORITHMS_HPP
#define GRAPH_ALGORITHMS_HPP

#include "Graph.h"
#include <vector>

std::vector<int> dijkstra(const Graph& graph, int start, int end);
std::vector<int> secondShortestPath(const Graph& graph, int start, int end);
double calculatePathDistance(const Graph& graph, const std::vector<int>& path);


#endif
