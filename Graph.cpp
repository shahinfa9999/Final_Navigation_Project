#include "Graph.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <queue>
#include <unordered_map>
#include <vector>
#include <algorithm>

const double EARTH_RADIUS_KM = 6371.0;

const std::unordered_map<long long, std::vector<Edge>>& Graph::getAdjacencyList() const {
    return adjacencyList;
}


double haversine(double lat1, double lon1, double lat2, double lon2) {
    double dlat = (lat2 - lat1) * M_PI / 180.0;
    double dlon = (lon2 - lon1) * M_PI / 180.0;

    lat1 = lat1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;

    double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
               std::sin(dlon / 2) * std::sin(dlon / 2) * std::cos(lat1) * std::cos(lat2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

    return EARTH_RADIUS_KM * c;
}

double Graph::calculateDistance(double lat1, double lon1, double lat2, double lon2) const {
    return haversine(lat1, lon1, lat2, lon2);
}

void Graph::addEdge(long long source, long long target, double weight) {
    adjacencyList[source].push_back({target, weight});
    adjacencyList[target].push_back({source, weight}); // Assuming undirected graph
    //std::cout << "Edge added: " << source << " -> " << target << " with weight " << weight << std::endl;
}

std::vector<Edge>& Graph::getNeighbors(long long node) {
    if (adjacencyList[node].empty())
    {
        static std::vector<Edge> empty;
        return empty;
    }
    
    return adjacencyList[node];
}
    //return adjacencyList[node];

const std::vector<Edge>& Graph::getNeighbors(long long node) const {
    auto it = adjacencyList.find(node);
    if (it == adjacencyList.end()) {
        throw std::out_of_range("Node not found in adjacency list");
    }
    return it->second;
}

const std::unordered_map<long long, Node>& Graph::getNodes() const {
    return nodes;
}

void Graph::loadFromOSMFile(const std::string& filename) {
    // Load the OSM file using libxml2
    xmlDocPtr doc = xmlReadFile(filename.c_str(), nullptr, 0);
    if (doc == nullptr) {
        throw std::runtime_error("Error loading OSM file");
    }

    // Parse the XML tree
    xmlNodePtr root = xmlDocGetRootElement(doc);
    for (xmlNodePtr node = root->children; node; node = node->next) {
        // Parse nodes and ways
        if (node->type == XML_ELEMENT_NODE) {
            // Parse nodes
            if (xmlStrcmp(node->name, BAD_CAST "node") == 0) {
                const char* idStr = (char*)xmlGetProp(node, BAD_CAST "id");
                const char* latStr = (char*)xmlGetProp(node, BAD_CAST "lat");
                const char* lonStr = (char*)xmlGetProp(node, BAD_CAST "lon");
    

                long long id = std::stoll(idStr);
                double lat = std::stod(latStr);
                double lon = std::stod(lonStr);
                nodes[id] = {lat, lon};
            } else if (xmlStrcmp(node->name, BAD_CAST "way") == 0) {
                std::vector<long long> wayNodes;
                for (xmlNodePtr nd = node->children; nd; nd = nd->next) {
                    if (nd->type == XML_ELEMENT_NODE && xmlStrcmp(nd->name, BAD_CAST "nd") == 0) {
                        const char* refStr = (char*)xmlGetProp(nd, BAD_CAST "ref");
                        long long ref = std::stoll(refStr);
                        wayNodes.push_back(ref);
                    }
                }
                // Add edges to the graph
                for (size_t i = 0; i < wayNodes.size() - 1; ++i) {
                    long long source = wayNodes[i];
                    long long target = wayNodes[i + 1];
                    double weight = calculateDistance(nodes[source].lat, nodes[source].lon, nodes[target].lat, nodes[target].lon);
                    addEdge(source, target, weight);
                }
            }
        }
    }

    xmlFreeDoc(doc);

    // Print adjacency list
    std::cout << "Adjacency list:" << std::endl;
    for (const auto& pair : adjacencyList) {
        std::cout << "Node " << pair.first << " has edges to: ";
        for (const auto& edge : pair.second) {
            std::cout << edge.target << " (weight " << edge.weight << "), ";
        }
        std::cout << std::endl;
    }
}

double Graph::calculatePathDistance(const std::vector<long long>& path) const {
    double totalDistance = 0.0;

    for (size_t i = 0; i < path.size() - 1; ++i) {
        long long current = path[i];
        long long next = path[i + 1];

        //std::cout << "Checking edge from " << current << " to " << next << std::endl;

        const auto& neighbors = getNeighbors(current);
        auto it = std::find_if(neighbors.begin(), neighbors.end(),
            [next](const Edge& edge) { return edge.target == next; });


        if (it == neighbors.end()) {
            throw std::runtime_error("Path is invalid: no edge between " +
                                     std::to_string(current) + " and " +
                                     std::to_string(next));
        }

        totalDistance += it->weight;
    }

    return totalDistance;
}

std::vector<long long> Graph::findShortestPath(long long startNode, long long endNode) const {
    std::unordered_map<long long, double> distances;
    std::unordered_map<long long, long long> previous;
    // min heap
    std::priority_queue<std::pair<double, long long>, std::vector<std::pair<double, long long>>, std::greater<std::pair<double, long long>>> pq;

    // Initialize distances
    for (const auto& node : nodes) {
        distances[node.first] = std::numeric_limits<double>::max();
    }
    distances[startNode] = 0.0;
    pq.push({0.0, startNode});

    // Dijkstra's algorithm loop
    while (!pq.empty()) {
        long long current = pq.top().second;
        pq.pop();

        if (current == endNode) {
            break;
        }

        for (const auto& edge : getNeighbors(current)) {
            long long neighbor = edge.target;
            double newDist = distances[current] + edge.weight;

            if (newDist < distances[neighbor]) {
                distances[neighbor] = newDist;
                previous[neighbor] = current;
                pq.push({newDist, neighbor});
                //std::cout << "Updating previous: " << neighbor << " <- " << current << std::endl;
            }
        }
    }

    // Handle invalid path
    if (previous.find(endNode) == previous.end()) {
        //std::cerr << "Error finding paths: No path found from " << startNode << " to " << endNode << std::endl;
        return {};
    }

    // Reconstruct the path
    std::vector<long long> path;
    long long at = endNode;
    while (at != startNode) {
        path.push_back(at);
        if (previous.find(at) == previous.end()) {
            //std::cerr << "Error: No predecessor found for node " << at << ", path reconstruction failed." << std::endl;
            return {};  // Return empty if no valid path
        }
        at = previous[at];
    }
    path.push_back(startNode);
    std::reverse(path.begin(), path.end());

    

    return path;
}


/*
std::vector<long long> Graph::findShortestPath(long long startNode, long long endNode) const {

    std::unordered_map<long long, double> distances;
    std::unordered_map<long long, long long> previous;
    std::priority_queue<std::pair<double, long long>, std::vector<std::pair<double, long long>>, std::greater<std::pair<double, long long>>> pq;

    for (const auto& node : nodes) {
        distances[node.first] = std::numeric_limits<double>::max();
    }
    distances[startNode] = 0.0;
    pq.push({0.0, startNode});

    while (!pq.empty()) {
        long long current = pq.top().second;
        pq.pop();

        if (current == endNode) {
            break;
        }

        for (const auto& edge : getNeighbors(current)) {
            long long neighbor = edge.target;
            double newDist = distances[current] + edge.weight;

            if (newDist < distances[neighbor]) {
                distances[neighbor] = newDist;
                previous[neighbor] = current;
                pq.push({newDist, neighbor});
                std::cout << "Updating previous: " << neighbor << " <- " << current << std::endl;
            }
        }
    }

    std::vector<long long> path;
    long long at = endNode;
    while (at != startNode) {
        if (previous.find(at) == previous.end()) {
            std::cerr << "Error finding paths: Path is invalid: no edge between " << at << " and " << previous[at] << std::endl;
            return {};
        }
        path.push_back(at);
        at = previous[at];
    }
    path.push_back(startNode);
    std::reverse(path.begin(), path.end());

    // Calculate path distance
    double pathDistance = calculatePathDistance(path);

    // Debugging output
    std::cout << "Path found: ";
    for (long long node : path) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
    std::cout << "Path distance: " << pathDistance << " units" << std::endl;

    return path;
}
*/
std::vector<long long> Graph::findSecondShortestPath(long long startNode, long long endNode) {
    std::vector<long long> shortestPath = findShortestPath(startNode, endNode);
    double shortestDistance = calculatePathDistance(shortestPath);

    std::vector<long long> secondShortestPath;
    double secondShortestDistance = std::numeric_limits<double>::max();

    // Backup the original adjacency list
    auto originalAdjacencyList = adjacencyList;

    for (size_t i = 0; i < shortestPath.size() - 1; ++i) {
        // Get the edge (u, v)
        long long u = shortestPath[i];
        long long v = shortestPath[i + 1];

        // Temporarily remove the edge (u, v)
        auto& neighbors = adjacencyList[u];
        auto it = std::remove_if(neighbors.begin(), neighbors.end(),
            [v](const Edge& edge) { return edge.target == v; });

        if (it != neighbors.end()) {
            neighbors.erase(it, neighbors.end());

            // Debugging output
            //std::cout << "Removed edge: " << u << " -> " << v << std::endl;

            // Find the shortest path without the edge (u, v)
            std::vector<long long> path = findShortestPath(startNode, endNode);

            // Restore the adjacency list from the backup
            adjacencyList = originalAdjacencyList;

            // Debugging output
            //std::cout << "Restored edge: " << u << " -> " << v << std::endl;

            if (!path.empty()) {
                double distance = calculatePathDistance(path);
                if (distance > shortestDistance && distance < secondShortestDistance) {
                    secondShortestDistance = distance;
                    secondShortestPath = path;
                }
            } else {
                std::cout << "Path invalid after removing edge: " << u << " -> " << v << std::endl;
            }
        }
    }

    

    return secondShortestPath;
}

long long Graph::findNearestNode(double lat, double lon) {
    long long nearestNode = -1;
    double minDistance = std::numeric_limits<double>::max();

    for (const auto& node : nodes) {
        double distance = std::sqrt(std::pow(node.second.lat - lat, 2) + std::pow(node.second.lon - lon, 2));
        if (distance < minDistance) {
            minDistance = distance;
            nearestNode = node.first;
        }
    }

    return nearestNode;
}