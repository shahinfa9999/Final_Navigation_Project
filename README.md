# Navigation Project

## Project Data Structure Implemented

This project implements a graph data structure to represent a road network. The graph is constructed from OpenStreetMap (OSM) data and supports finding the shortest and second shortest paths between nodes using Dijkstra's algorithm and a modified version of it.

## Explanation of the Data Structure

The graph is represented using an adjacency list, where each node (intersection) is connected to other nodes via edges (roads). Each edge has a weight representing the distance between the nodes. The nodes and edges are stored in the following structures:

- **Node**: Represents an intersection with latitude and longitude coordinates.
- **Edge**: Represents a road connecting two nodes with a weight (distance).
- **Graph**: Contains the adjacency list and methods to load data from an OSM file, add edges, and find shortest paths.

## How to Run the Project

### Prerequisites

- C++ compiler (e.g., g++)
- Python 3
- Folium library for Python (install using `pip install folium`)
- libxml2 library for parsing OSM files

### Instructions

1. **Compile the C++ Code**:

    Use `osmium-tool` to filter the OSM data for roads. You can install `osmium-tool` using Homebrew:
   ```sh
   brew install osmium-tool
   osmium tags-filter input.osm highway=* -o roads.osm

   ```sh
   g++ -std=c++11 -I/usr/local/include -I/opt/homebrew/include -o navigation_tool main.cpp Graph.cpp GraphAlgorithms.cpp -lexpat -lbz2 -lz -lxml2

2. **Run the C++ Program**:
    ./navigation_tool

This will generate the following files:
    - nodes.txt: Contains the nodes with their IDs and coordinates.
    - ways.txt: Contains the edges with their source and target node IDs.
    - shortest_path.txt: Contains the node IDs of the shortest path.
    - second_shortest_path.txt: Contains the node IDs of the second shortest path.
3. **Run the Python Script to Generate the Map**:
    python3 generate_osm.py

    This will create an HTML file named map.html with an interactive map showing the nodes, edges, shortest path, and second shortest path.

## Project Overview

The main goal of this project is to load a road network from an OSM file, construct a graph, and find the shortest and second shortest paths between two nodes. The project demonstrates the use of graph algorithms and data structures to solve real-world problems.

Screenshots
Graph Construction: Graph Construction Description: The graph is constructed from the OSM data, showing nodes and edges.

Shortest Path: Shortest Path Description: The shortest path between two nodes is highlighted in green.

Second Shortest Path: Second Shortest Path Description: The second shortest path between two nodes is highlighted in red.

Conclusion
This project demonstrates the implementation of a graph data structure to represent a road network and the use of graph algorithms to find the shortest and second shortest paths. The interactive map generated using Folium provides a visual representation of the paths, making it easier to understand the results.