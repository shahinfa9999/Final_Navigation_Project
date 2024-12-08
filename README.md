# Navigation Project

## Problem Summary:

The objective is to design and implement a map navigation tool that computes the shortest route between two locations in a city, as well as the second shortest route. The city map can be represented as a graph where intersections (nodes) are connected by streets (edges) with weights corresponding to the distance or travel time. Dijkstra’s algorithm will be used to find the shortest path, while Yin's algorithm will be used to calculate the second shortest path.

## Data Structure and Algorithm:

- **Graph Representation**: The city map will be represented using an adjacency list, where each node (intersection) stores a list of neighboring nodes (adjacent intersections) and the edge weights (distances or travel times).
- **Dijkstra’s Algorithm**: This well-known algorithm will be implemented using a priority queue (min-heap) to efficiently find the shortest path from the start node to all other nodes in the graph.
- **Yin's Algorithm (Second Shortest Path)**: After finding the shortest path with Dijkstra’s, Yin’s algorithm will be used to determine the second shortest path. This algorithm works by removing one edge at a time from the shortest path and recalculating the shortest path on the modified graph. The second shortest path will be found by this approach, ensuring it is distinct from the original shortest path.

## Project Data Structure Implemented

This project implements a graph data structure to represent a road network, a priority queue to build the graph, and Yin's and Dijkstra's algorithm to find the second shortest and shortest path. The graph is constructed from OpenStreetMap (OSM) data and supports finding the shortest and second shortest paths between nodes using Dijkstra's algorithm and a modified version of it.

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
   ```

   Then compile the code:
   ```sh
   g++ -std=c++11 -I/usr/local/include -I/opt/homebrew/include -o navigation_tool main.cpp Graph.cpp GraphAlgorithms.cpp -lexpat -lbz2 -lz -lxml2
   ```

2. **Run the C++ Program**:

   ```sh
   ./navigation_tool
   ```

   This will generate the following files:
   - `nodes.txt`: Contains the nodes with their IDs and coordinates.
   - `ways.txt`: Contains the edges with their source and target node IDs.
   - `shortest_path.txt`: Contains the node IDs of the shortest path.
   - `second_shortest_path.txt`: Contains the node IDs of the second shortest path.

3. **Run the Python Script to Generate the Map**:

   ```sh
   python3 generate_osm.py
   ```

   This will create an HTML file named `map.html` with an interactive map showing the nodes, edges, shortest path, and second shortest path.

## Project Overview

The main goal of this project is to load a road network from an OSM file, construct a graph, and find the shortest and second shortest paths between two nodes. The project demonstrates the use of graph algorithms and data structures to solve real-world problems.

## Example run Video Link

https://drive.google.com/file/d/1nOCDOY9D2e5pTDgGEPSr3e5pkv-gk1bk/view?usp=sharing

## Conclusion

This project demonstrates the implementation of a graph data structure to represent a road network and the use of graph algorithms to find the shortest and second shortest paths. The interactive map generated using Folium provides a visual representation of the paths, making it easier to understand the results.
