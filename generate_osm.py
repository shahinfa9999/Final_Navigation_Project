import folium

def read_nodes_from_file(filename):
    nodes = []
    node_id_map = {}
    with open(filename, 'r') as file:
        for line in file:
            node_id, lat, lon = line.strip().split()
            node_id = int(node_id)
            lat = float(lat)
            lon = float(lon)
            nodes.append((node_id, lat, lon))
            node_id_map[node_id] = len(nodes) - 1
    return nodes, node_id_map

def read_ways_from_file(filename):
    ways = []
    with open(filename, 'r') as file:
        for line in file:
            way = tuple(map(int, line.strip().split()))
            ways.append(way)
    return ways

def read_path_from_file(filename):
    with open(filename, 'r') as file:
        path = list(map(int, file.read().strip().split()))
    return path

def create_map(nodes, ways, shortest_path, second_shortest_path, node_id_map):
    # Create a base map
    start_coords = (nodes[0][1], nodes[0][2])
    m = folium.Map(location=start_coords, zoom_start=14)

    # Add nodes to the map
    for node_id, lat, lon in nodes:
        folium.CircleMarker(location=(lat, lon), radius=2, color='blue', fill=True).add_to(m)

    # Add ways to the map
    for way in ways:
        node1 = nodes[node_id_map[way[0]]]
        node2 = nodes[node_id_map[way[1]]]
        folium.PolyLine(locations=[(node1[1], node1[2]), (node2[1], node2[2])], color='gray').add_to(m)

    # Add shortest path to the map
    shortest_path_coords = [(nodes[node_id_map[node]][1], nodes[node_id_map[node]][2]) for node in shortest_path]
    folium.PolyLine(locations=shortest_path_coords, color='green', weight=5, opacity=0.7).add_to(m)

    # Add second shortest path to the map
    second_shortest_path_coords = [(nodes[node_id_map[node]][1], nodes[node_id_map[node]][2]) for node in second_shortest_path]
    folium.PolyLine(locations=second_shortest_path_coords, color='red', weight=2, opacity=0.7).add_to(m)

    # Create the legend HTML
    legend_html = '''
     <div style="
     position: fixed; 
     bottom: 50px; left: 50px; width: 150px; height: 130px; 
     border:2px solid grey; z-index:9999; font-size:14px;
     background-color:white;
     ">
     &nbsp; <b>Legend</b> <br>
     &nbsp; <i class="fa fa-circle" style="color:blue"></i>&nbsp; Nodes <br>
     &nbsp; <i class="fa fa-minus" style="color:gray"></i>&nbsp; Ways <br>
     &nbsp; <i class="fa fa-minus" style="color:green"></i>&nbsp; Shortest Path <br>
     &nbsp; <i class="fa fa-minus" style="color:red"></i>&nbsp; Second Shortest Path <br>
     </div>
     '''

    # Add the legend to the map
    m.get_root().html.add_child(folium.Element(legend_html))

    # Save the map to an HTML file
    m.save('map.html')

# Read the nodes, ways, and node IDs from files
nodes, node_id_map = read_nodes_from_file("nodes.txt")
ways = read_ways_from_file("ways.txt")

# Read the shortest and second shortest paths from files
shortest_path = read_path_from_file("shortest_path.txt")
second_shortest_path = read_path_from_file("second_shortest_path.txt")

# Create the map with the paths
create_map(nodes, ways, shortest_path, second_shortest_path, node_id_map)