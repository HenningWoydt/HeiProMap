import networkx as nx
import matplotlib.pyplot as plt


def read_metis(filename):
    G = nx.Graph()

    with open(filename, "r") as f:
        lines = f.readlines()

    # Skip comments
    lines = [l.strip() for l in lines if not l.startswith("%") and l.strip()]

    # Header
    header = lines[0].split()
    n = int(header[0])
    m = int(header[1])
    fmt = int(header[2]) if len(header) > 2 else 0

    has_edge_weights = fmt == 1 or fmt == 11

    # Add nodes
    for i in range(n):
        G.add_node(i)

    # Read adjacency
    for u in range(n):
        parts = list(map(int, lines[u + 1].split()))
        i = 0

        while i < len(parts):
            v = parts[i] - 1  # METIS is 1-based

            if has_edge_weights:
                w = parts[i + 1]
                i += 2
            else:
                w = 1
                i += 1

            # Avoid duplicate edges
            if not G.has_edge(u, v):
                G.add_edge(u, v, weight=w)

    return G


def plot_graph(G):
    pos = nx.spring_layout(
        G,
        seed=42,
        iterations=500,
        k=2.0,
        weight="weight"
    )

    # Draw nodes and edges
    nx.draw(G, pos, with_labels=True, node_size=500)

    # Draw edge weights
    edge_labels = nx.get_edge_attributes(G, "weight")
    if edge_labels:
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

    plt.show()


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python draw_metis.py <file.graph>")
        exit(1)

    filename = sys.argv[1]
    G = read_metis(filename)
    plot_graph(G)