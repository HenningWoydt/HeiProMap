import os
import argparse

def save_graph(filepath, n, edges, vertex_weights=None):
    """
    Saves a graph in the METIS / CSRGraph compatible format.
    - n: number of vertices
    - edges: list of (u, v, w) tuples for undirected edges (1-indexed)
    """
    adj = {i: [] for i in range(1, n + 1)}
    
    m = len(edges)
    
    for u, v, w in edges:
        adj[u].append((v, w))
        adj[v].append((u, w))
        
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w') as f:
        # Header: n m fmt (011 means vertex and edge weights are included)
        f.write(f"{n} {m} 011\n")
        
        for i in range(1, n + 1):
            vw = vertex_weights[i] if vertex_weights else 1
            line = [str(vw)]
            # Sort neighbors by vertex ID just to be clean
            for neighbor, weight in sorted(adj[i]):
                line.append(str(neighbor))
                line.append(str(weight))
            f.write(" ".join(line) + "\n")

def gen_grid(width, height):
    n = width * height
    edges = []
    for y in range(height):
        for x in range(width):
            u = y * width + x + 1
            if x < width - 1:
                v = y * width + (x + 1) + 1
                edges.append((u, v, 1))
            if y < height - 1:
                v = (y + 1) * width + x + 1
                edges.append((u, v, 1))
    return n, edges

def gen_ring(n):
    edges = [(i, i + 1, 1) for i in range(1, n)]
    if n > 2:
        edges.append((n, 1, 1))
    return n, edges

def gen_hierarchical_squares(square_dim, grid_dim_x, grid_dim_y, w_intra=2, w_inter=18):
    """
    Generates a 2D grid of smaller 2D squares (chiplets/clusters).
    """
    n_x = grid_dim_x * square_dim
    n_y = grid_dim_y * square_dim
    
    n = n_x * n_y
    edges = []
    
    def get_id(x, y):
        return y * n_x + x + 1
        
    for y in range(n_y):
        for x in range(n_x):
            u = get_id(x, y)
            
            # Edge in X
            if x < n_x - 1:
                v = get_id(x + 1, y)
                w = w_inter if (x + 1) % square_dim == 0 else w_intra
                edges.append((u, v, w))
                
            # Edge in Y
            if y < n_y - 1:
                v = get_id(x, y + 1)
                w = w_inter if (y + 1) % square_dim == 0 else w_intra
                edges.append((u, v, w))
                
    return n, edges

def main():
    parser = argparse.ArgumentParser(description="Generate arbitrary topologies for WaverMap.")
    parser.add_argument("--type", choices=['grid', 'ring', 'cubes'], required=True, help="Topology type")
    parser.add_argument("--out", default="data/topologies/topology.graph", help="Output filepath")
    
    # Generic arguments
    parser.add_argument("--width", type=int, default=4, help="Grid width")
    parser.add_argument("--height", type=int, default=4, help="Grid height")
    parser.add_argument("--nodes", type=int, default=16, help="Number of nodes for ring")
    
    # Cubes/Squares specific arguments
    parser.add_argument("--cube-dim", type=int, default=8, help="Size of the small 2D squares (e.g. 8 for 8x8)")
    parser.add_argument("--grid-dim-x", type=int, default=2, help="How many small squares in X direction")
    parser.add_argument("--grid-dim-y", type=int, default=2, help="How many small squares in Y direction")
    parser.add_argument("--w-intra", type=int, default=2, help="Cost of links inside the small squares")
    parser.add_argument("--w-inter", type=int, default=18, help="Cost of links between the small squares")
    
    args = parser.parse_args()
    
    if args.type == 'grid':
        n, edges = gen_grid(args.width, args.height)
        print(f"Generating {args.width}x{args.height} grid topology...")
    elif args.type == 'ring':
        n, edges = gen_ring(args.nodes)
        print(f"Generating {args.nodes}-node ring topology...")
    elif args.type == 'cubes':
        n, edges = gen_hierarchical_squares(args.cube_dim, args.grid_dim_x, args.grid_dim_y, args.w_intra, args.w_inter)
        print(f"Generating hierarchical 2D topology (Square Size: {args.cube_dim}^2, Grid of Squares: {args.grid_dim_x}x{args.grid_dim_y})...")
        print(f"Weights -> Intra-Square: {args.w_intra}, Inter-Square: {args.w_inter}")
        
    save_graph(args.out, n, edges)
    print(f"Successfully saved {args.type} topology to '{args.out}' ({n} vertices, {len(edges)} edges).")

if __name__ == "__main__":
    main()
