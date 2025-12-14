"""
proto_03_build_subgraphs.py

Purpose: Extract hub neighborhood subgraphs for visualization
- Load top hubs from Step 2
- Extract 1-hop ego-graphs (neighborhoods) for top 3 hubs
- Build separate graphs for tumor and normal conditions
- Save as NetworkX pickle files for efficient loading

Expected Runtime: ~10 seconds
Memory Usage: ~200 MB
"""

import numpy as np
import networkx as nx
import yaml
import json
import logging
import pickle
from pathlib import Path
from datetime import datetime
import sys

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('output/prototype_execution.log', mode='a'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


def load_config(config_path='config_proto.yaml'):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def load_correlation_matrix(npz_path):
    """Load correlation matrix from NPZ file."""
    data = np.load(npz_path, allow_pickle=True)
    return data['matrix'], data['genes']


def load_top_hubs(hubs_json_path, top_k=3):
    """
    Load top K hubs from JSON file.
    
    Args:
        hubs_json_path: Path to top_hubs.json from Step 2
        top_k: Number of top hubs to extract
    
    Returns:
        List of hub dictionaries
    """
    try:
        with open(hubs_json_path, 'r') as f:
            all_hubs = json.load(f)
        
        top_hubs = all_hubs[:top_k]
        logger.info(f"Loaded top {top_k} hubs:")
        for hub in top_hubs:
            logger.info(f"  {hub['rank']}. {hub['gene_symbol']} (Δ = {hub['delta_connectivity']:.1f})")
        
        return top_hubs
        
    except FileNotFoundError:
        logger.error(f"Hub file not found: {hubs_json_path}")
        logger.error("Please run proto_02_differential.py first")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading hubs: {e}")
        sys.exit(1)


def build_correlation_network(corr_matrix, genes, threshold=0.7):
    """
    Build weighted network from correlation matrix.
    
    Args:
        corr_matrix: Correlation matrix (genes × genes)
        genes: Gene identifiers
        threshold: Minimum |correlation| for edge inclusion
    
    Returns:
        NetworkX Graph object
    """
    logger.info(f"Building network with threshold |r| ≥ {threshold}...")
    
    G = nx.Graph()
    
    # Add all genes as nodes
    for gene in genes:
        G.add_node(gene)
    
    # Add edges for correlations above threshold
    n_genes = len(genes)
    edge_count = 0
    
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            corr_value = corr_matrix[i, j]
            if abs(corr_value) >= threshold:
                G.add_edge(genes[i], genes[j], weight=abs(corr_value))
                edge_count += 1
    
    logger.info(f"Network built: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    return G


def extract_hub_subgraph(G, hub_gene, max_neighbors=50):
    """
    Extract 1-hop ego-graph (neighborhood) centered on hub gene.
    
    Args:
        G: Full network (NetworkX Graph)
        hub_gene: Hub gene identifier
        max_neighbors: Maximum neighbors to include (for visual clarity)
    
    Returns:
        Subgraph containing hub and its neighbors
    """
    if hub_gene not in G:
        logger.warning(f"Hub gene {hub_gene} not found in network")
        return nx.Graph()
    
    # Get neighbors
    neighbors = list(G.neighbors(hub_gene))
    
    # If too many neighbors, select top K by edge weight
    if len(neighbors) > max_neighbors:
        logger.info(f"  Limiting {len(neighbors)} neighbors to top {max_neighbors} by edge weight")
        
        # Get edge weights
        neighbor_weights = [(n, G[hub_gene][n]['weight']) for n in neighbors]
        neighbor_weights.sort(key=lambda x: x[1], reverse=True)
        
        # Keep top K
        neighbors = [n for n, w in neighbor_weights[:max_neighbors]]
    
    # Extract subgraph (hub + neighbors)
    subgraph_nodes = [hub_gene] + neighbors
    subgraph = G.subgraph(subgraph_nodes).copy()
    
    logger.info(f"  Extracted subgraph: {subgraph.number_of_nodes()} nodes, {subgraph.number_of_edges()} edges")
    
    return subgraph


def save_subgraph(subgraph, output_path):
    """Save subgraph as NetworkX pickle file."""
    try:
        with open(output_path, 'wb') as f:
            pickle.dump(subgraph, f, protocol=pickle.HIGHEST_PROTOCOL)
        
        file_size_kb = output_path.stat().st_size / 1e3
        logger.info(f"Saved to: {output_path} ({file_size_kb:.1f} KB)")
        
    except Exception as e:
        logger.error(f"Error saving subgraph: {e}")
        sys.exit(1)


def main():
    """Main execution function."""
    logger.info("")
    logger.info("="*60)
    logger.info("PROTOTYPE STEP 3: BUILD HUB SUBGRAPHS")
    logger.info("="*60)
    logger.info("")
    
    start_time = datetime.now()
    
    # Load configuration
    config = load_config()
    
    # Extract parameters
    output_dir = Path(config['paths']['output_base'])
    corr_dir = output_dir / 'correlation_matrices'
    diff_dir = output_dir / 'differential_results'
    subgraph_dir = output_dir / 'hub_subgraphs'
    subgraph_dir.mkdir(parents=True, exist_ok=True)
    
    corr_threshold = config['network']['correlation_threshold']
    max_neighbors = config['subgraph']['max_neighbors']
    top_k_visualize = config['hubs']['top_k_visualize']
    
    # Step 1: Load top hubs
    logger.info("Step 1: Loading top hubs from Step 2")
    logger.info("-" * 40)
    hubs_path = diff_dir / 'top_hubs.json'
    top_hubs = load_top_hubs(hubs_path, top_k=top_k_visualize)
    logger.info("")
    
    # Step 2: Load correlation matrices
    logger.info("Step 2: Loading correlation matrices")
    logger.info("-" * 40)
    
    tumor_corr_path = corr_dir / 'tumor_corr_proto.npz'
    normal_corr_path = corr_dir / 'normal_corr_proto.npz'
    
    logger.info("Loading tumor correlation...")
    tumor_corr, tumor_genes = load_correlation_matrix(tumor_corr_path)
    logger.info(f"  Tumor: {tumor_corr.shape}")
    
    logger.info("Loading normal correlation...")
    normal_corr, normal_genes = load_correlation_matrix(normal_corr_path)
    logger.info(f"  Normal: {normal_corr.shape}")
    logger.info("")
    
    # Step 3: Build full networks
    logger.info("Step 3: Building full correlation networks")
    logger.info("-" * 40)
    
    logger.info("Building tumor network...")
    G_tumor = build_correlation_network(tumor_corr, tumor_genes, corr_threshold)
    
    logger.info("Building normal network...")
    G_normal = build_correlation_network(normal_corr, normal_genes, corr_threshold)
    logger.info("")
    
    # Step 4: Extract hub subgraphs
    logger.info("Step 4: Extracting hub subgraphs")
    logger.info("-" * 40)
    
    subgraph_count = 0
    
    for hub_info in top_hubs:
        hub_gene = hub_info['gene']
        hub_symbol = hub_info['gene_symbol']
        
        logger.info(f"Processing hub {hub_info['rank']}: {hub_symbol}")
        
        # Extract tumor subgraph
        logger.info("  Tumor condition:")
        subgraph_tumor = extract_hub_subgraph(G_tumor, hub_gene, max_neighbors)
        tumor_path = subgraph_dir / f"{hub_gene}_tumor.pkl"
        save_subgraph(subgraph_tumor, tumor_path)
        subgraph_count += 1
        
        # Extract normal subgraph
        logger.info("  Normal condition:")
        subgraph_normal = extract_hub_subgraph(G_normal, hub_gene, max_neighbors)
        normal_path = subgraph_dir / f"{hub_gene}_normal.pkl"
        save_subgraph(subgraph_normal, normal_path)
        subgraph_count += 1
        
        logger.info("")
    
    # Summary
    total_time = (datetime.now() - start_time).total_seconds()
    
    logger.info("="*60)
    logger.info("STEP 3 COMPLETE")
    logger.info("="*60)
    logger.info(f"Total execution time: {total_time:.1f} seconds")
    logger.info("")
    logger.info(f"Subgraphs created: {subgraph_count} files")
    logger.info(f"  • {top_k_visualize} hubs × 2 conditions (tumor/normal)")
    logger.info(f"  • Saved to: {subgraph_dir}")
    logger.info("")
    logger.info("Next step: python scripts/proto_04_visualize.py")
    logger.info("="*60)


if __name__ == "__main__":
    main()