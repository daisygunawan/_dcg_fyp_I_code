"""
proto_04_visualize.py

Purpose: Generate core visualizations for prototype validation
- Volcano plot (delta_r vs significance) using optimized dataset
- Hub ranking bar chart (top 20 hubs)
- Connectivity comparison scatter (tumor vs normal)
- Hub comparison grid (2×3 network visualization)

Expected Runtime: 3-4 minutes (with optimized volcano dataset)
Memory Usage: ~1 GB
"""

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import yaml
import json
import pickle
import logging
from pathlib import Path
from datetime import datetime
import sys

# Configure matplotlib
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'

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


def load_subgraph(pkl_path):
    """Load NetworkX subgraph from pickle file."""
    with open(pkl_path, 'rb') as f:
        return pickle.load(f)

def plot_volcano(df_volcano, output_path, fdr_threshold=0.05, min_effect_size=0.1):
    """
    Create TWO-PANEL volcano plot with optimized layout.
    
    Top panel: High significance (y >= 20) - shows strong rewiring signal
    Bottom panel: Low significance (y < 20) - shows stable relationships with gray points visible
    
    Args:
        df_volcano: DataFrame with volcano dataset
        output_path: Path to save figure
        fdr_threshold: FDR cutoff
        min_effect_size: Effect size cutoff
    """
    logger.info("Creating two-panel volcano plot...")
    logger.info(f"  Using {len(df_volcano):,} sampled points")
    
    # Prepare data
    df = df_volcano.copy()
    df['-log10_p_fdr'] = -np.log10(df['p_fdr'].clip(lower=1e-300))
    df['significant'] = df['is_significant']
    
    sig = df[df['significant']]
    non_sig = df[~df['significant']]
    
    logger.info(f"  Significant: {len(sig):,} ({100*len(sig)/len(df):.1f}%)")
    logger.info(f"  Non-significant: {len(non_sig):,} ({100*len(non_sig)/len(df):.1f}%)")
    
    # ==============================================================
    # CREATE TWO-PANEL FIGURE WITH OPTIMIZED LAYOUT
    # ==============================================================
    
    # Panel split threshold
    y_split = 20  # Split at -log10(p) = 20
    
    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1, 
        figsize=(12, 9),  # Reduced from 10 to 9 for less white space
        sharex=True,
        gridspec_kw={
            'height_ratios': [1, 2],  # REVERSED: bottom now bigger (was [2, 1])
            'hspace': 0.15  # Increased from 0.05 for more breathing room
        }
    )
    
    # ==============================================================
    # TOP PANEL: High significance region (y >= 20)
    # ==============================================================
    
    # Filter to high-significance region
    sig_top = sig[sig['-log10_p_fdr'] >= y_split]
    non_sig_top = non_sig[non_sig['-log10_p_fdr'] >= y_split]
    
    logger.info(f"  Top panel: {len(sig_top):,} sig + {len(non_sig_top):,} non-sig")
    
    # Plot points
    if len(sig_top) > 0:
        ax_top.scatter(sig_top['delta_r'], sig_top['-log10_p_fdr'],
                      s=3, alpha=0.7, c='crimson',  # Slightly larger for visibility
                      label=f'Significantly Rewired ({len(sig_top):,})',
                      rasterized=True, zorder=2)
    
    if len(non_sig_top) > 0:
        ax_top.scatter(non_sig_top['delta_r'], non_sig_top['-log10_p_fdr'],
                      s=4, alpha=0.8, c='dimgray',
                      label=f'Non-significant ({len(non_sig_top):,})',
                      rasterized=True, zorder=3)
    
    # Threshold lines (top panel)
    fdr_y = -np.log10(fdr_threshold)
    if fdr_y >= y_split:
        ax_top.axhline(fdr_y, color='blue', linestyle='--', linewidth=1.5, alpha=0.7, zorder=1)
    
    ax_top.axvline(min_effect_size, color='green', linestyle='--', linewidth=1.5, alpha=0.7, zorder=1)
    ax_top.axvline(-min_effect_size, color='green', linestyle='--', linewidth=1.5, alpha=0.7, zorder=1)
    
    # Labels (top panel)
    ax_top.set_ylabel(r'$-\log_{10}(\text{FDR p-value})$', fontsize=11, fontweight='bold')
    ax_top.set_title('High Significance Region (Strong Rewiring Signal)', 
                     fontsize=11, fontweight='bold', color='darkred', pad=10)
    
    # Legend in top-left to avoid crowding
    ax_top.legend(loc='upper left', framealpha=0.95, fontsize=9, markerscale=2)
    ax_top.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, zorder=0)
    
    # Set y-limits for top panel - TIGHTER to reduce white space
    y_max = df['-log10_p_fdr'].max()
    ax_top.set_ylim(y_split - 1, y_max + 5)  # Start just below y_split, less top padding
    
    # Arrow annotations for top panel
    arrow_y = y_max * 0.70  # Adjusted position
    ax_top.annotate('Loss of correlation\nin cancer', 
                   xy=(-0.5, arrow_y), xytext=(-0.8, y_max * 0.85),
                   ha='center', fontsize=9, color='darkred',
                   arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5),
                   zorder=4)
    
    ax_top.annotate('Gain of correlation\nin cancer', 
                   xy=(0.5, arrow_y), xytext=(0.8, y_max * 0.85),
                   ha='center', fontsize=9, color='darkred',
                   arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5),
                   zorder=4)
    
    # Remove x-axis tick labels from top panel (shared x-axis)
    ax_top.tick_params(labelbottom=False)
    
    # ==============================================================
    # BOTTOM PANEL: Low significance region (y < 20)
    # ==============================================================
    
    # Filter to low-significance region - INCLUDE boundary points
    sig_bottom = sig[sig['-log10_p_fdr'] < y_split]
    non_sig_bottom = non_sig[non_sig['-log10_p_fdr'] < y_split]
    
    logger.info(f"  Bottom panel: {len(sig_bottom):,} sig + {len(non_sig_bottom):,} non-sig")
    
    # Plot RED points FIRST (bottom layer)
    if len(sig_bottom) > 0:
        ax_bottom.scatter(sig_bottom['delta_r'], sig_bottom['-log10_p_fdr'],
                         s=1.5, alpha=0.4, c='crimson',
                         label=f'Marginally Sig. ({len(sig_bottom):,})',
                         rasterized=True, zorder=1)
    
    # Plot GRAY points SECOND (top layer) - NOW VISIBLE!
    if len(non_sig_bottom) > 0:
        ax_bottom.scatter(non_sig_bottom['delta_r'], non_sig_bottom['-log10_p_fdr'],
                         s=3, alpha=0.8, c='dimgray',
                         label=f'Non-significant ({len(non_sig_bottom):,})',
                         rasterized=True, zorder=2)
    
    # Threshold lines (bottom panel)
    if fdr_y < y_split:
        ax_bottom.axhline(fdr_y, color='blue', linestyle='--', 
                         linewidth=2, label=f'FDR = {fdr_threshold}', zorder=3)
    
    ax_bottom.axvline(min_effect_size, color='green', linestyle='--', 
                     linewidth=2, alpha=0.7, zorder=3)
    ax_bottom.axvline(-min_effect_size, color='green', linestyle='--', 
                     linewidth=2, alpha=0.7, label=rf'$|\Delta r| = {min_effect_size}$', zorder=3)
    
    # Labels (bottom panel)
    ax_bottom.set_xlabel(r'$\Delta r$ (Correlation Change: Tumor - Normal)', 
                        fontsize=12, fontweight='bold')
    ax_bottom.set_ylabel(r'$-\log_{10}(\text{FDR p-value})$', fontsize=11, fontweight='bold')
    ax_bottom.set_title('Low Significance Region (Stable Relationships)', 
                       fontsize=11, fontweight='bold', color='dimgray', pad=15)  # Increased pad
    ax_bottom.legend(loc='upper right', framealpha=0.95, fontsize=9, markerscale=2)
    ax_bottom.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, zorder=0)
    
    # Set y-limits for bottom panel - EXTEND slightly above y_split to show continuity
    ax_bottom.set_ylim(-0.5, y_split + 0.5)  # Extended upper limit to show transition
    
    # ==============================================================
    # Overall title and styling
    # ==============================================================
    
    fig.suptitle('Differential Co-Expression Volcano Plot (Two-Panel View)\n' +
                 f'Smart Sampling: {len(df):,} points ({len(sig):,} sig + {len(non_sig):,} non-sig) from 89.6M total',
                 fontsize=13, fontweight='bold', y=0.98)
    
    # Add interpretation note
    note_text = (
        f"Two-panel design: TOP = high-confidence rewiring (y≥{y_split}), "
        f"BOTTOM = low-significance region (y<{y_split}).\n"
        f"Gray points (non-significant) are visible in bottom panel. "
        f"Panels overlap slightly at y={y_split} to show continuity."
    )
    fig.text(0.5, 0.01, note_text, 
            transform=fig.transFigure, fontsize=8, 
            ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', 
                     alpha=0.8, edgecolor='blue', linewidth=1))
    
    # Set shared x-axis limits
    x_min, x_max = df['delta_r'].min(), df['delta_r'].max()
    ax_bottom.set_xlim(x_min - 0.05, x_max + 0.05)
    
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Two-panel volcano plot saved: {output_path}")
    logger.info(f"  - Top panel: {len(sig_top):,} points (y ≥ {y_split})")
    logger.info(f"  - Bottom panel: {len(sig_bottom) + len(non_sig_bottom):,} points (y < {y_split})")
    logger.info(f"  - Gray points in bottom: {len(non_sig_bottom):,} (FULLY VISIBLE!)")


def plot_hub_ranking(df_connectivity, output_path, top_k=20):
    """
    Create horizontal bar chart of top K hubs by ABSOLUTE delta_connectivity.
    Color by direction of change (gain/loss).
    
    Args:
        df_connectivity: DataFrame with connectivity values
        output_path: Path to save figure
        top_k: Number of top hubs to display
    """
    logger.info(f"Creating hub ranking plot (top {top_k})...")
    
    # Select top K hubs by absolute delta
    top_hubs = df_connectivity.nlargest(top_k, 'delta_connectivity').copy()
    
    # Extract gene symbols
    top_hubs['symbol'] = top_hubs['gene'].apply(
        lambda x: x.split('|')[1] if '|' in x else x
    )
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Color by direction of change
    colors = ['lightcoral' if d == 'loss' else 'lightsteelblue' for d in top_hubs['connectivity_direction']]
    
    # Horizontal bar chart (reversed order for top-to-bottom)
    y_pos = np.arange(len(top_hubs))
    
    # Plot bars with direction-based coloring
    bars = ax.barh(y_pos, top_hubs['delta_connectivity'], color=colors, alpha=0.8)
    
    # Add value labels on bars
    for i, (bar, delta) in enumerate(zip(bars, top_hubs['delta_connectivity'])):
        width = bar.get_width()
        direction = top_hubs['connectivity_direction'].iloc[i]
        symbol = top_hubs['symbol'].iloc[i]
        
        # Position text inside or outside bar based on width
        if width > max(top_hubs['delta_connectivity']) * 0.1:  # If bar is wide enough
            ax.text(width * 0.95, bar.get_y() + bar.get_height()/2, 
                   f'{symbol} ({direction})', 
                   ha='right', va='center', fontsize=9, fontweight='bold',
                   color='darkred' if direction == 'loss' else 'darkblue')
        else:
            ax.text(width + max(top_hubs['delta_connectivity']) * 0.01, 
                   bar.get_y() + bar.get_height()/2, 
                   f'{symbol} ({direction})', 
                   ha='left', va='center', fontsize=9, fontweight='bold',
                   color='darkred' if direction == 'loss' else 'darkblue')
    
    # Labels and styling
    ax.set_yticks(y_pos)
    ax.set_yticklabels([])  # Remove y-axis labels since we have inline labels
    ax.set_xlabel('Absolute Δ Connectivity |Tumor - Normal|', fontsize=12, fontweight='bold')
    ax.set_ylabel('Gene Symbol', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {top_k} Rewired Hub Genes\n(by Absolute Connectivity Change)', 
                 fontsize=14, fontweight='bold')
    
    # Add legend for direction
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='lightsteelblue', label='Connectivity Gain (Tumor > Normal)'),
        Patch(facecolor='lightcoral', label='Connectivity Loss (Tumor < Normal)')
    ]
    ax.legend(handles=legend_elements, loc='lower right', framealpha=0.9)
    
    ax.grid(True, alpha=0.3, linestyle=':', axis='x')
    
    # Invert y-axis so rank 1 is at top
    ax.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Hub ranking plot saved to: {output_path}")

def plot_connectivity_scatter(df_connectivity, output_path):
    """
    Create scatter plot: tumor connectivity vs normal connectivity.
    
    Args:
        df_connectivity: DataFrame with connectivity values
        output_path: Path to save figure
    """
    logger.info("Creating connectivity scatter plot...")
    
    df = df_connectivity.copy()
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Scatter plot
    ax.scatter(df['normal_connectivity'], df['tumor_connectivity'], 
               s=10, alpha=0.4, c='steelblue')
    
    # Add diagonal line (equal connectivity)
    max_val = max(df['normal_connectivity'].max(), df['tumor_connectivity'].max())
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='Equal connectivity')
    
    # Highlight top 3 hubs
    top3 = df.nlargest(3, 'delta_connectivity')
    ax.scatter(top3['normal_connectivity'], top3['tumor_connectivity'], 
               s=100, c='red', marker='*', edgecolors='black', linewidths=1,
               label='Top 3 hubs', zorder=5)
    
    # Annotate top 3
    for _, row in top3.iterrows():
        symbol = row['gene'].split('|')[1] if '|' in row['gene'] else row['gene']
        ax.annotate(symbol, 
                   (row['normal_connectivity'], row['tumor_connectivity']),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=10, fontweight='bold', color='darkred')
    
    # Labels and styling
    ax.set_xlabel('Normal Connectivity (# edges)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Tumor Connectivity (# edges)', fontsize=12, fontweight='bold')
    ax.set_title('Gene Connectivity: Tumor vs Normal', fontsize=14, fontweight='bold')
    ax.legend(loc='upper left', framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Connectivity scatter saved to: {output_path}")


def plot_hub_comparison_grid(top_hubs, subgraph_dir, output_path):
    """
    Create 2×3 grid comparing tumor vs normal subgraphs for top 3 hubs.
    THIS IS THE VISUAL PROOF - the most important figure!
    
    Args:
        top_hubs: List of top hub dictionaries
        subgraph_dir: Directory containing subgraph pickle files
        output_path: Path to save figure
    """
    logger.info("Creating hub comparison grid (2×3 layout)...")
    logger.info("  This is THE VISUAL PROOF of network rewiring!")
    
    n_hubs = len(top_hubs)
    fig, axes = plt.subplots(n_hubs, 2, figsize=(14, 5*n_hubs))
    
    if n_hubs == 1:
        axes = axes.reshape(1, -1)
    
    for idx, hub_info in enumerate(top_hubs):
        hub_gene = hub_info['gene']
        hub_symbol = hub_info['gene_symbol']
        
        logger.info(f"  Plotting hub {idx+1}: {hub_symbol}")
        
        # Load subgraphs
        tumor_path = subgraph_dir / f"{hub_gene}_tumor.pkl"
        normal_path = subgraph_dir / f"{hub_gene}_normal.pkl"
        
        G_tumor = load_subgraph(tumor_path)
        G_normal = load_subgraph(normal_path)
        
        # --- Tumor subgraph (left column) ---
        ax_tumor = axes[idx, 0]
        
        if G_tumor.number_of_nodes() > 0:
            # Compute layout
            pos_tumor = nx.spring_layout(G_tumor, k=0.5, iterations=50, seed=42)
            
            # Draw network
            nx.draw_networkx_edges(G_tumor, pos_tumor, alpha=0.2, width=0.5, 
                                  edge_color='gray', ax=ax_tumor)
            nx.draw_networkx_nodes(G_tumor, pos_tumor, node_size=50, 
                                  node_color='lightcoral', alpha=0.6, ax=ax_tumor)
            
            # Highlight hub
            nx.draw_networkx_nodes(G_tumor, pos_tumor, nodelist=[hub_gene],
                                  node_size=200, node_color='darkred', 
                                  node_shape='*', ax=ax_tumor)
            
            ax_tumor.set_title(f'{hub_symbol} - TUMOR\n'
                             f'{G_tumor.number_of_nodes()} nodes, '
                             f'{G_tumor.number_of_edges()} edges',
                             fontsize=12, fontweight='bold', color='darkred')
        else:
            ax_tumor.text(0.5, 0.5, 'No connections\nabove threshold', 
                        ha='center', va='center', fontsize=12, color='red')
            ax_tumor.set_title(f'{hub_symbol} - TUMOR (Disconnected)', 
                             fontsize=12, fontweight='bold', color='darkred')
        
        ax_tumor.axis('off')
        
        # --- Normal subgraph (right column) ---
        ax_normal = axes[idx, 1]
        
        if G_normal.number_of_nodes() > 0:
            # Compute layout
            pos_normal = nx.spring_layout(G_normal, k=0.5, iterations=50, seed=42)
            
            # Draw network
            nx.draw_networkx_edges(G_normal, pos_normal, alpha=0.2, width=0.5,
                                  edge_color='gray', ax=ax_normal)
            nx.draw_networkx_nodes(G_normal, pos_normal, node_size=50,
                                  node_color='lightblue', alpha=0.6, ax=ax_normal)
            
            # Highlight hub
            nx.draw_networkx_nodes(G_normal, pos_normal, nodelist=[hub_gene],
                                  node_size=200, node_color='darkblue',
                                  node_shape='*', ax=ax_normal)
            
            ax_normal.set_title(f'{hub_symbol} - NORMAL\n'
                              f'{G_normal.number_of_nodes()} nodes, '
                              f'{G_normal.number_of_edges()} edges',
                              fontsize=12, fontweight='bold', color='darkblue')
        else:
            ax_normal.text(0.5, 0.5, 'No connections\nabove threshold',
                         ha='center', va='center', fontsize=12, color='blue')
            ax_normal.set_title(f'{hub_symbol} - NORMAL (Disconnected)',
                              fontsize=12, fontweight='bold', color='darkblue')
        
        ax_normal.axis('off')
    
    # Overall title
    fig.suptitle('Hub Gene Network Comparison: Tumor vs Normal\n'
                 '(Red = Tumor Sparse, Blue = Normal Dense)',
                 fontsize=16, fontweight='bold', y=0.995)
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Hub comparison grid saved to: {output_path}")
    logger.info("  KEY FIGURE: This is the main visualization for demonstration!")


def main():
    """Main execution function."""
    logger.info("")
    logger.info("="*60)
    logger.info("PROTOTYPE STEP 4: GENERATE VISUALIZATIONS")
    logger.info("="*60)
    logger.info("")
    
    start_time = datetime.now()
    
    # Load configuration
    config = load_config()
    
    # Extract parameters
    output_dir = Path(config['paths']['output_base'])
    diff_dir = output_dir / 'differential_results'
    subgraph_dir = output_dir / 'hub_subgraphs'
    viz_dir = output_dir / 'visualizations'
    viz_dir.mkdir(parents=True, exist_ok=True)
    
    fdr_threshold = config['differential']['fdr_threshold']
    min_effect_size = config['differential']['min_effect_size']
    top_k_visualize = config['hubs']['top_k_visualize']
    
    # Step 1: Load data
    logger.info("Step 1: Loading results from Steps 2-3")
    logger.info("-" * 40)
    
    # Load OPTIMIZED VOLCANO dataset (not just significant pairs)
    volcano_path = diff_dir / 'differential_pairs_volcano.tsv'
    df_volcano = pd.read_csv(volcano_path, sep='\t')
    logger.info(f"Loaded optimized volcano dataset: {len(df_volcano):,} points")
    logger.info(f"  - Significant in dataset: {df_volcano['is_significant'].sum():,}")
    logger.info(f"  - Non-significant in dataset: {(~df_volcano['is_significant']).sum():,}")
    
    # Also load significant pairs (for reference)
    sig_path = diff_dir / 'differential_pairs_significant.tsv'
    df_sig = pd.read_csv(sig_path, sep='\t')
    logger.info(f"Loaded all significant pairs: {len(df_sig):,} (for reference)")
    
    # Load connectivity results
    connectivity_path = diff_dir / 'differential_connectivity.tsv'
    df_connectivity = pd.read_csv(connectivity_path, sep='\t')
    logger.info(f"Loaded connectivity for {len(df_connectivity):,} genes")
    
    # Load top hubs
    hubs_path = diff_dir / 'top_hubs.json'
    with open(hubs_path, 'r') as f:
        all_hubs = json.load(f)
    top_hubs = all_hubs[:top_k_visualize]
    logger.info(f"Loaded top {top_k_visualize} hubs")
    logger.info("")
    
    # Step 2: Volcano plot (using optimized dataset)
    logger.info("Step 2: Creating volcano plot (smart sampling)")
    logger.info("-" * 40)
    volcano_path_out = viz_dir / '01_volcano_plot.png'
    plot_volcano(df_volcano, volcano_path_out, fdr_threshold, min_effect_size)
    logger.info("")
    
    # Step 3: Hub ranking
    logger.info("Step 3: Creating hub ranking bar chart")
    logger.info("-" * 40)
    ranking_path = viz_dir / '02_hub_ranking_bar.png'
    plot_hub_ranking(df_connectivity, ranking_path, top_k=20)
    logger.info("")
    
    # Step 4: Connectivity scatter
    logger.info("Step 4: Creating connectivity scatter plot")
    logger.info("-" * 40)
    scatter_path = viz_dir / '03_connectivity_scatter.png'
    plot_connectivity_scatter(df_connectivity, scatter_path)
    logger.info("")
    
    # Step 5: Hub comparison grid (THE IMPACT VISUALIZATION)
    logger.info("Step 5: Creating hub comparison grid")
    logger.info("-" * 40)
    grid_path = viz_dir / '04_hub_comparison_grid.png'
    plot_hub_comparison_grid(top_hubs, subgraph_dir, grid_path)
    logger.info("")
    
    # Summary
    total_time = (datetime.now() - start_time).total_seconds()
    
    logger.info("="*60)
    logger.info("STEP 4 COMPLETE - ALL VISUALIZATIONS GENERATED")
    logger.info("="*60)
    logger.info(f"Total execution time: {total_time:.1f} seconds")
    logger.info("")
    logger.info("Visualizations saved:")
    logger.info(f"  1. {viz_dir / '01_volcano_plot.png'} (optimized smart sampling)")
    logger.info(f"  2. {viz_dir / '02_hub_ranking_bar.png'}")
    logger.info(f"  3. {viz_dir / '03_connectivity_scatter.png'}")
    logger.info(f"  4. {viz_dir / '04_hub_comparison_grid.png'}  KEY FIGURE")
    logger.info("")
    logger.info("Smart Sampling Statistics:")
    logger.info(f"  • Volcano dataset points: {len(df_volcano):,}")
    logger.info(f"  • Significant points in plot: {df_volcano['is_significant'].sum():,}")
    logger.info(f"  • Full analysis significant: {len(df_sig):,}")
    logger.info(f"  • Storage saved: ~{(89559036 - len(df_volcano)) * 120 / 1e6:.0f} MB")
    logger.info("")
    logger.info("PROTOTYPE VALIDATION COMPLETE!")
    logger.info("="*60)
    logger.info("")
    logger.info("Next steps:")
    logger.info("  • Review visualizations in output/visualizations/")
    logger.info("  • Check summary_stats.json for quantitative metrics")
    logger.info("  • Use 04_hub_comparison_grid.png for thesis defense")
    logger.info("  • Proceed to full pipeline if results are satisfactory")
    logger.info("")
    logger.info("="*60)


if __name__ == "__main__":
    main()