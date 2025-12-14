"""
proto_02_differential.py

Purpose: Differential co-expression analysis and hub identification
WITH EFFICIENT VOLCANO SAMPLING (no full 89.6M DataFrame)

Expected Runtime: 14-15 minutes
Memory Usage: ~4 GB peak
Storage Output: ~4.3 GB total (4.1 GB full sig + 180 MB volcano)
"""

import numpy as np
import pandas as pd
import yaml
import json
import logging
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests
from datetime import datetime
from tqdm import tqdm
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


def load_correlation_matrices(corr_dir):
    """
    Load correlation matrices from NPZ files.
    
    Returns:
        tumor_corr, normal_corr, genes, sample_info
    """
    logger.info("Loading correlation matrices from Step 1...")
    
    try:
        # Load tumor correlation
        tumor_data = np.load(corr_dir / 'tumor_corr_proto.npz', allow_pickle=True)
        tumor_corr = tumor_data['matrix']
        genes = tumor_data['genes']
        
        # Load normal correlation
        normal_data = np.load(corr_dir / 'normal_corr_proto.npz', allow_pickle=True)
        normal_corr = normal_data['matrix']
        
        # Load sample info
        with open(corr_dir / 'sample_info.json', 'r') as f:
            sample_info = json.load(f)
        
        logger.info(f"Tumor correlation: {tumor_corr.shape}")
        logger.info(f"Normal correlation: {normal_corr.shape}")
        logger.info(f"Genes: {len(genes)}")
        logger.info(f"Samples: {sample_info['n_tumor']} tumor + {sample_info['n_normal']} normal")
        
        return tumor_corr, normal_corr, genes, sample_info
        
    except FileNotFoundError as e:
        logger.error(f"Correlation files not found: {e}")
        logger.error("Please run proto_01_prepare_data.py first")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading correlations: {e}")
        sys.exit(1)


def fisher_z_transform(r, clip_value=0.9999):
    """
    Apply Fisher Z-transformation to correlation coefficients.
    
    Args:
        r: Correlation coefficient(s)
        clip_value: Maximum absolute correlation (to avoid inf)
    
    Returns:
        z: Fisher Z-transformed value(s)
    """
    # Clip to avoid infinities
    r_clipped = np.clip(r, -clip_value, clip_value)
    z = 0.5 * np.log((1 + r_clipped) / (1 - r_clipped))
    return z


def compute_fisher_z_pvalue(r1, r2, n1, n2):
    """
    Compute p-value for difference between two correlations using Fisher Z.
    
    Args:
        r1, r2: Correlation coefficients
        n1, n2: Sample sizes
    
    Returns:
        p_value: Two-tailed p-value
    """
    # Fisher Z-transformation
    z1 = fisher_z_transform(r1)
    z2 = fisher_z_transform(r2)
    
    # Standard error of difference
    se = np.sqrt(1/(n1-3) + 1/(n2-3))
    
    # Z-score for difference
    z_diff = (z1 - z2) / se
    
    # Two-tailed p-value
    p_value = 2 * (1 - stats.norm.cdf(np.abs(z_diff)))
    
    return p_value

def downsample_volcano_dataset(df_volcano, target_total=1500000, seed=42):
    """
    Downsample volcano dataset with SIMPLE random sampling.
    Keep all non-significant (gray) + sample significant (red) proportionally.
    
    Args:
        df_volcano: DataFrame with volcano points
        target_total: Target total points (~1.5M)
        seed: Random seed
    """
    logger.info("Downsampling volcano dataset (simple random sampling)...")
    
    df = df_volcano.copy()
    
    # Separate significant and non-significant
    df_sig = df[df['is_significant']].copy()
    df_non = df[~df['is_significant']].copy()
    
    logger.info(f"  Before: {len(df_sig):,} sig + {len(df_non):,} non = {len(df):,} total")
    
    # ==================================================================
    # STRATEGY: Keep ALL gray (500K), sample red to fill target
    # ==================================================================
    
    # Keep all non-significant (these are rare and important for null hypothesis)
    df_non_sample = df_non
    n_non_kept = len(df_non_sample)
    
    # Calculate how many significant points we can fit
    n_sig_target = target_total - n_non_kept
    n_sig_target = max(0, min(n_sig_target, len(df_sig)))  # Bounds check
    
    logger.info(f"  Target: {n_non_kept:,} non-sig (all) + {n_sig_target:,} sig (sampled)")
    
    # Simple random sampling of significant points
    if n_sig_target > 0 and len(df_sig) > n_sig_target:
        df_sig_sample = df_sig.sample(n=n_sig_target, random_state=seed)
        logger.info(f"  Sampled sig: {len(df_sig):,} -> {len(df_sig_sample):,} ({100*len(df_sig_sample)/len(df_sig):.1f}%)")
    else:
        df_sig_sample = df_sig
        logger.info(f"  Kept all sig: {len(df_sig):,}")
    
    # Combine and shuffle
    df_final = pd.concat([df_sig_sample, df_non_sample], ignore_index=True)
    df_final = df_final.sample(frac=1, random_state=seed).reset_index(drop=True)
    
    logger.info(f"  After sampling: {len(df_final):,} total")
    logger.info(f"  Final ratio: {len(df_sig_sample):,} sig ({100*len(df_sig_sample)/len(df_final):.1f}%) + {len(df_non_sample):,} non ({100*len(df_non_sample)/len(df_final):.1f}%)")
    
    return df_final

def compute_differential_coexpression_with_sampling(
    tumor_corr, normal_corr, n_tumor, n_normal, genes,
    fdr_threshold=0.05, min_effect_size=0.1, 
    volcano_max_points=1500000, seed=42
):
    """
    Compute differential co-expression with EFFICIENT VOLCANO SAMPLING.
    
    NEVER builds the full 89.6M DataFrame. Instead:
    1. Collects ALL significant pairs (34.6M)
    2. Samples non-significant pairs during iteration
    
    Args:
        tumor_corr, normal_corr: Correlation matrices
        n_tumor, n_normal: Sample sizes
        genes: Gene identifiers
        fdr_threshold: FDR cutoff for significance
        min_effect_size: Minimum |delta_r| for biological relevance
        volcano_max_points: Target points for volcano plot
        seed: Random seed
    
    Returns:
        df_sig_pairs: DataFrame with ALL significant pairs
        df_volcano: DataFrame with sampled volcano dataset
        rewiring_score: Mean |delta_r| for significant pairs
    """
    logger.info("Computing differential co-expression WITH EFFICIENT SAMPLING...")
    n_genes = tumor_corr.shape[0]
    
    # Get upper triangle indices
    rows, cols = np.triu_indices(n_genes, k=1)
    total_pairs = len(rows)
    
    logger.info(f"  Analyzing {total_pairs:,} unique gene pairs")
    logger.info(f"  FDR threshold: {fdr_threshold}")
    logger.info(f"  Min effect size: {min_effect_size}")
    logger.info(f"  Volcano target points: {volcano_max_points:,}")
    
    # Set random seed
    np.random.seed(seed)
    
    # Estimate sampling rates
    # We expect ~34.6M significant pairs (keep ALL)
    # Need ~500K non-significant for volcano (sampled from ~55M)
    expected_sig = int(total_pairs * 0.386)  # Based on previous runs
    expected_non = total_pairs - expected_sig
    
    # Calculate sampling probability for non-significant
    target_non_sig_points = volcano_max_points - expected_sig
    if target_non_sig_points > 0:
        non_sig_sampling_prob = target_non_sig_points / expected_non
    else:
        non_sig_sampling_prob = 0.01  # Default 1%
    
    logger.info(f"  Expected significant: ~{expected_sig:,}")
    logger.info(f"  Expected non-significant: ~{expected_non:,}")
    logger.info(f"  Non-significant sampling probability: {non_sig_sampling_prob:.4f}")
    
    # Initialize lists for collecting results
    sig_rows = []
    sig_cols = []
    sig_delta_r = []
    sig_p_values = []
    sig_p_fdr = []
    
    volcano_rows = []
    volcano_cols = []
    volcano_delta_r = []
    volcano_p_values = []
    volcano_p_fdr = []
    volcano_is_sig = []
    
    # Process in chunks to manage memory
    chunk_size = 500000
    total_chunks = (total_pairs + chunk_size - 1) // chunk_size
    
    # First pass: compute p-values
    logger.info("  Computing p-values (chunked)...")
    p_values_all = np.zeros(total_pairs)
    
    for chunk_idx in tqdm(range(total_chunks), desc="  Computing p-values"):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, total_pairs)
        
        chunk_rows = rows[start:end]
        chunk_cols = cols[start:end]
        
        r1 = tumor_corr[chunk_rows, chunk_cols]
        r2 = normal_corr[chunk_rows, chunk_cols]
        
        p_values_all[start:end] = compute_fisher_z_pvalue(r1, r2, n_tumor, n_normal)
    
    # Apply FDR correction
    logger.info("  Applying FDR correction...")
    reject, p_fdr_all, _, _ = multipletests(p_values_all, alpha=fdr_threshold, method='fdr_bh')
    
    # Second pass: collect results with sampling
    logger.info("  Collecting results with sampling...")
    
    # Track counters for progress
    sig_count = 0
    non_sig_sampled_count = 0
    
    for idx in tqdm(range(total_pairs), desc="  Sampling volcano data"):
        row = rows[idx]
        col = cols[idx]
        
        delta_r_val = np.abs(tumor_corr[row, col]) - np.abs(normal_corr[row, col])
        p_val = p_values_all[idx]
        p_fdr_val = p_fdr_all[idx]
        
        # Check if significant
        is_significant = reject[idx] and (np.abs(delta_r_val) >= min_effect_size)
        
        if is_significant:
            # ALWAYS keep significant pairs
            sig_rows.append(row)
            sig_cols.append(col)
            sig_delta_r.append(delta_r_val)
            sig_p_values.append(p_val)
            sig_p_fdr.append(p_fdr_val)
            sig_count += 1
            
            # Also add to volcano dataset
            volcano_rows.append(row)
            volcano_cols.append(col)
            volcano_delta_r.append(delta_r_val)
            volcano_p_values.append(p_val)
            volcano_p_fdr.append(p_fdr_val)
            volcano_is_sig.append(True)
            
        else:
            # Sample non-significant pairs for volcano
            if np.random.random() < non_sig_sampling_prob:
                volcano_rows.append(row)
                volcano_cols.append(col)
                volcano_delta_r.append(delta_r_val)
                volcano_p_values.append(p_val)
                volcano_p_fdr.append(p_fdr_val)
                volcano_is_sig.append(False)
                non_sig_sampled_count += 1
    
    # Log sampling statistics
    logger.info(f"  Collected {sig_count:,} significant pairs")
    logger.info(f"  Sampled {non_sig_sampled_count:,} non-significant pairs")
    logger.info(f"  Total volcano points: {sig_count + non_sig_sampled_count:,}")
    logger.info(f"  Actual sampling rate for non-sig: {non_sig_sampled_count/expected_non:.4f}")
    
    # Create significant pairs DataFrame
    logger.info("  Creating significant pairs DataFrame...")
    df_sig_pairs = pd.DataFrame({
        'gene1_idx': sig_rows,
        'gene2_idx': sig_cols,
        'r_tumor': tumor_corr[sig_rows, sig_cols],
        'r_normal': normal_corr[sig_rows, sig_cols],
        'delta_r': sig_delta_r,
        'p_value': sig_p_values,
        'p_fdr': sig_p_fdr
    })
    
    # Sort by effect size
    df_sig_pairs = df_sig_pairs.sort_values('delta_r', key=abs, ascending=False)
    
    # Compute rewiring score
    rewiring_score = np.mean(np.abs(df_sig_pairs['delta_r']))
    logger.info(f"  Rewiring score (mean |Δr|): {rewiring_score:.3f}")
    
    # Create volcano dataset DataFrame
    logger.info("  Creating volcano dataset DataFrame...")
    df_volcano = pd.DataFrame({
        'gene1_idx': volcano_rows,
        'gene2_idx': volcano_cols,
        'r_tumor': tumor_corr[volcano_rows, volcano_cols],
        'r_normal': normal_corr[volcano_rows, volcano_cols],
        'delta_r': volcano_delta_r,
        'p_value': volcano_p_values,
        'p_fdr': volcano_p_fdr,
        'is_significant': volcano_is_sig
    })
    
    # Shuffle volcano dataset
    df_volcano = df_volcano.sample(frac=1, random_state=seed).reset_index(drop=True)

    # Downsample to target sizes
    df_volcano = downsample_volcano_dataset(df_volcano, 
                                            target_total=1500000,
                                            seed=seed)
    
    return df_sig_pairs, df_volcano, rewiring_score


def compute_differential_connectivity(tumor_corr, normal_corr, genes, threshold=0.7):
    """
    Compute differential connectivity for each gene.
    
    Connectivity = number of edges above threshold
    Delta connectivity = ABSOLUTE VALUE of (tumor connectivity - normal connectivity)
    
    Args:
        tumor_corr, normal_corr: Correlation matrices
        genes: Gene identifiers
        threshold: Correlation threshold for edges
    
    Returns:
        DataFrame with genes sorted by absolute delta_connectivity
    """
    logger.info("Computing differential connectivity per gene...")
    logger.info(f"  Using threshold |r| ≥ {threshold}")
    
    # Threshold correlations
    tumor_adj = np.abs(tumor_corr) >= threshold
    normal_adj = np.abs(normal_corr) >= threshold
    
    # Compute connectivity (sum of edges, excluding self)
    tumor_connectivity = np.sum(tumor_adj, axis=1) - 1  # Subtract self-loop
    normal_connectivity = np.sum(normal_adj, axis=1) - 1
    
    # Use ABSOLUTE VALUE for delta connectivity to show magnitude of change
    delta_connectivity = np.abs(tumor_connectivity - normal_connectivity)
    
    # Also store the raw difference for reference
    raw_delta_connectivity = tumor_connectivity - normal_connectivity
    
    df_connectivity = pd.DataFrame({
        'gene': genes,
        'tumor_connectivity': tumor_connectivity,
        'normal_connectivity': normal_connectivity,
        'delta_connectivity': delta_connectivity,  # Absolute value for ranking
        'raw_delta_connectivity': raw_delta_connectivity,  # Raw difference for reference
        'connectivity_direction': np.where(raw_delta_connectivity > 0, 'gain', 'loss')  # Direction of change
    })
    
    # Sort by absolute delta (descending for largest changes)
    df_connectivity = df_connectivity.sort_values('delta_connectivity', ascending=False)
    
    logger.info(f"Connectivity computed for {len(genes)} genes")
    logger.info(f"  Max absolute delta: {df_connectivity['delta_connectivity'].max():.1f}")
    logger.info(f"  Min absolute delta: {df_connectivity['delta_connectivity'].min():.1f}")
    logger.info(f"  Mean absolute delta: {df_connectivity['delta_connectivity'].mean():.1f}")
    
    # Count genes with connectivity gains vs losses
    n_gains = (raw_delta_connectivity > 0).sum()
    n_losses = (raw_delta_connectivity < 0).sum()
    logger.info(f"  Genes with connectivity gain: {n_gains:,}")
    logger.info(f"  Genes with connectivity loss: {n_losses:,}")
    
    return df_connectivity


def identify_top_hubs(df_connectivity, top_k=50, selection_mode='top', config=None):
    """
    Identify top rewired hubs by ABSOLUTE connectivity change.
    
    Args:
        df_connectivity: DataFrame with connectivity values
        top_k: Number of top hubs to select
        selection_mode: 'top' (most extreme), 'moderate' (balanced), 'median' (middle), 'percentile' (custom range)
        config: Full config dict (needed for percentile mode)
    
    Returns:
        List of dicts with hub information
    """
    logger.info(f"Identifying top {top_k} rewired hubs (mode: {selection_mode})...")
    
    if selection_mode == 'moderate':
        # Select hubs with substantial connectivity (not too sparse/dense)
        # and meaningful rewiring (not too extreme)
        
        # Filter: require minimum connectivity in at least one condition
        min_connectivity = 50  # At least 50 edges in one condition
        max_connectivity = 800  # Not super-hubs (< 800 edges)
        
        # Also filter by reasonable delta range
        max_delta = df_connectivity['delta_connectivity'].quantile(0.90)  # Top 10% max
        
        filtered = df_connectivity[
            (
                ((df_connectivity['tumor_connectivity'] >= min_connectivity) & 
                 (df_connectivity['tumor_connectivity'] <= max_connectivity)) |
                ((df_connectivity['normal_connectivity'] >= min_connectivity) & 
                 (df_connectivity['normal_connectivity'] <= max_connectivity))
            ) &
            (df_connectivity['delta_connectivity'] <= max_delta)
        ].copy()
        
        # Sort by absolute delta
        filtered = filtered.sort_values('delta_connectivity', ascending=False)
        
        # Take top K from filtered set
        top_hubs_df = filtered.head(top_k)
        
        logger.info(f"  Filtered to {len(filtered)} moderate hubs (50-800 edges, delta < {max_delta:.0f})")
        logger.info(f"  Selected top {len(top_hubs_df)} from filtered set")
        
    elif selection_mode == 'median':
        # Select from middle range of delta values
        sorted_df = df_connectivity.sort_values('delta_connectivity', ascending=False)
        
        # Skip top 20% (too extreme) and take from 20%-40% range
        start_idx = int(len(sorted_df) * 0.20)
        end_idx = int(len(sorted_df) * 0.40)
        
        middle_range = sorted_df.iloc[start_idx:end_idx]
        top_hubs_df = middle_range.head(top_k)
        
        logger.info(f"  Selected from percentile range: 20%-40% (median zone)")
        logger.info(f"  Delta range: {top_hubs_df['delta_connectivity'].min():.0f} - {top_hubs_df['delta_connectivity'].max():.0f}")
        
    elif selection_mode == 'visual':
        # Select hubs optimized for visualization clarity
        # Want: decent tumor connectivity + large delta
        
        min_tumor_conn = 80   # Tumor must have visible structure
        max_tumor_conn = 500  # But not too dense
        min_normal_conn = 300 # Normal should be reasonably dense
        min_delta = df_connectivity['delta_connectivity'].quantile(0.50)  # At least median delta
        
        # Filter by these criteria
        filtered = df_connectivity[
            (df_connectivity['tumor_connectivity'] >= min_tumor_conn) &
            (df_connectivity['tumor_connectivity'] <= max_tumor_conn) &
            (df_connectivity['normal_connectivity'] >= min_normal_conn) &
            (df_connectivity['delta_connectivity'] >= min_delta)
        ].copy()
        
        if len(filtered) == 0:
            logger.warning("  No hubs met 'visual' criteria, relaxing constraints...")
            # Fallback: just require minimum tumor connectivity
            filtered = df_connectivity[
                df_connectivity['tumor_connectivity'] >= 50
            ].copy()
        
        # Rank by a composite score: tumor visibility + delta magnitude
        # This balances "can see tumor network" with "shows contrast"
        filtered['visual_score'] = (
            filtered['tumor_connectivity'] * 0.3 +  # 30% weight on tumor visibility
            filtered['delta_connectivity'] * 0.7     # 70% weight on contrast
        )
        
        filtered = filtered.sort_values('visual_score', ascending=False)
        top_hubs_df = filtered.head(top_k)
        
        logger.info(f"  Visual mode: filtered to {len(filtered)} hubs")
        logger.info(f"  Tumor connectivity range: {min_tumor_conn}-{max_tumor_conn}")
        logger.info(f"  Normal connectivity min: {min_normal_conn}")
        logger.info(f"  Selected {len(top_hubs_df)} hubs by visual score")
        if len(top_hubs_df) > 0:
            logger.info(f"  Tumor edges range: {top_hubs_df['tumor_connectivity'].min():.0f}-{top_hubs_df['tumor_connectivity'].max():.0f}")
            logger.info(f"  Normal edges range: {top_hubs_df['normal_connectivity'].min():.0f}-{top_hubs_df['normal_connectivity'].max():.0f}")
            logger.info(f"  Delta range: {top_hubs_df['delta_connectivity'].min():.0f}-{top_hubs_df['delta_connectivity'].max():.0f}")
        
    elif selection_mode == 'percentile':
        # Select from custom percentile range
        if config is None:
            logger.error("Config required for 'percentile' mode")
            raise ValueError("Config parameter required for percentile selection mode")
        
        # Get percentile range from config (e.g., [75, 90])
        p_range = config.get('hubs', {}).get('percentile_range', [75, 90])
        p_low, p_high = p_range
        
        # Calculate percentile bounds
        lower_bound = df_connectivity['delta_connectivity'].quantile(p_low / 100)
        upper_bound = df_connectivity['delta_connectivity'].quantile(p_high / 100)
        
        # Filter to this range
        filtered = df_connectivity[
            (df_connectivity['delta_connectivity'] >= lower_bound) &
            (df_connectivity['delta_connectivity'] <= upper_bound)
        ].copy()
        
        # Sort by delta and take top K from this range
        filtered = filtered.sort_values('delta_connectivity', ascending=False)
        top_hubs_df = filtered.head(top_k)
        
        logger.info(f"  Selected from percentile range: {p_low}%-{p_high}%")
        logger.info(f"  Delta bounds: {lower_bound:.0f} - {upper_bound:.0f}")
        logger.info(f"  Selected {len(top_hubs_df)} hubs from {len(filtered)} in range")
        logger.info(f"  Actual delta range: {top_hubs_df['delta_connectivity'].min():.0f} - {top_hubs_df['delta_connectivity'].max():.0f}")
        
    else:  # 'top' (original behavior)
        # Select top K by ABSOLUTE delta_connectivity (magnitude of change)
        top_hubs_df = df_connectivity.nlargest(top_k, 'delta_connectivity')
        logger.info(f"  Using top {top_k} by absolute delta (most extreme)")
    
    # Convert to list of dicts
    top_hubs = []
    for _, row in top_hubs_df.iterrows():
        hub_info = {
            'gene': row['gene'],
            'gene_symbol': row['gene'].split('|')[1] if '|' in row['gene'] else row['gene'],
            'tumor_connectivity': int(row['tumor_connectivity']),
            'normal_connectivity': int(row['normal_connectivity']),
            'delta_connectivity': float(row['delta_connectivity']),  # Absolute value
            'raw_delta_connectivity': float(row['raw_delta_connectivity']),  # Raw difference
            'connectivity_direction': row['connectivity_direction'],  # 'gain' or 'loss'
            'rank': len(top_hubs) + 1
        }
        top_hubs.append(hub_info)
    
    logger.info(f"Top {len(top_hubs)} hubs identified by {selection_mode} selection")
    logger.info(f"  Top hub: {top_hubs[0]['gene_symbol']} (|Δ| = {top_hubs[0]['delta_connectivity']:.1f}, direction: {top_hubs[0]['connectivity_direction']})")
    
    # Log breakdown of gains vs losses in top hubs
    n_gains = sum(1 for hub in top_hubs if hub['connectivity_direction'] == 'gain')
    n_losses = sum(1 for hub in top_hubs if hub['connectivity_direction'] == 'loss')
    logger.info(f"  Top hubs breakdown: {n_gains} gains, {n_losses} losses")
    
    # Log connectivity range
    conn_tumor_range = (min(h['tumor_connectivity'] for h in top_hubs), 
                        max(h['tumor_connectivity'] for h in top_hubs))
    conn_normal_range = (min(h['normal_connectivity'] for h in top_hubs), 
                         max(h['normal_connectivity'] for h in top_hubs))
    logger.info(f"  Tumor connectivity range: {conn_tumor_range[0]}-{conn_tumor_range[1]}")
    logger.info(f"  Normal connectivity range: {conn_normal_range[0]}-{conn_normal_range[1]}")
    
    return top_hubs


def main():
    """Main execution function."""
    logger.info("")
    logger.info("="*60)
    logger.info("PROTOTYPE STEP 2: DIFFERENTIAL ANALYSIS")
    logger.info("="*60)
    logger.info("")
    
    start_time = datetime.now()
    
    # Load configuration
    config = load_config()
    
    # Extract parameters
    output_dir = Path(config['paths']['output_base'])
    corr_dir = output_dir / 'correlation_matrices'
    diff_dir = output_dir / 'differential_results'
    diff_dir.mkdir(parents=True, exist_ok=True)
    
    fdr_threshold = config['differential']['fdr_threshold']
    min_effect_size = config['differential']['min_effect_size']
    corr_threshold = config['network']['correlation_threshold']
    top_k = config['hubs']['top_k']
    
    # Get volcano target points (new parameter)
    volcano_max_points = config.get('differential', {}).get('volcano_max_points', 1500000)
    
    # Get selection mode from config (default to 'top' if not specified)
    selection_mode = config.get('hubs', {}).get('selection_mode', 'top')
    
    # Step 1: Load correlation matrices
    logger.info("Step 1: Loading correlation matrices")
    logger.info("-" * 40)
    tumor_corr, normal_corr, genes, sample_info = load_correlation_matrices(corr_dir)
    logger.info("")
    
    # Step 2: Compute differential co-expression WITH EFFICIENT SAMPLING
    logger.info("Step 2: Differential co-expression analysis (with efficient sampling)")
    logger.info("-" * 40)
    df_sig_pairs, df_volcano, rewiring_score = compute_differential_coexpression_with_sampling(
        tumor_corr, normal_corr,
        sample_info['n_tumor'], sample_info['n_normal'],
        genes,
        fdr_threshold, min_effect_size,
        volcano_max_points, seed=42
    )
    
    # Step 3: Map gene indices to gene names
    logger.info("  Mapping gene indices to names...")
    
    # For significant pairs
    df_sig_pairs['gene1'] = df_sig_pairs['gene1_idx'].map(lambda i: genes[i])
    df_sig_pairs['gene2'] = df_sig_pairs['gene2_idx'].map(lambda i: genes[i])
    df_sig_pairs = df_sig_pairs.drop(columns=['gene1_idx', 'gene2_idx'])
    
    # For volcano dataset
    df_volcano['gene1'] = df_volcano['gene1_idx'].map(lambda i: genes[i])
    df_volcano['gene2'] = df_volcano['gene2_idx'].map(lambda i: genes[i])
    df_volcano = df_volcano.drop(columns=['gene1_idx', 'gene2_idx'])
    
    # Save volcano dataset (for visualization)
    volcano_path = diff_dir / 'differential_pairs_volcano.tsv'
    df_volcano.to_csv(volcano_path, sep='\t', index=False)
    logger.info(f"Saved volcano dataset to: {volcano_path}")
    logger.info(f"  Size: {len(df_volcano):,} points ({len(df_volcano[df_volcano['is_significant']]):,} sig + {len(df_volcano[~df_volcano['is_significant']]):,} non-sig)")
    logger.info(f"  Memory: {df_volcano.memory_usage(deep=True).sum()/1e6:.1f} MB")
    
    # Save ALL significant pairs (for downstream analysis)
    sig_path = diff_dir / 'differential_pairs_significant.tsv'
    df_sig_pairs.to_csv(sig_path, sep='\t', index=False)
    logger.info(f"Saved all significant pairs to: {sig_path}")
    logger.info(f"  Size: {len(df_sig_pairs):,} pairs ({df_sig_pairs.memory_usage(deep=True).sum()/1e6:.1f} MB)")
    logger.info("")
    
    # Step 4: Compute differential connectivity
    logger.info("Step 4: Differential connectivity analysis")
    logger.info("-" * 40)
    df_connectivity = compute_differential_connectivity(
        tumor_corr, normal_corr, genes, corr_threshold
    )
    
    # Save connectivity results
    connectivity_path = diff_dir / 'differential_connectivity.tsv'
    df_connectivity.to_csv(connectivity_path, sep='\t', index=False)
    logger.info(f"Saved connectivity to: {connectivity_path}")
    logger.info("")
    
    # Step 5: Identify top hubs with specified selection mode
    logger.info("Step 5: Hub identification")
    logger.info("-" * 40)
    top_hubs = identify_top_hubs(df_connectivity, top_k, selection_mode, config)
    
    # Save top hubs
    hubs_path = diff_dir / 'top_hubs.json'
    with open(hubs_path, 'w') as f:
        json.dump(top_hubs, f, indent=2)
    logger.info(f"Saved top hubs to: {hubs_path}")
    logger.info("")
    
    # Step 6: Create summary statistics
    logger.info("Step 6: Generating summary statistics")
    logger.info("-" * 40)
    
    # Calculate additional statistics
    n_connectivity_gains = (df_connectivity['raw_delta_connectivity'] > 0).sum()
    n_connectivity_losses = (df_connectivity['raw_delta_connectivity'] < 0).sum()
    
    # Calculate actual sampling statistics
    n_volcano_sig = df_volcano['is_significant'].sum()
    n_volcano_non = len(df_volcano) - n_volcano_sig
    n_total_non = (len(genes) * (len(genes) - 1) // 2) - len(df_sig_pairs)  # Estimated
    
    sampling_rate_non_sig = n_volcano_non / n_total_non if n_total_non > 0 else 0
    
    summary_stats = {
        'timestamp': datetime.now().isoformat(),
        'sample_counts': {
            'tumor': sample_info['n_tumor'],
            'normal': sample_info['n_normal']
        },
        'genes_analyzed': len(genes),
        'total_pairs_tested': len(genes) * (len(genes) - 1) // 2,
        'significant_pairs': len(df_sig_pairs),
        'rewiring_score': float(rewiring_score),
        'hub_selection_mode': selection_mode,
        'top_hub': top_hubs[0]['gene_symbol'],
        'top_hub_absolute_delta': float(top_hubs[0]['delta_connectivity']),
        'top_hub_direction': top_hubs[0]['connectivity_direction'],
        'volcano_sampling': {
            'target_points': volcano_max_points,
            'actual_points': len(df_volcano),
            'significant_points': int(n_volcano_sig),
            'non_significant_points': int(n_volcano_non),
            'sampling_rate_non_sig': float(sampling_rate_non_sig),
            'sampling_rate_sig': 1.0  # All significant included
        },
        'connectivity_changes': {
            'genes_with_gain': int(n_connectivity_gains),
            'genes_with_loss': int(n_connectivity_losses),
            'mean_absolute_delta': float(df_connectivity['delta_connectivity'].mean())
        },
        'parameters': {
            'fdr_threshold': fdr_threshold,
            'min_effect_size': min_effect_size,
            'correlation_threshold': corr_threshold,
            'volcano_max_points': volcano_max_points
        }
    }
    
    # Save summary
    summary_path = diff_dir / 'summary_stats.json'
    with open(summary_path, 'w') as f:
        json.dump(summary_stats, f, indent=2)
    logger.info(f"Saved summary to: {summary_path}")
    
    # Print key metrics
    logger.info("")
    logger.info("KEY METRICS:")
    logger.info(f"  • Significant pairs: {summary_stats['significant_pairs']:,}")
    logger.info(f"  • Rewiring score: {summary_stats['rewiring_score']:.3f}")
    logger.info(f"  • Hub selection mode: {selection_mode}")
    logger.info(f"  • Top hub: {summary_stats['top_hub']} (|Δ| = {summary_stats['top_hub_absolute_delta']:.1f}, {summary_stats['top_hub_direction']})")
    logger.info(f"  • Volcano dataset: {summary_stats['volcano_sampling']['actual_points']:,} points")
    logger.info(f"  • Sampling rate (non-sig): {summary_stats['volcano_sampling']['sampling_rate_non_sig']:.3%}")
    logger.info(f"  • Connectivity gains: {summary_stats['connectivity_changes']['genes_with_gain']:,} genes")
    logger.info(f"  • Connectivity losses: {summary_stats['connectivity_changes']['genes_with_loss']:,} genes")
    logger.info("")
    
    # Summary
    total_time = (datetime.now() - start_time).total_seconds()
    summary_stats['processing_time_seconds'] = int(total_time)
    
    # Re-save summary with timing
    with open(summary_path, 'w') as f:
        json.dump(summary_stats, f, indent=2)
    
    logger.info("="*60)
    logger.info("STEP 2 COMPLETE")
    logger.info("="*60)
    logger.info(f"Total execution time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    logger.info("")
    logger.info("Outputs generated:")
    logger.info(f"  • {sig_path} (all significant pairs, {len(df_sig_pairs):,} rows)")
    logger.info(f"  • {volcano_path} (optimized volcano dataset, {len(df_volcano):,} rows)")
    logger.info(f"  • {connectivity_path} (connectivity results)")
    logger.info(f"  • {hubs_path} (top hubs)")
    logger.info(f"  • {summary_path} (summary statistics)")
    logger.info("")
    logger.info("Storage estimate:")
    logger.info(f"  • Significant pairs: ~{len(df_sig_pairs) * 120 / 1e6:.1f} MB")
    logger.info(f"  • Volcano dataset: ~{len(df_volcano) * 120 / 1e6:.1f} MB")
    logger.info(f"  • Total Step 2 output: ~{(len(df_sig_pairs) + len(df_volcano)) * 120 / 1e6:.1f} MB")
    logger.info("")
    logger.info("Next step: python scripts/proto_03_build_subgraphs.py")
    logger.info("="*60)


if __name__ == "__main__":
    main()