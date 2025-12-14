"""
proto_01_prepare_data.py

Purpose: Sample selection and correlation matrix computation
- Load preprocessed expression matrices from 00_b
- Randomly select 100 tumor + 10 normal samples (stratified)
- Compute Spearman correlation matrices (tumor/normal separately)
- Save correlation matrices as compressed NPZ files

Expected Runtime: 1-2 minutes
Memory Usage: ~2.5 GB peak
"""

import numpy as np
import pandas as pd
import yaml
import json
import logging
from pathlib import Path
from scipy.stats import spearmanr
from datetime import datetime
import sys

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('output/prototype_execution.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


def load_config(config_path='config_proto.yaml'):
    """Load configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        logger.info(f"Configuration loaded from {config_path}")
        return config
    except FileNotFoundError:
        logger.error(f"Configuration file not found: {config_path}")
        logger.error("Please ensure config_proto.yaml exists in the prototype directory")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading configuration: {e}")
        sys.exit(1)


def load_expression_matrix(file_path):
    """
    Load expression matrix from TSV file.
    
    Expected format: genes in rows, samples in columns
    First column: gene_key (ENSEMBL_ID|GENE_SYMBOL)
    """
    try:
        logger.info(f"Loading matrix from: {file_path}")
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        logger.info(f"  Shape: {df.shape[0]} genes × {df.shape[1]} samples")
        logger.info(f"  Memory: {df.memory_usage(deep=True).sum() / 1e6:.1f} MB")
        return df
    except FileNotFoundError:
        logger.error(f"Matrix file not found: {file_path}")
        logger.error("Please ensure input data is linked/copied to input/preprocessed/matrices/")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading matrix: {e}")
        sys.exit(1)


def sample_data(tumor_df, normal_df, n_tumor, n_normal, random_seed=42):
    """
    Randomly select subset of samples (stratified by condition).
    
    Args:
        tumor_df: Full tumor expression matrix
        normal_df: Full normal expression matrix
        n_tumor: Number of tumor samples to select
        n_normal: Number of normal samples to select
        random_seed: Random seed for reproducibility
    
    Returns:
        tumor_subset, normal_subset, sample_info dict
    """
    logger.info(f"Selecting random subset: {n_tumor} tumor + {n_normal} normal samples")
    
    np.random.seed(random_seed)
    
    # Check if requested samples exceed available
    if n_tumor > tumor_df.shape[1]:
        logger.warning(f"Requested {n_tumor} tumor samples, but only {tumor_df.shape[1]} available")
        n_tumor = tumor_df.shape[1]
    
    if n_normal > normal_df.shape[1]:
        logger.warning(f"Requested {n_normal} normal samples, but only {normal_df.shape[1]} available")
        n_normal = normal_df.shape[1]
    
    # Random sampling
    tumor_sample_indices = np.random.choice(tumor_df.shape[1], size=n_tumor, replace=False)
    normal_sample_indices = np.random.choice(normal_df.shape[1], size=n_normal, replace=False)
    
    tumor_subset = tumor_df.iloc[:, tumor_sample_indices]
    normal_subset = normal_df.iloc[:, normal_sample_indices]
    
    # Store sample information
    sample_info = {
        'n_tumor': n_tumor,
        'n_normal': n_normal,
        'tumor_samples': tumor_subset.columns.tolist(),
        'normal_samples': normal_subset.columns.tolist(),
        'random_seed': random_seed,
        'timestamp': datetime.now().isoformat()
    }
    
    logger.info(f"Selected {n_tumor} tumor + {n_normal} normal samples")
    return tumor_subset, normal_subset, sample_info


def compute_correlation_matrix(expression_df, method='spearman'):
    """
    Compute correlation matrix using specified method.
    
    Args:
        expression_df: Expression matrix (genes × samples)
        method: Correlation method ('spearman' or 'pearson')
    
    Returns:
        corr_matrix: Correlation matrix (genes × genes)
        gene_names: List of gene identifiers
    """
    logger.info(f"Computing {method} correlation matrix...")
    logger.info(f"  Matrix shape: {expression_df.shape[0]} genes × {expression_df.shape[1]} samples")
    logger.info(f"  Total pairs to compute: {expression_df.shape[0] * (expression_df.shape[0] - 1) // 2:,}")
    
    start_time = datetime.now()
    
    try:
        if method.lower() == 'spearman':
            # Transpose: correlation needs samples in rows, genes in columns
            corr_matrix, _ = spearmanr(expression_df.T, axis=0, nan_policy='omit')
        elif method.lower() == 'pearson':
            corr_matrix = expression_df.T.corr(method='pearson').values
        else:
            raise ValueError(f"Unknown correlation method: {method}")
        
        # Handle any NaN values (replace with 0)
        corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)
        
        gene_names = expression_df.index.tolist()
        
        elapsed = (datetime.now() - start_time).total_seconds()
        logger.info(f"Correlation computed in {elapsed:.1f} seconds")
        logger.info(f"  Output shape: {corr_matrix.shape}")
        logger.info(f"  Correlation range: [{corr_matrix.min():.3f}, {corr_matrix.max():.3f}]")
        logger.info(f"  Mean absolute correlation: {np.abs(corr_matrix[np.triu_indices_from(corr_matrix, k=1)]).mean():.3f}")
        
        return corr_matrix, gene_names
        
    except MemoryError:
        logger.error("Memory error during correlation computation")
        logger.error("Try reducing sample size or number of genes in config_proto.yaml")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error computing correlation: {e}")
        sys.exit(1)


def save_correlation_matrix(corr_matrix, gene_names, output_path):
    """
    Save correlation matrix to compressed NPZ format.
    
    Args:
        corr_matrix: Numpy array (genes × genes)
        gene_names: List of gene identifiers
        output_path: Path to save NPZ file
    """
    try:
        np.savez_compressed(
            output_path,
            matrix=corr_matrix.astype(np.float32),  # Convert to float32 to save space
            genes=np.array(gene_names, dtype=object)
        )
        
        file_size_mb = Path(output_path).stat().st_size / 1e6
        logger.info(f"Saved to: {output_path}")
        logger.info(f"  File size: {file_size_mb:.1f} MB")
        
    except Exception as e:
        logger.error(f"Error saving correlation matrix: {e}")
        sys.exit(1)


def main():
    """Main execution function."""
    logger.info("="*60)
    logger.info("PROTOTYPE STEP 1: DATA PREPARATION & CORRELATION")
    logger.info("="*60)
    logger.info("")
    
    start_time = datetime.now()
    
    # Load configuration
    config = load_config()
    
    # Extract parameters
    input_dir = Path(config['paths']['input_matrices'])
    output_dir = Path(config['paths']['output_base'])
    n_tumor = config['sampling']['n_tumor']
    n_normal = config['sampling']['n_normal']
    random_seed = config['sampling']['random_seed']
    corr_method = config['network']['correlation_method']
    
    # Create output directories
    corr_output_dir = output_dir / 'correlation_matrices'
    corr_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Load expression matrices
    logger.info("Step 1: Loading expression matrices")
    logger.info("-" * 40)
    tumor_df = load_expression_matrix(input_dir / 'tumor_matrix.tsv')
    normal_df = load_expression_matrix(input_dir / 'normal_matrix.tsv')
    logger.info("")
    
    # Verify genes match between tumor and normal
    if not tumor_df.index.equals(normal_df.index):
        logger.error("Gene mismatch between tumor and normal matrices")
        logger.error("Matrices must have identical gene sets")
        sys.exit(1)
    
    logger.info(f"Gene sets match: {len(tumor_df.index)} genes")
    logger.info("")
    
    # Step 2: Sample subset
    logger.info("Step 2: Selecting sample subset")
    logger.info("-" * 40)
    tumor_subset, normal_subset, sample_info = sample_data(
        tumor_df, normal_df, n_tumor, n_normal, random_seed
    )
    
    # Save sample information
    sample_info_path = corr_output_dir / 'sample_info.json'
    with open(sample_info_path, 'w') as f:
        json.dump(sample_info, f, indent=2)
    logger.info(f"Sample info saved to: {sample_info_path}")
    logger.info("")
    
    # Step 3: Compute tumor correlation matrix
    logger.info("Step 3: Computing tumor correlation matrix")
    logger.info("-" * 40)
    tumor_corr, tumor_genes = compute_correlation_matrix(tumor_subset, method=corr_method)
    tumor_corr_path = corr_output_dir / 'tumor_corr_proto.npz'
    save_correlation_matrix(tumor_corr, tumor_genes, tumor_corr_path)
    logger.info("")
    
    # Step 4: Compute normal correlation matrix
    logger.info("Step 4: Computing normal correlation matrix")
    logger.info("-" * 40)
    normal_corr, normal_genes = compute_correlation_matrix(normal_subset, method=corr_method)
    normal_corr_path = corr_output_dir / 'normal_corr_proto.npz'
    save_correlation_matrix(normal_corr, normal_genes, normal_corr_path)
    logger.info("")
    
    # Summary
    total_time = (datetime.now() - start_time).total_seconds()
    logger.info("="*60)
    logger.info("STEP 1 COMPLETE")
    logger.info("="*60)
    logger.info(f"Total execution time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    logger.info("")
    logger.info("Outputs generated:")
    logger.info(f"  • {tumor_corr_path}")
    logger.info(f"  • {normal_corr_path}")
    logger.info(f"  • {sample_info_path}")
    logger.info("")
    logger.info("Next step: python scripts/proto_02_differential.py")
    logger.info("="*60)


if __name__ == "__main__":
    main()