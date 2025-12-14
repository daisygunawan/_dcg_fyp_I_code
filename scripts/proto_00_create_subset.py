"""
proto_00_create_subset.py

Purpose: Create small subset of expression matrices for prototype
- Load full preprocessed matrices from ../output/00_b_data_preprocess/matrices/
- Randomly select 25 tumor + 5 normal samples
- Save to prototype/input/preprocessed/matrices/
- Reduces file size from 300+ MB to ~10 MB for faster iteration

Expected Runtime: ~30 seconds
Memory Usage: ~1 GB peak (loading full matrices)
"""

import numpy as np
import pandas as pd
import yaml
import json
import logging
from pathlib import Path
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
        sys.exit(1)


def get_source_matrix_path(config, matrix_filename):
    """
    Get path to source matrix file from config.
    Tries multiple possible locations based on where script is run from.
    
    Args:
        config: Configuration dictionary
        matrix_filename: Name of matrix file (e.g., 'tumor_matrix.tsv')
    
    Returns:
        Path object to matrix file
    """
    base_dir = Path(config['source_data']['matrices_dir'])
    
    # Try different relative paths depending on execution context
    possible_base_paths = [
        base_dir,  # From prototype/ directory
        Path('..') / base_dir,  # From prototype/scripts/ directory
        Path('../..') / base_dir,  # From prototype/scripts/subfolder
    ]
    
    for base_path in possible_base_paths:
        full_path = base_path / matrix_filename
        if full_path.exists():
            logger.info(f"Found source matrix: {full_path.resolve()}")
            return full_path
    
    # If not found, report error with expected location
    expected_path = base_dir / matrix_filename
    logger.error(f"Source matrix not found: {matrix_filename}")
    logger.error(f"Expected location: {expected_path.resolve()}")
    logger.error(f"Configured path: {config['source_data']['matrices_dir']}")
    logger.error("")
    logger.error("Please ensure:")
    logger.error("  1. Full pipeline (00_b_data_preprocess) has been run")
    logger.error("  2. Source path in config_proto.yaml is correct")
    logger.error(f"  3. File exists: {expected_path}")
    sys.exit(1)


def load_full_matrix(file_path):
    """Load full expression matrix from TSV file."""
    try:
        logger.info(f"Loading full matrix: {file_path.name}")
        logger.info(f"  From: {file_path.parent}")
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        logger.info(f"  Shape: {df.shape[0]} genes × {df.shape[1]} samples")
        logger.info(f"  Memory: {df.memory_usage(deep=True).sum() / 1e6:.1f} MB")
        return df
    except FileNotFoundError:
        logger.error(f"Matrix file not found: {file_path}")
        logger.error("Please run the full pipeline (00_b_data_preprocess) first")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading matrix: {e}")
        sys.exit(1)


def create_subset(df, n_samples, condition_name, random_seed=42):
    """
    Randomly select subset of samples from matrix.
    
    Args:
        df: Full expression matrix (genes × samples)
        n_samples: Number of samples to select
        condition_name: 'tumor' or 'normal' for logging
        random_seed: Random seed for reproducibility
    
    Returns:
        subset_df: Sampled matrix
        selected_samples: List of selected sample IDs
    """
    np.random.seed(random_seed)
    
    if n_samples > df.shape[1]:
        logger.warning(f"Requested {n_samples} {condition_name} samples, but only {df.shape[1]} available")
        n_samples = df.shape[1]
    
    # Random sampling
    sample_indices = np.random.choice(df.shape[1], size=n_samples, replace=False)
    subset_df = df.iloc[:, sample_indices]
    selected_samples = subset_df.columns.tolist()
    
    logger.info(f"Selected {n_samples} {condition_name} samples")
    logger.info(f"  Subset size: {subset_df.memory_usage(deep=True).sum() / 1e6:.1f} MB")
    
    return subset_df, selected_samples


def save_subset_matrix(df, output_path):
    """Save subset matrix to TSV file."""
    try:
        df.to_csv(output_path, sep='\t')
        file_size_mb = output_path.stat().st_size / 1e6
        logger.info(f"Saved to: {output_path}")
        logger.info(f"  File size: {file_size_mb:.1f} MB")
    except Exception as e:
        logger.error(f"Error saving matrix: {e}")
        sys.exit(1)


def main():
    """Main execution function."""
    logger.info("="*60)
    logger.info("PROTOTYPE STEP 0: CREATE SUBSET MATRICES")
    logger.info("="*60)
    logger.info("")
    
    start_time = datetime.now()
    
    # Load configuration
    config = load_config()
    
    # Extract parameters
    n_tumor = config['sampling']['n_tumor']
    n_normal = config['sampling']['n_normal']
    random_seed = config['sampling']['random_seed']
    
    # Create output directory
    output_dir = Path('input/preprocessed/matrices')
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_dir.resolve()}")
    logger.info("")
    
    # Step 1: Locate source matrices using config
    logger.info("Step 1: Locating source matrices from config")
    logger.info("-" * 40)
    logger.info(f"Source directory: {config['source_data']['matrices_dir']}")
    logger.info(f"Tumor file: {config['source_data']['tumor_matrix_file']}")
    logger.info(f"Normal file: {config['source_data']['normal_matrix_file']}")
    logger.info("")
    
    # Get paths from config
    tumor_full_path = get_source_matrix_path(config, config['source_data']['tumor_matrix_file'])
    normal_full_path = get_source_matrix_path(config, config['source_data']['normal_matrix_file'])
    logger.info("")
    
    # Step 2: Load full tumor matrix
    logger.info("Step 2: Loading full tumor matrix")
    logger.info("-" * 40)
    tumor_full = load_full_matrix(tumor_full_path)
    logger.info("")
    
    # Step 3: Create tumor subset
    logger.info(f"Step 3: Creating tumor subset ({n_tumor} samples)")
    logger.info("-" * 40)
    tumor_subset, tumor_samples = create_subset(tumor_full, n_tumor, 'tumor', random_seed)
    tumor_output_path = output_dir / 'tumor_matrix.tsv'
    save_subset_matrix(tumor_subset, tumor_output_path)
    
    # Clear memory
    del tumor_full
    logger.info("")
    
    # Step 4: Load full normal matrix
    logger.info("Step 4: Loading full normal matrix")
    logger.info("-" * 40)
    normal_full = load_full_matrix(normal_full_path)
    logger.info("")
    
    # Step 5: Create normal subset
    logger.info(f"Step 5: Creating normal subset ({n_normal} samples)")
    logger.info("-" * 40)
    normal_subset, normal_samples = create_subset(normal_full, n_normal, 'normal', random_seed)
    normal_output_path = output_dir / 'normal_matrix.tsv'
    save_subset_matrix(normal_subset, normal_output_path)
    
    # Clear memory
    del normal_full
    logger.info("")
    
    # Step 6: Verify gene consistency
    logger.info("Step 6: Verifying gene consistency")
    logger.info("-" * 40)
    if not tumor_subset.index.equals(normal_subset.index):
        logger.error("Gene mismatch between tumor and normal subsets!")
        sys.exit(1)
    
    logger.info(f"Gene sets match: {len(tumor_subset.index)} genes")
    logger.info("")
    
    # Step 7: Save metadata
    logger.info("Step 7: Saving subset metadata")
    logger.info("-" * 40)
    
    metadata = {
        'timestamp': datetime.now().isoformat(),
        'random_seed': random_seed,
        'n_tumor_requested': n_tumor,
        'n_normal_requested': n_normal,
        'n_tumor_actual': len(tumor_samples),
        'n_normal_actual': len(normal_samples),
        'n_genes': len(tumor_subset.index),
        'tumor_samples': tumor_samples,
        'normal_samples': normal_samples,
        'source_config': {
            'matrices_dir': config['source_data']['matrices_dir'],
            'tumor_file': config['source_data']['tumor_matrix_file'],
            'normal_file': config['source_data']['normal_matrix_file']
        },
        'source_files_resolved': {
            'tumor': str(tumor_full_path.resolve()),
            'normal': str(normal_full_path.resolve())
        },
        'file_sizes_mb': {
            'tumor_subset': tumor_output_path.stat().st_size / 1e6,
            'normal_subset': normal_output_path.stat().st_size / 1e6
        }
    }
    
    metadata_path = output_dir / 'subset_metadata.json'
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    logger.info(f"Metadata saved to: {metadata_path}")
    logger.info("")
    
    # Summary
    total_time = (datetime.now() - start_time).total_seconds()
    
    logger.info("="*60)
    logger.info("STEP 0 COMPLETE")
    logger.info("="*60)
    logger.info(f"Total execution time: {total_time:.1f} seconds")
    logger.info("")
    logger.info("Subset matrices created:")
    logger.info(f"  • {tumor_output_path} ({metadata['file_sizes_mb']['tumor_subset']:.1f} MB)")
    logger.info(f"  • {normal_output_path} ({metadata['file_sizes_mb']['normal_subset']:.1f} MB)")
    logger.info(f"  • {metadata_path}")
    logger.info("")
    logger.info(f"Space saved: ~{300 - sum(metadata['file_sizes_mb'].values()):.0f} MB")
    logger.info("")
    logger.info("Next step: python scripts/proto_01_prepare_data.py")
    logger.info("="*60)


if __name__ == "__main__":
    main()