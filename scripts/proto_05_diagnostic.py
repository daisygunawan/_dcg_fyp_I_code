"""
proto_05_diagnostic.py

Purpose: Diagnostic analysis to understand connectivity distributions
- Analyze connectivity patterns across all genes
- Identify different "zones" of rewiring
- Help select appropriate visualization hubs

Expected Runtime: ~5 seconds
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
from pathlib import Path
import sys

def analyze_connectivity_distribution():
    """Analyze connectivity distribution and rewiring patterns."""
    
    # Load connectivity data
    df = pd.read_csv('output/differential_results/differential_connectivity.tsv', sep='\t')
    
    print("="*70)
    print("CONNECTIVITY DISTRIBUTION ANALYSIS")
    print("="*70)
    print()
    
    # Basic statistics
    print("CONNECTIVITY STATISTICS:")
    print("-" * 70)
    print(f"Total genes analyzed: {len(df):,}")
    print()
    
    print("TUMOR CONNECTIVITY:")
    print(f"  Mean: {df['tumor_connectivity'].mean():.1f}")
    print(f"  Median: {df['tumor_connectivity'].median():.1f}")
    print(f"  Range: {df['tumor_connectivity'].min():.0f} - {df['tumor_connectivity'].max():.0f}")
    print()
    
    print("NORMAL CONNECTIVITY:")
    print(f"  Mean: {df['normal_connectivity'].mean():.1f}")
    print(f"  Median: {df['normal_connectivity'].median():.1f}")
    print(f"  Range: {df['normal_connectivity'].min():.0f} - {df['normal_connectivity'].max():.0f}")
    print()
    
    print("DELTA CONNECTIVITY (ABSOLUTE):")
    print(f"  Mean: {df['delta_connectivity'].mean():.1f}")
    print(f"  Median: {df['delta_connectivity'].median():.1f}")
    print(f"  Range: {df['delta_connectivity'].min():.0f} - {df['delta_connectivity'].max():.0f}")
    print()
    
    # Percentile analysis
    print("ABSOLUTE DELTA PERCENTILES:")
    print("-" * 70)
    percentiles = [10, 25, 50, 75, 90, 95, 99]
    for p in percentiles:
        val = df['delta_connectivity'].quantile(p/100)
        print(f"  {p:2d}th percentile: {val:.1f}")
    print()
    
    # Classify genes into rewiring zones
    print("REWIRING ZONES:")
    print("-" * 70)
    
    # Define zones
    extreme_threshold = df['delta_connectivity'].quantile(0.95)  # Top 5%
    high_threshold = df['delta_connectivity'].quantile(0.75)     # Top 25%
    moderate_threshold = df['delta_connectivity'].quantile(0.50) # Top 50%
    
    extreme = df[df['delta_connectivity'] >= extreme_threshold]
    high = df[(df['delta_connectivity'] >= high_threshold) & 
              (df['delta_connectivity'] < extreme_threshold)]
    moderate = df[(df['delta_connectivity'] >= moderate_threshold) & 
                  (df['delta_connectivity'] < high_threshold)]
    low = df[df['delta_connectivity'] < moderate_threshold]
    
    print(f"1. EXTREME rewiring (top 5%, Δ ≥ {extreme_threshold:.0f}): {len(extreme):,} genes")
    print(f"   - These are the most dramatic changes")
    print(f"   - Best for showing biological significance")
    print(f"   - May be visually too extreme for presentations")
    print()
    
    print(f"2. HIGH rewiring (75-95%, Δ = {high_threshold:.0f}-{extreme_threshold:.0f}): {len(high):,} genes")
    print(f"   - Strong rewiring signals")
    print(f"   - Good balance of biological meaning + visual clarity")
    print()
    
    print(f"3. MODERATE rewiring (50-75%, Δ = {moderate_threshold:.0f}-{high_threshold:.0f}): {len(moderate):,} genes")
    print(f"   - Substantial but not extreme changes")
    print(f"   - Clearest visual contrast for demonstrations")
    print(f"   - RECOMMENDED for prototype visualization")
    print()
    
    print(f"4. LOW rewiring (bottom 50%, Δ < {moderate_threshold:.0f}): {len(low):,} genes")
    print(f"   - Subtle changes, may not show visually")
    print()
    
    # Sample genes from each zone
    print("SAMPLE GENES FROM EACH ZONE:")
    print("-" * 70)
    
    for zone_name, zone_df, n_samples in [
        ("EXTREME", extreme, 5),
        ("HIGH", high, 5),
        ("MODERATE", moderate, 5)
    ]:
        print(f"\n{zone_name} REWIRING (sample of {n_samples}):")
        sample = zone_df.head(n_samples)
        for _, row in sample.iterrows():
            symbol = row['gene'].split('|')[1] if '|' in row['gene'] else row['gene']
            print(f"  {symbol:15s} | Tumor: {row['tumor_connectivity']:4.0f} | "
                  f"Normal: {row['normal_connectivity']:4.0f} | "
                  f"Δ: {row['delta_connectivity']:4.0f} | "
                  f"Direction: {row['connectivity_direction']}")
    print()
    
    # Recommendations
    print("="*70)
    print("RECOMMENDATIONS FOR VISUALIZATION:")
    print("="*70)
    print()
    print("Based on your current results:")
    print()
    
    # Check if extreme genes are TOO extreme
    top3 = df.nlargest(3, 'delta_connectivity')
    
    for i, (_, row) in enumerate(top3.iterrows(), 1):
        symbol = row['gene'].split('|')[1] if '|' in row['gene'] else row['gene']
        t_conn = row['tumor_connectivity']
        n_conn = row['normal_connectivity']
        
        print(f"Current Hub #{i}: {symbol}")
        print(f"  Tumor: {t_conn:.0f} edges, Normal: {n_conn:.0f} edges")
        
        # Diagnose issues
        if t_conn < 10 and n_conn > 1000:
            print(f"  WARNING: EXTREME: Near-complete loss of connectivity")
            print(f"     Visual: Tumor will be nearly empty, Normal very dense")
            print(f"     Suggestion: Try 'moderate' mode for better visual balance")
        elif n_conn < 10 and t_conn > 1000:
            print(f"  WARNING: EXTREME: Near-complete gain of connectivity")
            print(f"     Visual: Normal will be nearly empty, Tumor very dense")
            print(f"     Suggestion: Try 'moderate' mode for better visual balance")
        elif abs(t_conn - n_conn) > 500:
            print(f"  WARNING: Very large difference (Δ > 500)")
            print(f"     Visual contrast will be very strong")
        else:
            print(f"  OK: Reasonable for visualization")
        print()
    
    print("SUGGESTED CONFIG CHANGES:")
    print("-" * 70)
    print()
    print("Option 1 (RECOMMENDED): Use 'moderate' mode")
    print("  hubs:")
    print("    selection_mode: 'moderate'  # Change from 'top'")
    print()
    print("Option 2: Use 'median' mode")
    print("  hubs:")
    print("    selection_mode: 'median'  # Select from middle range")
    print()
    print("Option 3: Increase sample size")
    print("  sampling:")
    print("    n_tumor: 100  # Increase from 25")
    print("    n_normal: 20  # Increase from 5")
    print("  (More samples -> more stable correlations -> less extreme patterns)")
    print()
    print("Option 4: Adjust correlation threshold")
    print("  network:")
    print("    correlation_threshold: 0.75  # Try different values (0.6-0.85)")
    print()
    
    # Create diagnostic plot
    create_diagnostic_plots(df, extreme_threshold, high_threshold, moderate_threshold)
    
    print("="*70)
    print("Diagnostic plots saved to: output/visualizations/")
    print("  - 05_connectivity_distribution.png")
    print("  - 06_rewiring_zones.png")
    print("="*70)


def create_diagnostic_plots(df, extreme_thresh, high_thresh, moderate_thresh):
    """Create diagnostic visualizations."""
    
    Path('output/visualizations').mkdir(parents=True, exist_ok=True)
    
    # Plot 1: Connectivity distribution
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Tumor connectivity histogram
    axes[0, 0].hist(df['tumor_connectivity'], bins=50, color='lightcoral', alpha=0.7, edgecolor='black')
    axes[0, 0].axvline(df['tumor_connectivity'].median(), color='red', linestyle='--', linewidth=2, label='Median')
    axes[0, 0].set_xlabel('Tumor Connectivity')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Tumor Connectivity Distribution')
    axes[0, 0].legend()
    axes[0, 0].grid(alpha=0.3)
    
    # Normal connectivity histogram
    axes[0, 1].hist(df['normal_connectivity'], bins=50, color='lightblue', alpha=0.7, edgecolor='black')
    axes[0, 1].axvline(df['normal_connectivity'].median(), color='blue', linestyle='--', linewidth=2, label='Median')
    axes[0, 1].set_xlabel('Normal Connectivity')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Normal Connectivity Distribution')
    axes[0, 1].legend()
    axes[0, 1].grid(alpha=0.3)
    
    # Absolute delta histogram with zones
    axes[1, 0].hist(df['delta_connectivity'], bins=50, color='gray', alpha=0.7, edgecolor='black')
    axes[1, 0].axvline(extreme_thresh, color='red', linestyle='--', linewidth=2, label='Extreme (95%)')
    axes[1, 0].axvline(high_thresh, color='orange', linestyle='--', linewidth=2, label='High (75%)')
    axes[1, 0].axvline(moderate_thresh, color='green', linestyle='--', linewidth=2, label='Moderate (50%)')
    axes[1, 0].set_xlabel('Absolute Δ Connectivity')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Rewiring Magnitude Distribution')
    axes[1, 0].legend()
    axes[1, 0].grid(alpha=0.3)
    
    # Tumor vs Normal scatter with zones
    colors = ['red' if d >= extreme_thresh else 
              'orange' if d >= high_thresh else 
              'yellow' if d >= moderate_thresh else 'lightgray' 
              for d in df['delta_connectivity']]
    axes[1, 1].scatter(df['normal_connectivity'], df['tumor_connectivity'], 
                      c=colors, s=5, alpha=0.5)
    
    max_val = max(df['normal_connectivity'].max(), df['tumor_connectivity'].max())
    axes[1, 1].plot([0, max_val], [0, max_val], 'k--', linewidth=1, alpha=0.5)
    axes[1, 1].set_xlabel('Normal Connectivity')
    axes[1, 1].set_ylabel('Tumor Connectivity')
    axes[1, 1].set_title('Connectivity Comparison (colored by rewiring zone)')
    axes[1, 1].grid(alpha=0.3)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', label='Extreme (top 5%)'),
        Patch(facecolor='orange', label='High (75-95%)'),
        Patch(facecolor='yellow', label='Moderate (50-75%)'),
        Patch(facecolor='lightgray', label='Low (bottom 50%)')
    ]
    axes[1, 1].legend(handles=legend_elements, loc='upper left')
    
    plt.tight_layout()
    plt.savefig('output/visualizations/05_connectivity_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Rewiring zones comparison
    fig, ax = plt.subplots(figsize=(12, 8))
    
    zones_data = df.copy()
    zones_data['zone'] = pd.cut(zones_data['delta_connectivity'], 
                                bins=[0, moderate_thresh, high_thresh, extreme_thresh, float('inf')],
                                labels=['Low', 'Moderate', 'High', 'Extreme'])
    
    # Box plot by zone
    sns.boxplot(data=zones_data, x='zone', y='delta_connectivity', 
               palette=['lightgray', 'yellow', 'orange', 'red'], ax=ax)
    ax.set_xlabel('Rewiring Zone', fontsize=12, fontweight='bold')
    ax.set_ylabel('Absolute Δ Connectivity', fontsize=12, fontweight='bold')
    ax.set_title('Connectivity Changes by Rewiring Zone', fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('output/visualizations/06_rewiring_zones.png', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    print("\nRunning connectivity diagnostic analysis...")
    print("This will help determine the best hub selection strategy.\n")
    
    try:
        analyze_connectivity_distribution()
    except FileNotFoundError:
        print("ERROR: Could not find differential_connectivity.tsv")
        print("Please run proto_02_differential.py first")
        sys.exit(1)