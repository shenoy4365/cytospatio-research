#!/usr/bin/env python3
"""
Generate 85th and 95th percentile input files for CytoSpatio.

This script creates new percentile datasets by:
1. Taking all original cells (types 0-4) from baseline
2. Sampling appropriate number of type 5 cells from existing percentile files
3. The number of type 5 cells follows the linear pattern observed in existing files

Methodology:
- Existing pattern: type5_count = -875.13 * percentile + 87537.10
- 85th percentile: ~13,151 type 5 cells → 57,051 total
- 95th percentile: ~4,400 type 5 cells → 48,300 total
"""

import pandas as pd
import numpy as np

def generate_percentile_file(percentile, output_path):
    """
    Generate a percentile input file.

    Parameters:
    -----------
    percentile : int
        Target percentile (e.g., 85 or 95)
    output_path : str
        Path to save the output CSV file
    """
    print(f"\nGenerating {percentile}th percentile file...")
    print("="*70)

    # Load baseline data (all original cells, types 0-4)
    baseline_path = "example/cell_data.csv"
    baseline = pd.read_csv(baseline_path)
    print(f"Loaded baseline: {len(baseline):,} cells (types 0-4)")

    # Calculate target number of type 5 cells using linear model
    # Based on pattern: type5_count = -875.13 * percentile + 87537.10
    type5_target = int(round(-875.13 * percentile + 87537.10))
    print(f"Target type 5 cells: {type5_target:,}")

    # Load type 5 cells from a source percentile file
    # Use 80th percentile as source (has 17,546 type 5 cells)
    # If we need more, we'll use 70th (26,269 type 5 cells)
    if type5_target > 17546:
        source_pct = 70
        source_path = f"example/cell_data_percentile_{source_pct}.csv"
    else:
        source_pct = 80
        source_path = f"example/cell_data_percentile_{source_pct}.csv"

    source_df = pd.read_csv(source_path)
    type5_cells = source_df[source_df['marks'] == 5].copy()
    print(f"Loaded {len(type5_cells):,} type 5 cells from {source_pct}th percentile")

    # Check if we have enough type 5 cells
    if len(type5_cells) < type5_target:
        print(f"WARNING: Not enough type 5 cells in source!")
        print(f"  Needed: {type5_target:,}")
        print(f"  Available: {len(type5_cells):,}")
        print(f"  Using all available type 5 cells")
        sampled_type5 = type5_cells
    else:
        # Randomly sample the target number of type 5 cells
        # Set seed for reproducibility
        np.random.seed(42 + percentile)
        sampled_indices = np.random.choice(len(type5_cells), size=type5_target, replace=False)
        sampled_type5 = type5_cells.iloc[sampled_indices].copy()
        print(f"Sampled {len(sampled_type5):,} type 5 cells")

    # Combine baseline cells with sampled type 5 cells
    result_df = pd.concat([baseline, sampled_type5], ignore_index=True)

    # Sort by coordinates for consistency
    result_df = result_df.sort_values(['x', 'y']).reset_index(drop=True)

    total_cells = len(result_df)
    type5_actual = len(sampled_type5)
    type5_percent = (type5_actual / total_cells) * 100

    print(f"\nResult:")
    print(f"  Total cells: {total_cells:,}")
    print(f"  Type 0-4: {len(baseline):,}")
    print(f"  Type 5: {type5_actual:,} ({type5_percent:.1f}%)")

    # Verify cell type distribution
    print(f"\nCell type distribution:")
    for cell_type in sorted(result_df['marks'].unique()):
        count = len(result_df[result_df['marks'] == cell_type])
        percent = (count / total_cells) * 100
        print(f"  Type {cell_type}: {count:,} ({percent:.1f}%)")

    # Save to CSV
    result_df.to_csv(output_path, index=False)
    print(f"\nSaved to: {output_path}")
    print(f"File size: {len(result_df):,} rows")

    return result_df

def verify_pattern():
    """
    Verify the linear pattern across all existing percentile files.
    """
    print("\nVerifying existing percentile pattern...")
    print("="*70)

    percentiles = [50, 60, 70, 80, 90]
    type5_counts = []

    for p in percentiles:
        df = pd.read_csv(f'example/cell_data_percentile_{p}.csv')
        type5_count = len(df[df['marks'] == 5])
        type5_counts.append(type5_count)
        print(f"{p}th percentile: {type5_count:,} type 5 cells")

    # Fit linear model
    coeffs = np.polyfit(percentiles, type5_counts, 1)
    print(f"\nLinear model: type5_count = {coeffs[0]:.2f} * percentile + {coeffs[1]:.2f}")

    # Show predictions for new percentiles
    for p in [85, 95]:
        predicted = int(round(coeffs[0] * p + coeffs[1]))
        print(f"  {p}th percentile: {predicted:,} type 5 cells (predicted)")

    return coeffs

def main():
    """
    Main function to generate both 85th and 95th percentile files.
    """
    print("="*70)
    print("GENERATING NEW PERCENTILE INPUT FILES")
    print("="*70)

    # Verify pattern
    verify_pattern()

    # Generate 85th percentile file
    file_85 = generate_percentile_file(85, "example/cell_data_percentile_85.csv")

    # Generate 95th percentile file
    file_95 = generate_percentile_file(95, "example/cell_data_percentile_95.csv")

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("\nGenerated files:")
    print("  1. example/cell_data_percentile_85.csv")
    print("  2. example/cell_data_percentile_95.csv")
    print("\nThese files follow the same pattern as the existing percentile files:")
    print("  - All 43,900 original cells (types 0-4)")
    print("  - Variable number of simulated type 5 cells")
    print("  - 85th: ~57,051 total cells (~13,151 type 5)")
    print("  - 95th: ~48,300 total cells (~4,400 type 5)")
    print("\nYour professor can now run CytoSpatio on these new percentile files")
    print("using the same parameters as the other percentile models.")
    print("="*70 + "\n")

if __name__ == "__main__":
    main()
