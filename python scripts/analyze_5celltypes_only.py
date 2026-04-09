#!/usr/bin/env python3
"""
Analysis script to calculate RMSE/MAE for ONLY the 5 original cell types (0-4),
excluding the 6th simulated cell type (5).

This answers the professor's question:
"Which model most accurately fits the 5 cell types? Need to calculate the RMSE/MAE
for just the residuals from the points of the 5 cell types, ignoring the added 6th type."

Key insight:
- Each model has 6 residuals per quadrature point (one for each cell type)
- Residuals are ordered: [type0, type1, type2, type3, type4, type5, type0, type1, ...]
- We extract only residuals for types 0-4 and calculate metrics
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Define paths
# set your own base directory here
base_dir = Path("")
outcomes_dir = base_dir / "cytospatio outcomes"
output_baseline_dir = base_dir / "output_baseline"
plots_dir = base_dir / "plots"
plots_dir.mkdir(exist_ok=True)

# Percentile datasets to analyze (plus baseline)
percentiles = [50, 60, 70, 80, 85, 90, 95]

def load_baseline_residuals():
    """
    Load residuals from the baseline model (5 cell types only).

    Returns:
    --------
    numpy array : All residuals from baseline model (already 5 types only)
    """
    resid_path = output_baseline_dir / "cell_data_resid_TR_500_IR_100_HR_1.csv"

    if not resid_path.exists():
        print(f"Warning: {resid_path} not found")
        return None

    # Load residuals
    df = pd.read_csv(resid_path)

    # Extract residual values
    if 'x' in df.columns:
        residuals = df['x'].values
    elif len(df.columns) == 2:
        residuals = df.iloc[:, 1].values
    else:
        residuals = df.values.flatten()

    return residuals

def load_residuals_filtered(percentile):
    """
    Load residuals for ONLY the 43,900 original cells and ONLY types 0-4.

    Parameters:
    -----------
    percentile : int or None
        If None: load baseline model (5 types only)
        If int: Percentile model (50, 60, 70, 80, or 90) - filter to 5 types

    Returns:
    --------
    numpy array : Residuals for 43,900 original cells, types 0-4 only (219,500 residuals)

    Methodology:
    ------------
    The residuals file is a linearized (43,900 + n) * 6 matrix where n is number of type 5 cells.
    We extract ONLY the 43,900 * 5 submatrix:
    - Only the first 43,900 cells (original cells, not type 5 cells)
    - Only their residuals for types 0-4 (not type 5)
    - Result: 43,900 * 5 = 219,500 residuals for ALL models (same as baseline!)
    """
    if percentile is None:
        # Baseline model - already has only 5 types
        return load_baseline_residuals()

    # Percentile model - needs filtering
    resid_path = outcomes_dir / f"outputpercentile{percentile}" / \
                 f"cell_data_percentile_{percentile}_resid_TR_500_IR_100_HR_1.csv"

    if not resid_path.exists():
        print(f"Warning: {resid_path} not found")
        return None

    # Load residuals
    df = pd.read_csv(resid_path)

    # Extract residual values (handle different CSV formats)
    if 'x' in df.columns:
        residuals_all = df['x'].values
    elif len(df.columns) == 2:
        residuals_all = df.iloc[:, 1].values
    else:
        residuals_all = df.values.flatten()

    # CRITICAL FIX: Only look at residuals for the 43,900 ORIGINAL cells
    # The percentile files have structure:
    #   - First 43,900 cells: original cells (types 0-4)
    #   - Remaining cells: type 5 cells
    # We only want the 43,900 * 5 = 219,500 submatrix (same as baseline!)

    n_original_cells = 43900
    n_residuals_per_cell = 6

    # Step 1: Extract only residuals for the first 43,900 cells
    n_residuals_original = n_original_cells * n_residuals_per_cell
    residuals_original_cells = residuals_all[:n_residuals_original]

    # Step 2: Filter to keep only types 0-4 (skip every 6th residual which is type 5)
    indices_to_keep = []
    for i in range(n_original_cells):
        base_idx = i * 6
        indices_to_keep.extend([base_idx + j for j in range(5)])  # Types 0-4 only

    residuals_5types = residuals_original_cells[indices_to_keep]

    return residuals_5types

def get_cell_counts():
    """
    Get the number of cells in each dataset.

    Returns:
    --------
    dict : Dictionary mapping percentile/baseline to cell count
    """
    counts = {}

    # Baseline dataset
    baseline_path = base_dir / "example" / "cell_data.csv"
    if baseline_path.exists():
        df = pd.read_csv(baseline_path)
        counts['baseline'] = len(df)

    # Percentile datasets
    for p in percentiles:
        pct_path = base_dir / "example" / f"cell_data_percentile_{p}.csv"
        if pct_path.exists():
            df = pd.read_csv(pct_path)
            counts[p] = len(df)
    return counts

def analyze_5celltypes_only():
    """
    Calculate RMSE/MAE for only the 5 original cell types (0-4).

    This provides a fair comparison across models by excluding the simulated
    cell type 5 which was artificially added.

    Returns:
    --------
    DataFrame with RMSE/MAE metrics for each model
    """

    print("="*80)
    print("ANALYSIS: RMSE/MAE FOR 5 ORIGINAL CELL TYPES ONLY")
    print("="*80)
    print("\nExcluding residuals from simulated cell type 5")
    print("Calculating metrics for cell types 0-4 only\n")

    cell_counts = get_cell_counts()
    results = []

    # First, analyze the BASELINE model (5 cell types only)
    print(f"\nAnalyzing BASELINE model (5 cell types: 0-4)...")
    residuals_baseline = load_residuals_filtered(None)

    if residuals_baseline is not None:
        n_cells = cell_counts.get('baseline', 0)
        n_residuals = len(residuals_baseline)

        print(f"  Number of cells: {n_cells:,}")
        print(f"  Number of residuals: {n_residuals:,}")
        print(f"  Residuals per cell: {n_residuals / n_cells:.1f}")

        # Calculate metrics
        mse = np.mean(residuals_baseline ** 2)
        rmse = np.sqrt(mse)
        mae = np.mean(np.abs(residuals_baseline))
        median_ae = np.median(np.abs(residuals_baseline))
        q95_ae = np.percentile(np.abs(residuals_baseline), 95)

        result = {
            'model': 'baseline',
            'n_cells': n_cells,
            'n_residuals_5types': n_residuals,
            'MSE_5types': mse,
            'RMSE_5types': rmse,
            'MAE_5types': mae,
            'Median_AE_5types': median_ae,
            'Q95_AE_5types': q95_ae
        }
        results.append(result)

        print(f"  RMSE (5 types): {rmse:.6f}")
        print(f"  MAE (5 types):  {mae:.6f}")
        print(f"  Median AE:      {median_ae:.6f}")

    # Now analyze percentile models
    for p in percentiles:
        print(f"\nAnalyzing {p}th percentile model...")

        # Load and filter residuals
        residuals_5types = load_residuals_filtered(p)

        if residuals_5types is None:
            continue

        # CRITICAL: We're analyzing the same 43,900 original cells across all models
        n_cells_analyzed = 43900
        n_residuals_analyzed = len(residuals_5types)
        n_expected = 43900 * 5  # Should be 219,500 for all models

        n_cells_total = cell_counts.get(p, 0)  # Total in dataset (including type 5)

        print(f"  Dataset has {n_cells_total:,} total cells (43,900 original + {n_cells_total - 43900:,} type 5)")
        print(f"  Analyzing: {n_cells_analyzed:,} original cells only")
        print(f"  Residuals analyzed: {n_residuals_analyzed:,} (43,900 cells * 5 types)")
        print(f"  Expected: {n_expected:,}")
        if n_residuals_analyzed == n_expected:
            print(f"  ✓ Correct! Same as baseline.")

        # Calculate metrics for 5 cell types only
        mse = np.mean(residuals_5types ** 2)
        rmse = np.sqrt(mse)
        mae = np.mean(np.abs(residuals_5types))
        median_ae = np.median(np.abs(residuals_5types))
        q95_ae = np.percentile(np.abs(residuals_5types), 95)

        result = {
            'model': f'{p}th',
            'n_cells': n_cells_analyzed,  # Always 43,900
            'n_residuals_5types': n_residuals_analyzed,  # Always 219,500
            'MSE_5types': mse,
            'RMSE_5types': rmse,
            'MAE_5types': mae,
            'Median_AE_5types': median_ae,
            'Q95_AE_5types': q95_ae
        }
        results.append(result)

        print(f"  RMSE (5 types): {rmse:.6f}")
        print(f"  MAE (5 types):  {mae:.6f}")
        print(f"  Median AE:      {median_ae:.6f}")

    results_df = pd.DataFrame(results)

    # Rank models by RMSE and MAE
    results_df['RMSE_rank'] = results_df['RMSE_5types'].rank()
    results_df['MAE_rank'] = results_df['MAE_5types'].rank()
    results_df['avg_rank'] = (results_df['RMSE_rank'] + results_df['MAE_rank']) / 2
    results_df = results_df.sort_values('avg_rank')

    # Save results
    output_path = plots_dir / "rmse_mae_5celltypes_only.csv"
    results_df.to_csv(output_path, index=False)

    print("\n" + "="*80)
    print("RESULTS: WHICH MODEL BEST FITS THE 5 ORIGINAL CELL TYPES?")
    print("="*80)
    print("\nRanked by combined RMSE/MAE performance (lower is better):\n")
    print(results_df[['model', 'RMSE_5types', 'MAE_5types', 'avg_rank']].to_string(index=False))

    best_model = results_df.iloc[0]
    print("\n" + "-"*80)
    print(f"ANSWER: The {best_model['model']} model most accurately")
    print(f"fits the 5 original cell types (0-4).")
    print(f"  - RMSE: {best_model['RMSE_5types']:.6f}")
    print(f"  - MAE:  {best_model['MAE_5types']:.6f}")
    print("-"*80)

    print(f"\nResults saved to: {output_path}")

    return results_df

def compare_6types_vs_5types():
    """
    Compare metrics when using all 6 types vs only 5 types.

    This shows how including the simulated type 5 affects the results.
    """

    print("\n" + "="*80)
    print("COMPARISON: 6 TYPES VS 5 TYPES ONLY")
    print("="*80)

    # Load previous analysis that used all 6 types
    prev_analysis_path = plots_dir / "residuals_analysis.csv"
    if not prev_analysis_path.exists():
        print("\nPrevious analysis (6 types) not found. Skipping comparison.")
        return

    df_6types = pd.read_csv(prev_analysis_path)

    # Load our new analysis (5 types only)
    df_5types = pd.read_csv(plots_dir / "rmse_mae_5celltypes_only.csv")

    # Merge for comparison (only percentile models have 6-type analysis)
    comparison = pd.merge(
        df_6types[['percentile', 'RMSE', 'MAE']],
        df_5types[df_5types['model'] != 'baseline'][['model', 'RMSE_5types', 'MAE_5types']],
        left_on='percentile',
        right_on='model'
    )

    comparison['RMSE_diff'] = comparison['RMSE'] - comparison['RMSE_5types']
    comparison['MAE_diff'] = comparison['MAE'] - comparison['MAE_5types']

    print("\nHow do metrics change when excluding type 5?\n")
    print(comparison.to_string(index=False))

    comparison.to_csv(plots_dir / "comparison_6types_vs_5types.csv", index=False)
    print(f"\nComparison saved to: {plots_dir / 'comparison_6types_vs_5types.csv'}")

def main():
    """
    Main analysis pipeline to answer the professor's question.
    """
    print("\n" + "="*80)
    print("ANSWERING PROFESSOR'S QUESTION")
    print("="*80)
    print("\nQuestion: Which model most accurately fits the 5 cell types?")
    print("Method: Calculate RMSE/MAE for residuals from cell types 0-4 only,")
    print("        ignoring the added 6th type (type 5).\n")

    # Run main analysis
    results_df = analyze_5celltypes_only()

    # Compare with previous analysis that used all 6 types
    compare_6types_vs_5types()

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()
