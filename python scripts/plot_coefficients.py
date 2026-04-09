#!/usr/bin/env python3
"""
Coefficient Analysis and Plotting Script

This script:
1. Reads coefficient CSV files from all CytoSpatio output folders
2. Extracts coefficient values for each threshold
3. Generates plots showing how each coefficient changes with threshold
4. Scales coefficients by their max value for readability
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from pathlib import Path


def find_coef_files(base_dir="."):
    """
    Find all coefficient CSV files in output directories.

    Returns:
        Dictionary mapping threshold values to coefficient file paths
    """
    coef_files = {}

    # Baseline (original data, no threshold)
    baseline_pattern = os.path.join(base_dir, "output_baseline", "*_coef_TR_500_IR_100_HR_1.csv")
    baseline_files = glob.glob(baseline_pattern)
    if baseline_files:
        coef_files[0] = baseline_files[0]  # Use 0 for baseline

    # Threshold files
    for threshold in [50, 60, 70, 80, 90]:
        pattern = os.path.join(base_dir, f"output_percentile_{threshold}", "*_coef_TR_500_IR_100_HR_1.csv")
        files = glob.glob(pattern)
        if files:
            coef_files[threshold] = files[0]

    return coef_files


def load_all_coefficients(coef_files):
    """
    Load coefficient data from all files.

    Args:
        coef_files: Dictionary mapping thresholds to file paths

    Returns:
        Dictionary mapping thresholds to DataFrames of coefficients
    """
    all_data = {}

    for threshold, filepath in sorted(coef_files.items()):
        print(f"Loading coefficients from threshold {threshold}: {filepath}")
        df = pd.read_csv(filepath)
        all_data[threshold] = df
        print(f"  Found {len(df)} coefficients")

    return all_data


def extract_coefficient_series(all_data):
    """
    Extract time series for each coefficient across thresholds.

    Args:
        all_data: Dictionary mapping thresholds to coefficient DataFrames

    Returns:
        DataFrame with columns = coefficient names, index = thresholds
    """
    # Get all unique coefficient names from the first dataset
    first_df = list(all_data.values())[0]

    # Assume first column is coefficient name, second is value
    # Check the actual column names
    print(f"\nColumn names in coefficient file: {list(first_df.columns)}")

    # Build a dataframe with thresholds as rows, coefficients as columns
    thresholds = sorted(all_data.keys())
    coef_names = first_df.iloc[:, 0].values  # First column = coefficient names

    # Initialize result dataframe
    result = pd.DataFrame(index=thresholds, columns=coef_names)

    for threshold in thresholds:
        df = all_data[threshold]
        for i, coef_name in enumerate(coef_names):
            # Find this coefficient in the current threshold's data
            matching_rows = df[df.iloc[:, 0] == coef_name]
            if len(matching_rows) > 0:
                # Get the coefficient value (typically second column)
                coef_value = matching_rows.iloc[0, 1]
                result.loc[threshold, coef_name] = coef_value
            else:
                result.loc[threshold, coef_name] = np.nan

    # Convert to numeric
    result = result.apply(pd.to_numeric, errors='coerce')

    return result


def plot_coefficients(coef_series, output_dir="plots"):
    """
    Generate plots of coefficients vs threshold.

    Args:
        coef_series: DataFrame with thresholds as index, coefficients as columns
        output_dir: Directory to save plots
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    thresholds = coef_series.index.values
    n_coefficients = len(coef_series.columns)

    print(f"\nGenerating plots for {n_coefficients} coefficients...")

    # Plot 1: All coefficients on one plot (scaled by max)
    fig, ax = plt.subplots(figsize=(12, 8))

    for coef_name in coef_series.columns:
        values = coef_series[coef_name].values

        # Scale by max absolute value
        max_val = np.nanmax(np.abs(values))
        if max_val > 0:
            scaled_values = values / max_val
        else:
            scaled_values = values

        ax.plot(thresholds, scaled_values, marker='o', label=coef_name, linewidth=2, markersize=6)

    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Threshold Percentile', fontsize=12)
    ax.set_ylabel('Coefficient Value (scaled by max)', fontsize=12)
    ax.set_title('All Coefficients vs Threshold (Scaled by Max)', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)

    # Set x-axis ticks
    ax.set_xticks([0] + list(range(50, 100, 10)))
    ax.set_xticklabels(['Baseline'] + [f'{t}%' for t in range(50, 100, 10)])

    plt.tight_layout()
    output_file = os.path.join(output_dir, "all_coefficients_scaled.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

    # Plot 2: Individual plots for each coefficient (unscaled)
    n_cols = 4
    n_rows = int(np.ceil(n_coefficients / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4 * n_rows))
    axes = axes.flatten() if n_coefficients > 1 else [axes]

    for i, coef_name in enumerate(coef_series.columns):
        ax = axes[i]
        values = coef_series[coef_name].values

        ax.plot(thresholds, values, marker='o', color='steelblue', linewidth=2, markersize=8)
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax.set_xlabel('Threshold Percentile', fontsize=10)
        ax.set_ylabel('Coefficient Value', fontsize=10)
        ax.set_title(coef_name, fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Set x-axis ticks
        ax.set_xticks([0] + list(range(50, 100, 10)))
        ax.set_xticklabels(['Base'] + [f'{t}%' for t in range(50, 100, 10)], fontsize=8)

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    plt.suptitle('Individual Coefficient Trends vs Threshold', fontsize=16, fontweight='bold', y=1.00)
    plt.tight_layout()

    output_file = os.path.join(output_dir, "individual_coefficients.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

    # Plot 3: Summary CSV
    summary_file = os.path.join(output_dir, "coefficient_summary.csv")
    coef_series.to_csv(summary_file)
    print(f"  Saved: {summary_file}")

    print(f"\nAll plots saved to: {output_dir}/")


def main():
    """
    Main function to orchestrate coefficient analysis.
    """
    print("=" * 60)
    print("CytoSpatio Coefficient Analysis")
    print("=" * 60)

    # Find all coefficient files
    print("\nSearching for coefficient files...")
    coef_files = find_coef_files()

    if not coef_files:
        print("Error: No coefficient files found!")
        print("Make sure CytoSpatio analyses have completed.")
        return

    print(f"Found {len(coef_files)} coefficient files:")
    for threshold, filepath in sorted(coef_files.items()):
        label = "Baseline" if threshold == 0 else f"{threshold}th percentile"
        print(f"  {label}: {filepath}")

    # Load all coefficients
    print("\nLoading coefficient data...")
    all_data = load_all_coefficients(coef_files)

    # Extract coefficient series
    print("\nExtracting coefficient time series...")
    coef_series = extract_coefficient_series(all_data)
    print(f"Extracted {len(coef_series.columns)} coefficients across {len(coef_series)} thresholds")

    # Generate plots
    print("\nGenerating plots...")
    plot_coefficients(coef_series)

    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
