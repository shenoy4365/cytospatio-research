#!/usr/bin/env python3
"""
Baseline Coefficient Visualization Script

This script visualizes the coefficients from the baseline CytoSpatio analysis.
Since the centroid-modified datasets exceed memory limits, we only have baseline data.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def load_baseline_coefficients(filepath="output_baseline/cell_data_coef_TR_500_IR_100_HR_1.csv"):
    """Load baseline coefficient data."""
    print(f"Loading baseline coefficients from: {filepath}")
    df = pd.read_csv(filepath)
    print(f"  Found {len(df)} coefficients")
    return df


def plot_coefficient_categories(df, output_dir="plots"):
    """
    Generate plots organized by coefficient type.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Extract coefficient names and values
    coef_names = df.iloc[:, 0].values  # "Unnamed: 0" column has the names
    coef_values = df['Coefficient'].values  # "Coefficient" column has the values

    # Separate into categories
    intercept = []
    marks = []
    interactions = []

    for i, name in enumerate(coef_names):
        if name == "(Intercept)":
            intercept.append((name, coef_values[i]))
        elif name.startswith("marks"):
            marks.append((name, coef_values[i]))
        elif name.startswith("Interaction"):
            interactions.append((name, coef_values[i]))

    print(f"\nCoefficient breakdown:")
    print(f"  Intercept: {len(intercept)}")
    print(f"  Marks: {len(marks)}")
    print(f"  Interactions: {len(interactions)}")

    # Plot 1: All coefficients in order
    fig, ax = plt.subplots(figsize=(20, 8))

    colors = []
    for name in coef_names:
        if name == "(Intercept)":
            colors.append('red')
        elif name.startswith("marks"):
            colors.append('blue')
        else:
            colors.append('green')

    ax.bar(range(len(coef_names)), coef_values, color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.set_xlabel('Coefficient Index', fontsize=12)
    ax.set_ylabel('Coefficient Value', fontsize=12)
    ax.set_title('Baseline CytoSpatio Coefficients (All)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', alpha=0.7, label='Intercept'),
        Patch(facecolor='blue', alpha=0.7, label='Marks'),
        Patch(facecolor='green', alpha=0.7, label='Interactions')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    output_file = os.path.join(output_dir, "baseline_all_coefficients.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n  Saved: {output_file}")
    plt.close()

    # Plot 2: Marks coefficients only
    fig, ax = plt.subplots(figsize=(10, 6))

    if marks:
        mark_names, mark_values = zip(*marks)
        x_pos = range(len(mark_names))
        ax.bar(x_pos, mark_values, color='steelblue', alpha=0.8, edgecolor='black', linewidth=1)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(mark_names, rotation=0, fontsize=11)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
        ax.set_ylabel('Coefficient Value', fontsize=12)
        ax.set_title('Baseline Mark Coefficients', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        output_file = os.path.join(output_dir, "baseline_marks_coefficients.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_file}")
        plt.close()

    # Plot 3: Interaction coefficients by range
    fig, axes = plt.subplots(5, 1, figsize=(16, 20))

    ranges = [100, 200, 300, 400, 500]
    for i, r in enumerate(ranges):
        ax = axes[i]

        # Filter interactions for this range
        range_interactions = [(name, val) for name, val in interactions if f"x{r}" in name]

        if range_interactions:
            int_names, int_values = zip(*range_interactions)
            x_pos = range(len(int_names))
            ax.bar(x_pos, int_values, color='forestgreen', alpha=0.8, edgecolor='black', linewidth=0.5)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(int_names, rotation=90, fontsize=8)
            ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
            ax.set_ylabel('Coefficient Value', fontsize=10)
            ax.set_title(f'Interaction Coefficients at Range {r}', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    output_file = os.path.join(output_dir, "baseline_interactions_by_range.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

    # Save summary CSV
    summary_file = os.path.join(output_dir, "baseline_coefficient_summary.csv")
    df.to_csv(summary_file, index=False)
    print(f"  Saved: {summary_file}")

    print(f"\nAll plots saved to: {output_dir}/")


def main():
    """Main function."""
    print("=" * 60)
    print("CytoSpatio Baseline Coefficient Visualization")
    print("=" * 60)

    # Load baseline coefficients
    df = load_baseline_coefficients()

    # Generate plots
    print("\nGenerating plots...")
    plot_coefficient_categories(df)

    print("\n" + "=" * 60)
    print("Visualization complete!")
    print("=" * 60)
    print("\nNote: Centroid-modified datasets could not be analyzed")
    print("due to memory limitations (see CENTROID_ANALYSIS_SUMMARY.md)")


if __name__ == "__main__":
    main()
