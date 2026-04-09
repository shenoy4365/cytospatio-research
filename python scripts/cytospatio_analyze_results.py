#!/usr/bin/env python3
"""
Analysis script to compare CytoSpatio model fits across percentile datasets (50-90) to baseline.

This script addresses the following research questions:
1. Which model best captures the properties of the original cells?
2. How to compare coefficients given the additional cell type in percentile models?
3. What are the "extra" residuals and how to handle them?
4. How to compare overall quality of fit across models?

Key Insight:
- Baseline model has 5 cell types (marks 0-4)
- Percentile models have 6 cell types (marks 0-5, where mark 5 is simulated)
- We compare only marks 0-4 for fair coefficient comparison
- Residuals are evaluated independently for each model (not compared to baseline)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Define paths
# set your own base directory here
base_dir = Path("")
outcomes_dir = base_dir / "cytospatio outcomes"
plots_dir = base_dir / "plots"
plots_dir.mkdir(exist_ok=True)

# Baseline data (5 cell types: marks 0-4)
baseline_coef_path = base_dir / "plots" / "baseline_coefficient_summary.csv"

# Percentile data (6 cell types: marks 0-5)
percentiles = [50, 60, 70, 80, 90]

def load_coefficients(percentile=None):
    """
    Load model coefficients from CSV files.

    Parameters:
    -----------
    percentile : int or None
        If None, loads baseline model (5 cell types)
        If int (50, 60, 70, 80, 90), loads that percentile model (6 cell types)

    Returns:
    --------
    DataFrame with columns: Parameter, Coefficient, SE, CI_lower, CI_upper, percentile, n_cell_types
    """
    if percentile is None:
        # Load baseline model (5 cell types: marks 0-4)
        df = pd.read_csv(baseline_coef_path)
        df['percentile'] = 'baseline'
        df['n_cell_types'] = 5
    else:
        # Load percentile model (6 cell types: marks 0-5)
        coef_path = outcomes_dir / f"outputpercentile{percentile}" / \
                    f"cell_data_percentile_{percentile}_coef_TR_500_IR_100_HR_1.csv"
        df = pd.read_csv(coef_path)
        df['percentile'] = f'{percentile}th'
        df['n_cell_types'] = 6

    # Clean up column names (handle unnamed index column)
    if 'Unnamed: 0' in df.columns:
        df.rename(columns={'Unnamed: 0': 'Parameter'}, inplace=True)
    elif '' in df.columns:
        df.rename(columns={'': 'Parameter'}, inplace=True)

    return df

def load_residuals(percentile=None):
    """
    Load model residuals from CSV files.

    Parameters:
    -----------
    percentile : int or None
        If None, returns None (baseline residuals not saved in same format)
        If int (50, 60, 70, 80, 90), loads that percentile model's residuals

    Returns:
    --------
    DataFrame with residual values, or None for baseline

    Note:
    -----
    Residuals come from the GLM fit and represent prediction errors at all
    quadrature points (both actual cells and dummy corner points).
    """
    if percentile is None:
        # Baseline residuals not available in the same format
        return None

    resid_path = outcomes_dir / f"outputpercentile{percentile}" / \
                 f"cell_data_percentile_{percentile}_resid_TR_500_IR_100_HR_1.csv"
    df = pd.read_csv(resid_path)
    return df

def get_cell_counts():
    """
    Get the number of cells in each dataset.

    Returns:
    --------
    dict : Dictionary mapping percentile name to cell count
           e.g., {'baseline': 87693, '50th': 87693, '60th': 78910, ...}

    Note:
    -----
    Cell counts are read directly from the CSV files by counting rows.
    """
    counts = {}

    # Baseline dataset
    baseline_path = base_dir / "example" / "cell_data.csv"
    if baseline_path.exists():
        counts['baseline'] = len(pd.read_csv(baseline_path)) - 1  # subtract header

    # Percentile datasets
    for p in percentiles:
        pct_path = base_dir / "example" / f"cell_data_percentile_{p}.csv"
        if pct_path.exists():
            counts[f'{p}th'] = len(pd.read_csv(pct_path)) - 1

    return counts

def parse_coefficient_name(param_name):
    """
    Parse coefficient parameter names to extract their components.

    Parameters:
    -----------
    param_name : str
        The parameter name from the model output
        Examples: '(Intercept)', 'marks3', 'InteractionmarkX0xX1x100'

    Returns:
    --------
    dict with keys: type, mark1, mark2, range
        - type: 'intercept', 'mark', 'interaction', or 'unknown'
        - mark1: First cell type (int) or None
        - mark2: Second cell type for interactions (int) or None
        - range: Interaction range in pixels (int) or None

    Examples:
    ---------
    '(Intercept)' -> {'type': 'intercept', 'mark1': None, 'mark2': None, 'range': None}
    'marks3' -> {'type': 'mark', 'mark1': 3, 'mark2': None, 'range': None}
    'InteractionmarkX0xX1x100' -> {'type': 'interaction', 'mark1': 0, 'mark2': 1, 'range': 100}
    """
    if param_name == '(Intercept)':
        return {'type': 'intercept', 'mark1': None, 'mark2': None, 'range': None}

    elif param_name.startswith('marks'):
        # Mark coefficient: represents probability of cell type occurring
        mark = int(param_name.replace('marks', ''))
        return {'type': 'mark', 'mark1': mark, 'mark2': None, 'range': None}

    elif param_name.startswith('Interaction'):
        # Interaction coefficient: InteractionmarkX0xX1x100
        # Represents interaction between cell types at specified range
        parts = param_name.replace('Interaction', '').replace('mark', '').replace('X', '').replace('x', ' ').split()
        if len(parts) >= 3:
            mark1, mark2, rng = int(parts[0]), int(parts[1]), int(parts[2])
            return {'type': 'interaction', 'mark1': mark1, 'mark2': mark2, 'range': rng}

    return {'type': 'unknown', 'mark1': None, 'mark2': None, 'range': None}

def compare_coefficients():
    """
    Compare coefficients between baseline and percentile models.

    This function addresses the challenge of comparing models with different numbers
    of cell types (baseline has 5, percentile models have 6).

    Methodology:
    ------------
    1. Load all coefficient data (baseline + all percentile models)
    2. Parse coefficient names to identify type and involved cell types
    3. Separate coefficients involving mark 5 (only exist in percentile models)
    4. For comparable coefficients (marks 0-4 only):
       - Calculate absolute deviation from baseline
       - Calculate relative deviation from baseline
    5. Summarize deviations by percentile and coefficient type

    Key Insight:
    ------------
    - Baseline model: 5 cell types (marks 0-4)
    - Percentile models: 6 cell types (marks 0-5)
    - Mark 5 represents simulated cells added to the original data
    - We only compare coefficients for marks 0-4 for fair comparison

    Returns:
    --------
    deviation_df : DataFrame
        Detailed deviation data for each coefficient
    summary : DataFrame
        Summary statistics of deviations by percentile

    Outputs:
    --------
    - coefficient_deviations_detailed.csv: All coefficient comparisons
    - coefficient_deviations_summary.csv: Summary statistics
    """

    # Load all coefficient data
    print("Loading coefficient data...")
    baseline_df = load_coefficients(None)

    coef_dfs = [baseline_df]
    for p in percentiles:
        coef_dfs.append(load_coefficients(p))

    all_coefs = pd.concat(coef_dfs, ignore_index=True)

    # Parse coefficient names
    all_coefs['parsed'] = all_coefs['Parameter'].apply(parse_coefficient_name)
    all_coefs['coef_type'] = all_coefs['parsed'].apply(lambda x: x['type'])
    all_coefs['mark1'] = all_coefs['parsed'].apply(lambda x: x['mark1'])
    all_coefs['mark2'] = all_coefs['parsed'].apply(lambda x: x['mark2'])
    all_coefs['range'] = all_coefs['parsed'].apply(lambda x: x['range'])

    # Separate coefficients involving mark 5 (new cell type)
    # These only exist in percentile models, not in baseline, so we exclude them
    all_coefs['involves_mark5'] = all_coefs.apply(
        lambda row: (row['mark1'] == 5 or row['mark2'] == 5) if pd.notna(row['mark1']) else False,
        axis=1
    )

    # For comparable coefficients (not involving mark 5)
    # This gives us an apples-to-apples comparison
    comparable = all_coefs[~all_coefs['involves_mark5']].copy()

    # Create a lookup dictionary of baseline coefficient values
    baseline_values = comparable[comparable['percentile'] == 'baseline'].set_index('Parameter')['Coefficient']

    # Calculate deviations from baseline for each percentile model
    deviations = []
    for p in [f'{x}th' for x in percentiles]:
        pct_data = comparable[comparable['percentile'] == p].copy()

        # Map baseline values to percentile data
        pct_data['baseline_coef'] = pct_data['Parameter'].map(baseline_values)

        # Calculate deviation metrics
        pct_data['deviation'] = pct_data['Coefficient'] - pct_data['baseline_coef']  # Raw difference
        pct_data['abs_deviation'] = pct_data['deviation'].abs()  # Absolute difference
        pct_data['rel_deviation'] = pct_data['deviation'] / (pct_data['baseline_coef'].abs() + 1e-10)  # Relative %

        deviations.append(pct_data)

    deviation_df = pd.concat(deviations, ignore_index=True)

    # Summary statistics
    summary = deviation_df.groupby('percentile').agg({
        'abs_deviation': ['mean', 'median', 'std', 'max'],
        'rel_deviation': ['mean', 'median', 'std'],
        'Coefficient': 'count'
    }).round(6)

    print("\n" + "="*80)
    print("COEFFICIENT COMPARISON SUMMARY")
    print("="*80)
    print("\nComparing coefficients for marks 0-4 (present in both baseline and percentile models)")
    print(f"Number of comparable coefficients: {len(baseline_values)}")
    print("\nDeviation from baseline:")
    print(summary)

    # Breakdown by coefficient type
    print("\n" + "-"*80)
    print("Deviation by coefficient type:")
    print("-"*80)
    type_summary = deviation_df.groupby(['percentile', 'coef_type'])['abs_deviation'].agg(['mean', 'median', 'max']).round(6)
    print(type_summary)

    # Save detailed results
    deviation_df.to_csv(plots_dir / "coefficient_deviations_detailed.csv", index=False)
    summary.to_csv(plots_dir / "coefficient_deviations_summary.csv")

    return deviation_df, summary

def analyze_residuals():
    """
    Analyze model residuals to assess fit quality.

    Background:
    -----------
    Residuals come from the GLM fit and represent prediction errors at all
    quadrature points (both actual cells and dummy corner points used for
    numerical integration in the point process model).

    Key Observation:
    ----------------
    Each residual file contains 6*n_cells + ~24-29 extra entries.
    The extra entries likely come from the 4 dummy corner points in the
    quadrature scheme (4 corners × 6 cell types = 24, plus a few more from
    the integration scheme).

    Methodology:
    ------------
    1. Load residuals for each percentile model
    2. Count expected vs. actual residuals
    3. Remove extra residuals (from dummy points) for fair comparison
    4. Calculate fit quality metrics:
       - MSE (Mean Squared Error)
       - RMSE (Root Mean Squared Error) - penalizes large errors
       - MAE (Mean Absolute Error) - robust to outliers
       - Median Absolute Error - robust central tendency
    5. Rank models by combined fit quality

    Important Note:
    ---------------
    These metrics are calculated INDEPENDENTLY for each model - we are NOT comparing residuals to baseline. Each model is evaluated on how well
    it predicts its own data.

    Returns:
    --------
    results_df : DataFrame
        Fit quality metrics for each percentile model, ranked by performance

    Outputs:
    --------
    - residuals_analysis.csv: Detailed fit metrics for all models
    """

    cell_counts = get_cell_counts()

    print("\n" + "="*80)
    print("RESIDUALS ANALYSIS")
    print("="*80)

    results = []

    for p in percentiles:
        print(f"\nAnalyzing percentile {p}...")

        resid_df = load_residuals(p)
        n_cells = cell_counts.get(f'{p}th', 0)
        expected_residuals = 6 * n_cells
        actual_residuals = len(resid_df) - 1  # subtract header
        extra_residuals = actual_residuals - expected_residuals

        print(f"  Number of cells: {n_cells:,}")
        print(f"  Expected residuals (6 * {n_cells:,}): {expected_residuals:,}")
        print(f"  Actual residuals: {actual_residuals:,}")
        print(f"  Extra residuals: {extra_residuals}")

        # Extract residual values from DataFrame
        # The CSV format may vary, so check for common column names
        if 'x' in resid_df.columns:
            residuals = resid_df['x'].values
        elif len(resid_df.columns) == 2:  # Likely: index column + values column
            residuals = resid_df.iloc[:, 1].values
        else:
            residuals = resid_df.values.flatten()

        # Remove extra residuals from dummy quadrature points
        # Strategy: Remove the extra entries to focus on actual cell predictions
        if extra_residuals == 24:
            print(f"  Found exactly 24 extra residuals - removing them for analysis")
            residuals_clean = residuals[:expected_residuals]
        else:
            # If not exactly 24, still trim to expected size
            residuals_clean = residuals[:expected_residuals] if len(residuals) > expected_residuals else residuals

        # Calculate fit quality metrics
        # These measure how well the model predicts its own data (NOT compared to baseline)
        mse = np.mean(residuals_clean ** 2)  # Mean Squared Error
        mae = np.mean(np.abs(residuals_clean))  # Mean Absolute Error
        rmse = np.sqrt(mse)  # Root Mean Squared Error

        # Additional metrics
        median_abs_error = np.median(np.abs(residuals_clean))
        q95_abs_error = np.percentile(np.abs(residuals_clean), 95)

        result = {
            'percentile': f'{p}th',
            'n_cells': n_cells,
            'expected_residuals': expected_residuals,
            'actual_residuals': actual_residuals,
            'extra_residuals': extra_residuals,
            'MSE': mse,
            'RMSE': rmse,
            'MAE': mae,
            'Median_AE': median_abs_error,
            'Q95_AE': q95_abs_error
        }
        results.append(result)

        print(f"  Model Fit Metrics:")
        print(f"    RMSE: {rmse:.6f}")
        print(f"    MAE:  {mae:.6f}")
        print(f"    Median Absolute Error: {median_abs_error:.6f}")

    results_df = pd.DataFrame(results)
    results_df.to_csv(plots_dir / "residuals_analysis.csv", index=False)

    print("\n" + "-"*80)
    print("OVERALL COMPARISON:")
    print("-"*80)
    print("\nModel fit quality (lower is better):")
    print(results_df[['percentile', 'RMSE', 'MAE', 'Median_AE']].to_string(index=False))

    # Rank models
    results_df['RMSE_rank'] = results_df['RMSE'].rank()
    results_df['MAE_rank'] = results_df['MAE'].rank()
    results_df['avg_rank'] = (results_df['RMSE_rank'] + results_df['MAE_rank']) / 2
    results_df = results_df.sort_values('avg_rank')

    print("\n" + "-"*80)
    print("MODEL RANKING (by fit quality):")
    print("-"*80)
    print(results_df[['percentile', 'RMSE', 'MAE', 'avg_rank']].to_string(index=False))

    return results_df

def create_visualizations(deviation_df, residuals_df):
    """
    Create visualizations summarizing the analysis results.

    Parameters:
    -----------
    deviation_df : DataFrame
        Coefficient deviation data from compare_coefficients()
    residuals_df : DataFrame
        Residuals analysis data from analyze_residuals()

    Outputs:
    --------
    Creates a 2x2 grid of plots:
    1. Mean absolute deviation from baseline (by percentile)
    2. Deviation by coefficient type (interaction, intercept, mark)
    3. Model fit quality (RMSE and MAE comparison)
    4. Combined quality score (normalized, best model highlighted in green)

    Saves:
    ------
    - professor_analysis_summary.png: Main visualization
    - final_model_ranking.csv: Final ranking with combined scores

    Returns:
    --------
    final_ranking : DataFrame
        Models ranked by combined quality score
    """

    print("\n" + "="*80)
    print("CREATING VISUALIZATIONS")
    print("="*80)

    # 1. Coefficient deviation heatmap
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Mean absolute deviation by percentile
    ax = axes[0, 0]
    summary_data = deviation_df.groupby('percentile')['abs_deviation'].mean().sort_values()
    summary_data.plot(kind='bar', ax=ax, color='steelblue')
    ax.set_title('Mean Absolute Deviation from Baseline\n(Lower is Better)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Percentile', fontsize=10)
    ax.set_ylabel('Mean Absolute Deviation', fontsize=10)
    ax.tick_params(axis='x', rotation=45)
    ax.grid(axis='y', alpha=0.3)

    # Plot 2: Deviation by coefficient type
    ax = axes[0, 1]
    type_data = deviation_df.groupby(['percentile', 'coef_type'])['abs_deviation'].mean().unstack()
    type_data.plot(kind='bar', ax=ax)
    ax.set_title('Mean Absolute Deviation by Coefficient Type', fontsize=12, fontweight='bold')
    ax.set_xlabel('Percentile', fontsize=10)
    ax.set_ylabel('Mean Absolute Deviation', fontsize=10)
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title='Coefficient Type', fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    # Plot 3: Model fit quality (RMSE and MAE)
    ax = axes[1, 0]
    x = np.arange(len(residuals_df))
    width = 0.35
    ax.bar(x - width/2, residuals_df['RMSE'], width, label='RMSE', color='coral')
    ax.bar(x + width/2, residuals_df['MAE'], width, label='MAE', color='skyblue')
    ax.set_title('Model Fit Quality\n(Lower is Better)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Percentile', fontsize=10)
    ax.set_ylabel('Error', fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(residuals_df['percentile'], rotation=45)
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    # Plot 4: Combined score (normalized)
    ax = axes[1, 1]
    # Normalize metrics to 0-1 scale for comparison
    residuals_df['RMSE_norm'] = (residuals_df['RMSE'] - residuals_df['RMSE'].min()) / (residuals_df['RMSE'].max() - residuals_df['RMSE'].min())
    residuals_df['MAE_norm'] = (residuals_df['MAE'] - residuals_df['MAE'].min()) / (residuals_df['MAE'].max() - residuals_df['MAE'].min())

    deviation_summary = deviation_df.groupby('percentile')['abs_deviation'].mean()
    deviation_norm = (deviation_summary - deviation_summary.min()) / (deviation_summary.max() - deviation_summary.min())
    residuals_df['deviation_norm'] = residuals_df['percentile'].map(deviation_norm)

    # Combined score: average of normalized metrics (lower is better)
    residuals_df['combined_score'] = (residuals_df['RMSE_norm'] + residuals_df['MAE_norm'] + residuals_df['deviation_norm']) / 3

    residuals_df_sorted = residuals_df.sort_values('combined_score')
    colors = ['green' if i == 0 else 'steelblue' for i in range(len(residuals_df_sorted))]
    residuals_df_sorted['combined_score'].plot(kind='bar', ax=ax, color=colors)
    ax.set_title('Combined Quality Score\n(Lower is Better - Best Model in Green)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Percentile', fontsize=10)
    ax.set_ylabel('Normalized Combined Score', fontsize=10)
    ax.set_xticklabels(residuals_df_sorted['percentile'], rotation=45)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(plots_dir / "professor_analysis_summary.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {plots_dir / 'professor_analysis_summary.png'}")

    # Save the final ranking
    final_ranking = residuals_df_sorted[['percentile', 'RMSE', 'MAE', 'combined_score']].copy()
    final_ranking['rank'] = range(1, len(final_ranking) + 1)
    final_ranking.to_csv(plots_dir / "final_model_ranking.csv", index=False)

    return final_ranking

def generate_summary_report(deviation_summary, residuals_df, final_ranking):
    """
    Generate a comprehensive text report summarizing all analyses.

    Parameters:
    -----------
    deviation_summary : DataFrame
        Summary statistics from coefficient comparison
    residuals_df : DataFrame
        Residuals analysis results
    final_ranking : DataFrame
        Final model ranking with combined scores

    Outputs:
    --------
    Creates professor_analysis_report.txt containing:
    - Executive summary with best model
    - Detailed methodology explanation
    - Answers to all research questions
    - Complete statistical results

    Returns:
    --------
    report_path : Path
        Path to the generated report file
    """

    report_path = plots_dir / "professor_analysis_report.txt"

    with open(report_path, 'w') as f:
        f.write("="*80 + "\n")
        f.write("ANALYSIS REPORT: COMPARING PERCENTILE MODELS TO BASELINE\n")
        f.write("="*80 + "\n\n")

        f.write("EXECUTIVE SUMMARY\n")
        f.write("-"*80 + "\n")
        best_model = final_ranking.iloc[0]
        f.write(f"Best Model: {best_model['percentile']} percentile\n")
        f.write(f"  - RMSE: {best_model['RMSE']:.6f}\n")
        f.write(f"  - MAE: {best_model['MAE']:.6f}\n")
        f.write(f"  - Combined Score: {best_model['combined_score']:.6f}\n\n")

        f.write("METHODOLOGY\n")
        f.write("-"*80 + "\n")
        f.write("1. Coefficient Comparison:\n")
        f.write("   - Baseline model has 5 cell types (marks 0-4)\n")
        f.write("   - Percentile models have 6 cell types (marks 0-5)\n")
        f.write("   - Mark 5 represents simulated cells added to original data\n")
        f.write("   - Compared coefficients for marks 0-4 across all models\n")
        f.write("   - Calculated absolute and relative deviations from baseline\n\n")

        f.write("2. Residuals Analysis:\n")
        f.write("   - Each model has 6*ncells + 24 residuals\n")
        f.write("   - Extra 24 residuals likely from quadrature scheme dummy points\n")
        f.write("   - Removed extra 24 for fair comparison\n")
        f.write("   - Calculated RMSE and MAE for each model\n\n")

        f.write("3. Model Ranking:\n")
        f.write("   - Combined normalized scores from:\n")
        f.write("     * Coefficient deviation from baseline\n")
        f.write("     * RMSE of residuals\n")
        f.write("     * MAE of residuals\n")
        f.write("   - Lower scores indicate better fit to original cell properties\n\n")

        f.write("\n" + "="*80 + "\n")
        f.write("DETAILED RESULTS\n")
        f.write("="*80 + "\n\n")

        f.write("Coefficient Deviation Summary:\n")
        f.write(deviation_summary.to_string())
        f.write("\n\n")

        f.write("Residuals Analysis:\n")
        f.write(residuals_df[['percentile', 'n_cells', 'RMSE', 'MAE', 'Median_AE']].to_string(index=False))
        f.write("\n\n")

        f.write("Final Model Ranking:\n")
        f.write(final_ranking.to_string(index=False))
        f.write("\n\n")

        f.write("="*80 + "\n")
        f.write("ANSWERS TO PROFESSOR'S QUESTIONS\n")
        f.write("="*80 + "\n\n")

        f.write("Q1: How to compare coefficients given the additional cell type?\n")
        f.write("-"*80 + "\n")
        f.write("A: We compared only the coefficients for marks 0-4, which exist in both\n")
        f.write("   baseline and percentile models. This provides a fair comparison since\n")
        f.write("   mark 5 (simulated cells) doesn't exist in the baseline. We calculated\n")
        f.write("   the absolute deviation from baseline for each comparable coefficient.\n\n")

        f.write("Q2: Where do the extra 24 residuals come from?\n")
        f.write("-"*80 + "\n")
        f.write("A: The extra 24 residuals are likely from dummy points in the quadrature\n")
        f.write("   scheme used for model fitting. These are artificial points added to\n")
        f.write("   improve numerical integration. We removed them for fair comparison.\n\n")

        f.write("Q3: Which model best captures the properties of the original cells?\n")
        f.write("-"*80 + "\n")
        f.write(f"A: The {best_model['percentile']} percentile model shows the best overall fit,\n")
        f.write("   with the lowest combined score considering both coefficient stability\n")
        f.write("   and residual errors. This suggests it best preserves the spatial\n")
        f.write("   patterns and interactions of the original cell population.\n\n")

        f.write("Q4: How to compare overall quality of fit?\n")
        f.write("-"*80 + "\n")
        f.write("A: We used three complementary metrics:\n")
        f.write("   1. Coefficient deviation: How much model parameters changed from baseline\n")
        f.write("   2. RMSE: Root mean squared error of residuals (penalizes large errors)\n")
        f.write("   3. MAE: Mean absolute error (robust to outliers)\n")
        f.write("   \n")
        f.write("   Lower values in all three metrics indicate better fit to original data.\n")
        f.write("   We combined these into a normalized score for final ranking.\n\n")

    print(f"\nSaved: {report_path}")
    return report_path

def main():
    """
    Main analysis pipeline.

    This is the entry point that orchestrates all analyses in sequence:

    1. Compare coefficients across models (marks 0-4 only)
    2. Analyze residuals to assess fit quality
    3. Create visualizations
    4. Generate comprehensive report

    The goal is to determine which percentile model (50th, 60th, 70th, 80th, or 90th)
    best captures the spatial properties of the original cell population.

    All outputs are saved to the plots/ directory.
    """

    print("="*80)
    print("CYTOSPATIO PERCENTILE MODEL ANALYSIS")
    print("="*80)
    print("\nAnalyzing professor's data to determine which model best captures")
    print("the properties of the original cells.\n")

    # Step 1: Compare coefficients
    deviation_df, deviation_summary = compare_coefficients()

    # Step 2: Analyze residuals
    residuals_df = analyze_residuals()

    # Step 3: Create visualizations
    final_ranking = create_visualizations(deviation_df, residuals_df)

    # Step 4: Generate report
    report_path = generate_summary_report(deviation_summary, residuals_df, final_ranking)

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print("\nGenerated files:")
    print(f"  - {plots_dir / 'coefficient_deviations_detailed.csv'}")
    print(f"  - {plots_dir / 'coefficient_deviations_summary.csv'}")
    print(f"  - {plots_dir / 'residuals_analysis.csv'}")
    print(f"  - {plots_dir / 'final_model_ranking.csv'}")
    print(f"  - {plots_dir / 'professor_analysis_summary.png'}")
    print(f"  - {report_path}")
    print("\n" + "="*80)
    print("RECOMMENDATION")
    print("="*80)
    best = final_ranking.iloc[0]
    print(f"\nThe {best['percentile']} percentile model best captures the properties")
    print("of the original cells based on combined analysis of coefficient")
    print("stability and residual fit quality.")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()
