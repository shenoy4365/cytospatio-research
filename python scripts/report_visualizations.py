#!/usr/bin/env python3
"""
Create visualizations for the updated research report.

This script generates:
1. RMSE vs Percentile plot
2. MAE vs Percentile plot
3. Combined error metrics plot
4. Model ranking visualization
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Data from the analysis
models = ['Baseline', '50th', '60th', '70th', '80th', '85th', '90th', '95th']
percentiles = [None, 50, 60, 70, 80, 85, 90, 95]
rmse = [0.000687, 0.001323, 0.001152, 0.000985, 0.000823, 0.000776, 0.000691, 0.000682]
mae = [0.000379, 0.000875, 0.000741, 0.000609, 0.000485, 0.000446, 0.000391, 0.000377]
centroids_added = [0, 43794, 35011, 26269, 17546, 13151, 8770, 4400]

# For percentile models only (excluding baseline)
percentile_models = ['50th', '60th', '70th', '80th', '85th', '90th', '95th']
percentile_values = [50, 60, 70, 80, 85, 90, 95]
rmse_percentile = [0.001323, 0.001152, 0.000985, 0.000823, 0.000776, 0.000691, 0.000682]
mae_percentile = [0.000875, 0.000741, 0.000609, 0.000485, 0.000446, 0.000391, 0.000377]
centroids_percentile = [43794, 35011, 26269, 17546, 13151, 8770, 4400]

# Create figure with subplots
fig = plt.figure(figsize=(16, 10))

# Plot 1: RMSE vs Percentile
ax1 = plt.subplot(2, 3, 1)
ax1.plot(percentile_values, rmse_percentile, 'o-', linewidth=2, markersize=8, color='#e74c3c', label='Percentile Models')
ax1.axhline(y=0.000687, color='#3498db', linestyle='--', linewidth=2, label='Baseline Model')
ax1.axhline(y=0.000682, color='#2ecc71', linestyle=':', linewidth=2, alpha=0.7, label='Best (95th)')
ax1.set_xlabel('Percentile Threshold', fontsize=12, fontweight='bold')
ax1.set_ylabel('RMSE', fontsize=12, fontweight='bold')
ax1.set_title('RMSE vs Percentile Threshold', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(45, 100)

# Plot 2: MAE vs Percentile
ax2 = plt.subplot(2, 3, 2)
ax2.plot(percentile_values, mae_percentile, 'o-', linewidth=2, markersize=8, color='#e74c3c', label='Percentile Models')
ax2.axhline(y=0.000379, color='#3498db', linestyle='--', linewidth=2, label='Baseline Model')
ax2.axhline(y=0.000377, color='#2ecc71', linestyle=':', linewidth=2, alpha=0.7, label='Best (95th)')
ax2.set_xlabel('Percentile Threshold', fontsize=12, fontweight='bold')
ax2.set_ylabel('MAE', fontsize=12, fontweight='bold')
ax2.set_title('MAE vs Percentile Threshold', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(45, 100)

# Plot 3: Centroids Added vs Percentile
ax3 = plt.subplot(2, 3, 3)
ax3.plot(percentile_values, centroids_percentile, 'o-', linewidth=2, markersize=8, color='#9b59b6')
ax3.set_xlabel('Percentile Threshold', fontsize=12, fontweight='bold')
ax3.set_ylabel('Number of Centroids Added', fontsize=12, fontweight='bold')
ax3.set_title('Centroid Count vs Percentile Threshold', fontsize=14, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.set_xlim(45, 100)

# Plot 4: Combined RMSE and MAE on same plot
ax4 = plt.subplot(2, 3, 4)
ax4_twin = ax4.twinx()
line1 = ax4.plot(percentile_values, rmse_percentile, 'o-', linewidth=2, markersize=8, color='#e74c3c', label='RMSE')
line2 = ax4_twin.plot(percentile_values, mae_percentile, 's-', linewidth=2, markersize=8, color='#3498db', label='MAE')
ax4.axhline(y=0.000687, color='#e74c3c', linestyle='--', linewidth=1, alpha=0.5)
ax4_twin.axhline(y=0.000379, color='#3498db', linestyle='--', linewidth=1, alpha=0.5)
ax4.set_xlabel('Percentile Threshold', fontsize=12, fontweight='bold')
ax4.set_ylabel('RMSE', fontsize=12, fontweight='bold', color='#e74c3c')
ax4_twin.set_ylabel('MAE', fontsize=12, fontweight='bold', color='#3498db')
ax4.set_title('Combined Error Metrics vs Percentile', fontsize=14, fontweight='bold')
ax4.tick_params(axis='y', labelcolor='#e74c3c')
ax4_twin.tick_params(axis='y', labelcolor='#3498db')
ax4.grid(True, alpha=0.3)
ax4.set_xlim(45, 100)
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax4.legend(lines, labels, loc='upper right', fontsize=10)

# Plot 5: Model Ranking Bar Chart
ax5 = plt.subplot(2, 3, 5)
ranking_order = [7, 0, 6, 5, 4, 3, 2, 1]  # 95th, Baseline, 90th, 85th, 80th, 70th, 60th, 50th
models_ordered = ['95th', 'Baseline', '90th', '85th', '80th', '70th', '60th', '50th']
rmse_ordered = [rmse[i] for i in ranking_order]
colors = ['#2ecc71', '#3498db', '#95a5a6', '#95a5a6', '#95a5a6', '#e67e22', '#e67e22', '#e74c3c']
bars = ax5.barh(models_ordered, rmse_ordered, color=colors, edgecolor='black', linewidth=1.5)
ax5.set_xlabel('RMSE', fontsize=12, fontweight='bold')
ax5.set_ylabel('Model', fontsize=12, fontweight='bold')
ax5.set_title('Model Ranking by RMSE (Lower is Better)', fontsize=14, fontweight='bold')
ax5.invert_yaxis()
ax5.grid(True, alpha=0.3, axis='x')
# Add value labels
for i, (bar, val) in enumerate(zip(bars, rmse_ordered)):
    ax5.text(val + 0.00005, i, f'{val:.6f}', va='center', fontsize=9, fontweight='bold')

# Plot 6: Scatter plot - Centroids vs Error
ax6 = plt.subplot(2, 3, 6)
scatter = ax6.scatter(centroids_added[1:], rmse[1:], s=200, c=percentile_values, cmap='RdYlGn_r',
                     edgecolor='black', linewidth=2, alpha=0.8)
ax6.scatter([0], [0.000687], s=300, c='blue', marker='*', edgecolor='black', linewidth=2,
           label='Baseline', zorder=5)
ax6.scatter([4400], [0.000682], s=300, c='green', marker='*', edgecolor='black', linewidth=2,
           label='95th (Best)', zorder=5)
ax6.set_xlabel('Number of Centroids Added', fontsize=12, fontweight='bold')
ax6.set_ylabel('RMSE', fontsize=12, fontweight='bold')
ax6.set_title('Error vs Augmentation Level', fontsize=14, fontweight='bold')
ax6.legend(fontsize=10, loc='upper left')
ax6.grid(True, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax6)
cbar.set_label('Percentile Threshold', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig('plots/research_report_visualizations.png', dpi=300, bbox_inches='tight')
print("✓ Saved: plots/research_report_visualizations.png")

# Create individual high-res plots for report

# Individual Plot 1: Clean RMSE vs Percentile for publication
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(percentile_values, rmse_percentile, 'o-', linewidth=3, markersize=10, color='#e74c3c', label='Percentile Models')
ax.axhline(y=0.000687, color='#3498db', linestyle='--', linewidth=2.5, label='Baseline Model', alpha=0.8)
ax.scatter([95], [0.000682], s=400, c='#2ecc71', marker='*', edgecolor='black', linewidth=2,
          label='Best Model (95th)', zorder=5)
ax.set_xlabel('Percentile Threshold', fontsize=14, fontweight='bold')
ax.set_ylabel('RMSE (Root Mean Squared Error)', fontsize=14, fontweight='bold')
ax.set_title('Model Performance: RMSE vs Percentile Threshold', fontsize=16, fontweight='bold', pad=20)
ax.legend(fontsize=12, loc='upper right')
ax.grid(True, alpha=0.3, linewidth=1)
ax.set_xlim(45, 100)
plt.tight_layout()
plt.savefig('plots/rmse_vs_percentile_clean.png', dpi=300, bbox_inches='tight')
print("✓ Saved: plots/rmse_vs_percentile_clean.png")

# Individual Plot 2: Summary table as image
fig, ax = plt.subplots(figsize=(12, 6))
ax.axis('tight')
ax.axis('off')

table_data = [
    ['Rank', 'Model', 'RMSE', 'MAE', 'Centroids Added', '% Worse than Best'],
    ['1', '95th percentile', '0.000682', '0.000377', '4,400', '—'],
    ['2', 'Baseline', '0.000687', '0.000379', '0', '+0.7%'],
    ['3', '90th percentile', '0.000691', '0.000391', '8,770', '+1.3%'],
    ['4', '85th percentile', '0.000776', '0.000446', '13,151', '+13.8%'],
    ['5', '80th percentile', '0.000823', '0.000485', '17,546', '+20.7%'],
    ['6', '70th percentile', '0.000985', '0.000609', '26,269', '+44.4%'],
    ['7', '60th percentile', '0.001152', '0.000741', '35,011', '+68.9%'],
    ['8', '50th percentile', '0.001323', '0.000875', '43,794', '+94.0%'],
]

table = ax.table(cellText=table_data, cellLoc='center', loc='center',
                colWidths=[0.1, 0.2, 0.15, 0.15, 0.2, 0.2])
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1, 2.5)

# Style header row
for i in range(6):
    cell = table[(0, i)]
    cell.set_facecolor('#3498db')
    cell.set_text_props(weight='bold', color='white', fontsize=12)

# Style best model row
for i in range(6):
    cell = table[(1, i)]
    cell.set_facecolor('#2ecc71')
    cell.set_text_props(weight='bold')

# Add horizontal lines
for i in range(9):
    for j in range(6):
        table[(i, j)].set_edgecolor('black')
        table[(i, j)].set_linewidth(1.5)

plt.title('Model Ranking: Complete Results', fontsize=16, fontweight='bold', pad=20)
plt.savefig('plots/model_ranking_table.png', dpi=300, bbox_inches='tight')
print("✓ Saved: plots/model_ranking_table.png")

print("\n✓ All visualizations created successfully!")
print("  - plots/research_report_visualizations.png (6-panel overview)")
print("  - plots/rmse_vs_percentile_clean.png (publication-ready)")
print("  - plots/model_ranking_table.png (formatted table)")
