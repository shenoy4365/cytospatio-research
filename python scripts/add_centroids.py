#!/usr/bin/env python3
"""
Delaunay Triangulation Centroid Insertion Script

This script:
1. Reads a CSV file with columns: x, y, marks
2. Generates Delaunay triangulation of x,y points
3. Filters triangles larger than a percentile threshold
4. Calculates centroids of selected triangles
5. Appends centroids with new cell type (max_marks + 1)
6. Writes output CSV file
"""

import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
import sys

def calculate_triangle_area(points, simplex):
    """
    Calculate the area of a triangle given three points.

    Args:
        points: Array of all point coordinates
        simplex: Indices of the three vertices forming the triangle

    Returns:
        Area of the triangle
    """
    # Get the three vertices
    p1, p2, p3 = points[simplex]

    # Calculate area using cross product formula
    # Area = 0.5 * |cross product of (p2-p1) and (p3-p1)| from calculus
    v1 = p2 - p1
    v2 = p3 - p1
    area = 0.5 * abs(np.cross(v1, v2))

    return area


def calculate_centroid(points, simplex):
    """
    Calculate the centroid of a triangle.

    Args:
        points: Array of all point coordinates
        simplex: Indices of the three vertices forming the triangle

    Returns:
        Centroid coordinates (x, y)
    """
    # Centroid is the average of the three vertices
    triangle_points = points[simplex]
    centroid = np.mean(triangle_points, axis=0)

    return centroid


def add_centroids_to_data(input_csv, threshold_percentile, output_csv):
    """
    Main function to process CSV data and add centroids.

    Args:
        input_csv: Path to input CSV file
        threshold_percentile: Percentile threshold (0-100) for filtering triangles
        output_csv: Path to output CSV file
    """
    # Step 1: Read CSV file
    print(f"Reading input file: {input_csv}")
    df = pd.read_csv(input_csv)

    # Extract columns
    x = df['x'].values
    y = df['y'].values
    marks = df['marks'].values

    print(f"  Loaded {len(df)} cells with marks: {sorted(set(marks))}")

    # Step 2: Generate Delaunay triangulation
    print(f"Generating Delaunay triangulation...")
    points = np.column_stack([x, y])
    tri = Delaunay(points)

    print(f"  Generated {len(tri.simplices)} triangles")

    # Step 3: Calculate areas and filter by threshold
    print(f"Calculating triangle areas...")
    areas = np.array([calculate_triangle_area(points, simplex) for simplex in tri.simplices])

    # Calculate the threshold area based on percentile
    threshold_area = np.percentile(areas, threshold_percentile)
    print(f"  Area range: [{areas.min():.2f}, {areas.max():.2f}]")
    print(f"  {threshold_percentile}th percentile threshold: {threshold_area:.2f}")

    # Select triangles larger than threshold
    selected_indices = areas > threshold_area
    selected_simplices = tri.simplices[selected_indices]
    selected_areas = areas[selected_indices]

    print(f"  Selected {len(selected_simplices)} triangles (>{threshold_percentile}th percentile)")

    # Step 4: Calculate centroids of selected triangles
    print(f"Calculating centroids...")
    centroids = np.array([calculate_centroid(points, simplex) for simplex in selected_simplices])

    # Step 5: Determine new cell type (max + 1)
    max_mark = int(marks.max())
    new_mark = max_mark + 1
    print(f"  New cell type for centroids: {new_mark} (max existing mark: {max_mark})")

    # Create array of marks for centroids
    centroid_marks = np.full(len(centroids), new_mark)

    # Step 6: Combine original data with centroids
    print(f"Appending centroids to data...")
    new_x = np.concatenate([x, centroids[:, 0]])
    new_y = np.concatenate([y, centroids[:, 1]])
    new_marks = np.concatenate([marks, centroid_marks])

    # Create output dataframe
    output_df = pd.DataFrame({
        'x': new_x,
        'y': new_y,
        'marks': new_marks.astype(int)
    })

    # Step 7: Write output CSV
    print(f"Writing output file: {output_csv}")
    output_df.to_csv(output_csv, index=False)

    print(f"  Total points: {len(output_df)} ({len(df)} original + {len(centroids)} centroids)")
    print(f"  Cell type distribution: {dict(output_df['marks'].value_counts().sort_index())}")
    print("Done!\n")


def main():
    """
    Command-line interface for the script.
    """
    if len(sys.argv) != 4:
        print("Usage: python add_centroids.py <input_csv> <threshold_percentile> <output_csv>")
        print("  input_csv: Path to input CSV file (with columns: x, y, marks)")
        print("  threshold_percentile: Percentile threshold (0-100) for filtering triangles")
        print("  output_csv: Path to output CSV file")
        print("\nExample:")
        print("  python add_centroids.py example/cell_data.csv 90 example/cell_data_percentile_90.csv")
        sys.exit(1)

    input_csv = sys.argv[1]
    threshold_percentile = float(sys.argv[2])
    output_csv = sys.argv[3]

    # Validate threshold
    if not (0 <= threshold_percentile <= 100):
        print(f"Error: threshold_percentile must be between 0 and 100 (got {threshold_percentile})")
        sys.exit(1)

    add_centroids_to_data(input_csv, threshold_percentile, output_csv)


if __name__ == "__main__":
    main()