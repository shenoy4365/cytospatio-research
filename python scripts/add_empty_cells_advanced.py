"""
Advanced Empty Cell Insertion via Triangulation

This script implements a sophisticated algorithm for adding empty cells to spatial cell data:
1. Generates Delaunay triangulation of all cell positions
2. Calculates average triangle side length (= average neighbor distance)
3. For each triangle:
   - Checks if area >= beta (minimum threshold)
   - Generates all integer pixel positions inside triangle
   - Removes pixels within average_radius from triangle corners
   - Randomly samples alpha × remaining_pixels positions
   - Adds sampled positions as empty cell type
4. Outputs combined CSV with original cells + empty cells

Parameters:
- alpha: Fraction of pixels to choose (suggested: 1/average_cell_area)
- beta: Minimum triangle area threshold (suggested: average_cell_area)
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
    # Area = 0.5 * |cross product of (p2-p1) and (p3-p1)|
    v1 = p2 - p1
    v2 = p3 - p1
    area = 0.5 * abs(np.cross(v1, v2))

    return area


def calculate_triangle_side_lengths(points, simplex):
    """
    Calculate the three side lengths of a triangle.

    Args:
        points: Array of all point coordinates
        simplex: Indices of the three vertices forming the triangle

    Returns:
        Array of three side lengths
    """
    # Get the three vertices
    p1, p2, p3 = points[simplex]

    # Calculate three side lengths
    side1 = np.linalg.norm(p2 - p1)
    side2 = np.linalg.norm(p3 - p2)
    side3 = np.linalg.norm(p1 - p3)

    return np.array([side1, side2, side3])


def point_in_triangle(p, triangle):
    """
    Check if point p is inside triangle using barycentric coordinates.

    Args:
        p: Point coordinates [x, y]
        triangle: Array of three triangle vertices [[x1,y1], [x2,y2], [x3,y3]]

    Returns:
        True if point is inside triangle, False otherwise
    """
    # Extract triangle vertices
    v0, v1, v2 = triangle

    # Compute vectors
    v0v1 = v1 - v0
    v0v2 = v2 - v0
    v0p = p - v0

    # Compute dot products
    dot00 = np.dot(v0v2, v0v2)
    dot01 = np.dot(v0v2, v0v1)
    dot02 = np.dot(v0v2, v0p)
    dot11 = np.dot(v0v1, v0v1)
    dot12 = np.dot(v0v1, v0p)

    # Compute barycentric coordinates
    inv_denom = 1 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom

    # Check if point is in triangle
    return (u >= 0) and (v >= 0) and (u + v <= 1)


def get_pixels_inside_triangle(points, simplex):
    """
    Get all integer pixel positions inside a triangle.

    Args:
        points: Array of all point coordinates
        simplex: Indices of the three vertices forming the triangle

    Returns:
        Array of pixel positions inside the triangle [[x1,y1], [x2,y2], ...]
    """
    # Get the three vertices
    triangle = points[simplex]

    # Get bounding box
    min_x = int(np.floor(triangle[:, 0].min()))
    max_x = int(np.ceil(triangle[:, 0].max()))
    min_y = int(np.floor(triangle[:, 1].min()))
    max_y = int(np.ceil(triangle[:, 1].max()))

    # Generate all integer coordinates in bounding box
    pixels = []
    for x in range(min_x, max_x + 1):
        for y in range(min_y, max_y + 1):
            p = np.array([x, y], dtype=float)
            if point_in_triangle(p, triangle):
                pixels.append([x, y])

    return np.array(pixels) if pixels else np.array([]).reshape(0, 2)


def remove_pixels_near_corners(pixels, triangle_vertices, radius):
    """
    Remove pixels that are within 'radius' distance from any triangle corner.

    Args:
        pixels: Array of pixel positions
        triangle_vertices: Array of three triangle vertices
        radius: Distance threshold

    Returns:
        Filtered array of pixels
    """
    if len(pixels) == 0:
        return pixels

    # For each corner, calculate distances to all pixels
    mask = np.ones(len(pixels), dtype=bool)

    for corner in triangle_vertices:
        distances = np.linalg.norm(pixels - corner, axis=1)
        mask &= (distances > radius)

    return pixels[mask]


def add_empty_cells(input_csv, alpha, beta, output_csv, verbose=True):
    """
    Main function to add empty cells using triangulation-based algorithm.

    Args:
        input_csv: Path to input CSV file (columns: x, y, marks)
        alpha: Fraction of pixels to choose per triangle
        beta: Minimum triangle area threshold
        output_csv: Path to output CSV file
        verbose: Print progress messages
    """
    # Step 1: Read CSV file
    if verbose:
        print(f"Reading input file: {input_csv}")
    df = pd.read_csv(input_csv)

    # Extract columns
    x = df['x'].values
    y = df['y'].values
    marks = df['marks'].values

    if verbose:
        print(f"  Loaded {len(df)} cells with marks: {sorted(set(marks))}")

    # Step 2: Generate Delaunay triangulation
    if verbose:
        print(f"Generating Delaunay triangulation...")
    points = np.column_stack([x, y])
    tri = Delaunay(points)

    if verbose:
        print(f"  Generated {len(tri.simplices)} triangles")

    # Step 3: Calculate average side length (average neighbor distance)
    if verbose:
        print(f"Calculating average triangle side length...")
    all_side_lengths = []
    for simplex in tri.simplices:
        side_lengths = calculate_triangle_side_lengths(points, simplex)
        all_side_lengths.extend(side_lengths)

    avg_side_length = np.mean(all_side_lengths)
    avg_radius = avg_side_length  # Use as radius for corner removal

    if verbose:
        print(f"  Average side length (neighbor distance): {avg_side_length:.2f} pixels")
        print(f"  Using corner removal radius: {avg_radius:.2f} pixels")

    # Step 4: Calculate triangle areas
    if verbose:
        print(f"Calculating triangle areas...")
    areas = np.array([calculate_triangle_area(points, simplex) for simplex in tri.simplices])

    # Calculate average cell area (for parameter suggestions)
    avg_area = np.mean(areas)

    if verbose:
        print(f"  Average triangle area: {avg_area:.2f} square pixels")
        print(f"  Area range: [{areas.min():.2f}, {areas.max():.2f}]")
        print(f"  Beta threshold: {beta:.2f}")
        print(f"  Alpha sampling fraction: {alpha:.6f}")

    # Step 5: Process each triangle and collect empty cell positions
    if verbose:
        print(f"Processing triangles to generate empty cells...")

    empty_cell_positions = []
    triangles_processed = 0
    triangles_skipped = 0

    for i, simplex in enumerate(tri.simplices):
        area = areas[i]

        # Skip if triangle area is below threshold
        if area < beta:
            triangles_skipped += 1
            continue

        # Get all pixels inside triangle
        pixels = get_pixels_inside_triangle(points, simplex)

        if len(pixels) == 0:
            triangles_skipped += 1
            continue

        # Remove pixels near corners
        triangle_vertices = points[simplex]
        filtered_pixels = remove_pixels_near_corners(pixels, triangle_vertices, avg_radius)

        if len(filtered_pixels) == 0:
            triangles_skipped += 1
            continue

        # Calculate number of pixels to sample
        num_to_sample = int(alpha * len(filtered_pixels))

        if num_to_sample == 0:
            triangles_skipped += 1
            continue

        # Randomly sample pixels
        if num_to_sample >= len(filtered_pixels):
            sampled_pixels = filtered_pixels
        else:
            indices = np.random.choice(len(filtered_pixels), num_to_sample, replace=False)
            sampled_pixels = filtered_pixels[indices]

        empty_cell_positions.extend(sampled_pixels)
        triangles_processed += 1

        # Progress update every 10000 triangles
        if verbose and (i + 1) % 10000 == 0:
            print(f"    Processed {i + 1}/{len(tri.simplices)} triangles...")

    if verbose:
        print(f"  Triangles processed: {triangles_processed}")
        print(f"  Triangles skipped: {triangles_skipped}")
        print(f"  Total empty cells generated: {len(empty_cell_positions)}")

    # Step 6: Create empty cell data
    if len(empty_cell_positions) == 0:
        if verbose:
            print(f"Warning: No empty cells generated. Try adjusting alpha/beta parameters.")
        empty_x = np.array([])
        empty_y = np.array([])
        empty_marks = np.array([])
    else:
        empty_cell_positions = np.array(empty_cell_positions)
        empty_x = empty_cell_positions[:, 0]
        empty_y = empty_cell_positions[:, 1]

        # Assign new mark type (max + 1)
        max_mark = int(marks.max())
        new_mark = max_mark + 1
        empty_marks = np.full(len(empty_cell_positions), new_mark)

        if verbose:
            print(f"  Empty cell type assigned: {new_mark} (max existing mark: {max_mark})")

    # Step 7: Combine original data with empty cells
    if verbose:
        print(f"Combining original cells with empty cells...")

    new_x = np.concatenate([x, empty_x])
    new_y = np.concatenate([y, empty_y])
    new_marks = np.concatenate([marks, empty_marks])

    # Create output dataframe
    output_df = pd.DataFrame({
        'x': new_x,
        'y': new_y,
        'marks': new_marks.astype(int)
    })

    # Step 8: Write output CSV
    if verbose:
        print(f"Writing output file: {output_csv}")
    output_df.to_csv(output_csv, index=False)

    if verbose:
        print(f"  Total points: {len(output_df)} ({len(df)} original + {len(empty_cell_positions)} empty)")
        print(f"  Cell type distribution: {dict(output_df['marks'].value_counts().sort_index())}")
        print("Done!\n")


def main():
    """
    Command-line interface for the script.
    """
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python add_empty_cells_advanced.py <input_csv> <alpha> <beta> [output_csv]")
        print("\nArguments:")
        print("  input_csv: Path to input CSV file (with columns: x, y, marks)")
        print("  alpha: Fraction of pixels to choose per triangle (suggested: 1/avg_cell_area)")
        print("  beta: Minimum triangle area threshold (suggested: avg_cell_area)")
        print("  output_csv: Path to output CSV file (optional, defaults to input_empty.csv)")
        print("\nSuggested Parameters:")
        print("  alpha = 1 / average_cell_area (typically 0.0001 - 0.001)")
        print("  beta = average_cell_area (typically 1000 - 5000)")
        print("\nExamples:")
        print("  python add_empty_cells_advanced.py example/cell_data.csv 0.0005 2000")
        print("  python add_empty_cells_advanced.py example/cell_data.csv 0.001 1500 example/cell_data_with_empties.csv")
        sys.exit(1)

    input_csv = sys.argv[1]
    alpha = float(sys.argv[2])
    beta = float(sys.argv[3])

    # Default output filename
    if len(sys.argv) == 5:
        output_csv = sys.argv[4]
    else:
        # Generate default output name
        import os
        base = os.path.splitext(input_csv)[0]
        output_csv = f"{base}_empty_alpha{alpha}_beta{beta}.csv"

    # Validate parameters
    if alpha <= 0 or alpha > 1:
        print(f"Error: alpha must be between 0 and 1 (got {alpha})")
        sys.exit(1)

    if beta <= 0:
        print(f"Error: beta must be positive (got {beta})")
        sys.exit(1)

    # Set random seed for reproducibility
    np.random.seed(42)

    add_empty_cells(input_csv, alpha, beta, output_csv, verbose=True)


if __name__ == "__main__":
    main()