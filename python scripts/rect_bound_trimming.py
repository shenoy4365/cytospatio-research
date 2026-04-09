import numpy as np

def divide_into_n_regions(pts, n_regions):
    """
    Divide point cloud into n equal regions (assumes square grid).
    n_regions should be a perfect square (4, 9, 16, etc.)
    """
    grid_size = int(np.sqrt(n_regions))

    x_min, x_max = pts['x'].min(), pts['x'].max()
    y_min, y_max = pts['y'].min(), pts['y'].max()

    x_bins = np.linspace(x_min, x_max, grid_size + 1)
    y_bins = np.linspace(y_min, y_max, grid_size + 1)

    regions = []
    for i in range(grid_size):
        for j in range(grid_size):
            mask = ((pts['x'] >= x_bins[i]) & (pts['x'] < x_bins[i+1]) &
                    (pts['y'] >= y_bins[j]) & (pts['y'] < y_bins[j+1]))
            regions.append(pts[mask])

    return regions

# Usage: divide into 4 equal regions
regions = divide_into_n_regions(df, n_regions=4)