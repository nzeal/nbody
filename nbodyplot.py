import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import glob
import os
import argparse

# --------------------------
# Parse command-line arguments
# --------------------------
parser = argparse.ArgumentParser(description='Generate 3D plots from N-body simulation data.')
parser.add_argument('--indir', type=str, default='./data', help='Input directory containing particles_*.txt files (default: ./data)')
parser.add_argument('--outdir', type=str, default='./plots', help='Output directory for PNG files (default: ./plots)')
args = parser.parse_args()

data_dir = os.path.abspath(args.indir)
output_dir = os.path.abspath(args.outdir)
file_pattern = 'particles_*.txt'
file_list = sorted(glob.glob(os.path.join(data_dir, file_pattern)))

if not os.path.exists(data_dir):
    raise FileNotFoundError(f"Input directory does not exist: {data_dir}")
if not file_list:
    raise FileNotFoundError(f"No matching files found in {data_dir}")

os.makedirs(output_dir, exist_ok=True)

print(f"Reading data from: {data_dir}")
print(f"Saving plots to: {output_dir}")

# --------------------------
# Read data files
# --------------------------
expected_columns = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'mass', 'fx', 'fy', 'fz']
data_frames = []

for file in file_list:
    try:
        df = pd.read_csv(file)
        if not all(col in df.columns for col in expected_columns):
            df = pd.read_csv(file, delim_whitespace=True)
        if all(col in df.columns for col in expected_columns):
            data_frames.append(df)
        else:
            print(f"Skipping file with missing columns: {file}")
    except Exception as e:
        print(f"Error reading file {file}: {e}")

if not data_frames:
    raise ValueError("No valid data frames could be read.")

# --------------------------
# Determine axis limits
# --------------------------
def add_padding(arr, padding=0.1):
    min_val, max_val = min(arr), max(arr)
    range_val = max_val - min_val if max_val > min_val else 1.0
    return min_val - padding * range_val, max_val + padding * range_val

all_x = np.concatenate([df['x'].values for df in data_frames])
all_y = np.concatenate([df['y'].values for df in data_frames])
all_z = np.concatenate([df['z'].values for df in data_frames])
xlim = add_padding(all_x)
ylim = add_padding(all_y)
zlim = add_padding(all_z)

# --------------------------
# Plot and save each frame
# --------------------------
for frame_idx, df in enumerate(data_frames):
    fig = plt.figure(figsize=(18, 5))

    # Position plot
    ax1 = fig.add_subplot(131, projection='3d')
    ax1.scatter(df['x'], df['y'], df['z'], c=df['mass'], cmap='viridis', s=50)
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.set_zlim(zlim)
    ax1.set_title(f'Particle Positions (Frame {frame_idx:04d})')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')

    # Velocity plot
    ax2 = fig.add_subplot(132, projection='3d')
    ax2.scatter(df['x'], df['y'], df['z'], c='blue', alpha=0.5, s=20)
    ax2.quiver(df['x'], df['y'], df['z'], df['vx'], df['vy'], df['vz'], length=0.5, normalize=True, color='red')
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    ax2.set_zlim(zlim)
    ax2.set_title('Velocity Field')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')

    # Force plot
    ax3 = fig.add_subplot(133, projection='3d')
    ax3.scatter(df['x'], df['y'], df['z'], c='blue', alpha=0.5, s=20)
    ax3.quiver(df['x'], df['y'], df['z'], df['fx'], df['fy'], df['fz'], length=0.5, normalize=True, color='green')
    ax3.set_xlim(xlim)
    ax3.set_ylim(ylim)
    ax3.set_zlim(zlim)
    ax3.set_title('Force Field')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('Z')

    plt.tight_layout()
    output_file = os.path.join(output_dir, f'frame_{frame_idx:04d}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close(fig)

# --------------------------
# Summary
# --------------------------
print(f"Processed {len(file_list)} files.")
print(f"Saved {len(data_frames)} plots to '{output_dir}'.")

