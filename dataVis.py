import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from pathlib import Path
import glob
import os
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap

# Set style for beautiful plots
plt.style.use('dark_background')
sns.set_palette("husl")

def setup_output_directories():
    """Create output directories for different plot types."""
    base_dir = Path("plots")
    dirs = ["3d_snapshots", "2d_projections", "velocity_analysis", "animations", "summary"]
    
    for dir_name in dirs:
        (base_dir / dir_name).mkdir(parents=True, exist_ok=True)
    
    return base_dir

def load_particle_data(data_dir="data"):
    """Load all particle data files and return sorted list of dataframes."""
    files = sorted(glob.glob(f"{data_dir}/particles_*.txt"))
    data_frames = []
    timesteps = []
    
    for file in files:
        # Extract timestep from filename
        timestep = int(Path(file).stem.split('_')[1])
        timesteps.append(timestep)
        
        # Load data
        df = pd.read_csv(file)
        df['timestep'] = timestep
        data_frames.append(df)
    
    return data_frames, timesteps

def create_custom_colormap():
    """Create beautiful cosmic-inspired colormap."""
    # Inspired by stellar/nebula imagery: dark -> red -> orange -> yellow -> white
    colors = ['#0a0a0a', '#1a0d0d', '#4d1a1a', '#800d00', '#b8001a', '#d43300', '#ff4500', '#ff6b35', '#ff8c42', '#ffb347', '#ffd700', '#fff8dc']
    return LinearSegmentedColormap.from_list("cosmic_cmap", colors)

def plot_3d_snapshot(df, timestep, output_dir, cmap):
    """Create beautiful 3D scatter plot for a single timestep."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Calculate velocity magnitude for coloring
    velocity_mag = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)
    
    # Create 3D scatter plot with cosmic styling
    scatter = ax.scatter(df['x'], df['y'], df['z'], 
                        c=velocity_mag, cmap=cmap, 
                        s=80, alpha=0.9, edgecolors='none')
    
    # Add glow effect by plotting larger, more transparent points behind
    ax.scatter(df['x'], df['y'], df['z'], 
              c=velocity_mag, cmap=cmap, 
              s=200, alpha=0.3, edgecolors='none')
    
    # Customize the plot with cosmic colors
    ax.set_xlabel('X Position', fontsize=12, color='#ffd700')
    ax.set_ylabel('Y Position', fontsize=12, color='#ffd700')
    ax.set_zlabel('Z Position', fontsize=12, color='#ffd700')
    ax.set_title(f'Particle Distribution (t={timestep:04d})', 
                fontsize=16, color='#ffb347', pad=20)
    
    # Add colorbar with cosmic styling
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.1)
    cbar.set_label('Velocity Magnitude', color='#ffd700', fontsize=12)
    cbar.ax.yaxis.set_tick_params(color='#ffd700')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='#ffd700')
    
    # Set background color
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "3d_snapshots" / f"particles_3d_{timestep:04d}.png", 
                dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()

def plot_2d_projections(df, timestep, output_dir, cmap):
    """Create 2D projection plots (XY, XZ, YZ)."""
    velocity_mag = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor('black')
    
    projections = [
        ('x', 'y', 'XY Projection'),
        ('x', 'z', 'XZ Projection'),
        ('y', 'z', 'YZ Projection')
    ]
    
    for i, (x_col, y_col, title) in enumerate(projections):
        # Main scatter plot
        scatter = axes[i].scatter(df[x_col], df[y_col], 
                                c=velocity_mag, cmap=cmap, 
                                s=60, alpha=0.9, edgecolors='none')
        
        # Add glow effect
        axes[i].scatter(df[x_col], df[y_col], 
                       c=velocity_mag, cmap=cmap, 
                       s=150, alpha=0.3, edgecolors='none')
        
        axes[i].set_xlabel(f'{x_col.upper()} Position', color='#ffd700', fontsize=12)
        axes[i].set_ylabel(f'{y_col.upper()} Position', color='#ffd700', fontsize=12)
        axes[i].set_title(f'{title} (t={timestep:04d})', color='#ffb347', fontsize=14)
        axes[i].grid(True, alpha=0.2, color='#4d1a1a')
        axes[i].set_facecolor('black')
        
        # Color tick labels
        axes[i].tick_params(colors='#ffd700')
    
    # Add colorbar to the last subplot with cosmic styling
    cbar = plt.colorbar(scatter, ax=axes[2], pad=0.1)
    cbar.set_label('Velocity Magnitude', color='#ffd700', fontsize=12)
    cbar.ax.yaxis.set_tick_params(color='#ffd700')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='#ffd700')
    
    plt.tight_layout()
    plt.savefig(output_dir / "2d_projections" / f"projections_{timestep:04d}.png", 
                dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()

def plot_velocity_analysis(df, timestep, output_dir):
    """Create velocity distribution and analysis plots."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.patch.set_facecolor('black')
    
    # Velocity magnitude histogram with cosmic colors
    velocity_mag = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)
    axes[0, 0].hist(velocity_mag, bins=20, alpha=0.8, color='#ff6b35', edgecolor='#ffd700', linewidth=1)
    axes[0, 0].set_xlabel('Velocity Magnitude', color='#ffd700')
    axes[0, 0].set_ylabel('Count', color='#ffd700')
    axes[0, 0].set_title(f'Velocity Distribution (t={timestep:04d})', color='#ffb347')
    axes[0, 0].grid(True, alpha=0.2, color='#4d1a1a')
    axes[0, 0].set_facecolor('black')
    
    # Velocity components with fiery colors
    vel_components = ['vx', 'vy', 'vz']
    colors = ['#ff4500', '#ff6b35', '#ff8c42']
    for i, (component, color) in enumerate(zip(vel_components, colors)):
        axes[0, 1].hist(df[component], bins=15, alpha=0.7, 
                       label=component, color=color, edgecolor='#ffd700', linewidth=0.5)
    axes[0, 1].set_xlabel('Velocity Component', color='#ffd700')
    axes[0, 1].set_ylabel('Count', color='#ffd700')
    axes[0, 1].set_title('Velocity Components', color='#ffb347')
    legend = axes[0, 1].legend()
    legend.get_frame().set_facecolor('black')
    legend.get_frame().set_edgecolor('#ffd700')
    for text in legend.get_texts():
        text.set_color('#ffd700')
    axes[0, 1].grid(True, alpha=0.2, color='#4d1a1a')
    axes[0, 1].set_facecolor('black')
    
    # Speed vs position scatter with glow
    axes[1, 0].scatter(df['x'], velocity_mag, alpha=0.8, c='#ff6b35', s=40, edgecolors='none')
    axes[1, 0].scatter(df['x'], velocity_mag, alpha=0.3, c='#ff8c42', s=100, edgecolors='none')
    axes[1, 0].set_xlabel('X Position', color='#ffd700')
    axes[1, 0].set_ylabel('Speed', color='#ffd700')
    axes[1, 0].set_title('Speed vs X Position', color='#ffb347')
    axes[1, 0].grid(True, alpha=0.2, color='#4d1a1a')
    axes[1, 0].set_facecolor('black')
    
    # Kinetic energy with stellar colors
    kinetic_energy = 0.5 * df['mass'] * velocity_mag**2
    axes[1, 1].scatter(range(len(df)), kinetic_energy, alpha=0.8, c='#ffd700', s=40, edgecolors='none')
    axes[1, 1].scatter(range(len(df)), kinetic_energy, alpha=0.3, c='#ffb347', s=100, edgecolors='none')
    axes[1, 1].set_xlabel('Particle Index', color='#ffd700')
    axes[1, 1].set_ylabel('Kinetic Energy', color='#ffd700')
    axes[1, 1].set_title('Kinetic Energy Distribution', color='#ffb347')
    axes[1, 1].grid(True, alpha=0.2, color='#4d1a1a')
    axes[1, 1].set_facecolor('black')
    
    # Style all axes
    for ax in axes.flat:
        ax.tick_params(colors='#ffd700')
    
    plt.tight_layout()
    plt.savefig(output_dir / "velocity_analysis" / f"velocity_analysis_{timestep:04d}.png", 
                dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()

def create_evolution_plots(data_frames, timesteps, output_dir):
    """Create plots showing evolution over time."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.patch.set_facecolor('black')
    
    # Calculate metrics over time
    total_kinetic_energy = []
    avg_speed = []
    position_variance = []
    particle_count = []
    
    for df in data_frames:
        velocity_mag = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)
        ke = 0.5 * df['mass'] * velocity_mag**2
        
        total_kinetic_energy.append(ke.sum())
        avg_speed.append(velocity_mag.mean())
        position_variance.append(df[['x', 'y', 'z']].var().mean())
        particle_count.append(len(df))
    
    # Plot evolution with cosmic/stellar colors
    axes[0, 0].plot(timesteps, total_kinetic_energy, 'o-', color='#ff6b35', linewidth=3, markersize=8, markerfacecolor='#ffd700', markeredgecolor='#ff4500', markeredgewidth=2)
    axes[0, 0].set_xlabel('Time Step', color='#ffd700')
    axes[0, 0].set_ylabel('Total Kinetic Energy', color='#ffd700')
    axes[0, 0].set_title('Energy Evolution', color='#ffb347')
    axes[0, 0].grid(True, alpha=0.2, color='#4d1a1a')
    axes[0, 0].set_facecolor('black')
    
    axes[0, 1].plot(timesteps, avg_speed, 'o-', color='#ff8c42', linewidth=3, markersize=8, markerfacecolor='#ffd700', markeredgecolor='#ff6b35', markeredgewidth=2)
    axes[0, 1].set_xlabel('Time Step', color='#ffd700')
    axes[0, 1].set_ylabel('Average Speed', color='#ffd700')
    axes[0, 1].set_title('Speed Evolution', color='#ffb347')
    axes[0, 1].grid(True, alpha=0.2, color='#4d1a1a')
    axes[0, 1].set_facecolor('black')
    
    axes[1, 0].plot(timesteps, position_variance, 'o-', color='#ffb347', linewidth=3, markersize=8, markerfacecolor='#ffd700', markeredgecolor='#ff8c42', markeredgewidth=2)
    axes[1, 0].set_xlabel('Time Step', color='#ffd700')
    axes[1, 0].set_ylabel('Position Variance', color='#ffd700')
    axes[1, 0].set_title('Spatial Distribution Evolution', color='#ffb347')
    axes[1, 0].grid(True, alpha=0.2, color='#4d1a1a')
    axes[1, 0].set_facecolor('black')
    
    axes[1, 1].plot(timesteps, particle_count, 'o-', color='#ffd700', linewidth=3, markersize=8, markerfacecolor='#fff8dc', markeredgecolor='#ffb347', markeredgewidth=2)
    axes[1, 1].set_xlabel('Time Step', color='#ffd700')
    axes[1, 1].set_ylabel('Particle Count', color='#ffd700')
    axes[1, 1].set_title('Particle Count Evolution', color='#ffb347')
    axes[1, 1].grid(True, alpha=0.2, color='#4d1a1a')
    axes[1, 1].set_facecolor('black')
    
    # Style all axes
    for ax in axes.flat:
        ax.tick_params(colors='#ffd700')
    
    plt.tight_layout()
    plt.savefig(output_dir / "summary" / "evolution_analysis.png", 
                dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()

def create_summary_dashboard(data_frames, timesteps, output_dir):
    """Create a comprehensive summary dashboard."""
    fig = plt.figure(figsize=(20, 16))
    fig.patch.set_facecolor('black')
    
    # Use first and last timesteps for comparison
    df_first = data_frames[0]
    df_last = data_frames[-1]
    
    # 3D comparison
    ax1 = plt.subplot(3, 3, 1, projection='3d')
    velocity_mag = np.sqrt(df_first['vx']**2 + df_first['vy']**2 + df_first['vz']**2)
    ax1.scatter(df_first['x'], df_first['y'], df_first['z'], 
               c=velocity_mag, s=30, alpha=0.7)
    ax1.set_title(f'Initial State (t={timesteps[0]:04d})', color='white')
    
    ax2 = plt.subplot(3, 3, 2, projection='3d')
    velocity_mag = np.sqrt(df_last['vx']**2 + df_last['vy']**2 + df_last['vz']**2)
    ax2.scatter(df_last['x'], df_last['y'], df_last['z'], 
               c=velocity_mag, s=30, alpha=0.7)
    ax2.set_title(f'Final State (t={timesteps[-1]:04d})', color='white')
    
    # Evolution plots
    ax3 = plt.subplot(3, 3, 3)
    total_ke = []
    for df in data_frames:
        v_mag = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)
        ke = 0.5 * df['mass'] * v_mag**2
        total_ke.append(ke.sum())
    ax3.plot(timesteps, total_ke, 'o-', color='cyan')
    ax3.set_title('Total Kinetic Energy', color='white')
    ax3.set_facecolor('black')
    ax3.grid(True, alpha=0.3)
    
    # Add more summary plots...
    # Position histograms
    ax4 = plt.subplot(3, 3, 4)
    ax4.hist(df_first['x'], alpha=0.7, label='Initial', color='blue', edgecolor='white')
    ax4.hist(df_last['x'], alpha=0.7, label='Final', color='red', edgecolor='white')
    ax4.set_title('X Position Distribution', color='white')
    ax4.legend()
    ax4.set_facecolor('black')
    ax4.grid(True, alpha=0.3)
    
    # Velocity distribution comparison
    ax5 = plt.subplot(3, 3, 5)
    v_mag_first = np.sqrt(df_first['vx']**2 + df_first['vy']**2 + df_first['vz']**2)
    v_mag_last = np.sqrt(df_last['vx']**2 + df_last['vy']**2 + df_last['vz']**2)
    ax5.hist(v_mag_first, alpha=0.7, label='Initial', color='blue', edgecolor='white')
    ax5.hist(v_mag_last, alpha=0.7, label='Final', color='red', edgecolor='white')
    ax5.set_title('Speed Distribution', color='white')
    ax5.legend()
    ax5.set_facecolor('black')
    ax5.grid(True, alpha=0.3)
    
    # Phase space plot
    ax6 = plt.subplot(3, 3, 6)
    ax6.scatter(df_last['x'], df_last['vx'], alpha=0.6, c='orange', s=20)
    ax6.set_xlabel('X Position', color='white')
    ax6.set_ylabel('X Velocity', color='white')
    ax6.set_title('Phase Space (X)', color='white')
    ax6.set_facecolor('black')
    ax6.grid(True, alpha=0.3)
    
    # Style all axes
    for ax in fig.get_axes():
        if hasattr(ax, 'tick_params'):
            ax.tick_params(colors='white')
        if hasattr(ax, 'set_xlabel'):
            ax.xaxis.label.set_color('white')
        if hasattr(ax, 'set_ylabel'):
            ax.yaxis.label.set_color('white')
    
    plt.suptitle('Particle Simulation Analysis Dashboard', 
                color='white', fontsize=20, y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / "summary" / "dashboard.png", 
                dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()

def main():
    """Main function to generate all plots."""
    print("🎨 Starting particle data visualization...")
    
    # Setup
    output_dir = setup_output_directories()
    cmap = create_custom_colormap()
    
    # Load data
    print("📊 Loading particle data...")
    data_frames, timesteps = load_particle_data()
    print(f"   Found {len(data_frames)} timesteps with {len(data_frames[0])} particles each")
    
    # Generate plots for each timestep
    print("🎯 Generating individual timestep plots...")
    for i, (df, timestep) in enumerate(zip(data_frames, timesteps)):
        print(f"   Processing timestep {timestep} ({i+1}/{len(data_frames)})")
        
        # 3D snapshot
        plot_3d_snapshot(df, timestep, output_dir, cmap)
        
        # 2D projections
        plot_2d_projections(df, timestep, output_dir, cmap)
        
        # Velocity analysis
        plot_velocity_analysis(df, timestep, output_dir)
    
    # Generate evolution plots
    print("📈 Creating evolution analysis...")
    create_evolution_plots(data_frames, timesteps, output_dir)
    
    # Generate summary dashboard
    print("📋 Creating summary dashboard...")
    create_summary_dashboard(data_frames, timesteps, output_dir)
    
    print(f"✅ All plots saved to '{output_dir}' directory!")
    print("\nGenerated plot categories:")
    print("  📦 3d_snapshots/     - 3D particle distributions")
    print("  📦 2d_projections/   - XY, XZ, YZ projections")
    print("  📦 velocity_analysis/ - Speed and energy analysis")
    print("  📦 summary/          - Evolution plots and dashboard")

if __name__ == "__main__":
    main()
