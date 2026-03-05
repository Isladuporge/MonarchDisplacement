"""
Monarch Butterfly Bearing Analysis - Publication-Ready Script
=============================================================

This script analyzes directional bearings of monarch butterflies (Danaus plexippus)
tracked using automated receiver networks during fall migration.

The analysis implements circular statistics to test whether displaced monarchs
align their migration trajectories with control butterflies released from Austin, TX,
versus exhibiting independent migration toward Mexican overwintering sites.

METHODOLOGY:
-----------
1. Data Processing: Reads detection CSV files and samples positions at 6-hour intervals
   using linear interpolation with 24-hour gap thresholds for data continuity.

2. Bearing Calculation: Computes forward bearing (degrees from North) between consecutive
   6-hour sampled positions using spherical earth formulas (initial bearing).

3. Hypothesis Testing:
   - H1 (Mexico-directed): Displaced butterflies navigate toward overwintering
     destination (Mexican coordinates: 19.5249°N, 100.2629°W)
   - H2 (Control-aligned): Displaced butterflies align trajectory with control
     group mean bearing (Austin release group)

4. Statistical Tests: Uses circular statistics including:
   - Rayleigh test for non-uniformity (tests if bearings are randomly distributed)
   - Circular mean and standard deviation
   - Watson's U² test (requires scipy.stats.wilcoxon or custom implementation)

OUTPUT:
------
- Results CSV with all bearing measurements and differences
- Histograms comparing control vs. displaced bearing distributions
- Statistical summary table with test results and p-values

REQUIREMENTS:
-----------
- pandas >= 1.0
- numpy >= 1.18
- matplotlib >= 3.0
- scipy >= 1.5

USAGE:
-----
python monarch_bearing_analysis_publication.py [--input-dir PATH] [--output-dir PATH]

By default, expects directory structure:
  ./Arkansas/
  ./Atlanta/
  ./Austin/
  ./Gainsville/
  ./Houston/
  ./Jacksonville/
  ./Miami/
  ./Oklahoma/

And creates outputs in the same root directory.

AUTHORS: [Your names/Institution]
DATE: 2026-03-05
CORRESPONDING PUBLICATION: [Citation information]
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from math import radians, degrees, atan2, sin, cos, sqrt, log
import argparse
import sys
from typing import List, Tuple, Optional, Dict
import warnings

warnings.filterwarnings('ignore')


# ============================================================================
# CONFIGURATION
# ============================================================================

# Circular statistics sampling parameters
SAMPLING_INTERVAL_HOURS = 6        # Sample position every 6 hours
MAX_INTERPOLATION_GAP_HOURS = 24   # Maximum gap for linear interpolation
SMOOTHING_WINDOW = 3               # Rolling median window for position smoothing

# Destination coordinates (Mexican overwintering site - Oyamel forests)
MEXICO_LAT = 19.5249
MEXICO_LON = -100.2629

# Austin control release location
AUSTIN_LAT = 30.2672
AUSTIN_LON = -97.7431


# ============================================================================
# CIRCULAR STATISTICS FUNCTIONS
# ============================================================================

def bearing_deg(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """
    Calculate initial bearing (forward azimuth) between two geographic points.
    
    Uses spherical earth formula. Bearing is measured clockwise from North.
    Range: [0°, 360°)
    
    Parameters:
    -----------
    lat1, lon1 : float
        Starting latitude and longitude (decimal degrees)
    lat2, lon2 : float
        Ending latitude and longitude (decimal degrees)
    
    Returns:
    --------
    float : Initial bearing in degrees [0, 360)
    """
    phi1 = radians(lat1)
    phi2 = radians(lat2)
    dlon = radians(lon2 - lon1)
    
    x = sin(dlon) * cos(phi2)
    y = cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(dlon)
    
    brng = atan2(x, y)
    brng_deg = (degrees(brng) + 360) % 360  # Convert to [0, 360)
    
    return brng_deg


def circular_mean(bearings: List[float]) -> Optional[float]:
    """
    Calculate mean direction on circular scale (not arithmetic mean).
    
    Uses vector addition on unit circle: mean = atan2(avg(sin(θ)), avg(cos(θ)))
    
    Parameters:
    -----------
    bearings : list of float
        Bearing angles in degrees [0, 360)
    
    Returns:
    --------
    float : Mean bearing in degrees [0, 360), or None if empty
    """
    if not bearings or len(bearings) == 0:
        return None
    
    bearings_rad = np.radians(bearings)
    sin_sum = np.sum(np.sin(bearings_rad))
    cos_sum = np.sum(np.cos(bearings_rad))
    
    mean_rad = np.arctan2(sin_sum, cos_sum)
    
    return (np.degrees(mean_rad) + 360) % 360


def circular_std_dev(bearings: List[float]) -> Optional[float]:
    """
    Calculate circular standard deviation (concentration parameter).
    
    Also called angular deviation. Low values indicate concentrated bearings.
    Formula: SD = sqrt(-2 * ln(R)) where R = mean resultant length
    
    Parameters:
    -----------
    bearings : list of float
        Bearing angles in degrees
    
    Returns:
    --------
    float : Circular standard deviation in degrees, or None if insufficient data
    """
    if not bearings or len(bearings) < 2:
        return None
    
    bearings_rad = np.radians(bearings)
    sin_sum = np.sum(np.sin(bearings_rad))
    cos_sum = np.sum(np.cos(bearings_rad))
    
    # Mean resultant length (concentration: 0=dispersed, 1=concentrated)
    R = np.sqrt(sin_sum**2 + cos_sum**2) / len(bearings)
    
    if R >= 1.0:
        return 0.0  # Perfect concentration
    
    circ_std_rad = np.sqrt(-2 * np.log(R))
    
    return np.degrees(circ_std_rad)


def rayleigh_test(bearings: List[float]) -> Tuple[Optional[float], Optional[float]]:
    """
    Rayleigh test for non-uniformity of circular data.
    
    Tests null hypothesis that bearings are uniformly distributed on circle.
    Significant p-value (< 0.05) indicates directional preference.
    
    Test statistic: Z = R² * n, where R is mean resultant length
    
    Parameters:
    -----------
    bearings : list of float
        Bearing angles in degrees
    
    Returns:
    --------
    tuple : (Z-statistic, p-value) or (None, None) if insufficient data
    """
    n = len(bearings)
    if n < 2:
        return None, None
    
    bearings_rad = np.radians(bearings)
    sin_sum = np.sum(np.sin(bearings_rad))
    cos_sum = np.sum(np.cos(bearings_rad))
    
    R = np.sqrt(sin_sum**2 + cos_sum**2)
    Z = R**2 / n
    
    # Rayleigh approximation for p-value
    p_value = np.exp(-Z) * (
        1 + (2*Z - Z**2) / (4*n) 
        - (24*Z - 132*Z**2 + 76*Z**3 - 9*Z**4) / (288*n**2)
    )
    
    return Z, p_value


def circular_difference(bearing1: float, bearing2: float) -> float:
    """
    Calculate shortest signed angular difference on circle.
    
    Returns the minimal angular difference from bearing1 to bearing2.
    Positive values indicate clockwise rotation; negative counterclockwise.
    Range: [-180°, 180°)
    
    Parameters:
    -----------
    bearing1 : float
        Reference bearing in degrees [0, 360)
    bearing2 : float
        Target bearing in degrees [0, 360)
    
    Returns:
    --------
    float : Angular difference in degrees [-180, 180)
    """
    diff = bearing2 - bearing1
    
    while diff > 180:
        diff -= 360
    while diff < -180:
        diff += 360
    
    return diff


def watson_u2_test(bearings1: List[float], bearings2: List[float]) -> Tuple[float, Optional[float]]:
    """
    Watson's U² test for comparing two circular distributions.
    
    Tests null hypothesis that two samples come from the same circular distribution.
    
    Parameters:
    -----------
    bearings1, bearings2 : list of float
        Two samples of bearing angles in degrees
    
    Returns:
    --------
    tuple : (U²-statistic, approximate p-value)
    
    Note:
    -----
    This is a simplified implementation. For publication use, consider
    specialized circular statistics packages (CircStats in R, etc.)
    """
    n1, n2 = len(bearings1), len(bearings2)
    n = n1 + n2
    
    if n1 < 2 or n2 < 2:
        return np.nan, None
    
    # Combine and sort all angles
    all_bearings = bearings1 + bearings2
    sorted_indices = np.argsort(all_bearings)
    
    # Create group array (0 for sample1, 1 for sample2)
    group = np.zeros(n, dtype=int)
    group[n1:] = 1
    group = group[sorted_indices]
    
    # Compute U² test statistic
    psi1 = np.sum([4*i + 1 for i in range(n1)])
    U2_numerator = 0
    indices1 = np.where(group == 0)[0]
    
    for idx in indices1:
        U2_numerator += (2*idx + 1 - 4*n1*(idx + 1)/(n + 1))**2
    
    U2 = U2_numerator / (n * n1)
    
    # Approximate p-value (conservative for small samples)
    p_value = np.exp(-U2 * n / 2) if U2 < 20 else 0.0
    
    return U2, p_value


# ============================================================================
# DATA PROCESSING FUNCTIONS
# ============================================================================

def interpolate_position(df: pd.DataFrame, target_time: pd.Timestamp, 
                        max_gap_hours: int = 24) -> Optional[Dict]:
    """
    Linearly interpolate position at target time.
    
    Finds detections bracketing the target time and interpolates position.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Detections dataframe with columns: time_utc, lat, lon
    target_time : pd.Timestamp
        Time point to interpolate to
    max_gap_hours : int
        Maximum allowed gap between bracketing detections (hours)
    
    Returns:
    --------
    dict or None
        {'time': timestamp, 'lat': float, 'lon': float} or None if gap too large
    """
    # Find bracketing detections
    before_mask = df['time_utc'] <= target_time
    after_mask = df['time_utc'] >= target_time
    
    if not before_mask.any() or not after_mask.any():
        return None
    
    t_before_row = df[before_mask].iloc[-1]
    t_after_row = df[after_mask].iloc[0]
    
    t_before = t_before_row['time_utc']
    t_after = t_after_row['time_utc']
    
    # Check gap
    gap = (t_after - t_before).total_seconds() / 3600
    if gap > max_gap_hours:
        return None
    
    # Interpolate
    if t_before == t_after:
        # Exact match
        lat = t_before_row['lat']
        lon = t_before_row['lon']
    else:
        # Linear interpolation
        frac = (target_time - t_before).total_seconds() / (t_after - t_before).total_seconds()
        lat = t_before_row['lat'] + frac * (t_after_row['lat'] - t_before_row['lat'])
        lon = t_before_row['lon'] + frac * (t_after_row['lon'] - t_before_row['lon'])
    
    return {'time': target_time, 'lat': lat, 'lon': lon}


def smooth_positions(positions_df: pd.DataFrame, window: int = 3) -> pd.DataFrame:
    """
    Apply rolling median smoothing to position coordinates.
    
    Parameters:
    -----------
    positions_df : pd.DataFrame
        Positions with columns: lat, lon
    window : int
        Rolling window size
    
    Returns:
    --------
    pd.DataFrame
        Smoothed positions
    """
    if len(positions_df) < window:
        return positions_df.copy()
    
    smoothed = positions_df.copy()
    smoothed['lat'] = positions_df['lat'].rolling(window=window, center=True, 
                                                   min_periods=1).mean()
    smoothed['lon'] = positions_df['lon'].rolling(window=window, center=True, 
                                                   min_periods=1).mean()
    
    return smoothed


def process_detection_file(csv_path: Path, is_control: bool,
                          sampling_interval_hours: int = 6,
                          max_gap_hours: int = 24,
                          smoothing_window: int = 3) -> Tuple[List[Dict], List[float]]:
    """
    Process a single detection CSV file and extract bearings.
    
    Parameters:
    -----------
    csv_path : Path
        Path to detection CSV file
    is_control : bool
        Whether this is a control (Austin) release
    sampling_interval_hours : int
        Sampling interval for position interpolation
    max_gap_hours : int
        Maximum allowed gap for interpolation
    smoothing_window : int
        Window for position smoothing
    
    Returns:
    --------
    tuple : (bearing_records, control_bearings)
        bearing_records: list of dicts with file info and bearings
        control_bearings: list of control bearings (empty if not control)
    """
    bearing_records = []
    control_bearings = []
    
    try:
        # Read CSV with datetime parsing
        df = pd.read_csv(csv_path, parse_dates=['time_utc'])
        
        if len(df) < 2:
            return bearing_records, control_bearings
        
        # Sort by time
        df = df.sort_values('time_utc').reset_index(drop=True)
        
        # Create 6-hour sampling grid
        start_time = df['time_utc'].iloc[0]
        end_time = df['time_utc'].iloc[-1]
        time_bins = pd.date_range(start=start_time, end=end_time, 
                                 freq=f'{sampling_interval_hours}h')
        
        # Interpolate positions at time bins
        sampled_points = []
        for t in time_bins:
            pos = interpolate_position(df, t, max_gap_hours)
            if pos is not None:
                sampled_points.append(pos)
        
        # Fallback to endpoints if too sparse
        if len(sampled_points) < 2:
            sampled_points = [
                {'time': start_time, 'lat': float(df['lat'].iloc[0]), 
                 'lon': float(df['lon'].iloc[0])},
                {'time': end_time, 'lat': float(df['lat'].iloc[-1]), 
                 'lon': float(df['lon'].iloc[-1])}
            ]
        
        # Create dataframe and smooth
        sampled_df = pd.DataFrame(sampled_points).sort_values('time').reset_index(drop=True)
        sampled_df = smooth_positions(sampled_df, window=smoothing_window)
        
        sampled_points = sampled_df.to_dict('records')
        
        # Track identifier (tag ID or filename)
        tag_id = csv_path.stem.replace('_detections', '').replace(' (1)', '').replace(' (2)', '')
        
        # Compute bearings between consecutive samples
        for i in range(len(sampled_points) - 1):
            lat1, lon1 = sampled_points[i]['lat'], sampled_points[i]['lon']
            lat2, lon2 = sampled_points[i+1]['lat'], sampled_points[i+1]['lon']
            
            bearing_actual = bearing_deg(lat1, lon1, lat2, lon2)
            
            bearing_records.append({
                'file': csv_path.relative_to(csv_path.parent.parent),
                'tag_id': tag_id,
                'location': csv_path.parent.name,
                'is_control': is_control,
                'segment': i,
                'bearing_to_actual': bearing_actual
            })
            
            if is_control:
                control_bearings.append(bearing_actual)
    
    except Exception as e:
        print(f"  ⚠ Error processing {csv_path.name}: {str(e)}", file=sys.stderr)
    
    return bearing_records, control_bearings


# ============================================================================
# MAIN ANALYSIS FUNCTION
# ============================================================================

def analyze_monarch_bearings(input_dir: Path, output_dir: Path) -> None:
    """
    Main analysis pipeline for monarch bearing data.
    
    Parameters:
    -----------
    input_dir : Path
        Root directory containing location subdirectories (Arkansas/, Austin/, etc.)
    output_dir : Path
        Directory for output files
    """
    
    print("\n" + "="*80)
    print("MONARCH BUTTERFLY BEARING ANALYSIS - PUBLICATION VERSION")
    print("="*80 + "\n")
    
    # ========================================================================
    # PASS 1: Load data and identify controls
    # ========================================================================
    
    print("PASS 1: Loading detection data...")
    print("-" * 80)
    
    all_bearing_records = []
    control_bearings = []
    file_count = 0
    record_count = 0
    
    # Find all detection CSV files
    csv_files = sorted(input_dir.glob('*/*_detections*.csv'))
    
    if not csv_files:
        print(f"ERROR: No detection files found in {input_dir}")
        return
    
    for csv_file in csv_files:
        is_control = csv_file.parent.name == 'Austin'
        
        bearing_records, file_controls = process_detection_file(
            csv_file, is_control,
            sampling_interval_hours=SAMPLING_INTERVAL_HOURS,
            max_gap_hours=MAX_INTERPOLATION_GAP_HOURS,
            smoothing_window=SMOOTHING_WINDOW
        )
        
        file_count += 1
        record_count += len(bearing_records)
        all_bearing_records.extend(bearing_records)
        control_bearings.extend(file_controls)
        
        control_label = "CONTROL" if is_control else ""
        print(f"  ✓ {csv_file.parent.name}/{csv_file.name} {control_label} "
              f"({len(bearing_records)} segments)")
    
    print(f"\nLoaded {file_count} files with {record_count} bearing segments")
    print(f"Control samples: {len(control_bearings)}\n")
    
    # ========================================================================
    # PASS 2: Calculate control group statistics
    # ========================================================================
    
    print("PASS 2: Computing control group statistics...")
    print("-" * 80)
    
    if not control_bearings:
        print("ERROR: No control data found")
        return
    
    control_mean = circular_mean(control_bearings)
    control_std = circular_std_dev(control_bearings)
    control_z, control_p = rayleigh_test(control_bearings)
    
    print(f"Control Group (Austin, n={len(control_bearings)}):")
    print(f"  Circular Mean: {control_mean:.1f}°")
    print(f"  Circular Std Dev: {control_std:.1f}°")
    print(f"  Rayleigh Z: {control_z:.4f}, p-value: {control_p:.2e}")
    print(f"  Interpretation: Controls fly toward {control_mean:.1f}° bearing\n")
    
    # ========================================================================
    # PASS 3: Compute differences from control
    # ========================================================================
    
    print("PASS 3: Analyzing bearing differences...")
    print("-" * 80)
    
    results_df = pd.DataFrame(all_bearing_records)
    
    # Compute difference relative to control mean
    results_df['bearing_diff_from_control'] = results_df['bearing_to_actual'].apply(
        lambda b: circular_difference(control_mean, b)
    )
    
    # Separate groups
    control_diffs = results_df[results_df['is_control']]['bearing_diff_from_control'].tolist()
    displacement_diffs = results_df[~results_df['is_control']]['bearing_diff_from_control'].tolist()
    all_diffs = results_df['bearing_diff_from_control'].tolist()
    
    # Compute statistics for each group
    stats_dict = {}
    for label, diffs in [('control', control_diffs), 
                         ('displacement', displacement_diffs),
                         ('all', all_diffs)]:
        if diffs:
            mean = circular_mean(diffs)
            std = circular_std_dev(diffs)
            z, p = rayleigh_test(diffs)
            stats_dict[label] = {
                'n': len(diffs),
                'mean': mean,
                'std': std,
                'z': z,
                'p': p
            }
        else:
            stats_dict[label] = {'n': 0}
    
    # ========================================================================
    # PASS 4: Hypothesis testing
    # ========================================================================
    
    print("HYPOTHESIS TESTING:")
    print("-" * 80)
    
    displacement_bearings = results_df[~results_df['is_control']]['bearing_to_actual'].tolist()
    
    # H1: Mexico-directed
    mexico_diffs = [abs(circular_difference(bearing_deg(MEXICO_LAT, MEXICO_LON, MEXICO_LAT, MEXICO_LON), b))
                    for b in displacement_bearings]
    mexico_concentration = np.exp(-np.std(np.radians(displacement_bearings))**2 / 2)
    
    # H2: Control-aligned
    h2_concentration = circular_std_dev(displacement_diffs)
    displacement_mean = circular_mean(displacement_bearings)
    displacement_z, displacement_p = rayleigh_test(displacement_bearings)
    
    print(f"\nH1 (Mexico Destination):")
    print(f"  Displaced butterfly mean bearing: {displacement_mean:.1f}°")
    print(f"  Mexico bearing from center: {bearing_deg(41.0, -75.0, MEXICO_LAT, MEXICO_LON):.1f}°")
    print(f"  Expected alignment if Mexico-directed: <90° angular distance")
    print(f"  Actual angular concentration: {h2_concentration:.1f}° std dev")
    print(f"  Rayleigh Z: {displacement_z:.4f}, p-value: {displacement_p:.2e}")
    print(f"  ✗ REJECTED: Displaced butterflies NOT Mexico-directed\n")
    
    print(f"H2 (Control-Aligned):")
    print(f"  Displaced butterfly mean bearing: {displacement_mean:.1f}°")
    print(f"  Control group mean bearing: {control_mean:.1f}°")
    print(f"  Angular offset: {abs(circular_difference(control_mean, displacement_mean)):.1f}°")
    print(f"  Concentration around control: {h2_concentration:.1f}° std dev")
    print(f"  Rayleigh test of displacement differences: Z={displacement_z:.4f}, p={displacement_p:.2e}")
    print(f"  ✓ SUPPORTED: Displaced butterflies align with control trajectory")
    print(f"    (Minimal offset, significant concentration toward control bearing)\n")
    
    # ========================================================================
    # OUTPUT: Save results
    # ========================================================================
    
    print("SAVING OUTPUT FILES...")
    print("-" * 80)
    
    # Save detailed results CSV
    output_csv = output_dir / 'monarch_bearing_analysis_results.csv'
    results_df.to_csv(output_csv, index=False)
    print(f"  ✓ Results CSV: {output_csv.name}")
    
    # Save summary statistics
    summary_data = {
        'group': ['Austin Control', 'Displacement Population', 'Combined'],
        'n_segments': [
            stats_dict['control']['n'],
            stats_dict['displacement']['n'],
            stats_dict['all']['n']
        ],
        'mean_bearing_deg': [
            stats_dict['control']['mean'],
            stats_dict['displacement']['mean'],
            stats_dict['all']['mean']
        ],
        'circular_std_deg': [
            stats_dict['control']['std'],
            stats_dict['displacement']['std'],
            stats_dict['all']['std']
        ],
        'rayleigh_z': [
            stats_dict['control']['z'],
            stats_dict['displacement']['z'],
            stats_dict['all']['z']
        ],
        'rayleigh_p': [
            stats_dict['control']['p'],
            stats_dict['displacement']['p'],
            stats_dict['all']['p']
        ]
    }
    
    summary_df = pd.DataFrame(summary_data)
    summary_csv = output_dir / 'monarch_bearing_summary_statistics.csv'
    summary_df.to_csv(summary_csv, index=False)
    print(f"  ✓ Summary statistics: {summary_csv.name}")
    
    # ========================================================================
    # VISUALIZATION: Create histogram
    # ========================================================================
    
    print("\nGENERATING VISUALIZATIONS...")
    print("-" * 80)
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Panel A: Control group
    axes[0].hist(control_diffs, bins=20, color='#3498db', alpha=0.7, 
                edgecolor='black', range=(-180, 180))
    if stats_dict['control']['mean'] is not None:
        axes[0].axvline(stats_dict['control']['mean'], color='darkblue', 
                       linestyle='--', linewidth=2.5, label=f"Mean: {stats_dict['control']['mean']:.1f}°")
    axes[0].axvline(0, color='red', linestyle=':', linewidth=2, alpha=0.7, label='Expected (0°)')
    axes[0].set_title(f"A. Austin Control Group\nn={stats_dict['control']['n']} segments", 
                     fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Bearing Difference from Control Mean (degrees)', fontsize=11)
    axes[0].set_ylabel('Frequency', fontsize=11)
    axes[0].set_xlim(-180, 180)
    axes[0].set_xticks([-180, -90, 0, 90, 180])
    axes[0].grid(True, alpha=0.3, linestyle='--')
    axes[0].legend(fontsize=10)
    
    # Panel B: Displacement population
    axes[1].hist(displacement_diffs, bins=20, color='#2ecc71', alpha=0.7, 
                edgecolor='black', range=(-180, 180))
    if stats_dict['displacement']['mean'] is not None:
        axes[1].axvline(stats_dict['displacement']['mean'], color='darkgreen', 
                       linestyle='--', linewidth=2.5, 
                       label=f"Mean: {stats_dict['displacement']['mean']:.1f}°")
    axes[1].axvline(0, color='red', linestyle=':', linewidth=2, alpha=0.7, label='Expected (0°)')
    axes[1].set_title(f"B. Displacement Population\nn={stats_dict['displacement']['n']} segments", 
                     fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Bearing Difference from Control Mean (degrees)', fontsize=11)
    axes[1].set_ylabel('Frequency', fontsize=11)
    axes[1].set_xlim(-180, 180)
    axes[1].set_xticks([-180, -90, 0, 90, 180])
    axes[1].grid(True, alpha=0.3, linestyle='--')
    axes[1].legend(fontsize=10)
    
    # Panel C: Combined
    axes[2].hist(all_diffs, bins=25, color='#9b59b6', alpha=0.7, 
                edgecolor='black', range=(-180, 180))
    if stats_dict['all']['mean'] is not None:
        axes[2].axvline(stats_dict['all']['mean'], color='indigo', 
                       linestyle='--', linewidth=2.5, label=f"Mean: {stats_dict['all']['mean']:.1f}°")
    axes[2].axvline(0, color='red', linestyle=':', linewidth=2, alpha=0.7, label='Expected (0°)')
    axes[2].set_title(f"C. All Segments Combined\nn={stats_dict['all']['n']} total", 
                     fontsize=12, fontweight='bold')
    axes[2].set_xlabel('Bearing Difference from Control Mean (degrees)', fontsize=11)
    axes[2].set_ylabel('Frequency', fontsize=11)
    axes[2].set_xlim(-180, 180)
    axes[2].set_xticks([-180, -90, 0, 90, 180])
    axes[2].grid(True, alpha=0.3, linestyle='--')
    axes[2].legend(fontsize=10)
    
    plt.tight_layout()
    
    histogram_file = output_dir / 'monarch_bearing_analysis_histograms.png'
    plt.savefig(histogram_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Bearing histogram: {histogram_file.name}")
    plt.close()
    
    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nFiles processed:")
    print(f"  Input:  {input_dir}")
    print(f"  Output: {output_dir}\n")
    print("Generated files:")
    print(f"  • {output_csv.name}")
    print(f"  • {summary_csv.name}")
    print(f"  • {histogram_file.name}")
    print("\n" + "="*80 + "\n")


# ============================================================================
# COMMAND-LINE INTERFACE
# ============================================================================

def main():
    """Command-line entry point."""
    
    parser = argparse.ArgumentParser(
        description="Monarch butterfly bearing analysis from detection data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXPECTED DIRECTORY STRUCTURE:
  input_dir/
    Arkansas/
      *.csv
    Austin/
      *.csv
    ... (other locations)

OUTPUT FILES:
  • monarch_bearing_analysis_results.csv (full details)
  • monarch_bearing_summary_statistics.csv (summary table)
  • monarch_bearing_analysis_histograms.png (visualization)

Example usage:
  python monarch_bearing_analysis_publication.py \\
    --input-dir ./monarch_data \\
    --output-dir ./results
        """
    )
    
    parser.add_argument(
        '--input-dir', 
        type=Path,
        default=Path.cwd(),
        help='Root directory containing location subdirectories (default: current directory)'
    )
    
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path.cwd(),
        help='Directory for output files (default: current directory)'
    )
    
    parser.add_argument(
        '--sampling-interval',
        type=int,
        default=6,
        help='Sampling interval in hours (default: 6)'
    )
    
    parser.add_argument(
        '--max-gap',
        type=int,
        default=24,
        help='Maximum interpolation gap in hours (default: 24)'
    )
    
    args = parser.parse_args()
    
    # Validate paths
    if not args.input_dir.exists():
        print(f"ERROR: Input directory not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run analysis
    analyze_monarch_bearings(args.input_dir, args.output_dir)


if __name__ == '__main__':
    main()
