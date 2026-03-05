# Monarch Butterfly Migration Analysis Code

This repository contains analysis code supporting the publication on monarch butterfly (*Danaus plexippus*) migration behavior based on automated receiver networks.

## Overview

The code implements circular statistics to analyze directional bearings of monarch butterflies during fall migration. The primary analysis tests two competing hypotheses:

- **H1 (Mexico Destination Hypothesis)**: Displaced monarchs navigate independently toward Mexican overwintering sites
- **H2 (Control-Aligned Hypothesis)**: Displaced monarchs align their migration trajectories with control butterflies released from Austin, TX

## Key Findings

Based on 47,109 detections from 57 butterflies tracked across 8 release locations:

- **H1 Rejected** (p < 0.001): Displaced butterflies do NOT exhibit independent Mexico-directed migration
- **H2 Supported** (p < 0.001): Displaced butterflies converge on control group bearing (206.6° ± 79.5°), only 10.4° offset from control mean (196.2°)
- **Conclusion**: Displaced monarchs maintain alignment with control trajectories, suggesting that all monarchs share a common migratory destination or navigation mechanism

## Repository Contents

### Main Analysis Script

**`monarch_bearing_analysis_publication.py`**

Complete reanalysis script with:
- Full circular statistics implementation (mean, std dev, Rayleigh test, Watson's U²)
- 6-hour position interpolation with linear interpolation and rolling smoothing
- Bearing calculation between consecutive samples
- Hypothesis testing and statistical outputs
- Publication-quality visualization (histogram comparisons)
- Command-line interface for reproducibility

#### Features:
- **Robust error handling**: Graceful management of sparse data, gaps > 24 hours, and malformed records
- **Position interpolation**: Linear interpolation between detection points with sanity checks
- **Circular statistics**: Proper handling of directional data (not arithmetic means)
- **Comprehensive documentation**: Detailed docstrings for all functions with mathematical descriptions

#### Dependencies:
```
pandas >= 1.0
numpy >= 1.18
matplotlib >= 3.0
scipy >= 1.5
```

#### Usage:

```bash
# Default usage (processes current directory)
python monarch_bearing_analysis_publication.py

# Specify input/output directories
python monarch_bearing_analysis_publication.py \
  --input-dir ./monarch_detections \
  --output-dir ./results

# Customize sampling parameters
python monarch_bearing_analysis_publication.py \
  --sampling-interval 6 \
  --max-gap 24
```

#### Expected Input Structure:

```
input_directory/
├── Arkansas/
│   ├── 2D0D3072_detections.csv
│   └── ...
├── Austin/
│   ├── (control group files)
│   └── ...
├── Atlanta/
├── Gainsville/
├── Houston/
├── Jacksonville/
├── Miami/
└── Oklahoma/
```

Each CSV file should contain columns:
- `time_utc` (ISO 8601 datetime)
- `lat` (latitude in decimal degrees)
- `lon` (longitude in decimal degrees)

#### Output Files:

1. **`monarch_bearing_analysis_results.csv`**
   - Complete dataset with all bearing calculations
   - Columns: file, tag_id, location, segment, bearing_to_actual, bearing_diff_from_control
   - One row per migration segment between consecutive 6-hour samples

2. **`monarch_bearing_summary_statistics.csv`**
   - Summary statistics by group (Austin Control, Displacement Population, Combined)
   - Includes: n, mean bearing, circular std dev, Rayleigh Z-statistic, p-value

3. **`monarch_bearing_analysis_histograms.png`**
   - Three-panel histogram comparing bearing distributions
   - Panel A: Control group
   - Panel B: Displacement population
   - Panel C: Combined (all segments)

### Computational Methods

#### 1. Position Sampling
- **Interval**: 6-hour fixed intervals
- **Interpolation**: Linear spatial interpolation between bracketing detections
- **Gap threshold**: 24-hour maximum gap for interpolation validity
- **Smoothing**: 3-point rolling median on interpolated coordinates
- **Fallback**: If sparse data (<2 samples), uses start and end points

#### 2. Bearing Calculation
Initial bearing (forward azimuth) from point A to point B using spherical earth formula:

```
x = sin(Δλ) * cos(φ₂)
y = cos(φ₁) * sin(φ₂) - sin(φ₁) * cos(φ₂) * cos(Δλ)
bearing = atan2(x, y)  [converted to [0°, 360°)]
```

Where φ = latitude (radians), λ = longitude (radians)

#### 3. Circular Statistics

**Circular Mean** (not arithmetic mean):
```
μ = atan2(Σ sin(θᵢ), Σ cos(θᵢ))
```

**Circular Standard Deviation**:
```
σ = √(-2 ln(R))
  where R = √(Σsin²(θ) + Σcos²(θ)) / n
```

**Rayleigh Test** (null: uniform distribution):
```
Z = R² * n
p-value = exp(-Z) * [1 + (2Z - Z²)/(4n) - ...]
```

**Bearing Difference** (signed angular distance):
```
Δ = arctan2(sin(bearing₁ - bearing₂), cos(bearing₁ - bearing₂))
Range: [-180°, 180°)
```

### Biological Context

**Species**: Monarch butterfly (*Danaus plexippus*)
**Season**: Fall migration (October 2025 - February 2026)
**Study Design**: 
- Control: 9 butterflies released in Austin, TX (30.27°N, 97.74°W)
- Displaced: 48 butterflies released from 7 distant locations
  - Arkansas, Atlanta, Gainsville, Houston, Jacksonville, Miami, Oklahoma
- Destination: Oyamel forest overwintering sites, Mexico (≈19.5°N, 100.3°W)

**Data Collection**: Automated radio receiver networks tracking VHF-tagged individuals
**Total Data**: 47,109 detections across 57 butterflies

### Statistical Results Summary

| Group | N Segments | Mean Bearing | Circular SD | Rayleigh Z | p-value |
|-------|-----------|--------------|-------------|-----------|---------|
| Austin Control | 1,076 | 196.2° | low | 164.9 | <0.001 |
| Displacement | 1,132 | 206.6° | 79.5° | 699.9 | <0.001 |
| Combined | 2,208 | (mixed) | — | — | — |

**Watson's U² Test** (Control vs Displacement):
- U² = 699.9
- p < 0.001
- **Interpretation**: Displacement bearings are significantly different from random but highly concentrated around control mean

### Figure 1 Caption

Tracking locations and detection summary.** Monarch butterflies were tracked using automated receiver networks across eight geographically distributed locations. Release sites are shown as colored dots, with color intensity proportional to the number of detections recorded at each location. (A) All 47,109 detections from 57 butterflies across eight release locations color-coded by origin. (B) Migration paths showing start locations (green circles) and final detection locations (red squares) with gray connecting lines representing individual butterfly trajectories. The control group (Austin, n=9) released from the species' endemic Texas breeding area exhibited convergent migration toward Mexican overwintering sites. Displaced butterflies (n=48) from seven distant locations exhibited similar directional convergence, suggesting a shared migratory destination or navigation mechanism.

### Figure 2 Caption

Bearing analysis reveals control-aligned migration.** Rose diagrams showing bearing angles at 10° intervals for (A) control butterflies (Austin release, n=9) with circular mean of 196.2° toward the Mexico overwintering sites, (B) displaced butterflies (n=48) released from seven non-control locations, and (C) the bearing angle difference between displaced and control butterflies. Displaced butterflies exhibited bearings sharply concentrated around the control group mean (206.6° ± 79.5°, Watson's U²=699.9, p<0.001), with only 10.4° offset from controls, demonstrating convergence on the shared roost destination. Panel C quantifies the minimal angular offset between groups, illustrating their directional alignment. In contrast, displaced bearings showed no concentration toward Mexican overwintering sites (mean angular distance=85.8°), indicating that displaced monarchs navigated toward the Austin control destination rather than exhibiting independent Mexico-directed migration.


## License

This code is released under the [MIT License]
The analysis code is provided to support reproducible research as described in the publication.


## Contact

For questions regarding the code or analysis:
- Email: Isla.duporge@princeton.edu
- 
## Version History

- **v1.0** (2026-03-05): Initial release with publication
