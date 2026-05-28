## Color Image Slope Gauge (CISG) - Modular Implementation

### Overview
This code designs a Color Image Slope Gauge system and processes images to measure water surface slopes via refraction.

### File Structure
```
cisg_main.m              - Main script (run this)
cisg_design_pattern.m    - Design system and generate pattern
cisg_generate_synthetic.m - Create synthetic test data
cisg_calculate_slopes.m  - Calculate slopes from images
cisg_plot_results.m      - Visualization
```

### Quick Start

1. **Run the main script:**
   ```matlab
   cisg_main
   ```
   This will:
   - Design the CISG system
   - Generate a pattern (CISG_pattern.png)
   - Create synthetic test data
   - Calculate slopes
   - Validate and visualize results

2. **Adjust parameters in cisg_main.m:**
   ```matlab
   params.camera_height = 100;  % cm
   params.water_depth = 30;     % cm
   params.focal_length_mm = 16; % mm
   ```

### Using with Real Data

After collecting experimental images:

```matlab
% Load design parameters
load('CISG_design.mat');

% Load your images
ref_image = imread('your_reference.png');
obs_image = imread('your_observed.png');

% Calculate slopes
[Sx, Sy] = cisg_calculate_slopes(ref_image, obs_image, design);

% Visualize (without ground truth)
x = linspace(0, design.FOV_width, size(Sx,2));
y = linspace(0, design.FOV_height, size(Sx,1));
cisg_plot_results([], [], Sx, Sy, [], x, y);
```

### Functions

**cisg_design_pattern(params)**
- Input: camera and optical parameters
- Output: design struct, RGB pattern image
- Saves: CISG_pattern.png, CISG_design.mat

**cisg_generate_synthetic(design, pattern, wave_params)**
- Input: design, pattern, wave characteristics
- Output: reference image, observed image, ground truth slopes

**cisg_calculate_slopes(ref_image, obs_image, design)**
- Input: reference and observed RGB images, design parameters
- Output: Sx, Sy slope arrays

**cisg_plot_results(Sx_true, Sy_true, Sx_meas, Sy_meas, eta, x, y)**
- Input: true and measured slopes, elevation, coordinates
- Creates comparison figures

### Workflow

1. Design: Determine pattern size and resolution
2. Print: Print CISG_pattern.png at recommended DPI
3. Setup: Install pattern under water, position camera
4. Acquire: Capture reference (flat) and observed (wavy) images
5. Process: Run cisg_calculate_slopes() on your images
6. Analyze: Visualize and export slope fields
