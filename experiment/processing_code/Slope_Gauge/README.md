## Color Image Slope Gauge (CISG) - Processing Pipeline

### Overview
Processes Color Image Slope Gauge (CISG) data to measure water-surface
slopes via refraction of an RGB pattern viewed through the water. The
pipeline reads CoreView `.raw` frames, builds a per-frame slope cache,
reconstructs surface elevation `eta`, and renders a Veron-style plot
suite.

### Pipeline

```
cisg_main.m              -> per-frame slope cache (Sx, Sy)
make_all_plots_csig.m    -> eta reconstruction + plot suite
```

Run `cisg_main.m` first to populate the cache, then
`make_all_plots_csig.m` to make figures.

### File Structure

```
cisg_main.m                  - Main processing script (cache builder)
make_all_plots_csig.m        - Plotting driver (reads the cache)

cisg_load_coreview.m         - Read CoreView .raw frame -> RGB
cisg_spatial_calibration.m   - Build cm/pixel calibration (run once)
cisg_calculate_slopes.m      - RGB ref + obs -> (Sx, Sy)
cisg_design_pattern.m        - Design the RGB pattern (one-time)
cisg_generate_synthetic.m    - Synthetic test images
cisg_plot_results.m          - Single-frame slope visualization
make_slope_video.m           - Render Sx/Sy MP4 from a slope cube

util/
  load_slope_frame.m         - Read one cached frame
  load_slope_cube.m          - Read a frame range, crop to inner ROI
  cisg_reconstruct_eta.m     - (Sx, Sy) -> eta via spectral integration
  cisg_local_wavenumber.m    - Local k via Hilbert / spectral methods
  dispersion_curves.m        - Deep-water dispersion overlays
  pick_roi_interactive.m     - Interactive ROI picker
  save_figure.m, inferno.m   - Helpers

plotting/
  plot_A1_snapshots_csig.m         plot_A2_spectrograms_csig.m
  plot_A3_kspectra_csig.m          plot_A4_slope_variance.m
  plot_A6_directional_B.m          plot_B1_kx_omega_csig.m
  plot_B3_phase_speed_dev.m        plot_B4_angular_dist.m
  plot_C1_local_k_map.m            plot_C4_cross_stream_coh.m
  plot_C5_heterogeneity_idx.m      plot_C6_local_stokes_map.m
  plot_E3_bandwidth_correction.m
  plot_kx_t_slopes_csig.m          plot_f_t_slopes_csig.m
  plot_kx_omega_Sx_csig.m
```

### cisg_main.m — build the per-frame slope cache

Reads a flat-water reference frame and a sequence of observed `.raw`
frames, computes `(Sx, Sy)` for each, crops to an outer ROI, and writes
one `frame_NNNN.mat` per processed frame plus a `metadata.mat`.

**Setup (top of script):**
```matlab
coreview_num     = 22;                       % observed run
coreview_ref_num = 19;                       % flat-water reference run
raw_root         = '\\Airseaserver28\D\HLAB_2026\SlopeGauge';
cache_root       = 'D:\HLAB_2026\SlopeGauge';

setup.camera_height = 41.7;   % cm above water surface
setup.water_depth   = 49;     % cm from surface to pattern
setup.n_water = 1.333; setup.n_air = 1.000;
setup.frame_rate = 50;        % Hz
```

A spatial calibration produced by `cisg_spatial_calibration.m` is loaded
from `<raw_root>\Processed\CISG_Calibration_CoreView_*.mat` and supplies
`cm_per_pixel_x`, `cm_per_pixel_y`.

**Workflow inside the script:**
1. Load the reference `.raw` and display its R/G/B channels.
2. Pick an **outer ROI** interactively (baked into the cache — pick
   generously, it cannot be expanded without re-running).
3. List all observed `.raw` frames, write `metadata.mat`
   (`setup`, `time`, `coreview_num`, image size).
4. Loop over `frames_to_process = start_frame:10:end_frame`, calling
   `cisg_calculate_slopes(ref_image, obs_image, setup)`, cropping to the
   outer ROI, saving `Sx`, `Sy`, `n`, `t` as `single` to
   `frame_<n>.mat`. Existing frames are skipped, so the loop is
   **resumable**.

**Output:**
```
<cache_root>\CISG_slopes_CoreView_<NN>\
    metadata.mat
    frame_0001.mat, frame_0011.mat, ...   (Sx, Sy, n, t)
```

The script ends at `return;` after the cache is built. The cells below
that line are legacy snippets (statistics, time series, video) kept for
reference and run manually.

### make_all_plots_csig.m — reconstruct eta and plot

Reads the cache built above, reconstructs `eta` from `(Sx, Sy)`, and
runs the plotting suite in `plotting/`.

**Setup (top of script):**
```matlab
coreview_num = 22;
cache_root   = 'D:\HLAB_2026\SlopeGauge';
cache_dir    = fullfile(cache_root, sprintf('CISG_slopes_CoreView_%d', coreview_num));
```

**Workflow:**
1. Load `metadata.mat` -> `setup`, `time`; derive `dx`, `dy` (m), `dt`.
2. Load a saved **inner ROI** (in original full-image coords) from
   `<raw_root>\Processed\CISG_ROI_CoreView_<NN>.mat`, translate to
   cached-frame coords, and clamp.
3. Preview one cached frame and (optionally) recompute the same frame
   from the raw `.raw` files as a sanity check.
4. Preview `eta` for the single frame using
   `util/cisg_reconstruct_eta.m`.
5. Materialize the inner-ROI slope cube:
   ```matlab
   frame_subset = 1111:10:length(time);
   [Sx, Sy] = load_slope_cube(cache_dir, frame_subset, inner_roi);
   ```
6. Subtract a **per-pixel rest-period offset** (mean of `rest_frames =
   1:50`) from `Sx`, `Sy` to remove time-invariant artifacts.
7. Reconstruct `eta = cisg_reconstruct_eta(Sx, Sy, dx, dy)`.
8. Pack into `data = struct('Sx', 'Sy', 'eta', 'time', 'setup', 'dx',
   'dy')` and call the `plotting/plot_*_csig.m` functions. The A1
   snapshot panel is driven by `target_times` (seconds); other plots
   take `data` and optional index sets.

**Output:** figures (PDF/PNG via `util/save_figure.m`) into a local
`figures/` folder.

### Mean square slope (ak)^2 time series

Both scripts include a cell computing `mean(Sx.^2)` and `mean(Sy.^2)`
over the inner ROI per frame, plotted vs `time`.

### Typical end-to-end run

```matlab
% one-time, per camera position:
cisg_spatial_calibration      % -> CISG_Calibration_CoreView_<N>.mat

% per CoreView dataset:
cisg_main                     % pick outer ROI, build the cache
make_all_plots_csig           % reconstruct eta, render figures
```
