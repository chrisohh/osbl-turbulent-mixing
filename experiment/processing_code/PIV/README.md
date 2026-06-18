# Wave–Turbulence Decomposition (Method A: Linear Wave Theory)

This README documents the two top-level scripts that perform Method A
(linear-wave-theory, LWT) wave–turbulence decomposition on longitudinal
PIV data, plus the helper functions they depend on.

## Pipeline overview

For every image pair the workflow is:

1. Load raw PIV image pair (`*_Piv_NNN_a.mat`, `*_Piv_NNN_b.mat`) and the
   matching surface-camera frames (`*_Pivsurf_NNN_*.mat`).
2. Detect the water surface and build a PIV mask (`FindSurfaceCapillary`),
   then extend the mask downward by `GLINT_BUFFER` pixels to exclude the
   surface glint/foam band (`apply_glint_buffer`).
3. Compute Cartesian velocity `(u, w)` in m/s from the image pair with
   multi-pass PIV (`ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt`,
   including inter-pass universal outlier detection), validate and clean it
   (`validatePIV` — UOD + group removal + gap fill), then smooth with
   `smoothn`.
4. Build the wave-following coordinate transform and linear-wave-theory
   orbital velocity field (`generateTransfo_LC_noLFV_2023`).
5. Interpolate the Cartesian velocity onto the wave-following grid
   (`transformVelField_decay_forFab`).
6. Subtract the LWT orbital field from the wave-following velocity to get
   the residual (mean + turbulence).
7. Transform the residual, orbital, and mean fields back to the lab frame
   (`reverseTransformVelField_decay_forFab`).

## Scripts

### `run_decomposition_linearWaveTheory.m`

Single-frame, step-by-step run that plots every stage of the pipeline for
sanity checking. **Run this first on a representative frame before
batching.**

Edit at the top of the script:

- `get_LCL_params('ExpLCL_2_01')` — change the experiment name string to select the experiment (see `get_LCL_params.m`)
- `rootpath`    — where `GC-Wave-Gen` lives locally (Fabrice's code is added to the path from here)
- `image_pair_number` — frame to inspect
- `rerun`       — `0` loads cached `compVel`; `1` recomputes PIV from raw images

Mid-script (STEP 2):

- `SU_OFFSET`   — pixel offset that maps surface-camera coords to PIV-image coords (experiment specific; default `0`)

Behavior:

- Loads `compVel` from `Chris_recompute/PIVMat/` if cached (and `rerun = 0`),
  otherwise computes it from raw images using a glint-buffered mask,
  inter-pass UOD, and `validatePIV`.
- Produces a sequence of figures: surface-camera image, raw PIV pair with
  detected surface, (when computing) the glint-buffer mask overlay,
  Cartesian `u`/`w`, constant-ζ lines overlaid on `u`, wave-following
  `u`/`w`, LWT orbital `u`/`w`, the residual decomposition in
  wave-following coords (`figure(6)`), and the same decomposition
  reverse-transformed to the lab frame (`figure(7)`).
- **Does not save** anything — purely for inspection.

### `run_decomposition_loop.m`

Batch version of the above. Loops over every `*_Piv_*_a.mat` frame in the
experiment, performs the full decomposition, and writes results to disk.
Change the experiment name in `get_LCL_params(...)` at the top to select a different experiment.

Other top-of-script knobs (keep PIV settings in sync with `run_decomposition_linearWaveTheory.m`):

- `recompute_piv` — `true` reruns `ComputeVelocities` even when a `PIVMat`
  file exists (reusing saved surfaces); `false` uses cached `compVel`.
- `SU_OFFSET`     — surface-image → PIV-image pixel offset (default `0`).
- `GLINT_BUFFER`  — pixels below the detected surface to exclude.
- `IntrWndw` / `GrdSpc` — multi-pass interrogation-window and grid-spacing pyramid.
- `iuod`          — inter-pass universal-outlier-detection settings.
- `val_opts`      — `validatePIV` options (UOD / group removal / gap fill).

Outputs (per frame) to `<experiment>/Chris_recompute/`:

- `PIVMat/<exp>_compVel_NNN.mat`
  Raw PIV velocity + surface structs. Written only when newly computed
  (cache miss).
  Variables: `compVel`, `imSurfa`, `imSurfb`.

- `PIVMat_TURB/<exp>_compVel_NNN.mat`
  Full decomposition. Variables: `decomposedVel`, `pivRes`.

  Fields in `decomposedVel.compVel` (velocities stored as `single`):

  | field                     | frame          | description                                   |
  |---------------------------|----------------|-----------------------------------------------|
  | `u_raw`, `w_raw`          | lab, Cartesian | smoothed measured PIV velocity (m/s)          |
  | `intrp_u_raw`, `intrp_w_raw` | wave-following | measured velocity in wave-following coords |
  | `ORBX_ms`, `ORBZ_ms`      | wave-following | LWT orbital velocity (m/s)                    |
  | `intrp_u_res`, `intrp_w_res` | wave-following | residual mean+turb = `intrp_raw − ORB_ms` (m/s) |
  | `u_res_lab`, `w_res_lab`  | lab            | residual reverse-transformed to lab (m/s)     |
  | `u_orb_lab`, `w_orb_lab`  | lab            | LWT orbital reverse-transformed to lab (m/s)  |
  | `SU`                      | —              | constant-ζ pixel rows on PIV image            |
  | `pf_surf`                 | —              | detected surface (`Surface_PIV`)              |
  | `dcor`                    | lab, Cartesian | PIV correlation quality (NaN=air; mask `dcor<0.4` before computing Reynolds stresses) |

The ensemble-average block at the bottom is commented out; uncomment it to
accumulate a 2D time mean in wave-following coordinates.

### `make_decomposition_video.m`

Visualizes the saved `PIVMat_TURB/` results as a 3×2 lab-frame panel
(raw `u`/`w`, LWT orbital `u`/`w`, and residual `u_raw−u_orb` / `w_raw−w_orb`).
**Run `run_decomposition_loop.m` first** to populate `PIVMat_TURB/`.

Set the experiment via `get_LCL_params(...)` at the top, then choose a mode:

- `makeVideo = 1` — render every frame to an MP4 in `D:\DelawareDataResult\`.
- `makeVideo = 0` — display a **single frame** only (set `frame` to the frame
  number shown in the panel title). Use this to inspect one frame with and
  without the wave removed: the top row is `w_raw` (wave present) and the
  bottom row is `w_raw − w_orb` (wave removed).

## Functions used (in this folder or in `GC-Wave-Gen`)

External (in `GC-Wave-Gen/M-Files_FabMarcNovDec2014/`, added to path at the top of each script):

- `FindSurfaceCapillary` — detects the air–water interface from the
  surface-camera frame; returns `surfacePIVImg` (surface row in PIV-image
  pixel coords) and `mask`.
- `ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt` — multi-pass PIV
  cross-correlation with image deformation and inter-pass UOD, producing
  `compVel` (with `delta_x`, `delta_z`, `xPIV`, `zPIV`, `GS`, `Mask`,
  `dcor`).
- `generateTransfo_LC_noLFV_2023` — builds the wave-following transform
  and LWT orbital velocities; returns `SU` (constant-ζ surfaces) and
  `ORBX`/`ORBZ` (orbitals in pixels/frame). Optional `opt.U_drift` adds a
  uniform wind-drift offset.
- `transformVelField_decay_forFab` — interpolates a Cartesian field onto
  the wave-following grid defined by `SU`.
- `reverseTransformVelField_decay_forFab` — inverse of the above; maps a
  wave-following field back to the Cartesian PIV grid.
- `smoothn` — robust spline smoothing of the velocity components.

Local helpers in this folder (`PIV/`):

- `get_LCL_params.m` — returns a struct of all timing, path, and ramp parameters for a given longitudinal experiment name (e.g. `'ExpLCL_2_01'`). All scripts call this instead of hardcoding `ii`, `DX`, `DT`, paths, or ramp delays.
- `apply_glint_buffer.m` — extends the surface mask downward by `GLINT_BUFFER` pixels to exclude the surface glint/foam band before PIV.
- `validatePIV.m` — post-PIV cleanup: universal outlier detection, group removal, gap filling (controlled by `val_opts`).
- `brewermap.m` — colormap used for all plots.
- `get_server_root.m` — resolves the data-server root.

## Conventions

### Timing

All cameras and lasers were synchronized. Acquisition rates:

| quantity | value |
|---|---|
| Camera image rate (`Fs_cam`) | 14.4 Hz |
| PIV velocity map rate (`Fs_PIV`) | 7.2 Hz (one pair = two images) |
| Time between consecutive pairs (`DT_pair`) | 1/7.2 ≈ 0.139 s — **use for time axis** |
| Inter-frame dt LCL (`DT`) | 10 ms — **use for velocity scaling only** (`u = delta_x * DX / DT`) |

**Time axis convention:** t = 0 is when the wind starts ramping up.

Wind protocol for ramp N:
- t = 0 s — wind trigger (stays at 0 V)
- t = 40 s — wind begins ramping
- t = 5 + IID s — imaging starts (IID = imaging initial delay, ramp-dependent)

So imaging starts at `t_imaging_wrt_ramp = (5 + IID) - 40` seconds after wind ramp-up.
`get_LCL_params` returns this as `p.t_imaging_wrt_ramp`; the time axis is:

```matlab
t = p.t_imaging_wrt_ramp + (0:N_frames-1) * DT_pair;
```

Ramp parameters:

| ramp | IID (s) | imaging start re wind ramp (s) | ramp-up duration (s) | voltage (V) |
|------|---------|-------------------------------|----------------------|-------------|
| 1 | 57 | 22 | 30 | 7 |
| 2 | 73 | 38 | 90 | 7 |
| 3 | 88 | 53 | 120 | 7 |
| 4 | 111 | 76 | 120 | 5 |

### Velocity and coordinates

- `DX` is m/pixel, `DT` is s (inter-frame, A→B within a pair).
- ORBX/ORBZ from `generateTransfo_LC_noLFV_2023` are **pixels/frame** and must be multiplied by `DX/DT` to get m/s.
- Velocities in `decomposedVel.compVel` are already in m/s.
- `intrp_u_res`, `intrp_w_res` are in the **curvilinear (wave-following)** coordinate system: ζ = 0 at the surface, increasing downward. Row i corresponds to i × GS pixels below the surface.
- `SU(1,:)` is the surface row and is removed by `SU = transfo.SU(2:end,:)`. The same `2:end` slice is applied to ORBX/ORBZ.
- `SU_OFFSET` shifts surface-image pixel rows into PIV-image pixel rows; the default `0` assumes the two cameras are already co-registered.
