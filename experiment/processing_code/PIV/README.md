# Wave–Turbulence Decomposition (Method A: Linear Wave Theory)

This README documents the two top-level scripts that perform Method A
(linear-wave-theory, LWT) wave–turbulence decomposition on longitudinal
PIV data, plus the helper functions they depend on.

## Pipeline overview

For every image pair the workflow is:

1. Load raw PIV image pair (`*_Piv_NNN_a.mat`, `*_Piv_NNN_b.mat`) and the
   matching surface-camera frames (`*_Pivsurf_NNN_*.mat`).
2. Detect the water surface and build a PIV mask (`FindSurfaceCapillary`).
3. Compute Cartesian velocity `(u, w)` in m/s from the image pair
   (`ComputeVelocities_Quick_NoFilt_Deform_Water`), then clean it with
   `removeOutliers` + `smoothn`.
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

- `LONG`        — root path to the longitudinal PIV experiments
- `rootpath`    — where `GC-Wave-Gen` lives locally (Fabrice's code is added to the path from here)
- `ii`          — experiment index into `dir(LONG)`
- `image_pair_number` — frame to inspect
- `DX`, `DT`    — m/pixel and s/frame calibration
- `SU_OFFSET`   — pixel offset that maps surface-camera coords to PIV-image coords (experiment specific)

Behavior:

- Loads `compVel` from `Chris_recompute/PIVMat/` if cached, otherwise
  computes it from raw images.
- Produces seven figures: surface-camera image, raw PIV pair with
  detected surface, Cartesian `u`/`w`, constant-ζ lines overlaid on `u`,
  wave-following `u`/`w`, LWT orbital `u`/`w`, the residual decomposition
  in wave-following coords, and the same decomposition reverse-transformed
  to the lab frame.
- **Does not save** anything — purely for inspection.

### `run_decomposition_loop.m`

Batch version of the above. Loops over every `*_Piv_*_a.mat` frame in the
experiment, performs the full decomposition, and writes results to disk.
Parameters at the top must be kept in sync with
`run_decomposition_linearWaveTheory.m`.

Outputs (per frame) to `<experiment>/Chris_recompute/`:

- `PIVMat/<exp>_compVel_NNN.mat`
  Raw PIV velocity + surface structs. Written only when newly computed
  (cache miss).
  Variables: `compVel`, `imSurfa`, `imSurfb`.

- `PIVMat_TURB/<exp>_compVel_NNN.mat`
  Full decomposition. Variables: `decomposedVel`, `pivRes`.

  Fields in `decomposedVel.compVel` (all `single`):

  | field          | frame          | description                          |
  |----------------|----------------|--------------------------------------|
  | `u`, `w`       | lab, Cartesian | smoothed PIV velocity (m/s)          |
  | `intrp_u/w`    | wave-following | velocity in wave-following coords    |
  | `ORBX`, `ORBZ` | wave-following | LWT orbitals, **pixels/frame**       |
  | `SU`           | —              | constant-ζ pixel rows on PIV image   |
  | `pf_surf`      | —              | detected surface (`Surface_PIV`)     |
  | `intrp_u/w_res`| wave-following | residual = `intrp − ORBX_ms` (m/s)   |
  | `u_res`, `w_res`     | lab      | residual reverse-transformed (m/s)   |
  | `u/w_mean_lab`       | lab      | wave-following mean back in lab (m/s)|
  | `u/w_orb_lab`        | lab      | LWT orbital back in lab (m/s)        |

The ensemble-average block at the bottom is commented out; uncomment it to
accumulate a 2D time mean in wave-following coordinates.

## Functions used (in this folder or in `GC-Wave-Gen`)

External (in `GC-Wave-Gen/M-Files_FabMarcNovDec2014/`, added to path at the top of each script):

- `FindSurfaceCapillary` — detects the air–water interface from the
  surface-camera frame; returns `surfacePIVImg` (surface row in PIV-image
  pixel coords) and `mask`.
- `ComputeVelocities_Quick_NoFilt_Deform_Water` — multi-pass PIV
  cross-correlation producing `compVel` (with `delta_x`, `delta_z`,
  `xPIV`, `zPIV`, `GS`, `Mask`, `dcor`).
- `generateTransfo_LC_noLFV_2023` — builds the wave-following transform
  and LWT orbital velocities; returns `SU` (constant-ζ surfaces) and
  `ORBX`/`ORBZ` (orbitals in pixels/frame). Optional `opt.U_drift` adds a
  uniform wind-drift offset.
- `transformVelField_decay_forFab` — interpolates a Cartesian field onto
  the wave-following grid defined by `SU`.
- `reverseTransformVelField_decay_forFab` — inverse of the above; maps a
  wave-following field back to the Cartesian PIV grid.
- `removeOutliers` — masks PIV vectors whose correlation peak (`dcor`) is
  below threshold.
- `smoothn` — robust spline smoothing of the velocity components.

Local helpers in this folder (`PIV/`):

- `brewermap.m` — colormap used for all plots.
- `get_server_root.m` — resolves the data-server root.

## Conventions

- Units: `DX` is m/pixel, `DT` is s/frame. ORBX/ORBZ from
  `generateTransfo_LC_noLFV_2023` are **pixels/frame** and must be
  multiplied by `DX/DT` to get m/s.
- `SU(1,:)` is the row corresponding to the surface itself and is removed
  by `SU = transfo.SU(2:end,:)`. The same `2:end` slice is applied to
  ORBX/ORBZ.
- `SU_OFFSET` shifts surface-image pixel rows into PIV-image pixel rows;
  the default `0` assumes the two cameras are already co-registered.
