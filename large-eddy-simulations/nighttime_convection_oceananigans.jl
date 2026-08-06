# # Nighttime convection LES with pure Oceananigans
#
# Large-eddy simulation of nighttime convective mixing in the presence of
# light winds and surface waves. Prescribes surface buoyancy and momentum
# fluxes directly without using bulk formula.
#
# This simulation is designed to compare with laboratory experiments of
# wind-driven turbulence over a freshwater tank.

using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie
using JLD2

# ## Physical parameters
#
# Prescribed surface fluxes for nighttime cooling with light winds.
# Values are chosen to be representative of typical conditions.

const Qᵇ = 1e-7              # Surface buoyancy flux (m²/s³) - cooling
const τˣ = 0.02              # Wind stress x-component (N/m²) - light wind ~3 m/s  
const τʸ = 0.0               # Wind stress y-component (N/m²)
const ρₒ = 1000.0            # Reference density (kg/m³) - freshwater
const g = 9.81               # Gravitational acceleration (m/s²)
const α = 2e-4               # Thermal expansion coefficient (1/K) for freshwater

# Convert buoyancy flux to heat flux for reference: Q = Qᵇ * ρₒ * cₚ / (g * α)
# With cₚ ≈ 4000 J/(kg·K), Q ≈ 200 W/m² cooling

# Wave parameters for Stokes drift (simple monochromatic approximation)
const wave_amplitude = 0.02  # Wave amplitude (m)
const wave_wavenumber = 2.0  # Wavenumber (rad/m)
const wave_frequency = sqrt(g * wave_wavenumber)  # Deep water dispersion

# Stokes drift: uˢ = a² k ω exp(2kz)
const uˢ_surface = wave_amplitude^2 * wave_wavenumber * wave_frequency

# ## Grid configuration
#
# Small 3D domain for initial testing (32³ at ~1m resolution)

arch = CPU()
Lx, Ly, Lz = 32, 32, 32  # meters
Nx, Ny, Nz = 32, 32, 32  # grid points

grid = RectilinearGrid(arch; size=(Nx, Ny, Nz), halo=(5, 5, 5),
                       topology=(Periodic, Periodic, Bounded),
                       x=(0, Lx), y=(0, Ly), z=(-Lz, 0))

# ## Boundary conditions
#
# Top: prescribed momentum and buoyancy fluxes
# Bottom: no-flux (free-slip)

u_top_bc = FluxBoundaryCondition(τˣ / ρₒ)  # Kinematic stress (m²/s²)
v_top_bc = FluxBoundaryCondition(τʸ / ρₒ)
b_top_bc = FluxBoundaryCondition(Qᵇ)       # Buoyancy flux (m²/s³)

u_bcs = FieldBoundaryConditions(top=u_top_bc)
v_bcs = FieldBoundaryConditions(top=v_top_bc)
b_bcs = FieldBoundaryConditions(top=b_top_bc)

# ## Stokes drift

stokes_drift = UniformStokesDrift(;
    ∂z_uˢ = (z, t) -> 2 * wave_wavenumber * uˢ_surface * exp(2 * wave_wavenumber * z)
)

# ## Model setup

coriolis = FPlane(latitude=33)

model = NonhydrostaticModel(;
    grid,
    advection = WENO(order=9),
    coriolis,
    stokes_drift,
    buoyancy = BuoyancyTracer(),
    tracers = :b,
    boundary_conditions = (u=u_bcs, v=v_bcs, b=b_bcs)
)

# Initial conditions: uniform buoyancy with small perturbations
bᵢ(x, y, z) = 0.01 * randn()
set!(model, b=bᵢ)

# ## Simulation

simulation = Simulation(model; Δt=1.0, stop_time=1hour)

# Adaptive time stepping
wizard = TimeStepWizard(cfl=0.7, max_Δt=10.0)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# ## Output

u, v, w = model.velocities
b = model.tracers.b

# Surface slices and profiles
uˢ = view(u, :, :, Nz)
wˢ = view(w, :, :, Nz)
bˢ = view(b, :, :, Nz)
u_avg = Average(u, dims=(1, 2))
b_avg = Average(b, dims=(1, 2))
w²_avg = Average(w^2, dims=(1, 2))  # Vertical velocity variance

outputs = (; uˢ, wˢ, bˢ, u_avg, b_avg, w²_avg)

simulation.output_writers[:snapshots] = JLD2Writer(model, outputs;
                                                          filename="nighttime_convection_oceananigans",
                                                          schedule=TimeInterval(5minutes),
                                                          overwrite_existing=true)

# ## Progress logging

function progress(sim)
    t = time(sim)
    min = t / 60
    u, v, w = sim.model.velocities
    b = sim.model.tracers.b
    
    bmax, bmin = maximum(interior(b)), minimum(interior(b))
    umax, wmax = maximum(abs, u), maximum(abs, w)
    
    # Estimate convective velocity scale: w* = (Qᵇ * Lz)^(1/3)
    w_star = (abs(Qᵇ) * Lz)^(1/3)
    
    @info @sprintf("t = %5.1f min | Δt: %s | b: %.2e–%.2e | max|u|: %.2e (u*/w*=%.1f), max|w|: %.2e (w/w*=%.1f)",
                   min, prettytime(sim.Δt), bmin, bmax, umax, umax/w_star, wmax, wmax/w_star)
end

add_callback!(simulation, progress, IterationInterval(50))

# ## Run

@info "Running nighttime convection simulation with pure Oceananigans..."
@info @sprintf("  Buoyancy flux: Qᵇ = %.2e m²/s³ (cooling)", Qᵇ)
@info @sprintf("  Wind stress: τ = %.3f N/m²", τˣ)
@info @sprintf("  Convective velocity scale: w* = %.3f m/s", (abs(Qᵇ) * Lz)^(1/3))
@info @sprintf("  Surface Stokes drift: uˢ = %.4f m/s", uˢ_surface)

run!(simulation)
@info "Simulation complete!"

# ## Visualization

@info "Creating visualization..."

file = jldopen("nighttime_convection_oceananigans.jld2")
iterations = keys(file["timeseries/t"])
times = [file["timeseries/t/$i"] for i in iterations]
close(file)

b_ts = FieldTimeSeries("nighttime_convection_oceananigans.jld2", "bˢ")
u_ts = FieldTimeSeries("nighttime_convection_oceananigans.jld2", "uˢ")
b_avg_ts = FieldTimeSeries("nighttime_convection_oceananigans.jld2", "b_avg")
w²_avg_ts = FieldTimeSeries("nighttime_convection_oceananigans.jld2", "w²_avg")

Nt = length(b_ts.times)
xc = xnodes(b_ts)
yc = ynodes(b_ts)
xf = xnodes(u_ts)
zc = znodes(b_avg_ts)
zf = znodes(w²_avg_ts)

fig = Figure(size=(1200, 800))

# Surface buoyancy
ax1 = Axis(fig[1, 1], xlabel="x (m)", ylabel="y (m)", title="Surface Buoyancy (m/s²)")
n = Observable(Nt)
b_slice = @lift interior(b_ts[$n], :, :, 1)
hm = heatmap!(ax1, xc, yc, b_slice, colormap=:thermal)
Colorbar(fig[1, 2], hm)

# Surface velocity
ax2 = Axis(fig[1, 3], xlabel="x (m)", ylabel="y (m)", title="Surface u (m/s)")
u_slice = @lift interior(u_ts[$n], :, :, 1)
hm2 = heatmap!(ax2, xf, yc, u_slice, colormap=:balance)
Colorbar(fig[1, 4], hm2)

# Buoyancy profile
ax3 = Axis(fig[2, 1:2], xlabel="b (m/s²)", ylabel="z (m)", title="Mean Buoyancy Profile")
b_profile = @lift interior(b_avg_ts[$n], 1, 1, :)
lines!(ax3, b_profile, zc, linewidth=2)

# Vertical velocity variance profile (on Face z-nodes)
ax4 = Axis(fig[2, 3:4], xlabel="⟨w²⟩ (m²/s²)", ylabel="z (m)", title="Vertical Velocity Variance")
w²_profile = @lift interior(w²_avg_ts[$n], 1, 1, :)
lines!(ax4, w²_profile, zf, linewidth=2, color=:blue)

save("nighttime_convection_oceananigans.png", fig)
@info "Saved figure to nighttime_convection_oceananigans.png"

