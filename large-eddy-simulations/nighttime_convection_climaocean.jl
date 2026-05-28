# # Nighttime convection LES with ClimaOcean
#
# Large-eddy simulation of nighttime convective mixing in the presence of
# light winds and surface waves. Uses ClimaOcean's bulk formula approach
# to compute surface fluxes from a prescribed constant atmosphere state.
#
# This simulation is designed to compare with laboratory experiments of
# wind-driven turbulence over a freshwater tank.

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Printf
using Statistics: mean
using CairoMakie
using JLD2

# ## Physical parameters
#
# Light wind, surface cooling conditions typical of nighttime.
# Freshwater (S=0) to match laboratory conditions.

const Tₐ = 288.0            # Atmospheric temperature (K) ~ 15°C
const Tₒ = 293.0            # Initial ocean temperature (K) ~ 20°C  
const qₐ = 0.008            # Atmospheric specific humidity (kg/kg)
const Uₐ = 3.0              # Wind speed (m/s) - light winds
const Qˡʷ = 300.0           # Downwelling longwave (W/m²) - clear night sky
const Qˢʷ = 0.0             # Shortwave radiation (W/m²) - nighttime

# Wave parameters for Stokes drift (simple monochromatic approximation)
const wave_amplitude = 0.02  # Wave amplitude (m)
const wave_wavenumber = 2.0  # Wavenumber (rad/m)
const wave_frequency = sqrt(9.81 * wave_wavenumber)  # Deep water dispersion

# Stokes drift vertical shear: ∂uˢ/∂z = 2 k² a² ω exp(2kz)
const uˢ_surface = wave_amplitude^2 * wave_wavenumber * wave_frequency  # Surface Stokes drift

# ## Grid configuration
#
# Small 3D domain for initial testing (32³ at ~1m resolution)

arch = CPU()
Lx, Ly, Lz = 32, 32, 32  # meters
Nx, Ny, Nz = 32, 32, 32  # grid points

grid = RectilinearGrid(arch; size=(Nx, Ny, Nz), halo=(5, 5, 5),
                       topology=(Periodic, Periodic, Bounded),
                       x=(0, Lx), y=(0, Ly), z=(-Lz, 0))

# Stokes drift with exponential profile
stokes_drift = UniformStokesDrift(; ∂z_uˢ = (z, t) -> 2 * wave_wavenumber * uˢ_surface * exp(2 * wave_wavenumber * z))

# Set up ocean simulation with Stokes drift
coriolis = FPlane(latitude=33)
ocean = nonhydrostatic_ocean_simulation(grid; coriolis, stokes_drift)
conjure_time_step_wizard!(ocean, cfl=0.7)

# Initial conditions: uniform temperature freshwater with small perturbations
Tᵢ(x, y, z) = 20 + 0.01 * randn()  # °C
Sᵢ(x, y, z) = 0.0                   # Freshwater
set!(ocean.model, T=Tᵢ, S=Sᵢ)

# ## Constant atmosphere
#
# Create a simple atmosphere with constant properties (no time variation)

atmosphere_grid = RectilinearGrid(arch; size=(1, 1),
                                  topology=(Periodic, Periodic, Flat),
                                  x=(0, Lx), y=(0, Ly))

atmosphere_times = [0.0, 2hours]  # Two time points for constant atmosphere
atmosphere = PrescribedAtmosphere(atmosphere_grid, atmosphere_times)

# Set constant atmospheric state
set!(atmosphere;
     u = (x, y, t) -> Uₐ,
     v = (x, y, t) -> 0.0,
     T = (x, y, t) -> Tₐ,
     q = (x, y, t) -> qₐ,
     shortwave = (x, y, t) -> Qˢʷ,
     longwave = (x, y, t) -> Qˡʷ)

# ## Coupled model

coupled_model = OceanSeaIceModel(ocean; atmosphere)
simulation = Simulation(coupled_model; Δt=1.0, stop_time=1hour)

add_callback!(simulation, ClimaOcean.OceanSeaIceModels.align_component_steps!)

# ## Output

u, v, w = ocean.model.velocities
T = ocean.model.tracers.T
Jᵀ = coupled_model.interfaces.net_fluxes.ocean.T
τx = coupled_model.interfaces.net_fluxes.ocean.u

ρₒ = coupled_model.interfaces.ocean_properties.reference_density
cₒ = coupled_model.interfaces.ocean_properties.heat_capacity
Q = ρₒ * cₒ * Jᵀ  # Heat flux (W/m²)

# Surface slices and profiles
Nz = grid.Nz
uˢ = view(u, :, :, Nz)
wˢ = view(w, :, :, Nz)
Tˢ = view(T, :, :, Nz)
u_avg = Average(u, dims=(1, 2))
T_avg = Average(T, dims=(1, 2))
Q_avg = Average(Q, dims=(1, 2))

outputs = (; uˢ, wˢ, Tˢ, Q, u_avg, T_avg, Q_avg)

simulation.output_writers[:snapshots] = JLD2Writer(ocean.model, outputs;
                                                   filename="nighttime_convection_climaocean",
                                                   schedule=TimeInterval(5minutes),
                                                   overwrite_existing=true)

# ## Progress logging

function progress(sim)
    t = time(sim)
    min = t / 60
    ocean_Δt = sim.model.ocean.Δt
    u, v, w = sim.model.ocean.model.velocities
    T = sim.model.ocean.model.tracers.T
    
    Tmax, Tmin = maximum(interior(T)), minimum(interior(T))
    umax, wmax = maximum(abs, u), maximum(abs, w)
    
    Jᵀ = sim.model.interfaces.net_fluxes.ocean.T
    ρₒ = sim.model.interfaces.ocean_properties.reference_density
    cₒ = sim.model.interfaces.ocean_properties.heat_capacity
    Q_mean = ρₒ * cₒ * mean(Jᵀ)
    
    @info @sprintf("t = %5.1f min | Δt: %s | SST: %.3f–%.3f°C | Q: %.1f W/m² | max|u|: %.2e, max|w|: %.2e",
                   min, prettytime(ocean_Δt), Tmin, Tmax, Q_mean, umax, wmax)
end

add_callback!(simulation, progress, IterationInterval(50))

# ## Run

@info "Running nighttime convection simulation with ClimaOcean..."
run!(simulation)
@info "Simulation complete!"

# ## Visualization

@info "Creating visualization..."

file = jldopen("nighttime_convection_climaocean.jld2")
iterations = keys(file["timeseries/t"])
times = [file["timeseries/t/$i"] for i in iterations]
close(file)

T_ts = FieldTimeSeries("nighttime_convection_climaocean.jld2", "Tˢ")
u_ts = FieldTimeSeries("nighttime_convection_climaocean.jld2", "uˢ")
Q_avg_ts = FieldTimeSeries("nighttime_convection_climaocean.jld2", "Q_avg")
T_avg_ts = FieldTimeSeries("nighttime_convection_climaocean.jld2", "T_avg")

Nt = length(T_ts.times)
xc = xnodes(T_ts)
yc = ynodes(T_ts)
zc = znodes(T_avg_ts)

fig = Figure(size=(1200, 800))

# Surface temperature
ax1 = Axis(fig[1, 1], xlabel="x (m)", ylabel="y (m)", title="Surface Temperature (°C)")
n = Observable(Nt)
T_slice = @lift interior(T_ts[$n], :, :, 1)
hm = heatmap!(ax1, xc, yc, T_slice, colormap=:thermal)
Colorbar(fig[1, 2], hm)

# Surface velocity
ax2 = Axis(fig[1, 3], xlabel="x (m)", ylabel="y (m)", title="Surface u (m/s)")
u_slice = @lift interior(u_ts[$n], :, :, 1)
hm2 = heatmap!(ax2, xc, yc, u_slice, colormap=:balance)
Colorbar(fig[1, 4], hm2)

# Temperature profile
ax3 = Axis(fig[2, 1:2], xlabel="T (°C)", ylabel="z (m)", title="Mean Temperature Profile")
T_profile = @lift collect(interior(T_avg_ts[$n], 1, 1, :))
zc_vec = collect(zc)
lines!(ax3, T_profile, zc_vec, linewidth=2)

# Heat flux time series
ax4 = Axis(fig[2, 3:4], xlabel="Time (min)", ylabel="Q (W/m²)", title="Surface Heat Flux")
Q_series = Float64[interior(Q_avg_ts[i], 1, 1, 1)[] for i in 1:Nt]
t_min = collect(T_ts.times ./ 60)
lines!(ax4, t_min, Q_series, linewidth=2, color=:red)

save("nighttime_convection_climaocean.png", fig)
@info "Saved figure to nighttime_convection_climaocean.png"

