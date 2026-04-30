function [v_turb, w_turb, phase_avg_info] = remove_wave_phase_avg(v_meas, w_meas, phi_t, Nbins)
%REMOVE_WAVE_PHASE_AVG  Phase-average wave removal for transverse PIV.
%
%   Bins transverse PIV frames by wave phase (from longitudinal eta),
%   computes the phase-averaged velocity in each bin, and subtracts it.
%
%   Inputs:
%     v_meas  [Nt x Ny x Nz]  cross-wind velocity
%     w_meas  [Nt x Ny x Nz]  vertical velocity
%     phi_t   [Nt x 1]        wave phase at PIV location for each frame
%     Nbins   scalar           number of phase bins (default 16)
%
%   Outputs:
%     v_turb, w_turb  [Nt x Ny x Nz]  turbulent fluctuations
%     phase_avg_info   struct with diagnostics

if nargin < 4 || isempty(Nbins)
    Nbins = 16;
end

[Nt, Ny, Nz] = size(w_meas);

% --- Wrap phase to [0, 2*pi) ---
phi_wrapped = mod(phi_t(:), 2*pi);

% --- Assign each frame to a phase bin ---
bin_edges = linspace(0, 2*pi, Nbins + 1);
bin_idx   = discretize(phi_wrapped, bin_edges);

% --- Compute phase-averaged velocity in each bin ---
v_phavg = zeros(Nbins, Ny, Nz);
w_phavg = zeros(Nbins, Ny, Nz);
counts  = zeros(Nbins, 1);

for ib = 1:Nbins
    mask = (bin_idx == ib);
    counts(ib) = sum(mask);
    if counts(ib) > 0
        v_phavg(ib,:,:) = mean(v_meas(mask,:,:), 1);
        w_phavg(ib,:,:) = mean(w_meas(mask,:,:), 1);
    end
end

% --- Subtract phase average from each frame ---
v_turb = zeros(size(v_meas));
w_turb = zeros(size(w_meas));

for it = 1:Nt
    ib = bin_idx(it);
    v_turb(it,:,:) = v_meas(it,:,:) - v_phavg(ib,:,:);
    w_turb(it,:,:) = w_meas(it,:,:) - w_phavg(ib,:,:);
end

% --- Pack diagnostics ---
phase_avg_info.Nbins     = Nbins;
phase_avg_info.bin_edges = bin_edges;
phase_avg_info.bin_idx   = bin_idx;
phase_avg_info.counts    = counts;
phase_avg_info.v_phavg   = v_phavg;
phase_avg_info.w_phavg   = w_phavg;

end
