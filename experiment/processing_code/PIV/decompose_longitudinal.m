function [D, pivRes] = decompose_longitudinal(u_raw, w_raw, compVel, Surface_PIV, SU_OFFSET, opt)
%DECOMPOSE_LONGITUDINAL  Method A (linear wave theory) wave-turbulence
% decomposition for LONGITUDINAL PIV.
%
% Valid only when the image x-axis is the wave-propagation direction: the FFT
% of Surface_PIV gives the horizontal wavenumber modes, exp(-h*k) decay and the
% sqrt(g*k) dispersion relation in generateTransfo_LC_noLFV_2023 then apply
% directly. For the transverse plane use decompose_transverse instead.
%
% Inputs
%   u_raw, w_raw  measured velocity, lab frame, m/s  [Nz x Nx]
%   compVel       PIV result (needs .zPIV .xPIV .GS .Mask .DX .DT)
%   Surface_PIV   detected surface row coords, PIV image px  [1 x w]
%   SU_OFFSET     surface-image -> PIV-image pixel offset (scalar)
%   opt           optional; opt.U_drift (m/s) for Doppler-shifted dispersion
%
% Output D (all fields m/s unless noted):
%   u_raw, w_raw                 measured, lab frame (pass-through)
%   intrp_u_raw, intrp_w_raw     measured, wave-following
%   ORBX_ms, ORBZ_ms             orbital, wave-following
%   intrp_u_res, intrp_w_res     mean+turb residual, wave-following
%   u_res_lab, w_res_lab         mean+turb, lab frame
%   u_orb_lab, w_orb_lab         orbital, lab frame
%   u_raw_lab, w_raw_lab         measured round-tripped to lab (sanity check)
%   SU                           curvilinear (constant-zeta) grid, PIV px
% pivRes is returned with .pf_surf set, for saving / reverse transforms.

if nargin < 6, opt = struct(); end
DX = compVel.DX;  DT = compVel.DT;

pivRes.zPIV = compVel.zPIV;
pivRes.xPIV = compVel.xPIV;
pivRes.GS   = compVel.GS;
pivRes.mask = compVel.Mask;

% --- linear wave theory transform + orbitals ---
if isfield(opt, 'U_drift')
    transfo = generateTransfo_LC_noLFV_2023(compVel, Surface_PIV, pivRes, opt);
else
    transfo = generateTransfo_LC_noLFV_2023(compVel, Surface_PIV, pivRes);
end
SU      = transfo.SU(2:end,:) + SU_OFFSET;   % drop zeta=0 (surface) row
ORBX    = transfo.ORBX(2:end,:);             % pixels/frame
ORBZ    = transfo.ORBZ(2:end,:);
ORBX_ms = ORBX * DX/DT;                       % m/s
ORBZ_ms = ORBZ * DX/DT;
pivRes.pf_surf = SU(1,:);

% --- transform raw velocity to wave-following frame ---
intrp_u_raw = transformVelField_decay_forFab(u_raw, pivRes, SU);
intrp_w_raw = transformVelField_decay_forFab(w_raw, pivRes, SU);

% --- subtract orbitals -> mean+turb residual (wave-following) ---
intrp_u_res = intrp_u_raw - ORBX_ms;
intrp_w_res = intrp_w_raw - ORBZ_ms;

% --- reverse-transform to lab frame ---
u_res_lab = reverseTransformVelField_decay_forFab(intrp_u_res, pivRes, SU);
w_res_lab = reverseTransformVelField_decay_forFab(intrp_w_res, pivRes, SU);
u_orb_lab = reverseTransformVelField_decay_forFab(ORBX_ms,     pivRes, SU);
w_orb_lab = reverseTransformVelField_decay_forFab(ORBZ_ms,     pivRes, SU);
u_raw_lab = reverseTransformVelField_decay_forFab(intrp_u_raw, pivRes, SU);
w_raw_lab = reverseTransformVelField_decay_forFab(intrp_w_raw, pivRes, SU);

% --- pack ---
D.u_raw       = u_raw;        D.w_raw       = w_raw;
D.intrp_u_raw = intrp_u_raw;  D.intrp_w_raw = intrp_w_raw;
D.ORBX_ms     = ORBX_ms;      D.ORBZ_ms     = ORBZ_ms;
D.intrp_u_res = intrp_u_res;  D.intrp_w_res = intrp_w_res;
D.u_res_lab   = u_res_lab;    D.w_res_lab   = w_res_lab;
D.u_orb_lab   = u_orb_lab;    D.w_orb_lab   = w_orb_lab;
D.u_raw_lab   = u_raw_lab;    D.w_raw_lab   = w_raw_lab;
D.SU          = SU;
end
