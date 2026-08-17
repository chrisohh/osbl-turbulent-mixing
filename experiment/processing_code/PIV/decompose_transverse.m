function [D, pivRes] = decompose_transverse(u_raw, w_raw, compVel, Surface_PIV, SU_OFFSET, opt)
%DECOMPOSE_TRANSVERSE  Wave removal for TRANSVERSE (cross-wind) PIV.
%
% Implements OSBL-notes/wave_removal.tex. The PIV plane is cross-wind (y,z) and
% waves propagate along-wind (x), so:
%   * the wave orbital velocity is UNIFORM across y;
%   * there is NO cross-wave (in-plane horizontal) orbital  -> u has no wave part;
%   * the only wave signal in-plane is the vertical heave w_wave(z,t);
%   * there is NO wave-following / curvilinear transform — everything stays in
%     the lab frame (the "intrp_*" fields are returned equal to the lab fields,
%     and SU is NaN, purely to honour the decompose_longitudinal contract so the
%     caller's save block is unchanged).
%
% Two methods (opt.method):
%   'ymean' (default) — y-mean subtraction. w_wave is uniform in y, so the
%       instantaneous cross-wind mean captures wave + mean flow:
%           u' = u_raw - <u_raw>_y ,   w' = w_raw - <w_raw>_y .
%       Self-contained (needs no eta data). Note: this removes the MEAN as well
%       as the wave, so the residual is the turbulent fluctuation directly
%       (ensemble-averaging it gives ~0, not <u>). This is what the y-averaged
%       second moments <w'w'>, <v'v'>, <v'w'> need.
%   'eta'  — Hilbert/eta method (wave_removal.tex §"Compute w_wave"). Removes
%       only the wave, leaving mean+turb (like the longitudinal convention).
%       Requires the longitudinal surface line at the matched time:
%           opt.eta     eta(x) at this frame's time, metres (vector along x)
%           opt.eta_dx  x spacing of opt.eta, metres
%       Cross-wave orbital is 0; vertical orbital is wave_w_from_eta(...).
%
% I/O contract is identical to decompose_longitudinal (see that file's header).

if nargin < 6, opt = struct(); end
if ~isfield(opt, 'method'), opt.method = 'ymean'; end

DX = compVel.DX;  DT = compVel.DT;   %#ok<NASGU>  % DT kept for parity / eta scaling
[Nz, Nx] = size(w_raw);

pivRes.zPIV = compVel.zPIV;
pivRes.xPIV = compVel.xPIV;
pivRes.GS   = compVel.GS;
pivRes.mask = compVel.Mask;
pivRes.pf_surf = Surface_PIV;        % surface line (px), kept for reference

switch lower(opt.method)
    case 'ymean'
        % y-uniform component = instantaneous cross-wind mean (wave + mean flow)
        ubar = mean(u_raw, 2, 'omitnan');     % [Nz x 1]
        wbar = mean(w_raw, 2, 'omitnan');
        ORBX_ms = repmat(ubar, 1, Nx);        % no wave in u; this is mean cross-wind flow
        ORBZ_ms = repmat(wbar, 1, Nx);        % wave heave + mean

    case 'eta'
        if ~isfield(opt, 'eta') || ~isfield(opt, 'eta_dx')
            error('decompose_transverse:missingEta', ...
                  ['method=''eta'' needs opt.eta (eta(x) at this frame''s time, m) ' ...
                   'and opt.eta_dx (x spacing, m). See wave_w_from_eta.m.']);
        end
        % depth of each PIV row below the mean surface (m, <= 0 in water)
        surf_row = mean(Surface_PIV, 'omitnan');
        z_depth  = -(compVel.zPIV(:) - surf_row) * DX;   % [Nz x 1]
        w_wave   = wave_w_from_eta(opt.eta, opt.eta_dx, z_depth, opt);  % [Nz x 1]
        ORBX_ms  = zeros(Nz, Nx);             % no cross-wave orbital
        ORBZ_ms  = repmat(w_wave, 1, Nx);     % vertical orbital, uniform in y

    otherwise
        error('decompose_transverse:badMethod', ...
              'Unknown opt.method ''%s'' (use ''ymean'' or ''eta'').', opt.method);
end

% mask to the water region so subtraction inherits the velocity NaNs
ORBX_ms = ORBX_ms .* (compVel.Mask ./ compVel.Mask);   % NaN where Mask is NaN
ORBZ_ms = ORBZ_ms .* (compVel.Mask ./ compVel.Mask);

% residual (see method notes above for what "residual" means in each case)
u_res = u_raw - ORBX_ms;
w_res = w_raw - ORBZ_ms;

% --- pack (no wave-following frame: intrp_* and *_lab all equal the lab fields) ---
D.u_raw       = u_raw;        D.w_raw       = w_raw;
D.intrp_u_raw = u_raw;        D.intrp_w_raw = w_raw;       % no transform applied
D.ORBX_ms     = ORBX_ms;      D.ORBZ_ms     = ORBZ_ms;
D.intrp_u_res = u_res;        D.intrp_w_res = w_res;
D.u_res_lab   = u_res;        D.w_res_lab   = w_res;
D.u_orb_lab   = ORBX_ms;      D.w_orb_lab   = ORBZ_ms;
D.u_raw_lab   = u_raw;        D.w_raw_lab   = w_raw;
D.SU          = nan(Nz, Nx);  % no curvilinear grid in the transverse plane
end
