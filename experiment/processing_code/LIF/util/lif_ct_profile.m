function P = lif_ct_profile(density_file, pyc_frac)
%LIF_CT_PROFILE Pre-run CT profile and its pycnocline, on the LIF depth axis.
%
%   P = LIF_CT_PROFILE(density_file)          pyc_frac 0.5
%   P = LIF_CT_PROFILE(density_file, pyc_frac)
%
% density_file - densityProfileData_*.mat (depth cm positive down, voltage);
%                [] or '' -> P = [] (unstratified run)
% P.z   - height, cm, negative below the surface (sorted ascending)
% P.v   - CT conductivity voltage (V) at P.z
% P.pyc - [bottom top] of the pycnocline, cm: where |dV/dz| is above
%         pyc_frac of its maximum (the density gradient refracts the
%         backlight there, so it is marked on the depth plots)
% Shared by compare_depth_runs.m and the side-view video (surf_theta_video).
    P = [];
    if isempty(density_file), return; end
    if nargin < 2 || isempty(pyc_frac), pyc_frac = 0.5; end
    S = load(density_file, 'depth', 'voltage');
    [zp, o] = sort(-S.depth(:));                       % -> negative below the surface
    vp = S.voltage(o);  vp = vp(:);
    g = abs(gradient(vp, zp));
    pyc = zp(g >= pyc_frac * max(g));
    P = struct('z', zp, 'v', vp, 'pyc', [min(pyc) max(pyc)], 'file', density_file);
end
