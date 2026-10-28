function [snr, th_ml, sig] = lif_dye_snr(theta_img, z_cm, z_surf, zeta_ml)
%LIF_DYE_SNR Per-pixel signal-to-noise of the dye layer in one theta frame.
%
%   [snr, th_ml, sig] = LIF_DYE_SNR(theta_img, z_cm, z_surf, zeta_ml)
%
% th_ml - median theta over zeta_ml (cm below the local surface, default
%         [-5 -2]: under the surface band, inside the mixed layer late in
%         the run)
% sig   - pixel noise there, from neighbour differences along x
%         (1.4826 * MAD / sqrt(2)): measured inside the dye itself, so it
%         needs no clear water and the dye's own level does not bias it
% snr   = th_ml / sig
%
% 118: ~14 (z99(x) line clean with a 1 cm window); 120: ~4 (needed 3 cm).

    if nargin < 4, zeta_ml = [-5 -2]; end
    z = z_cm(:);
    v = [];  d = [];
    for j = find(~isnan(z_surf(:).'))
        r = z - z_surf(j) >= zeta_ml(1) & z - z_surf(j) <= zeta_ml(2);
        v = [v; theta_img(r, j)]; %#ok<AGROW>
    end
    for j = find(~isnan(z_surf(:).') & [~isnan(z_surf(2:end)) false])
        r = z - z_surf(j) >= zeta_ml(1) & z - z_surf(j) <= zeta_ml(2);
        d = [d; theta_img(r, j+1) - theta_img(r, j)]; %#ok<AGROW>
    end
    th_ml = median(v, 'omitnan');
    d   = d(~isnan(d));
    sig = 1.4826 * median(abs(d - median(d))) / sqrt(2);
    snr = th_ml / sig;
end
