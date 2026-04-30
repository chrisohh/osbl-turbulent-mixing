function stats = compute_horizontal_avg_turbstats(v_turb, w_turb, z, t)

[Nt, Ny, Nz] = size(v_turb);

% Preallocate
v_mean = zeros(Nz, Nt);
w_mean = zeros(Nz, Nt);
vv     = zeros(Nz, Nt);
ww     = zeros(Nz, Nt);
vw     = zeros(Nz, Nt);
www    = zeros(Nz, Nt);

for it = 1:Nt
    for iz = 1:Nz

        v_slice = squeeze(v_turb(it,:,iz));
        w_slice = squeeze(w_turb(it,:,iz));

        % Means
        vbar = mean(v_slice);
        wbar = mean(w_slice);

        v_mean(iz,it) = vbar;
        w_mean(iz,it) = wbar;

        % Fluctuations
        vp = v_slice - vbar;
        wp = w_slice - wbar;

        % Moments
        vv(iz,it)  = mean(vp.^2);
        ww(iz,it)  = mean(wp.^2);
        vw(iz,it)  = mean(vp .* wp);
        www(iz,it) = mean(wp.^3);

    end
end

stats.v_mean = v_mean;
stats.w_mean = w_mean;
stats.vv     = vv;
stats.ww     = ww;   % ← KEY
stats.vw     = vw;
stats.www    = www;
stats.tke    = 0.5 * (vv + ww);

stats.z = z;
stats.t = t;

end