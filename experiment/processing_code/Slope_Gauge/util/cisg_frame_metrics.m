function m = cisg_frame_metrics(n, ref_b, st, frame_file, bin, full_size, ycm_rows, yq, do_joint)
%CISG_FRAME_METRICS Stripe / profile statistics of one frame (no maps kept), for run comparisons.
%
%   m = cisg_frame_metrics(n, ref_b, st, frame_file, bin, full_size, ycm_rows, yq, do_joint)
%
% Slopes with the 2-channel colour inversion (and, if do_joint, the joint dye solve, whose dye map c is
% also summarised).  x = image columns (streamwise), y = rows.  Per method:
%   rms Sx, Sy over the valid pixels;
%   row-mean profiles of Sx, Sy over x (streamwise average) interpolated on the y grid yq (cm);
%   row variance profiles (variance along x per row);
%   stripe strength = variance over y of the row-mean profile / total variance of the field.
% ycm_rows: y (cm) of the ROI rows.  m.ok = false if the raw file is unreadable.

    m = struct('n', n, 'ok', false);
    try
        obs = cisg_load_coreview(frame_file(n), full_size);
    catch
        return;
    end
    if bin > 1, obs = uint8(imresize(obs, 1 / bin, 'box')); end
    r = st.roi;  rg = @(A) A(r.row_min:r.row_max, r.col_min:r.col_max);
    px2Sx = st.pattern_scale * st.cm_per_pixel_x / ((st.n_water - st.n_air) * st.L_eff_cm);
    px2Sy = st.pattern_scale * st.cm_per_pixel_y / ((st.n_water - st.n_air) * st.L_eff_cm);
    [dx, dy] = cisg_colour_displacement(ref_b, obs, st);
    m = stats(m, '', rg(dx) * px2Sx, rg(dy) * px2Sy, ycm_rows, yq);
    m.masked = mean(isnan(rg(dx)), 'all');
    if do_joint
        [dxj, dyj, cm] = cisg_dye_unmix(ref_b, obs, st, dx, dy);
        m = stats(m, 'j', rg(dxj) * px2Sx, rg(dyj) * px2Sy, ycm_rows, yq);
        c = rg(cm);  rc = mean(c, 2, 'omitnan');
        m.c_mean = mean(c(:), 'omitnan');  m.c_p99 = prctile(c(isfinite(c)), 99);
        m.c_row = interp1(ycm_rows, rc, yq, 'linear', NaN);
        ok = isfinite(rc) & isfinite(mean(rg(dx) * px2Sx, 2, 'omitnan')) & isfinite(mean(rg(dy) * px2Sy, 2, 'omitnan'));
        m.corr_c_Sx = corr(rc(ok), mean(rg(dx) * px2Sx, 2, 'omitnan') * 0 + mean(rg(dx(ok, :)) * 0, 2) + mean(rg(dx) * px2Sx, 2, 'omitnan') .* 0 + rowmean_ok(rg(dx) * px2Sx, ok));
        m.corr_c_Sy = corr(rc(ok), rowmean_ok(rg(dy) * px2Sy, ok));
        m.corr_c_Sxj = corr(rc(ok), rowmean_ok(rg(dxj) * px2Sx, ok));
        m.corr_c_Syj = corr(rc(ok), rowmean_ok(rg(dyj) * px2Sy, ok));
    end
    m.ok = true;
end

function v = rowmean_ok(A, ok)
    v = mean(A, 2, 'omitnan');  v = v(ok);
end

function m = stats(m, tag, Sx, Sy, ycm_rows, yq)
    names = {'Sx', 'Sy'};  F = {Sx, Sy};
    for k = 1:2
        A = F{k};
        m.(['rms' names{k} tag]) = rms(A(:), 'omitnan');
        rmn = mean(A, 2, 'omitnan');  rvr = var(A, 0, 2, 'omitnan');
        v = isfinite(rmn);
        m.(['str' names{k} tag]) = var(rmn(v)) / var(A(:), 'omitnan');
        m.(['rm' names{k} tag]) = interp1(ycm_rows, rmn, yq, 'linear', NaN);
        m.(['rv' names{k} tag]) = interp1(ycm_rows, rvr, yq, 'linear', NaN);
    end
end
