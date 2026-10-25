function [setup_b, B] = cisg_flat_bias(R, ref_b, setup_b, bin, frames)
%CISG_FLAT_BIAS Fixed slope bias from flat-water frames, subtracted from every frame by cisg_cache_frame.
%
%   [setup_b, B] = cisg_flat_bias(R, ref_b, setup_b, bin, frames)
%
% R, ref_b, setup_b, bin as in cisg_slope_video (cisg_setup_run + cisg_bin_setup).  The slopes of the
% flat frames (water still flat; default use 60:10:200, after the reference frames 2:4:58) against the
% reference are averaged per pixel: what is left is the fixed-pattern error of the reference / colour map,
% not waves.  setup_b.bias_Sx, bias_Sy (ROI-sized single) are set; cisg_cache_frame subtracts them.
% Pixels valid in fewer than half of the flat frames get no correction (0).  frames = [] does nothing.
%
% Cached in Slope_Gauge/calibration/flat_bias_CoreView_<N>_bin<B>.mat, reused while the frames, the
% reference frames and the ROI are the same.  B: bias_Sx, bias_Sy, frames, n_valid, rms (before / after).

    B = [];
    if isfield(setup_b, 'bias_Sx'), setup_b = rmfield(setup_b, {'bias_Sx', 'bias_Sy'}); end   % measure without a previous bias
    if isempty(frames), return; end
    sg = fileparts(fileparts(mfilename('fullpath')));
    run = sscanf(regexp(R.obs_folder, 'CoreView_\d+', 'match', 'once'), 'CoreView_%d');
    f = fullfile(sg, 'calibration', sprintf('flat_bias_CoreView_%d_bin%d.mat', run, bin));
    key = struct('frames', frames(:).', 'ref_frames', R.ref_frames(:).', 'roi', setup_b.roi, ...
                 'dye_model', setup_b.dye_model, 'L_eff_cm', setup_b.L_eff_cm, 'drift_corr', setup_b.drift_corr, ...
                 'sy_up', true);                                    % Sy sign fixed 2026-10-25: older bias files are recomputed
    if isfield(setup_b, 'colour_channels') && ~strcmpi(setup_b.colour_channels, 'RG')   % 'RGB' / 'HYBRID' get their own key
        key.colour_channels = upper(setup_b.colour_channels);
    end
    if isfile(f)
        S = load(f, 'B', 'key');
        if isequaln(S.key, key)
            B = S.B;
            fprintf('flat-water bias from %s (frames %s)\n', f, mat2str(frames));
        end
    end
    if isempty(B)
        fprintf('flat-water bias: %d frames %s ... ', numel(frames), mat2str(frames));  t0 = tic;
        r = setup_b.roi;  sz = [r.row_max - r.row_min + 1, r.col_max - r.col_min + 1];
        sx = zeros(sz);  sy = zeros(sz);  nv = zeros(sz);  rms0 = nan(numel(frames), 2);  k = 0;
        for n = frames
            try
                [a, b] = cisg_cache_frame(n, ref_b, setup_b, R.frame_file, '', bin, R.dt, R.img_size);
            catch ME
                if ~strcmp(ME.identifier, 'cisg:badFrame'), rethrow(ME); end
                continue;
            end
            k = k + 1;  rms0(k, :) = [rms(a, 'all', 'omitnan'), rms(b, 'all', 'omitnan')];
            v = isfinite(a) & isfinite(b);
            sx(v) = sx(v) + double(a(v));  sy(v) = sy(v) + double(b(v));  nv = nv + v;
        end
        if k == 0, error('cisg_flat_bias:none', 'No readable flat frame in %s', mat2str(frames)); end
        use = nv >= 0.5 * k;
        bx = zeros(sz, 'single');  by = bx;
        bx(use) = sx(use) ./ nv(use);  by(use) = sy(use) ./ nv(use);
        B = struct('bias_Sx', bx, 'bias_Sy', by, 'frames', frames, 'n_used', k, 'n_valid', uint16(nv), ...
                   'rms_flat', mean(rms0(1:k, :), 1), 'rms_bias', [rms(bx(use)), rms(by(use))], ...
                   'mean_bias', [mean(bx(use)), mean(by(use))]);
        save(f, 'B', 'key');
        fprintf('%.0f s\n', toc(t0));
    end
    fprintf(['flat frames vs reference: rms Sx %.4f, Sy %.4f; fixed bias: mean %+.4f / %+.4f, rms %.4f / %.4f ' ...
             '(subtracted from every frame)\n'], B.rms_flat, B.mean_bias, B.rms_bias);
    setup_b.bias_Sx = B.bias_Sx;  setup_b.bias_Sy = B.bias_Sy;
end
