function CompVel = validatePIV(CompVel, varargin)
% validatePIV  Post-process a CompVel struct from ComputeVelocities_*_dcorFilt.
%
%   CompVel = validatePIV(CompVel)          uses default thresholds
%   CompVel = validatePIV(CompVel, opts)    override any field of opts
%
%   Pipeline:
%     1. dcor quality gate  (proxy for peak-ratio Q < 1.5)
%     2. Universal Outlier Detection (Westerweel & Scarano 2005)
%           5x5 neighbourhood, remove if normalised residual > 2,
%           reinsert if < 3, minimum 4 valid neighbours, 2 passes
%     3. Remove isolated groups  (< 5 vectors, 4-connectivity)
%     4. Denoising: 5x5 Gaussian smooth
%     5. Interpolate remaining gaps (water region only)
%
%   opts fields and defaults:
%     .do_dcor       false  enable dcor quality gate
%     .dcor_min      []     threshold for dcor gate ([] = skip even if do_dcor true)
%     .do_uod        true   enable Universal Outlier Detection
%     .uod_remove    2.0    UOD: remove  if normalised residual > this
%     .uod_reinsert  3.0    UOD: reinsert if normalised residual < this
%     .uod_minvec    4      UOD: skip if fewer valid neighbours than this
%     .uod_eps       0.1    UOD: noise floor (px) prevents div-by-zero
%     .uod_passes    2      number of UOD passes
%     .do_groups     true   remove isolated groups smaller than min_group
%     .min_group     5      minimum connected-component size (vectors)
%     .denoise       false  apply 5x5 Gaussian smooth (off: biases near-surface vectors)
%     .fill_gaps     true   interpolate NaN gaps after validation

opts.do_dcor      = false;
opts.dcor_min     = [];
opts.do_uod       = true;
opts.uod_remove   = 2.0;
opts.uod_reinsert = 3.0;
opts.uod_minvec   = 4;
opts.uod_eps      = 0.1;
opts.uod_passes   = 2;
opts.do_groups    = true;
opts.min_group    = 5;
opts.denoise      = false;
opts.fill_gaps    = true;

if nargin > 1 && isstruct(varargin{1})
    in = varargin{1};
    for f = fieldnames(in)'
        opts.(f{1}) = in.(f{1});
    end
end

u  = CompVel.delta_x;
w  = CompVel.delta_z;
qc = CompVel.dcor;
mk = CompVel.Mask;      % NaN = air/invalid, 1 = water

% --- 1. dcor gate ---
if opts.do_dcor && ~isempty(opts.dcor_min)
    bad = (qc < opts.dcor_min) & (mk == 1);
    u(bad) = NaN;
    w(bad) = NaN;
end

% --- 2. Universal Outlier Detection ---
if opts.do_uod
    for pass = 1:opts.uod_passes
        [u, w] = uod(u, w, mk, opts);
    end
end

% --- 3. Remove small isolated groups ---
if opts.do_groups
    valid = ~isnan(u) & (mk == 1);
    CC    = bwconncomp(valid, 4);
    for k = 1:CC.NumObjects
        if numel(CC.PixelIdxList{k}) < opts.min_group
            u(CC.PixelIdxList{k}) = NaN;
            w(CC.PixelIdxList{k}) = NaN;
        end
    end
end

% --- 4. Denoising: 5x5 Gaussian smooth ---
if opts.denoise
    gk = gaussianKernel5();
    u  = smoothField(u, mk, gk);
    w  = smoothField(w, mk, gk);
end

% --- 5. Interpolate gaps (water region only) ---
if opts.fill_gaps
    [nr, nc] = size(u);
    [xi, yi] = meshgrid(1:nc, 1:nr);
    water = (mk == 1);
    good  = water & ~isnan(u);
    if sum(good(:)) > 4
        u(water) = griddata(xi(good), yi(good), u(good), xi(water), yi(water), 'linear');
        w(water) = griddata(xi(good), yi(good), w(good), xi(water), yi(water), 'linear');
    end
end

CompVel.delta_x = u;
CompVel.delta_z = w;
CompVel.delx    = u;
CompVel.dely    = w;

end

% =========================================================================
function [u, w] = uod(u, w, mk, opts)
% Universal Outlier Detection — Westerweel & Scarano (2005).
u0   = u;
w0   = w;
half = 2;   % 5x5 neighbourhood
[nr, nc] = size(u);

for r = 1+half : nr-half
    for c = 1+half : nc-half
        if mk(r,c) ~= 1, continue, end

        nb_u = u(r-half:r+half, c-half:c+half);
        nb_w = w(r-half:r+half, c-half:c+half);
        nb_u(half+1, half+1) = NaN;
        nb_w(half+1, half+1) = NaN;

        vals_u = nb_u(~isnan(nb_u));
        vals_w = nb_w(~isnan(nb_w));

        if numel(vals_u) < opts.uod_minvec, continue, end

        med_u = median(vals_u);
        med_w = median(vals_w);
        rms_u = median(abs(vals_u - med_u)) + opts.uod_eps;
        rms_w = median(abs(vals_w - med_w)) + opts.uod_eps;

        if ~isnan(u(r,c))
            if abs(u(r,c) - med_u)/rms_u > opts.uod_remove || ...
               abs(w(r,c) - med_w)/rms_w > opts.uod_remove
                u(r,c) = NaN;
                w(r,c) = NaN;
            end
        elseif ~isnan(u0(r,c))
            if abs(u0(r,c) - med_u)/rms_u < opts.uod_reinsert && ...
               abs(w0(r,c) - med_w)/rms_w < opts.uod_reinsert
                u(r,c) = u0(r,c);
                w(r,c) = w0(r,c);
            end
        end
    end
end
end

% =========================================================================
function gk = gaussianKernel5()
[x, y] = meshgrid(-2:2, -2:2);
gk = exp(-(x.^2 + y.^2) / 2);
gk = gk / sum(gk(:));
end

% =========================================================================
function field = smoothField(field, mk, kernel)
water    = (mk == 1);
num      = conv2(field .* water, kernel, 'same');
den      = conv2(double(water & ~isnan(field)), kernel, 'same');
smoothed = num ./ max(den, 1e-6);
field(water & ~isnan(field)) = smoothed(water & ~isnan(field));
end
