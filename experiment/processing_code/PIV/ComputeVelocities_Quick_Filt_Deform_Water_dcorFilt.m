function [CompVel] = ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt(PIV1, PIV2, Mask1, Mask2, IntrWndw, GrdSpc, keep_best, inter_uod_in)
% ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt
%
%   Water-side PIV with image deformation plus four improvements from the
%   FAB6 algorithm family and optional RECTANGULAR interrogation windows.
%
%   RECTANGULAR WINDOWS  (new)
%     Pass IntrWndw and GrdSpc as Nx2 matrices, one row per pyramid level:
%         IntrWndw = [IW_x_1  IW_z_1;       e.g.  [128 32;
%                     IW_x_2  IW_z_2;               64 16;
%                        ...    ...  ]               32 16;
%                                                    16  8;
%                                                     8  4]
%     Wide x-windows capture large horizontal displacements; shallow
%     z-windows minimise air/interface inclusion near the surface.
%     For square windows (backward compatible) pass Nx1 vectors as before.
%
%   [1] Cascade keep-best
%         When a new level's peak is no better than the prior, the prior
%         displacement is carried forward instead of storing NaN.
%         Prevents a bad near-surface window from collapsing the chain.
%
%   [2] Incremental dcor gate
%         Accept new estimate only when:  dcor_new >= dcor_prior - 0.5/lvl
%         Tolerance tightens toward the finest level (lvl=4: tol=0.125).
%         Makes going to IW=4 safe — empty windows lose the gate.
%
%   [3] Decaying smoothing  S = 1/lvl/10  at global levels
%         Level 1: S=0.1 (heavy gap-fill for deformation).
%         Level N-1: S=1/(N-1)/10 (light, preserves structure).
%         Final level: dcor-weighted smoothn with S=0 (see [4]).
%
%   [4] dcor-weighted final smoothn
%         smoothn(delx, dcor_w, 0, 'robust') — high-dcor vectors
%         preserved; low-dcor vectors filled from good neighbours.
%
%   Mask convention: 1 = water (valid), NaN = air.
%   Accepts both NaN/1 and 0/1 input masks — normalised internally.
%   Output CompVel.Mask is always NaN/1.
%
%   Drop-in replacement for ComputeVelocities_Quick_NoFilt_Deform_Water.
%   Identical output struct fields (plus IW_x, IW_z, GS_x, GS_z).
%
%   Square example:
%       IntrWndw = [64 32 16 8];   GrdSpc = [32 16 8 4];
%   Rectangular example:
%       IntrWndw = [128 64 32 16 8; 32 16 16 8 4]';   % Nx2
%       GrdSpc   = [64  32 16  8 4; 16  8  8 4 2]';
%
%   Original algorithm: Fabrice Veron / Marc Buckley, Nov 2013.
%   dcorFilt + rectangular modifications: 2025.

% keep_best: true  = carry prior displacement when new dcor fails gate (original)
%            false = set to NaN when gate fails; let validatePIV fill gaps (default)
if nargin < 7, keep_best = false; end

% inter_uod: optional struct to apply UOD between pyramid passes before deformation
%   .enabled     false   set true to activate
%   .remove      2.0     normalised residual threshold for removal
%   .reinsert    3.0     normalised residual threshold for reinsertion
%   .minvec      5       minimum valid neighbours required (conservative near surface)
%   .eps         0.1     noise floor to prevent div-by-zero
inter_uod.enabled  = false;
inter_uod.remove   = 2.0;
inter_uod.reinsert = 3.0;
inter_uod.minvec   = 5;
inter_uod.eps      = 0.1;
if nargin >= 8 && isstruct(inter_uod_in)
    for f = fieldnames(inter_uod_in)'
        inter_uod.(f{1}) = inter_uod_in.(f{1});
    end
end

%% ---- Parse pyramid — support square (Nx1) and rectangular (Nx2) ----
IntrWndw = IntrWndw(:);   % ensure column if square
GrdSpc   = GrdSpc(:);

if size(IntrWndw, 2) == 2
    IWx = IntrWndw(:, 1);   % horizontal window sizes per level
    IWz = IntrWndw(:, 2);   % vertical   window sizes per level
    GSx = GrdSpc(:, 1);
    GSz = GrdSpc(:, 2);
else
    IWx = IntrWndw;          % square: same in both dimensions
    IWz = IntrWndw;
    GSx = GrdSpc;
    GSz = GrdSpc;
end
nLevels = length(IWx);

%% ---- Preprocessing ----
% Normalise mask to NaN/1 (handles 0/1 or NaN/1 input).
Mask1 = double(Mask1);  Mask1(Mask1 == 0) = NaN;
Mask2 = double(Mask2);  Mask2(Mask2 == 0) = NaN;

PIV1 = double(PIV1) .* Mask1;
PIV2 = double(PIV2) .* Mask2;
PIV1(isnan(PIV1)) = nanmean(PIV1(:));
PIV2(isnan(PIV2)) = nanmean(PIV2(:));
IM1_D = PIV1;
IM2_D = PIV2;

[h, w]   = size(PIV1);
[X1, Y1] = meshgrid(1:w, 1:h);

%%
% ======================= GLOBAL LEVELS =======================
for lvl = 1:nLevels-1

    iw_x = IWx(lvl);   % horizontal window size this level
    iw_z = IWz(lvl);   % vertical   window size this level
    gs_x = GSx(lvl);
    gs_z = GSz(lvl);

    x = iw_x/2 : gs_x : (w - iw_x/2);   % grid centre coords (horizontal)
    y = iw_z/2 : gs_z : (h - iw_z/2);   % grid centre coords (vertical)

    bxsNh = floor(1 + (h - iw_z) / gs_z);
    bxsNw = floor(1 + (w - iw_x) / gs_x);

    delx = NaN(bxsNh, bxsNw);
    dely = NaN(bxsNh, bxsNw);
    dcor = NaN(bxsNh, bxsNw);

    if lvl == 1
        pdelx = zeros(bxsNh, bxsNw);
        pdely = zeros(bxsNh, bxsNw);
        pdcor = zeros(bxsNh, bxsNw);   % prior dcor = 0 → gate always passes at lvl 1
    end

    dcor_tol = 0.5 / lvl;   % [2] loosest at coarse, tightest at fine

    bxCNTc = 1;
    for c = x

        bxCNTr = 1;
        bdryL = c - iw_x/2 + 1;
        bdryR = c + iw_x/2;

        for r = y

            bdryT = r - iw_z/2 + 1;
            bdryB = r + iw_z/2;

            VALID_BLOCK = ~isnan(Mask1(r, c));

            if VALID_BLOCK
                try
                    bxA   = IM1_D(bdryT:bdryB, bdryL:bdryR);   % iw_z × iw_x
                    bxB   = IM2_D(bdryT:bdryB, bdryL:bdryR);
                    bxAmm = bxA - mean(bxA(:));
                    bxBmm = bxB - mean(bxB(:));

                    fftCorr = fft2(bxBmm) .* conj(fft2(bxAmm));
                    Xcorr   = fftshift(real(ifft2(fftCorr))) ...
                              ./ sqrt(sum(bxAmm(:).^2)) ...
                              ./ sqrt(sum(bxBmm(:).^2));
                    % Xcorr size is iw_z × iw_x; zero-displacement at centre.
                    % For rectangular FFT: peak row→z-disp, peak col→x-disp.

                    [Xpky, Xpkx] = find(Xcorr == max(Xcorr(:)));

                    ldelx = Xpkx - iw_x/2 - 1;
                    ldely = Xpky - iw_z/2 - 1;

                    % 3-point Gaussian subpixel (3×3 patch, same for any IW shape)
                    T = log(Xcorr(Xpky-1:Xpky+1, Xpkx-1:Xpkx+1));
                    t = T(:, 2);
                    SubpixelY = (1/2)*(t(3)-t(1)) / (2*t(2)-t(1)-t(3));
                    t = T(2, :);
                    SubpixelX = (1/2)*(t(3)-t(1)) / (2*t(2)-t(1)-t(3));

                    new_dcor  = max(Xcorr(:));
                    dcor_prev = pdcor(bxCNTr, bxCNTc);

                    subpix_ok = isreal([SubpixelY SubpixelX])  && ...
                                abs(SubpixelY) < 1             && ...
                                abs(SubpixelX) < 1             && ...
                                abs(ldelx)     < iw_x/2        && ...
                                abs(ldely)     < iw_z/2;

                    if subpix_ok && (new_dcor >= dcor_prev - dcor_tol)
                        delx(bxCNTr, bxCNTc) = pdelx(bxCNTr, bxCNTc) + ldelx + SubpixelX;
                        dely(bxCNTr, bxCNTc) = pdely(bxCNTr, bxCNTc) + ldely + SubpixelY;
                        dcor(bxCNTr, bxCNTc) = new_dcor;
                    elseif keep_best
                        delx(bxCNTr, bxCNTc) = pdelx(bxCNTr, bxCNTc);
                        dely(bxCNTr, bxCNTc) = pdely(bxCNTr, bxCNTc);
                        dcor(bxCNTr, bxCNTc) = dcor_prev;
                    else
                        delx(bxCNTr, bxCNTc) = NaN;
                        dely(bxCNTr, bxCNTc) = NaN;
                        dcor(bxCNTr, bxCNTc) = NaN;
                    end

                catch
                    if keep_best
                        delx(bxCNTr, bxCNTc) = pdelx(bxCNTr, bxCNTc);
                        dely(bxCNTr, bxCNTc) = pdely(bxCNTr, bxCNTc);
                        dcor(bxCNTr, bxCNTc) = pdcor(bxCNTr, bxCNTc);
                    else
                        delx(bxCNTr, bxCNTc) = NaN;
                        dely(bxCNTr, bxCNTc) = NaN;
                        dcor(bxCNTr, bxCNTc) = NaN;
                    end
                end
            end

            bxCNTr = bxCNTr + 1;
        end

        bxCNTc = bxCNTc + 1;
    end

    % Inter-pass UOD (optional) — clean displacement field before deformation
    if inter_uod.enabled
        mk_lvl = double(~isnan(delx));   % 1=valid, NaN=air at this grid
        mk_lvl(mk_lvl == 0) = NaN;
        [delx, dely] = interpass_uod(delx, dely, mk_lvl, inter_uod);
    end

    % [3] Decaying smoothing
    S_lvl   = 1 / lvl / 10;
    delx_sm = smoothn(delx, S_lvl, 'robust');
    dely_sm = smoothn(dely, S_lvl, 'robust');

    [X, Y] = meshgrid(x, y);

    % Upsample smoothed field to image resolution for image deformation
    U1 = interp2(X, Y, delx_sm, X1, Y1, '*spline');
    V1 = interp2(X, Y, dely_sm, X1, Y1, '*spline');

    % Warp second image (one-sided, same as _NoFilt_Deform_Water)
    IM1_D = interp2(1:w, (1:h)', PIV1, X1,      Y1,      '*linear');
    IM2_D = interp2(1:w, (1:h)', PIV2, X1 + U1, Y1 + V1, '*linear');

    % Downsample to next level's grid
    iw2_x = IWx(lvl + 1);  iw2_z = IWz(lvl + 1);
    gs2_x = GSx(lvl + 1);  gs2_z = GSz(lvl + 1);
    x2    = iw2_x/2 : gs2_x : (w - iw2_x/2);
    y2    = iw2_z/2 : gs2_z : (h - iw2_z/2);
    pdelx = U1(y2, x2);
    pdely = V1(y2, x2);

    % Interpolate raw dcor to next level's grid (NaN→0: no prior quality)
    dcor_fill = dcor;
    dcor_fill(isnan(dcor_fill)) = 0;
    [X2, Y2] = meshgrid(x2, y2);
    pdcor = interp2(X, Y, dcor_fill, X2, Y2, '*linear', 0);

end
% ======================= END GLOBAL LEVELS =======================

%%
% ======================= FINAL (LOCAL) LEVEL =======================
lvl   = nLevels;
iw_x  = IWx(lvl);
iw_z  = IWz(lvl);
gs_x  = GSx(lvl);
gs_z  = GSz(lvl);

x = iw_x/2 : gs_x : (w - iw_x/2);
y = iw_z/2 : gs_z : (h - iw_z/2);

bxsNh = floor(1 + (h - iw_z) / gs_z);
bxsNw = floor(1 + (w - iw_x) / gs_x);

delx = pdelx;                      % start from best prior estimate
dely = pdely;
dcor = NaN(bxsNh, bxsNw);
MASK = NaN(bxsNh, bxsNw);

dcor_tol = 0.5 / lvl;

bxCNTc = 1;
for c = x

    bxCNTr = 1;
    bdryL = c - iw_x/2 + 1;
    bdryR = c + iw_x/2;

    for r = y

        bdryT = r - iw_z/2 + 1;
        bdryB = r + iw_z/2;

        VALID_BLOCK = ~isnan(Mask1(r, c));

        if VALID_BLOCK
            MASK(bxCNTr, bxCNTc) = 1;
            try
                bxA   = IM1_D(bdryT:bdryB, bdryL:bdryR);
                bxB   = IM2_D(bdryT:bdryB, bdryL:bdryR);
                bxAmm = bxA - mean(bxA(:));
                bxBmm = bxB - mean(bxB(:));

                fftCorr = fft2(bxBmm) .* conj(fft2(bxAmm));
                Xcorr   = fftshift(real(ifft2(fftCorr))) ...
                          ./ sqrt(sum(bxAmm(:).^2)) ...
                          ./ sqrt(sum(bxBmm(:).^2));

                [Xpky, Xpkx] = find(Xcorr == max(Xcorr(:)));

                ldelx = Xpkx - iw_x/2 - 1;
                ldely = Xpky - iw_z/2 - 1;

                T = log(Xcorr(Xpky-1:Xpky+1, Xpkx-1:Xpkx+1));
                t = T(:, 2);
                SubpixelY = (1/2)*(t(3)-t(1)) / (2*t(2)-t(1)-t(3));
                t = T(2, :);
                SubpixelX = (1/2)*(t(3)-t(1)) / (2*t(2)-t(1)-t(3));

                new_dcor  = max(Xcorr(:));
                dcor_prev = pdcor(bxCNTr, bxCNTc);

                subpix_ok = isreal([SubpixelY SubpixelX]) && ...
                            abs(SubpixelY) < 1            && ...
                            abs(SubpixelX) < 1;

                if subpix_ok && (new_dcor >= dcor_prev - dcor_tol)
                    delx(bxCNTr, bxCNTc) = pdelx(bxCNTr, bxCNTc) + ldelx + SubpixelX;
                    dely(bxCNTr, bxCNTc) = pdely(bxCNTr, bxCNTc) + ldely + SubpixelY;
                    dcor(bxCNTr, bxCNTc) = new_dcor;
                else
                    % [1] Cascade keep-best
                    delx(bxCNTr, bxCNTc) = pdelx(bxCNTr, bxCNTc);
                    dely(bxCNTr, bxCNTc) = pdely(bxCNTr, bxCNTc);
                    dcor(bxCNTr, bxCNTc) = dcor_prev;
                end

            catch
                delx(bxCNTr, bxCNTc) = pdelx(bxCNTr, bxCNTc);
                dely(bxCNTr, bxCNTc) = pdely(bxCNTr, bxCNTc);
                dcor(bxCNTr, bxCNTc) = pdcor(bxCNTr, bxCNTc);
            end
        end

        bxCNTr = bxCNTr + 1;
    end

    bxCNTc = bxCNTc + 1;
end

% [4] dcor-weighted final smoothn (S=0: gap-fill only, no blurring)
[X, Y] = meshgrid(x, y);
dcor_w  = dcor;
dcor_w(isnan(dcor_w)) = 0;
delx_sm = smoothn(delx, dcor_w, 0.4, 'robust');
dely_sm = smoothn(dely, dcor_w, 0.4, 'robust');

U1 = interp2(X, Y, delx_sm, X1, Y1, '*spline');
V1 = interp2(X, Y, dely_sm, X1, Y1, '*spline');
% ======================= END FINAL LEVEL =======================

%%
% Output — same fields as _Quick_NoFilt_Deform_Water (drop-in compatible)
% plus IW_x / IW_z / GS_x / GS_z for rectangular-window callers.
CompVel.INTdelx = pdelx;
CompVel.INTdelz = -pdely;
CompVel.dcor    = dcor;    % NaN=air; measured or best-prior where water
CompVel.Mask    = MASK;    % NaN/1  (NaN=air/uncomputed, 1=valid water)
CompVel.mask    = MASK;

CompVel.delta_x  = delx;   % positive = downwind
CompVel.delta_z  = -dely;  % positive = upward (away from surface)
CompVel.delx     = delx;
CompVel.dely     = -dely;

CompVel.delta_x1 = U1;     % dcor-smoothed at original image resolution
CompVel.delta_z1 = -V1;

CompVel.xPIV = x;
CompVel.zPIV = y;

% Scalar IW/GS kept for backward compatibility (x-dimension value)
CompVel.IW   = iw_x;
CompVel.GS   = gs_x;

% Explicit per-dimension values for rectangular-window callers
CompVel.IW_x = iw_x;
CompVel.IW_z = iw_z;
CompVel.GS_x = gs_x;
CompVel.GS_z = gs_z;

end

% =========================================================================
function [u, w] = interpass_uod(u, w, mk, p)
% UOD between pyramid passes — cleans displacement field before image deformation.
% Conservative: skips vectors with fewer than p.minvec valid neighbours.
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

        if numel(vals_u) < p.minvec, continue, end

        med_u = median(vals_u);
        med_w = median(vals_w);
        rms_u = median(abs(vals_u - med_u)) + p.eps;
        rms_w = median(abs(vals_w - med_w)) + p.eps;

        if ~isnan(u(r,c))
            if abs(u(r,c) - med_u)/rms_u > p.remove || ...
               abs(w(r,c) - med_w)/rms_w > p.remove
                u(r,c) = NaN;
                w(r,c) = NaN;
            end
        elseif ~isnan(u0(r,c))
            if abs(u0(r,c) - med_u)/rms_u < p.reinsert && ...
               abs(w0(r,c) - med_w)/rms_w < p.reinsert
                u(r,c) = u0(r,c);
                w(r,c) = w0(r,c);
            end
        end
    end
end
end
