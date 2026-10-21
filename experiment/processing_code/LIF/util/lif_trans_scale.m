function [mmpp, c0, info] = lif_trans_scale(fplate, opts)
%LIF_TRANS_SCALE mm-per-pixel from a plate image that is NOT parallel to the sensor.
%
%   [mmpp, c0, info] = LIF_TRANS_SCALE(fplate, opts)
%
% lif_plate_scale fits a UNIFORM scale (nonreflectivesimilarity), which is
% only right when the plate faces the camera.  The Water_PIV plate
% (CoreView_91) does not: it is a vertical strip rotated about a VERTICAL
% axis, so its squares image as tall rectangles --
%     vertical step  173.1 px,  horizontal step 88.0 px  (ratio 0.51),
% i.e. ~59 deg of yaw.  A uniform fit would average the two and be ~40%
% wrong.
%
% Rotation about a vertical axis leaves VERTICAL lengths unforeshortened,
% so the vertical corner step is the one usable scale, and it is the one
% that matters here: depths (z99, the concentration profile) are vertical.
% The horizontal (cross-tank) axis of the transverse view carries the same
% mm/px ONLY if the laser sheet is normal to the optical axis; it is drawn
% with this scale and should be read as approximate -- see the warning
% printed when info.foreshorten is far from 1.
%
% As in lif_plate_scale, the plate may sit opts.offset_mm TOWARD the
% camera from the measurement plane, which is then farther away:
%     mmpp = mmpp_plate * (D + d) / D,  D = f_px * mmpp_plate (pinhole),
%     d = offset_mm / n_medium (apparent shift seen through water).
%
% opts fields (defaults):
%   square_mm 11, offset_mm 0, n_medium 1.33, lens_f_mm 35,
%   pixel_um 5.5, verbose true
%
% Outputs:
%   mmpp - mm per pixel along the VERTICAL image axis, on the measurement
%          plane (full-resolution pixels)
%   c0   - [col row] of the detected board centre
%   info - struct: square_px_v, square_px_h, foreshorten (h/v), yaw_deg,
%          mmpp_plate, D_mm, boardSize, imgPts

    def = struct('square_mm', 11, 'offset_mm', 0, 'n_medium', 1.33, ...
                 'lens_f_mm', 35, 'pixel_um', 5.5, 'verbose', true);
    if nargin < 2, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end

    P   = lif_load_raw(fplate);
    lim = prctile(P(:), [0.5 99.5]);
    P8  = uint8(255 * min(max((P - lim(1)) / diff(lim), 0), 1));

    % Full (non-partial) detection only: it returns a complete rectangular
    % grid, which is what the step statistics below need.  The board is
    % only partly in frame, so the detected size is smaller than the
    % printed one -- that is fine, the SQUARE is what sets the scale.
    [imgPts, boardSize] = detectCheckerboardPoints(P8, 'PartialDetections', false);
    if isempty(imgPts) || any(isnan(imgPts(:))) || any(boardSize < 3)
        error('lif_trans_scale:detect', 'Checkerboard not detected in %s', fplate);
    end

    g  = reshape(imgPts, boardSize(1)-1, boardSize(2)-1, 2);   % [nr x nc x (x,y)]
    s1 = squeeze(mean(diff(g, 1, 1), [1 2]));                  % mean step along dim 1
    s2 = squeeze(mean(diff(g, 1, 2), [1 2]));                  % mean step along dim 2
    if abs(s1(2)) > abs(s2(2))                                 % the more vertical one
        sv = norm(s1);  sh = norm(s2);
    else
        sv = norm(s2);  sh = norm(s1);
    end

    mmpp_plate = opts.square_mm / sv;
    D_mm = (opts.lens_f_mm / (opts.pixel_um/1000)) * mmpp_plate;
    d_mm = opts.offset_mm / opts.n_medium;
    mmpp = mmpp_plate * (D_mm + d_mm) / D_mm;
    c0   = mean(imgPts, 1);

    info = struct('square_px_v', sv, 'square_px_h', sh, 'foreshorten', sh/sv, ...
                  'yaw_deg', acosd(min(sh/sv, 1)), 'mmpp_plate', mmpp_plate, ...
                  'D_mm', D_mm, 'boardSize', boardSize, 'imgPts', imgPts);

    if opts.verbose
        fprintf('Plate %s: %d corners, %dx%d grid\n', fplate, size(imgPts,1), boardSize(1), boardSize(2));
        fprintf('  square %g mm: %.2f px vertical, %.2f px horizontal (ratio %.2f -> yaw ~%.0f deg)\n', ...
                opts.square_mm, sv, sh, info.foreshorten, info.yaw_deg);
        fprintf('  vertical scale %.5f mm/px; camera-plate ~%.0f mm; offset %.1f mm (apparent %.1f) -> %.5f mm/px (+%.2f%%)\n', ...
                mmpp_plate, D_mm, opts.offset_mm, d_mm, mmpp, 100*d_mm/D_mm);
        if info.foreshorten < 0.9
            warning(['lif_trans_scale: the plate is yawed ~%.0f deg, so only the VERTICAL scale ' ...
                     'is trustworthy.  Depths are calibrated; the cross-tank (y) axis uses the ' ...
                     'same mm/px and is approximate until a plate is imaged IN the laser sheet.'], ...
                    info.yaw_deg);
        end
    end
end
