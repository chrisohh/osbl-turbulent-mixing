function [mmpp, c0, info] = lif_plate_scale(fplate, opts)
%LIF_PLATE_SCALE mm-per-pixel of a camera from one checkerboard-plate image.
%
%   [mmpp, c0, info] = LIF_PLATE_SCALE(fplate, opts)
%
% No rectification: fits a similarity transform (uniform scale + rotation
% + shift) to the detected corners, so the image can be shown as recorded
% with physical axes.  The plate is assumed to sit opts.offset_mm TOWARD
% the camera from the measurement plane, which is therefore farther away:
%     mmpp = mmpp_plate * (D + d) / D,   D = f_px * mmpp_plate (pinhole),
%     d = offset_mm / n_medium (apparent shift seen through water).
%
% Board ('coarse' 40 mm 11x10 or 'fine' 12 mm 11x8, from
% IR Camera/calibration/calib_target_params.m) is picked from the detected
% square count.  The plate was printed scaled down, so the true square is
% design * opts.plate_scale, unless opts.square_mm (measured) is given.
%
% opts fields (defaults):
%   plate_scale 1, square_mm [], offset_mm 0, n_medium 1.33,
%   lens_f_mm 35, pixel_um 5.5, verbose true
%
% Outputs:
%   mmpp - mm per pixel on the measurement plane (full-resolution pixels)
%   c0   - [col row] of the board centre in the image (axis origin)
%   info - struct: square_mm, board, mmpp_plate, D_mm, rms_px, imgPts

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '..', 'IR Camera', 'calibration'));   % calib_target_params

    def = struct('plate_scale', 1, 'square_mm', [], 'offset_mm', 0, ...
                 'n_medium', 1.33, 'lens_f_mm', 35, 'pixel_um', 5.5, 'verbose', true);
    if nargin < 2, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end

    % Detect corners
    P   = lif_load_raw(fplate);
    lim = prctile(P(:), [0.5 99.5]);
    P8  = uint8(255 * min(max((P - lim(1)) / diff(lim), 0), 1));
    [imgPts, boardSize] = detectCheckerboardPoints(P8, 'PartialDetections', false);
    if isempty(imgPts) || any(isnan(imgPts(:)))
        error('lif_plate_scale:detect', 'Checkerboard not detected in %s', fplate);
    end

    board = [];
    for nm = {'fine', 'coarse'}
        b = calib_target_params(nm{1});
        if isequal(sort(boardSize), sort([b.n_rows b.n_cols])), board = b; end
    end
    if isempty(board)
        error('lif_plate_scale:board', ...
              'Detected %dx%d squares: matches neither the fine (11x8) nor coarse (11x10) board.', ...
              boardSize(1), boardSize(2));
    end
    square_mm = opts.square_mm;
    if isempty(square_mm), square_mm = board.square_mm * opts.plate_scale; end
    wPts = generateCheckerboardPoints(boardSize, square_mm);   % plate mm

    % Uniform scale
    tf = fitgeotrans(wPts, imgPts, 'nonreflectivesimilarity');   % plate mm -> px
    mmpp_plate = 1 / norm(tf.T(1:2, 1));
    res_px = sqrt(sum((transformPointsForward(tf, wPts) - imgPts).^2, 2));

    % Plate -> measurement-plane offset
    D_mm = (opts.lens_f_mm / (opts.pixel_um/1000)) * mmpp_plate;
    d_mm = opts.offset_mm / opts.n_medium;
    mmpp = mmpp_plate * (D_mm + d_mm) / D_mm;
    c0   = mean(imgPts, 1);

    info = struct('square_mm', square_mm, 'board', board.name, ...
                  'mmpp_plate', mmpp_plate, 'D_mm', D_mm, ...
                  'rms_px', sqrt(mean(res_px.^2)), 'imgPts', imgPts);

    if opts.verbose
        fprintf('Plate %s: %d corners, ''%s'' board, square %.3f mm (design %g mm)\n', ...
                fplate, size(imgPts,1), board.name, square_mm, board.square_mm);
        fprintf('  plate %.4f mm/px (uniform-scale RMS %.2f px, max %.2f px)\n', ...
                mmpp_plate, info.rms_px, max(res_px));
        fprintf('  camera-plate ~%.0f mm; offset %.1f mm (apparent %.1f) -> %.4f mm/px (+%.2f%%)\n', ...
                D_mm, opts.offset_mm, d_mm, mmpp, 100*d_mm/D_mm);
    end
end
