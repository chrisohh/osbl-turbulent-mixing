function [S, n_removed] = cisg_despeckle(S, k)
%CISG_DESPECKLE Remove isolated speckle: |S - local median| > k * local sigma  ->  NaN.
%
%   [S, n_removed] = cisg_despeckle(S, k)
%
% Local median over 5 x 5 (separable median of medians, NaN ignored); local sigma = 1.4826 * median of
% |S - median5| over 15 x 15, floored at 0.005 (slope units; smooth fields have a MAD of ~0 and would be stripped).
% k = 0 or [] disables.  Applies after the colour inversion in every ROI mode (cisg_calculate_slopes).

    n_removed = 0;
    if isempty(k) || k <= 0, return; end
    S = double(S);
    m = movmedian(movmedian(S, 5, 1, 'omitnan'), 5, 2, 'omitnan');
    a = abs(S - m);
    sg = 1.4826 * movmedian(movmedian(a, 15, 1, 'omitnan'), 15, 2, 'omitnan');
    sg = max(sg, 0.005);                                   % absolute floor (slope units): smooth fields have MAD ~ 0
    out = isfinite(S) & a > k * sg;
    S(out) = NaN;
    n_removed = nnz(out);
end
