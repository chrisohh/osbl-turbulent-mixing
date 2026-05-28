function cmap = inferno(n)
% INFERNO  Matplotlib-style perceptually uniform colormap.
%   cmap = inferno(n)   -> n x 3 RGB array (default n = 256).
%
% Order (low -> high): near-black -> dark purple -> red-purple -> red ->
% orange -> yellow -> very pale yellow. Use flipud(inferno) to get the
% "white-yellow-red-purple-dark-purple" sequence (high -> low light).

if nargin < 1 || isempty(n), n = 256; end

% Key control points sampled from the Matplotlib inferno LUT.
pos = [0.00 0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95 1.00];
key = [...
    0.001462 0.000466 0.013866; ...   % 0.00 near-black
    0.046386 0.030374 0.150328; ...   % 0.05 deep purple
    0.137510 0.067880 0.293385; ...   % 0.15 dark purple
    0.243113 0.057105 0.402796; ...   % 0.25 purple
    0.358078 0.071729 0.432984; ...   % 0.35 red-purple
    0.473911 0.087502 0.450787; ...   % 0.45 magenta-purple
    0.597055 0.121733 0.434744; ...   % 0.55 magenta
    0.717520 0.182919 0.380284; ...   % 0.65 red-magenta
    0.833861 0.273809 0.298802; ...   % 0.75 red-orange
    0.927499 0.408479 0.179212; ...   % 0.85 orange
    0.987622 0.645320 0.039886; ...   % 0.95 yellow-orange
    0.988362 0.998364 0.644924];      % 1.00 light yellow

t    = linspace(0, 1, n).';
cmap = interp1(pos, key, t, 'pchip');
cmap = min(1, max(0, cmap));
end
