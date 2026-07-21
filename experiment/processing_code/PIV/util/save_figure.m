function save_figure(fig, name, varargin)
% SAVE_FIGURE  Save current figure as PNG (300 dpi) and .fig under figures/.
%
% Duplicate of Slope_Gauge/util/save_figure.m.

p = inputParser;
addParameter(p, 'dir', 'figures', @ischar);
parse(p, varargin{:});
outdir = p.Results.dir;

if ~exist(outdir, 'dir'), mkdir(outdir); end
png_path = fullfile(outdir, [name '.png']);
fig_path = fullfile(outdir, [name '.fig']);
exportgraphics(fig, png_path, 'Resolution', 300);
savefig(fig, fig_path);
fprintf('  saved %s\n', png_path);
end
