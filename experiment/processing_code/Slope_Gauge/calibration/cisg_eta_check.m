%% cisg_eta_check.m
% Sign / scale check of util/cisg_slope_to_eta.m with an analytic field,
%   eta = A sin(2 pi x / lam) + B cos(2 pi y / mu)   (mm, x right, y UP),
% slopes sampled on the pixel grid (rows run DOWN), integrated and compared
% with the mean-removed truth.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'util'));

mmpx = 0.0703;  nr = 1843;  nc = 2460;  q = 16;
x = (0:nc-1) * mmpx;                                    % image right
y = ((nr-1:-1:0)') * mmpx;                              % image UP: row 1 = largest y
A = 2;  lam = 40;  B = 1.5;  mu = 60;                   % mm
Sx = ones(nr, 1) * (A * 2*pi/lam * cos(2*pi*x/lam));
Sy = -B * 2*pi/mu * sin(2*pi*y/mu) * ones(1, nc);
Sx(300:500, 400:800) = NaN;  Sy(300:500, 400:800) = NaN;            % a gap
[eta, info] = cisg_slope_to_eta(Sx, Sy, mmpx, mmpx, q);

xc = ((0:floor(nc/q)-1) * q + (q-1)/2) * mmpx;
yc = (nr - 1 - ((0:floor(nr/q)-1)' * q + (q-1)/2)) * mmpx;
tru = A * sin(2*pi*xc/lam) + B * cos(2*pi*yc/mu);
tru = tru - mean(tru(:));
p = polyfit(tru(:), eta(:), 1);
fprintf('synthetic eta: rms truth %.3f mm, rms error %.3f mm (%.1f%%), gain %.3f, corr %.4f, filled cells %.1f%%\n', ...
        rms(tru(:)), rms(eta(:) - tru(:)), 100*rms(eta(:) - tru(:))/rms(tru(:)), p(1), corr(tru(:), eta(:)), 100*info.nan_fraction);
% round trip: gradient of eta vs the (coarse) input slopes
[gx, gr] = gradient(eta, info.dx_mm, info.dy_mm);
ok = ~isnan(info.Sx_c);
fprintf('round trip: rms(d eta/dx - Sx) %.4f vs rms Sx %.4f; rms(-d eta/drow - Sy) %.4f vs rms Sy %.4f\n', ...
        rms(gx(:) - info.Sx_c(:)), rms(info.Sx_c(:)), rms(-gr(:) - info.Sy_c(:)), rms(info.Sy_c(:)));
