% visualize_triggers.m
% Timeline diagram of experiment triggers / acquisition for one ramp.
%
% Layout:
%   Top axes  — wind voltage vs time (trigger, ramp-up, target hold)
%   Bottom axes — multi-lane timeline, one lane per instrument:
%                 PIV (Longitudinal or Transverse B), LIF, IR + heat markers,
%                 wave gauge. Acquisition windows drawn as bands; pulsed
%                 channels get tick marks (if not too dense) + a zoom inset
%                 showing the inter-frame dt within one PIV pair.
%
% t = 0 is the WIND TRIGGER. Wind ramps at t = 40 s; imaging starts at
% t = 5 + IID. All cameras/lasers were synchronized.

clear; clc;

%% =========================================================================
%% PARAMETERS
%% =========================================================================
exp_name = 'ExpLCL_2_01';        % which experiment (any setup: ExpLCL_/ExpLCTA_/ExpLCTB_)

p = get_piv_params(exp_name);    % ramp timing, DT, Fs, heat-marker rate for this setup

% --- known timing (from get_piv_params) ---
t_trigger   = 0;                 % s, wind trigger
t_wind_ramp = p.t_wind_ramp;     % s, wind starts ramping (40)
rampup      = p.rampup;          % s, ramp-up duration
voltage     = p.voltage;         % V, target voltage
t_img_start = p.t_imaging;       % s, imaging start (= 5 + IID), re wind trigger
img_dur     = 60;                % s, total imaging duration (always 60 s)
t_img_end   = t_img_start + img_dur;

% --- acquisition rates / inter-frame dt (per setup) ---
Fs_PIV  = p.Fs_PIV;              % velocity-map (pair) rate
dt_pair = p.DT;                  % inter-frame dt (A→B within a pair)
Fs_IR   = p.Fs_IR;              % Hz, IR frame rate
Fs_CO2  = p.Fs_CO2;             % Hz, heat-marker lay-down rate (NaN if N/A)

setup_label = struct('LCL','Longitudinal','LCTA','Transverse A','LCTB','Transverse B');
piv_label   = ['PIV (' setup_label.(p.setup) ')'];

% --- placeholders: EDIT with real values (timing not in my notes) -------
lif_on    = [t_img_start, t_img_end];   % LIF acquisition window [start end] (s)  <-- VERIFY
Fs_LIF    = NaN;                         % Hz, LIF frame rate                      <-- FILL IN
wave_on   = [t_trigger,   t_img_end];   % wave gauge window [start end] (s)        <-- VERIFY
Fs_wave   = NaN;                         % Hz, wave gauge sample rate              <-- FILL IN

%% =========================================================================
%% FIGURE
%% =========================================================================
fig = figure('Color','white','Position',[120 80 1100 640]);
tlo = tiledlayout(fig, 4, 1, 'TileSpacing','compact','Padding','compact');

t_max = t_img_end + 5;           % x-limit with a little margin

%% ---- TOP: wind voltage ----------------------------------------------------
ax1 = nexttile([1 1]);
tw = [t_trigger, t_wind_ramp, t_wind_ramp + rampup, t_max];
vw = [0,         0,           voltage,              voltage];
plot(ax1, tw, vw, '-', 'Color',[0.85 0.33 0.10], 'LineWidth', 2); hold(ax1,'on');
xline(ax1, t_trigger,   ':', 'trigger',  'Color',[0.4 0.4 0.4], 'LabelOrientation','horizontal');
xline(ax1, t_wind_ramp, ':', 'ramp-up',  'Color',[0.4 0.4 0.4], 'LabelOrientation','horizontal');
ylabel(ax1, 'wind (V)'); ylim(ax1, [-0.3, voltage+0.7]);
title(ax1, sprintf('Trigger timeline — %s  (ramp %d)', exp_name, p.ramp), 'Interpreter','none');
grid(ax1,'on'); set(ax1,'XTickLabel',[]); xlim(ax1,[t_trigger, t_max]);

%% ---- BOTTOM: instrument lanes (3 tiles tall) -----------------------------
ax2 = nexttile([3 1]); hold(ax2,'on'); box(ax2,'on');

% lane definitions: name, [start end], rate(Hz), color, pulsed?
lanes = {
    piv_label,        [t_img_start t_img_end], Fs_PIV, [0.00 0.45 0.74], false
    'LIF',            lif_on,                  Fs_LIF, [0.47 0.67 0.19], false
    'IR camera',      [t_img_start t_img_end], Fs_IR,  [0.64 0.08 0.18], false
    'CO_2 heat marks',[t_img_start t_img_end], Fs_CO2, [0.93 0.69 0.13], true
    'Wave gauge',     wave_on,                 Fs_wave,[0.49 0.18 0.56], false
};
nL  = size(lanes,1);
h   = 0.6;                       % band half-height fraction

for k = 1:nL
    yk   = nL - k + 1;           % top lane = first row
    win  = lanes{k,2};
    col  = lanes{k,4};
    rate = lanes{k,3};
    % acquisition band
    fill(ax2, [win(1) win(2) win(2) win(1)], yk + h*[-1 -1 1 1], col, ...
        'FaceAlpha',0.25, 'EdgeColor',col, 'LineWidth',1);
    % pulsed channels: tick marks if sparse enough
    if lanes{k,5} && isfinite(rate)
        tp = win(1):1/rate:win(2);
        if numel(tp) <= 200
            plot(ax2, [tp;tp], yk + h*[-1;1], '-', 'Color',col, 'LineWidth',0.8);
        end
    end
    % rate annotation
    if isfinite(rate)
        rate_str = sprintf('%.3g Hz', rate);
    else
        rate_str = 'rate: EDIT';
    end
    text(ax2, win(2)+0.5, yk, rate_str, 'FontSize',8, 'Color',col, ...
        'VerticalAlignment','middle');
end

% reference lines aligned with the wind axes
xline(ax2, t_trigger,   ':', 'Color',[0.4 0.4 0.4]);
xline(ax2, t_wind_ramp, ':', 'Color',[0.4 0.4 0.4]);
xline(ax2, t_img_start, '--','imaging start', 'Color',[0.2 0.2 0.2], 'LabelOrientation','horizontal');

set(ax2, 'YTick', 1:nL, 'YTickLabel', flipud(lanes(:,1)));
ylim(ax2, [0.3, nL+0.7]); xlim(ax2, [t_trigger, t_max]);
xlabel(ax2, 't (s)   —   t = 0 at wind trigger'); grid(ax2,'on');
linkaxes([ax1 ax2], 'x');

%% ---- ZOOM INSET: pulse pattern + inter-frame dt within one PIV pair -------
axz = axes(fig, 'Position',[0.62 0.13 0.30 0.12]); hold(axz,'on'); box(axz,'on');
n_pairs = 3;                                   % show first few pairs
for ip = 0:n_pairs-1
    tp0 = t_img_start + ip/Fs_PIV;             % start of pair ip
    plot(axz, [tp0 tp0],         [0 1], '-',  'Color',[0.00 0.45 0.74], 'LineWidth',1.5); % frame A
    plot(axz, [tp0 tp0]+dt_pair, [0 1], '--', 'Color',[0.00 0.45 0.74], 'LineWidth',1.5); % frame B
    if ip==0
        text(axz, tp0+dt_pair/2, 1.08, sprintf('dt = %.3g ms', dt_pair*1e3), ...
            'HorizontalAlignment','center','FontSize',8);
    end
end
xlim(axz, [t_img_start - 0.02, t_img_start + (n_pairs-1)/Fs_PIV + dt_pair + 0.03]);
ylim(axz, [0 1.25]); set(axz,'YTick',[]);
xlabel(axz, 't (s)','FontSize',8);
title(axz, sprintf('PIV pulse pattern (pairs @ %.3g Hz, A/B split by dt)', Fs_PIV), 'FontSize',8);

fprintf('Camera: %s | pair rate %.3g Hz | dt %.3g ms | imaging %.0f–%.0f s\n', ...
    camera, Fs_PIV, dt_pair*1e3, t_img_start, t_img_end);
