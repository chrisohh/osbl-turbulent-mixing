function p = get_piv_params(varargin)
% get_piv_params  Timing, path, and ramp parameters for a PIV experiment.
%
%   p = get_piv_params('ExpLCL_2_01')        % parse setup/ramp/run from name
%   p = get_piv_params('longitudinal', 2, 1) % explicit setup, ramp, run
%
% Setup may be given via the name prefix (ExpLCL / ExpLCTA / ExpLCTB) or as
% the first argument: 'longitudinal' | 'transverseA' | 'transverseB'
% (aliases: 'LCL' | 'LCTA' | 'LCTB').
%
% Fields returned:
%   p.exp_name            — e.g. 'ExpLCL_2_01'
%   p.setup               — 'LCL' | 'LCTA' | 'LCTB'
%   p.ramp, p.run         — ramp number (1-4), run number
%   p.load_path           — full path to experiment folder
%   p.turb_save, p.piv_save — PIVMat_TURB/ and PIVMat/ paths
%   p.DX                  — m/pixel
%   p.DT                  — s, inter-frame (A→B within pair), for velocity scaling
%   p.Fs_cam              — Hz, camera image acquisition rate
%   p.Fs_PIV              — Hz, velocity map rate (= Fs_cam/2)
%   p.DT_pair             — s, time between consecutive pairs (= 1/Fs_PIV)
%   p.Fs_IR               — Hz, IR camera frame rate
%   p.Fs_CO2              — Hz, CO2 heat-marker lay-down rate (NaN if N/A)
%   p.IID                 — s, imaging initial delay
%   p.t_imaging           — s, imaging start relative to wind trigger (= 5 + IID)
%   p.t_wind_ramp         — s, wind starts ramping (40 from trigger)
%   p.t_imaging_wrt_ramp  — s, imaging start relative to wind ramp-up
%   p.rampup              — s, wind ramp-up duration
%   p.voltage             — V, target wind voltage

%% --- resolve setup, ramp, run, exp_name from arguments ---
if nargin == 1
    exp_name = varargin{1};
    tok = regexp(exp_name, 'Exp(LCL|LCTA|LCTB)_(\d+)_(\d+)', 'tokens');
    if isempty(tok)
        error(['exp_name must match Exp<LCL|LCTA|LCTB>_<ramp>_<run>, got: %s'], exp_name);
    end
    setup = tok{1}{1};
    ramp  = str2double(tok{1}{2});
    run   = str2double(tok{1}{3});
elseif nargin == 3
    setup = canonical_setup(varargin{1});
    ramp  = varargin{2};
    run   = varargin{3};
    exp_name = sprintf('Exp%s_%d_%02d', setup, ramp, run);
else
    error('Call as get_piv_params(name) or get_piv_params(setup, ramp, run).');
end

p.exp_name = exp_name;
p.setup    = setup;
p.ramp     = ramp;
p.run      = run;

%% --- per-setup constants (paths, DX, dt, heat-marker rate) ---
switch setup
    case 'LCL'    % Longitudinal
        root    = 'D:\DelawareDataBackup\Longitudinal\PIV\';
        p.DX    = 1/17697.69;   % m/pixel
        p.DT    = 10e-3;        % s, inter-frame (A→B)
        p.Fs_CO2 = NaN;         % heat markers not the primary LCL channel
        p.piv_folder     = 'PIV';        % PIVRaw subfolder, PIV camera
        p.pivsurf_folder = 'PIVSURF';    % PIVRaw subfolder, surface camera
    case 'LCTA'   % Transverse A
        root    = 'D:\DelawareDataBackup\Transverse\PIV\';
        p.DX    = 6.7861e-05;   % m/pixel
        p.DT    = 12.5e-3;      % s, inter-frame (A→B)
        p.Fs_CO2 = 0.48;        % Hz
        p.piv_folder     = 'PIVCC';      % transverse uses the CC folders
        p.pivsurf_folder = 'PIVSURFCC';
    case 'LCTB'   % Transverse B
        root    = 'D:\DelawareDataBackup\Transverse\PIV\';
        p.DX    = 6.7861e-05;   % m/pixel
        p.DT    = 8e-3;         % s, inter-frame (A→B)
        p.Fs_CO2 = 1.8;         % Hz
        p.piv_folder     = 'PIVCC';      % transverse uses the CC folders
        p.pivsurf_folder = 'PIVSURFCC';
    otherwise
        error('Unknown setup: %s', setup);
end

p.load_path = fullfile(root, exp_name);
p.turb_save = fullfile(p.load_path, 'Chris_recompute', 'PIVMat_TURB');
p.piv_save  = fullfile(p.load_path, 'Chris_recompute', 'PIVMat');

%% --- shared camera / PIV timing ---
p.Fs_cam  = 14.4;         % Hz — camera acquisition rate (images)
p.Fs_PIV  = p.Fs_cam / 2; % Hz — velocity map rate (one pair = two images)
p.DT_pair = 1 / p.Fs_PIV; % s  — time between consecutive pairs (for time axis)
p.Fs_IR   = 43.2;         % Hz — IR camera frame rate

%% --- shared ramp table (same definitions for all setups) ---
% ramp N: wind triggered t=0, stays 0V until t=40s, ramps up, imaging at t=5+IID
ramp_table = [
%  ramp  IID   rampup  voltage
    1     57    30      7
    2     73    90      7
    3     88    120     7
    4     111   120     5
];
row = ramp_table(ramp_table(:,1) == ramp, :);
if isempty(row)
    error('Unknown ramp number: %d', ramp);
end
p.IID                = row(2);                       % s, imaging initial delay
p.t_imaging          = 5 + p.IID;                    % s, imaging start re wind trigger (t=0)
p.t_wind_ramp        = 40;                           % s, wind starts ramping from trigger
p.t_imaging_wrt_ramp = p.t_imaging - p.t_wind_ramp;  % s, imaging start re wind ramp-up
p.rampup             = row(3);                       % s, wind ramp-up duration
p.voltage            = row(4);                       % V, target wind voltage

end

%% =========================================================================
function s = canonical_setup(x)
% Map setup aliases to the canonical 'LCL' | 'LCTA' | 'LCTB' token.
switch lower(x)
    case {'lcl', 'longitudinal', 'long'},          s = 'LCL';
    case {'lcta', 'transversea', 'transverse_a'},  s = 'LCTA';
    case {'lctb', 'transverseb', 'transverse_b'},  s = 'LCTB';
    otherwise, error('Unknown setup alias: %s', x);
end
end
