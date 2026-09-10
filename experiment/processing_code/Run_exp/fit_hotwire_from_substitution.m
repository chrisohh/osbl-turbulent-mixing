%% fit_hotwire_from_substitution.m
% Build a hot-wire velocity calibration by SUBSTITUTION: the reference probe
% and the hot-wire cannot occupy the tunnel at the same time, so they are run
% separately and paired through the fan command voltage.
%
%   Run A (ref probe)  -> T_ref : FanV, U_mean          (from analyze_steps_from_trace)
%   Run B (hot-wire)   -> T_hw  : FanV, mean_Probe1..3  (stepTable from run_hotwire_calibration)
%
% Provide those two tables in the workspace, or set the two paths below.
%
% WHY PAIRING BY FAN VOLTAGE IS THE WEAK LINK: with simultaneous probes each
% (E,U) pair is measured at one instant and is valid whatever the flow was
% doing. Here the pair is stitched from two different runs, so it is only as
% good as the repeatability of fan voltage -> velocity between them. That is
% why both runs must use the UP-LEG ONLY (wave history is reproducible only
% from rest), the same settle protocol, and ideally the same day and ambient
% temperature. Treat the run-to-run spread as a floor on the result.

POLY_ORDER = 3;    % 3 keeps a 7-point fit honest; convert_E2U_fn reads C0..C4 so C4 = 0
REF_MAT = '';      % optional: path to a .mat holding T_ref
HW_MAT  = '';      % optional: path to a .mat holding T_hw (or stepTable)

CTA_DIR = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\CTA';
addpath(CTA_DIR);

%% 1. Get the two tables -----------------------------------------------------
if ~exist('T_ref','var')
    if isempty(REF_MAT), error('Need T_ref in the workspace, or set REF_MAT.'); end
    S = load(REF_MAT);
    if isfield(S,'T_ref'), T_ref = S.T_ref; elseif isfield(S,'T'), T_ref = S.T;
    else, error('%s contains neither T_ref nor T.', REF_MAT); end
end
if ~exist('T_hw','var')
    if isempty(HW_MAT), error('Need T_hw in the workspace, or set HW_MAT.'); end
    S = load(HW_MAT);
    if isfield(S,'T_hw'), T_hw = S.T_hw; elseif isfield(S,'stepTable'), T_hw = S.stepTable;
    else, error('%s contains neither T_hw nor stepTable.', HW_MAT); end
end

sensors = {'Probe1','Probe2','Probe3'};
hwCols  = "mean_" + string(sensors);
missing = hwCols(~ismember(hwCols, string(T_hw.Properties.VariableNames)));
if ~isempty(missing)
    error('Hot-wire table is missing column(s): %s', strjoin(missing, ', '));
end

%% 2. Pair on fan voltage ----------------------------------------------------
% Exact matches only. An interpolated pair would hide the fact that the two
% runs did not actually visit the same operating point, so anything unmatched
% is reported and dropped rather than filled in.
Vref = T_ref.FanV(:);  Uref = T_ref.U_mean(:);
Vhw  = T_hw.FanV(:);
tol  = 1e-6;

pairV = []; pairU = []; pairE = [];
for i = 1:numel(Vhw)
    j = find(abs(Vref - Vhw(i)) < tol, 1);
    if isempty(j) || isnan(Uref(j)), continue; end
    e = arrayfun(@(c) T_hw.(char(c))(i), hwCols);
    if any(isnan(e)), continue; end
    pairV(end+1,1) = Vhw(i);  %#ok<SAGROW>
    pairU(end+1,1) = Uref(j); %#ok<SAGROW>
    pairE(end+1,:) = e(:)';   %#ok<SAGROW>
end

unmatched = setdiff(round(Vhw,4), round(pairV,4));
if ~isempty(unmatched)
    fprintf(['Unmatched hot-wire levels (no ref point at the same fan voltage): %s V\n' ...
             '  -> the two runs did not visit the same operating points.\n'], ...
        strjoin(compose('%.2f', unmatched(:)'), ' '));
end
n = numel(pairV);
fprintf('Paired %d levels: %s V\n', n, strjoin(compose('%.2f', pairV'), ' '));
fprintf('  U range %.3f .. %.3f m/s\n', min(pairU), max(pairU));

if n <= POLY_ORDER + 1
    error(['%d paired points cannot support an order-%d fit (%d free parameters). ' ...
           'Lower POLY_ORDER or add levels to both runs.'], n, POLY_ORDER, POLY_ORDER+1);
elseif n < POLY_ORDER + 3
    warning(['Only %d points for an order-%d fit (%d dof) -- the fit will partly ' ...
             'trace noise. Check the residuals below before trusting it.'], ...
        n, POLY_ORDER, n-POLY_ORDER-1);
end

%% 3. Fit each sensor --------------------------------------------------------
calFit = struct('sensor',{},'p',{},'C',{},'rms',{},'maxres',{},'E',{},'U',{});
fprintf('\nOrder-%d fits, U = C0 + C1*E + ...\n', POLY_ORDER);
for i = 1:numel(sensors)
    E = pairE(:,i);
    p = polyfit(E, pairU, POLY_ORDER);
    r = pairU - polyval(p, E);
    C = zeros(1,5);                      % pad to C0..C4 for convert_E2U_fn
    Cl = fliplr(p);  C(1:numel(Cl)) = Cl;
    calFit(i) = struct('sensor',sensors{i},'p',p,'C',C, ...
        'rms',sqrt(mean(r.^2)),'maxres',max(abs(r)),'E',E,'U',pairU); %#ok<SAGROW>
    fprintf('  %s: %s\n           RMS %.4f m/s (%.2f%%), max %.4f m/s\n', sensors{i}, ...
        strjoin(compose('C%d=%.6g', (0:4)', C(:)), ' '), ...
        calFit(i).rms, 100*calFit(i).rms/mean(pairU), calFit(i).maxres);
end

%% 4. Plot -------------------------------------------------------------------
figure('Position',[80 80 1150 420],'Color','w');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i = 1:numel(sensors)
    nexttile; hold on;
    plot(pairE(:,i), pairU, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.5 0.8], ...
        'MarkerEdgeColor','none', 'DisplayName','paired points');
    Ef = linspace(min(pairE(:,i)), max(pairE(:,i)), 200);
    plot(Ef, polyval(calFit(i).p, Ef), 'k-', 'LineWidth', 1.4, ...
        'DisplayName', sprintf('order-%d fit', POLY_ORDER));
    hold off; grid on;
    xlabel(sprintf('E_%d (V)', i)); ylabel('U_{ref} (m/s)');
    title(sprintf('%s  (RMS %.3f m/s)', sensors{i}, calFit(i).rms));
    legend('Location','northwest');
end

fprintf(['\nCoefficients are in calFit(i).C as [C0 C1 C2 C3 C4] -- the order\n' ...
         'convert_E2U_fn expects. They are NOT written to the calibration file;\n' ...
         'copy them in once the residuals look acceptable.\n']);
