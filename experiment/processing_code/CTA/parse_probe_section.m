function cal = parse_probe_section(filename, probeName)
% PARSE_PROBE_SECTION  Read one probe's block out of a StreamWare cal/header file.
%
%   cal = parse_probe_section(filename, probeName)
%
%   A StreamWare export can hold several probes in one file -- probe4.txt holds
%   the tri-axial hot-wire (55P95) and the velocity reference probe (T29). Each
%   block starts at a "Probe name:" line. This returns the polynomial and the
%   calibrated range for the block whose name matches probeName, so the T29
%   reference probe is reduced with the coefficients from THIS file rather than
%   the certificate table hardcoded in convert_Eref2Uref.m.
%
%   probeName - name as it appears after "Probe name:" (e.g. 'T29', '55P95').
%               Matched case-insensitively, leading/trailing space ignored.
%
%   cal fields:
%     name              probe name as found in the file
%     nSensors          number of "Probe sensor no.:" entries in the block
%     C                 nSensors x 5, columns C0..C4 (lowest power first)
%     U_floor, E_floor  min calibration point per sensor (m/s, V)
%     U_ceil,  E_ceil   max calibration point per sensor (m/s, V)
%     T_ref             calibration reference temperature (degC), NaN if absent
%
%   For the full tri-axial reduction (k, h, Mp) use parse_calibration.m, which
%   reads the first probe block in the file.

    txt = fileread(filename);
    txt = regexprep(txt, '\r\n?', '\n');

    % Split the file into per-probe blocks at each "Probe name:" line
    starts = regexp(txt, 'Probe name:', 'start');
    if isempty(starts)
        error('parse_probe_section:noProbes', ...
            'No "Probe name:" line found in %s', filename);
    end
    bounds = [starts, numel(txt) + 1];

    block = '';
    for k = 1:numel(starts)
        b    = txt(bounds(k) : bounds(k+1) - 1);
        nameTok = regexp(b, 'Probe name:\s*([^\n]*)', 'tokens', 'once');
        thisName = strtrim(nameTok{1});
        if strcmpi(thisName, strtrim(probeName))
            block = b;
            cal.name = thisName;
            break;
        end
    end
    if isempty(block)
        allNames = regexp(txt, 'Probe name:\s*([^\n]*)', 'tokens');
        error('parse_probe_section:notFound', ...
            'Probe "%s" not found in %s (file holds: %s)', probeName, filename, ...
            strjoin(cellfun(@(c) strtrim(c{1}), allNames, 'UniformOutput', false), ', '));
    end

    % Per-sensor polynomial and calibrated range
    sensorStarts = regexp(block, 'Probe sensor no\.:', 'start');
    cal.nSensors = numel(sensorStarts);
    if cal.nSensors == 0
        error('parse_probe_section:noSensors', ...
            'Probe "%s" has no "Probe sensor no.:" entries.', cal.name);
    end
    sBounds = [sensorStarts, numel(block) + 1];

    cal.C       = nan(cal.nSensors, 5);
    cal.U_floor = nan(1, cal.nSensors);
    cal.E_floor = nan(1, cal.nSensors);
    cal.U_ceil  = nan(1, cal.nSensors);
    cal.E_ceil  = nan(1, cal.nSensors);

    for s = 1:cal.nSensors
        sb = block(sBounds(s) : sBounds(s+1) - 1);
        for c = 0:4
            tok = regexp(sb, sprintf('C%d:\\s*(-?[\\d\\.eE+-]+)', c), 'tokens', 'once');
            if ~isempty(tok)
                cal.C(s, c+1) = str2double(tok{1});
            end
        end
        mn = regexp(sb, 'Min\. calibration point[^:]*:\s*(-?[\d\.]+)\s+(-?[\d\.]+)', ...
                    'tokens', 'once');
        if ~isempty(mn)
            cal.U_floor(s) = str2double(mn{1});
            cal.E_floor(s) = str2double(mn{2});
        end
        mx = regexp(sb, 'Max\. calibration point[^:]*:\s*(-?[\d\.]+)\s+(-?[\d\.]+)', ...
                    'tokens', 'once');
        if ~isempty(mx)
            cal.U_ceil(s) = str2double(mx{1});
            cal.E_ceil(s) = str2double(mx{2});
        end
    end

    % Reference temperature. The value sits after the "[deg.C]:" label -- the
    % pattern has to skip that label, or it matches the dot inside "[deg.C]".
    tTok = regexp(block, 'Cal\. ref\. temp\.[^:]*:\s*(-?\d+\.?\d*)', 'tokens', 'once');
    if ~isempty(tTok)
        cal.T_ref = str2double(tTok{1});
    else
        cal.T_ref = NaN;
    end
end
