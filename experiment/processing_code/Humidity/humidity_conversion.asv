function [AH] = humidity_conversion(RH,T)
% HUMIDITY_CONVERSION Convert between relative and absolute humidity
%
% Syntax:
%   AH = humidity_conversion(T, RH)           % RH to AH
%   RH = humidity_conversion(T, AH, 'AH')     % AH to RH
%   [AH, RH] = humidity_conversion(T, RH)     % Both outputs
%
% Inputs:
%   T  - Temperature in degrees Celsius (scalar or array)
%   RH - Relative humidity in percent [0-100] (for RH to AH conversion)
%   AH - Absolute humidity in g/m^3 (for AH to RH conversion)
%
% Outputs:
%   AH - Absolute humidity in g/m^3
%   RH - Relative humidity in percent
%
% Example:
%   AH = humidity_conversion(20, 60)          % T=20°C, RH=60% -> AH
%   RH = humidity_conversion(20, 10.4, 'AH')  % T=20°C, AH=10.4 g/m³ -> RH


    % RH to AH conversion (default)
    RH_input = humidity_value;


% Calculate saturation vapor pressure using Magnus formula
% e_s in hPa (or mb)
e_s = 6.112 * exp((17.67 * T) ./ (T + 243.5));

% Perform conversion
        % RH to AH
        if any(RH_input < 0) || any(RH_input > 100)
            warning('Relative humidity should be between 0 and 100%');
        end

        % Actual vapor pressure
        e = (RH_input / 100) .* e_s;

        % Absolute humidity (g/m³)
        T_K = T + 273.15;
        AH = (2.16679 * e) ./ T_K;
end