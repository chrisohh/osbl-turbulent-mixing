function RT = thermistor_resistance(T_celsius, R0, beta, T0)
    % THERMISTOR_RESISTANCE Calculate thermistor resistance at given temperature
    %
    % RT = thermistor_resistance(T_celsius, R0, beta, T0)
    %
    % Inputs:
    %   T_celsius - Temperature in Celsius (scalar or array)
    %   R0        - Resistance at reference temperature (Ohms) [default: 10000]
    %   beta      - Material constant (Kelvin) [default: 3550]
    %   T0        - Reference temperature (Celsius) [default: 25]
    %
    % Output:
    %   RT        - Resistance at T_celsius (Ohms)
    %
    % Example:
    %   RT = thermistor_resistance(18, 10000, 3550, 25);
    %
    % Based on Steinhart-Hart equation (simple form):
    %   RT/R0 = exp(beta * (1/T - 1/T0))
    
    % Set default values if not provided
    if nargin < 4
        T0 = 25;        % Default reference temperature (°C)
    end
    if nargin < 3
        beta = 3521;    % Default beta for Material System A, Curve 5
    end
    if nargin < 2
        R0 = 10000;     % Default 10kΩ at 25°C
    end
    
    % Convert to Kelvin
    T_kelvin = T_celsius + 273.15;
    T0_kelvin = T0 + 273.15;
    
    % Calculate resistance using Steinhart-Hart equation
    RT = R0 .* exp(beta .* (1./T_kelvin - 1./T0_kelvin));
    
end