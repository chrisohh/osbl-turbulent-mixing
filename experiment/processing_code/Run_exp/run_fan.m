clear;clc

DEV_ID = "Dev4";   % <-- confirm this matches NI MAX for the USB-6451
d = daq("ni");
addoutput(d, DEV_ID, "ao0", "Voltage");
vStart = 2.1;
vEnd   = 2.1;
rampUpTime  = 0;   % seconds, ramp up duration
rampDownTime = 0;   % seconds, ramp up duration
holdTime  = 120;    % seconds, hold at V_END
dt     = 0.5;  % seconds, step interval
fan_control(d, vStart, vEnd, rampUpTime, holdTime, rampDownTime, dt)