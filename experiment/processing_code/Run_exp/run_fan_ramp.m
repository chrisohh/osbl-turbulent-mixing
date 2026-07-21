function run_fan_ramp(d, vStart, vEnd, rampUpTime, holdTime, rampDownTime, dt)
    try
        nStepsUp = round(rampUpTime / dt) + 1;
        nStepsDown = round(rampDownTime / dt) + 1;

        disp('Ramping up...');
        for v = linspace(vStart, vEnd, nStepsUp)
            write(d, v);
            fprintf('%.2f V\n', v);
            pause(dt);
        end

        fprintf('Holding at %.1f V...\n', vEnd);
        pause(holdTime);

        disp('Ramping down...');
        for v = linspace(vEnd, vStart, nStepsDown)
            write(d, v);
            fprintf('%.2f V\n', v);
            pause(dt);
        end

        disp('Fan ramp done.');
    catch ME
        disp('Interrupted -- setting fan to 0V.');
        write(d, 0);
        rethrow(ME);
    end
    endv