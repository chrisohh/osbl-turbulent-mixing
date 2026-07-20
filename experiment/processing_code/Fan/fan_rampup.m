d = daq("ni");
ch = addoutput(d, "Dev1", "ao0", "Voltage");

v_start = 1.5;
v_end = 5;
t_total = 30;
dt = 0.5;

try
    % Ramp up
    disp('Ramping up...');
    for v = linspace(v_start, v_end, t_total/dt + 1)
        write(d, v);
        fprintf('%.2f V\n', v);
        pause(dt);
    end

    disp('Holding at max...');
    pause(30);

    %stop fan
write(d, 0);
    % % Ramp down
    % disp('Ramping down...');
    % for v = linspace(v_end, v_start, t_total/dt + 1)
    %     write(d, v);
    %     fprintf('%.2f V\n', v);
    %     pause(dt);
    % end

    disp('Done.');

catch
    disp('Interrupted — setting to 0V.');
    write(d, 0);
end

%% Stop the fan
write(d, 0);