%% Function to read a specific frame
function frame = readFrame(v, frameNumber)
reset(v);
for i = 1:(frameNumber-1)
    step(v);
end
frame = step(v);
end

% % Usage:
% frame500 = readFrame(v, 500);  % Get frame 500
% frame1000 = readFrame(v, 1000); % Get frame 1000