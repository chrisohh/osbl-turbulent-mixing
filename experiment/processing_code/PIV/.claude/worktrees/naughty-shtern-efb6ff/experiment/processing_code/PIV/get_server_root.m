function root = get_server_root()
%GET_SERVER_ROOT  Returns the root path to DelawareDataBackup,
%   auto-detecting whether you are running on the server, locally with
%   F:\ mapped to the server share, or remotely over the network.
%
%   On the server  (D:\DelawareDataBackup\ exists):
%       root = 'D:\DelawareDataBackup\'
%
%   Local machine  (F:\ is mapped directly to D:\DelawareDataBackup\
%       so F:\Transverse\... == D:\DelawareDataBackup\Transverse\...):
%       root = 'F:\'
%
%   Remote / UNC   (neither local drive exists):
%       root = '\\Airseaserver41\D\DelawareDataBackup\'
%
%   Usage:
%     root = get_server_root();
%     load([root 'Longitudinal\PIV\ExpLCL_2_01\...']);

SERVER_LOCAL  = 'D:\DelawareDataBackup\';
MAPPED_DRIVE  = 'F:\';                             % F:\ == DelawareDataBackup\
SERVER_REMOTE = '\\Airseaserver41\D\DelawareDataBackup\';

if exist(SERVER_LOCAL, 'dir')
    root = SERVER_LOCAL;
elseif exist(MAPPED_DRIVE, 'dir')
    root = MAPPED_DRIVE;
else
    root = SERVER_REMOTE;
end

end
