function root = get_server_root()
%GET_SERVER_ROOT  Returns the root path to DelawareDataBackup,
%   switching automatically between running ON the server vs. remotely.
%
%   On the server:    root = 'D:\DelawareDataBackup\'
%   Off the server:   root = '\\Airseaserver41\D\DelawareDataBackup\'
%
%   Usage:
%     root = get_server_root();
%     load([root 'Longitudinal\PIV\ExpLCL_2_01\...']);

SERVER_LOCAL = 'D:\DelawareDataBackup\';
SERVER_REMOTE = '\\Airseaserver41\D\DelawareDataBackup\';

if exist(SERVER_LOCAL, 'dir')
    root = SERVER_LOCAL;
else
    root = SERVER_REMOTE;
end

end
