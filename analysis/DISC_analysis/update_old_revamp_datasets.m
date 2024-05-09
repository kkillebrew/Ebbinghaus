% update two old revamp datasets so that variables fit conventions of revamp_HORRM
% N.B. no values are changed from their original values.  only adding variables that weren't originally
%      explicitly defined.  so we aren't doing anything fishy here.

% BKN_1_revamp
clear all; clear all; % make sure there are no variables in the current workspace
load ../data/BKN_1_revamp_full.mat
% update or convert
condition_list = 1:5; % [-9 10 11 12 13] -> [1 2 3 4 5]
rawdata(:,1) = Remap(rawdata(:,1),[-9 10 11 12 13],1:5); % update conditions in rawdata
% add to match newdataset
duration_sec = duration/refresh; % duration of animation in seconds
diagrow = sgrowth; % old datasets using diagonal growth that was based on sgrowth in x and y direction (not an angled vector of sgrowth size)
pursuit_speed = (sgrowth * refresh) / ppd;    % (pix/frame * frame/s ) / (pix/deg) = deg/s
pursuit_dist  = pursuit_speed * duration_sec; % deg/s * s = deg
% add new variables for special cases
diag_pursuit_speed = (diagrow * refresh) / ppd;    % (pix/frame * frame/s ) / (pix/deg) = deg/s
diag_pursuit_dist  = diag_pursuit_speed * duration_sec; % deg/s * s = deg
% save, but DO NOT OVERWRITE
save('../data/BKN_1_revamp_UPDATED_full')
save('../data/BKN_1_revamp_UPDATED','rawdata','reversal_list','Pos_v_Trace');


% DNS_1_revamp
clear all; clear all; % make sure there are no variables in the current workspace
load ../data/DNS_1_revamp_full.mat
% update or convert
condition_list = 1:5; % [-9 10 11 12 13] -> [1 2 3 4 5]
rawdata(:,1) = Remap(rawdata(:,1),[-9 10 11 12 13],1:5); % update conditions in rawdata
% add to match newdataset
duration_sec = duration/refresh; % duration of animation in seconds
diagrow = sgrowth; % old datasets using diagonal growth that was based on sgrowth in x and y direction (not an angled vector of sgrowth size)
pursuit_speed = (sgrowth * refresh) / ppd;    % (pix/frame * frame/s ) / (pix/deg) = deg/s
pursuit_dist  = pursuit_speed * duration_sec; % deg/s * s = deg
% add new variables for special cases
diag_pursuit_speed = (diagrow * refresh) / ppd;    % (pix/frame * frame/s ) / (pix/deg) = deg/s
diag_pursuit_dist  = diag_pursuit_speed * duration_sec; % deg/s * s = deg
% save, but DO NOT OVERWRITE
save('../data/DNS_1_revamp_UPDATED_full')
save('../data/DNS_1_revamp_UPDATED','rawdata','reversal_list','Pos_v_Trace');