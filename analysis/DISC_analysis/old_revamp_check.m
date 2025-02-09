% script to check variables across datasets

old = load('../data/BKN_1_revamp_UPDATED_full.mat'); % dataset from before fixed QuickTime Vector
new = load('../data/KEN_1_revampHORRM_full.mat'); % dataset from before after QuickTime Vector

for fncell = {'ppd' 'condition_list' 'duration' 'duration_sec' 'refresh' 'fixduration' 'fixpause' 'sgrowth' 'diagrow' 'smin' 'rmin' 'pursuit_speed' 'pursuit_dist' 'stair_repetitions' 'trials_per_stair' 'nstairs' 'growth_list'} % list of variables to check
    fn = fncell{1}; % convert to string 
    fprintf('%s\n',fn)
    kickout = 0;
    
    % check that variable exists in both datasets
    if ~isfield(old,fn)
        fprintf('\tdoes not exist for old dataset\n')
        kickout = 1;
    end
    if ~isfield(new,fn)
        fprintf('\tdoes not exist for new dataset\n')
        kickout = 1;
    end
    if kickout; continue; end
    
    
    % same value in each dataset?
    if ~issame(old.(fn),new.(fn))
        fprintf('\tvalue not matched across datasets\n')
        disp(old.(fn))
        disp(new.(fn))
    end
        
end