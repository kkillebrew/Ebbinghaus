% script-based expeirmental code for Dynamic Ebbinghaus vs. Dynamic
% Corridor Illusions.

% This code, including some basic stimulus parameters, is largely based off
% of the Dynamic Ebbinghaus code (ebbinghaus_0) and the DISC experiment.

% P2 - using adaptive staircase and forced-choice (grow/shrink) responses to avoid issues with unmotivated subjects


% Experiment Planning
% P1 Dynamic vs. Static Ebbinghaus and Corridor (2 uncertaintly levels) - method of adjustment
% P2 Dynamic vs. Static Corridor, peripheral fixation - method of constant stimuli
% P3 Dynamic vs. Static Ebbinghaus and Corridor (2 uncertaintly levels) - method of constant stimuli
%
%    SEE THE CODE DEFINING THE CONDS STRUCTURE BELOW FOR MORE DETAILS...


clear conds data fix params staircase stim this timing tmp trials; %params conds trials timing data stim this tmp; % clear some variables/structures that are important

AssertOpenGL; % make sure Psychtoolbox-3 is installed properly on this system
rand('twister',sum(100.*clock)); % shuffle the random number generator

addpath('./shared/'); % for shared helper functions

%% parameters
% all necessary parameters should be stored in the params structure for
% organization and easy params.datafile saving/loading.
params.experiment = mfilename; % brief string to describe experiment.  will be part of datafile.  ideally short and only characters (no spaces, underscores, dashes, etc.). 

params.setup = 'eyetracker_emulation';
[kb_id, kb_names] = GetKeyboardIndices;
params.which_screen = max(Screen('Screens')); % use second monitor, if available
switch params.setup
    % N.B. when defining your setup, if you'd like a coding/params.debugging setup
    %      with a less-than-full-screen window, you can either make this
    %      small window show the actual stimulus sizes, or be a scaled down
    %      version of what would appear on a full screen.  For scaled-down
    %      version, set the mon_width_cm based on the FULL monitor, and set
    %      win_size smaller (e.g., [0 0 960 600].  if you want the small
    %      window to show things at the correct size (for the params
    %      defined), set mon_width_cm to reflect the width of the PTB
    %      portion of the monitor (i.e., the width of the small screen as
    %      projected on your monitor).
    case 'ryan_macair'
        % monitor dimensions
        params.mon_width_cm = 28.5; % height = 18 cm
        params.mon_dist_cm = 66;
        
        % which monitor, and what size display?
        params.win_size = [];%[0 0 600 400]; % leave as empty matrix [] for full screen
        %params.which_screen = 0; % override default, if desired
        
        % what keyboard to listen to
        params.kb_dev_id = kb_id(strcmp(kb_names,'Apple Keyboard')); % use -1 for all devices, or specify a specific device

        params.debugging = 1; % controls some features (e.g., HideCursor, dispaly messages, etc.)
    case {'eyetracker' 'eyetracker_emulation'}
        % monitor dimensions
        params.mon_width_cm = 40; % height = X cm
        params.mon_dist_cm = 73;
        params.mon_width_deg = 2*(180/pi)*atan((params.mon_width_cm/2)/params.mon_dist_cm);
        
        % which monitor, and what size display?
        params.win_size = []; % leave as empty matrix [] for full screen
        %params.which_screen = 0; % override default, if desired
        
        % what keyboard to listen to
        params.kb_dev_id = kb_id(strcmp(kb_names,'Apple Keyboard')); % use -1 for all devices, or specify a specific device

        params.debugging = 0; % controls some features (e.g., HideCursor, dispaly messages, etc.)
        
        % N.B. above is the default setting for the eyetracker.  the
        %      following if statement overrides some of these for params.debugging
        %      from other setups (e.g., ryan's laptop), but allowing for
        %      matched stimulus size with eyetracker.
        if strcmp(params.setup,'eyetracker_emulation')
            params.win_size = [0 0 1024 768]; % leave as empty matrix [] for full screen
            params.debugging = 0;%2;
        end
    case 'eeg'
        % monitor dimensions
        params.mon_width_cm = 40; % height 18 cm
        params.mon_dist_cm = 44;
        
        % which monitor, and what size display?
        params.win_size = []; % leave as empty matrix [] for full screen
        %params.which_screen = 0; % override default, if desired
        
        % what keyboard to listen to
        params.kb_dev_id = kb_id(strcmp(kb_names,'Apple Keyboard')); % use -1 for all devices, or specify a specific device

        params.debugging = 0; % controls some features (e.g., HideCursor, dispaly messages, etc.)
    case 'chris_imac'
        % monitor dimensions
        params.mon_width_cm = 59.5; % height = X cm
        params.mon_dist_cm = 57;
        
        % which monitor, and what size display?
        params.win_size = [0 0 1440 900]; % leave as empty matrix [] for full screen
        %params.which_screen = 0; % override default, if desired

        % what keyboard to listen to
        params.kb_dev_id = kb_id(strcmp(kb_names,'Apple Keyboard')); % use -1 for all devices, or specify a specific device

        params.debugging = 1; % controls some features (e.g., HideCursor, dispaly messages, etc.)
    otherwise
        error('unrecognized setup %s',params.setup)
end
params.mon_width_deg = 2*(180/pi)*atan((params.mon_width_cm/2)/params.mon_dist_cm);


% params that determine number of trials
% trials.n = length(params.included_conds) * length(params.standard.width_deg) * length(params.init_anchor_bases) * length(params.targ_context_sizes) * params.repetitions
params.included_conds = [103 113 303 313];%[101 111 201 211 301 311];%[103 113 303 313];%[101 111 201 211 301 311]; % which conditions should we include? defined in conds structure by conds.id (does not direclty match row index)
%%%params.init_anchor_bases = [0 1]; % possible initial anchor points (+/- params.init_anchor_jitter).  0 to start at anchor-small; 1 to start at anchor-large
params.targ_context_sizes = 2;%[1 2]; % possible sizes of inducers surrounding the "target" (opposite of those surrounding "standard").  1 for targ-in-small/standard-in-large[dynamic=small-to-large_growth, like reverse DISC];  2 for targ-in-large/standard-in-small/[dynamic=large-to-small_growth, like DISC]
% N.B. (for targ_context_sizes).  in this experiemnt, the standard is
% always on top.  to match the corridor illusion, the ebbinghaus should be
% presented with small inducers on top (perceived as larger).  this is
% equivilent to "far" condition on top in the corridor illusion (perceived
% as larger).
% N.B. also, the responses for this experiment make it such that there is
% no need for response reversals for dynamic vs. static conditions (as
% we've utilized in past experiments)
%%params.repetitions = 1; % how many repetitions of each condition? THIS SHOULD BE ACCOUNTED FOR IN STAIRCASE
params.break_trials = .10:.10:.90; % proportion of trials after which to allow subject to take a self-timed break.


% fixation spot
params.fix.width_deg = 0.1; % width of circular fixation spot in degrees
params.fix.color = [0 255 0]; % RGB value (0-255) of inner fix spot
% params.fix.innerwidth_deg = 0.1; % width of inner portion of fixation spot in degrees
% params.fix.inner_color = [255 255 0]; % RGB value (0-255) of inner fix spot
% params.fix.outerwidth_deg = 0.25; % width of outer portion of fixation spot in degrees
% params.fix.outer_color = [255 255 255]; % RGB value (0-255) of outer fix spot

params.fix.placement = {'target' 'periphery'}; % 'target', 'periphery' ('none' is implicit if conds.fix_placement_idx = 0)
params.fix.placement_args = {[] -3}; % for 'periphery' eccentricity of target from the halfway point between target and standard  along line perpendicular to trajectory; ignored for 'target' (and the implicit 'none')


% misc
params.instruction_text_size = 24; % size of instructions on screen
params.text_color = [0 0 0]; % RGB value (0-255) of instructional text
params.bg_color   = [128 128 128]; % RGB value (0-255) of background


% response keys
KbName('UnifyKeyNames'); % cross-platform defined key names
params.abort_key  = KbName('Escape'); % pressing escape kills the code
params.return_keys = KbName('Return'); % continuing after a break, or starting experiment, requires pressing the Return key
%%params.v_key = KbName('v'); % for forcing veridical anchor point, if params.debugging
%params.trigger_key = KbName('1!'); % key that serves as scanner trigger
params.resp_keys  = KbName({'O' 'M'}); % subject response key(s) index for KbCheck, etc.  For this experiment, maps to target bigger in the [start end] position (which corresponds to [top bottom] of screen)


% stimulus timing
params.fix_dur = .500; % duraction of fixation-only period at start of trial (only fixation spot visible)
params.cue_dur = .200; % duration of pre-animcation cue period in seconds (stimulus displayed at starting position, but not moving yet)
params.full_cycle_dur = 2;%2.8; % duration of a full dynamic stimulus cycle in seconds (from start position, trans to end position, trans back to start position)
params.half_cycle_dur = params.full_cycle_dur/2; % duration of a half dynamic sitmulus cycle in seconds (from start position, trans to end position)
params.animation_dur = params.half_cycle_dur;%Inf; % duration of the animation period.  If Inf, animation will continue until there is a response
params.resp_dur = Inf; % duration of response window in s.  Set to Inf to always wait until subject responds before continuing.
params.iti = 1;%.500; % inter-trial interval in sec


% some response and feedback options
params.show_targ_during_static_cue  = 1; % should we show target circle during cue period for static trials?  context and standard are always shown.
params.allow_resp_during_static_cue = params.show_targ_during_static_cue; % boolean, should we count responses made during the cue period?  they don't make sense for Dynamic conditions, but are potentially valid for Static conditions, unless we aren't showing the target until after the cue period
params.show_feedback_during_iti     = 1; % boolean, show feedback (correct or incorrect) during ITI?  only makes sense for a set of pre-experiment practice trials.
params.show_fixspot_during_iti      = 1; % boolean, show fixation spot during iti?  works well if fixation spot is static throughout experiment and doesn't change across trial types


% standard/target
params.scale_from_dynebb = 0.75; % a scaling factor applied to all stimulus sizes so that proportions match those used for DynEbb paper, but we can make the stimuli, overall, slightly smaller and avoid overlap of inducers in the static condition
params.standard.width_deg = params.scale_from_dynebb * 2.0; % [min max] possible sizes of the standard circle, chosen from a uniform distribution on a trial-by-trial basis. also the starting size for dynamic conditions.
params.standard.width_jitter = 0.02 * params.standard.width_deg; % maximum absolute size jitter applied on each trial (deg)
params.standard.color = [0 0 0]; % RGB value (0-255) of default (or normtarg) standard (and target)


% standard/target-related uncertainty manipulations
% % params.uncertainty.lowcontrast = params.bg_color + 0.05 * (params.standard.color-params.bg_color); % RGB value (0-255 of low contrast target.  may be a scaled version of params.standard.color vs. params.bg_color
% % params.uncertainty.blurred     = [1.4 0];%[1.6 .25 1]; % params to determine how to generate blurred target.  [sigma1 sigma2].  mask is difference of 2 gaussians (1-2), with center portion adjusted to be a flat plateau
% % params.uncertainty.isoluminant = [255 0 0]; % RGB value (0-255) of isoluminant standard/target (should come from a pre-saved file after a flicker-fusion self-adjustment paradigm)
params.uncertainty.pos_jitter      = params.scale_from_dynebb * .06; % positional jitter in deg, applied randomly frame-by-frame, for the target position


% Ebbinghaus inducers
params.ebbinghaus.color = [0 0 0]; % RGB value (0-255) of all inducers
params.ebbinghaus.small_width_deg = params.scale_from_dynebb * 0.7; % size of 'small' inducers
params.ebbinghaus.large_width_deg = params.scale_from_dynebb * 4.25; % size of 'large' inducers
params.ebbinghaus.n = 6; % how many (equally-spaced, starting at right horizontal position) inducers?
% for inducer eccen: 
%     large_eccen==small_eccen for constant center eccen (i.e., shrinking edge eccen)
%     large_eccen==small_eccen+(large_width-small_width)/2 for constant inner edge eccen
%     large_eccen > small_eccen+(large_width/2-small_width)/2 for growing edge eccen
params.ebbinghaus.small_eccen_deg = params.scale_from_dynebb * 1.5; %0.75 * min(params.standard.width_deg); % eccentricity of inducers when they are small
params.ebbinghaus.large_eccen_deg = params.scale_from_dynebb * 5; %params.ebbinghaus.small_eccen_deg + (params.ebbinghaus.large_width_deg-params.ebbinghaus.small_width_deg)/2;% + mean(params.standard.width_deg);% eccentricity of inducers when they are large


% Corridor background
switch params.setup
    % file used for corridor needs to be the appropriate dimensions for the
    % setup/monitor (based on its resolution)
    case 'ryan_macair'
        % 1440 x 900
        params.corridor.file = 'hallway_1440x900.png'; % image file for corridor background
    case {'eyetracker' 'eyetracker_emulation'}
        % 1024 x 768
        params.corridor.file = 'hallway_1024x768.png'; % image file for corridor background
        %params.corridor.file = 'hallway_1440x900.png'; % image file for corridor background

    otherwise
        error('an appropriate corridor background image has not been defined for setup %s',params.setup)
end
params.corridor.width_deg ='fullscreen'; %24.3673; %'fullscreen'; % width of corridor background in deg, or 'fullscreen'


% stimulus position and dynamics
%params.start_pos_deg = [4.85 4]; % relative to screen center (for ryan_macair 1440x900) PILOTING ONLY, MAYBE WSU
params.start_pos_deg = [4.7 4.5]; % relative to screen center (for UNR eyetracker setup 1024x768)
params.trans_ang_rad  = -(90+20.1) * pi/180;  % angle of full stimulus translation (for dynamic version) in radians (matlab: 0 right, + counter-clockwise)
params.trans_dist_deg = 8;%7;%params.ebbinghaus.large_eccen_deg-params.ebbinghaus.small_eccen_deg;    % distance of full stimulus translation (for dynamic version) in degrees (relative to start_pos_deg)
% N.B. you can set trans_ang_rad and trans_dist_deg (above) such that one
%      inducer will be stationary (yet growing).  could be good for an
%      inducer-fixation condition.
% e.g., params.trans_ang_rad = -(2*pi/params.ebbinghaus.n);
%       params.trans_dist_deg = params.ebbinghaus.large_eccen_deg-params.ebbinghaus.small_eccen_deg;
[tmp.x, tmp.y] = pol2cart(params.trans_ang_rad,params.trans_dist_deg); % match end position (for static conditions) to translation in dynamic conditions
params.end_pos_deg = params.start_pos_deg + [tmp.x tmp.y]; % position of isolated standard (when present) relative to screen center
params.global_pos_jitter_deg = .2; % maximum radius of circular position jitter for stimulus placement on each trial *** does not include a shift in the corridor image, so as to avoid adaptation of lines of the corridor


% staircase params
params.staircase.length      = 15;%25; % length, in trials, of each interleaved staircase
params.staircase.repetitions = 3; % how many repetitions of each staircase (condition X targ_context_sizes X anchor)
params.staircase.targ_change_props = linspace(-1.80,1.80,40); % actual target changes (as a proportion of the standard) to select from during staircase procedure
params.staircase.cap_targ_change_props = 1; % boolean, should we limit the possible targ_change_prop values so that the target never touches the inducers, or the fidation spot?
if params.staircase.cap_targ_change_props
    % determine a good range of potential targ_change props, such that target
    % should NEVER overlap with the inducers
    tmp.min_small_inducer_edge_eccen = params.ebbinghaus.small_eccen_deg - 0.5*params.ebbinghaus.small_width_deg; % the minimum eccentricity to touch an inducer's edge (i.e., the largetst that a target should EVER be)
    tmp.min_large_inducer_edge_eccen = params.ebbinghaus.large_eccen_deg - 0.5*params.ebbinghaus.large_width_deg; % the minimum eccentricity to touch an inducer's edge (i.e., the largetst that a target should EVER be)
    if length(params.targ_context_sizes) > 1
        % targ-in-small AND targ-in-large conditions interleaved - find min
        % position
        tmp.min_inducer_edge_eccen = min(tmp.min_small_inducer_edge_eccen,tmp.min_large_inducer_edge_eccen);
    elseif params.targ_context_sizes == 1
        % targ-in-small condition only
        tmp.min_inducer_edge_eccen =  tmp.min_small_inducer_edge_eccen;
    elseif params.targ_context_sizes == 2
        % targ-in-large condition only
        tmp.min_inducer_edge_eccen =  tmp.min_large_inducer_edge_eccen;
    else
        error('umm, not sure how, but cannot determine min_inducer_edge_eccen')
    end
    
    tmp.max_standard_halfwidth = (params.standard.width_deg + params.standard.width_jitter)/2; % largest possible halfwidth (radius) for a standard (accounting for max jitter)
    params.staircase.targ_change_prop_upper_limit = (tmp.min_inducer_edge_eccen - tmp.max_standard_halfwidth) / tmp.max_standard_halfwidth; % largest target halfwidth cannot exceed the inner edge of inducers
    
    %%this.targ_halfwidth_pix = this.standard_halfwidth_pix + this.cur_targ_change_prop*this.standard_halfwidth_pix; % size of target (at end pos) in pix
    
    tmp.min_standard_halfwidth = (params.standard.width_deg - params.standard.width_jitter)/2; % smallest possible halfwidth (radius) for a standard (accounting for max jitter)
    params.staircase.targ_change_prop_lower_limit = (params.fix.width_deg - tmp.min_standard_halfwidth) / tmp.min_standard_halfwidth; % smallest target halfwidth cannot be less than twice the fixation spot (min_halfwidth = fix_fullwidth)
else
    % not capping targ_change_prop values 
    params.staircase.targ_change_prop_upper_limit = Inf;
    params.staircase.targ_change_prop_lower_limit = -Inf;
end
tmp.ref = and(lt(params.staircase.targ_change_props,params.staircase.targ_change_prop_upper_limit),gt(params.staircase.targ_change_props,params.staircase.targ_change_prop_lower_limit));
params.staircase.targ_change_props = params.staircase.targ_change_props(tmp.ref);

params.staircase.n_targ_change_props = length(params.staircase.targ_change_props); % total number of different target change levels (i.e., resolution of psychometric curves)
params.staircase.anchors_idx = [1 params.staircase.n_targ_change_props]; % each staircase starts from one of the extremes
% N.B. params.staircase.steps_by_iter.iters/steps allows the experiment to
% use larger steps earlier in the staircase, and then smaller steps later
% in the staircase when the participant is presumably close to their PSE.
% Thus, we can have a large number params.staircase.n_targ_change_props
% and not require huge staircase lengths.
% params.staircase.steps_by_iter.iters is the max iterations at which
% params.staircase.steps_by_iter.steps is the number of steps moved along
% the adaptive staircase.
params.staircase.steps_by_iter.iters = [5  10 params.staircase.length]; %[5 15 params.staircase.length]; % values must increase with the last value should be params.staircase.length
params.staircase.steps_by_iter.steps = [5  2  1]; % values must decrease and the last value should be 1
if any(diff(params.staircase.steps_by_iter.iters) <= 0)
   error('params.staircase.steps_by_iter.iters must be a list of increasing values, without repeats')
end
if params.staircase.steps_by_iter.iters(end) ~= params.staircase.length
    error('the last value of params.staircase.steps_by_iter.iters must be params.staircase.length')
end
if any(diff(params.staircase.steps_by_iter.steps) >= 0)
    error('params.staircase.steps_by_iter.steps must be a list of decreasing values, without repeats')
end
if params.staircase.steps_by_iter.steps(end) ~= 1
    error('the last value of params.staircase.steps_by_iter.steps should be 1 for the best final resolution.  if this is undesired, you will need to adjust the code to remove this error message.')
end
if length(params.staircase.steps_by_iter.iters) ~= length(params.staircase.steps_by_iter.steps)
    error('params.staircase.steps_by_iter.steps and params.staircase.steps_by_iter.steps must be the same length')
end

%% validate params
% i.e., check for invalid stimulus configurations based on current
% parameters and throw errors if necessary

% verify that inducers are appropriately proportioned relative to standard
if max(params.standard.width_deg) >= params.ebbinghaus.large_width_deg
   error('''large'' inducers (%.02f) are smaller than largest possible standard (%.02f)',params.ebbinghaus.large_width_deg,max(params.standard.width_deg));
end
if min(params.standard.width_deg) <= params.ebbinghaus.small_width_deg
   error('''small'' inducers (%.02f) are larger than smallest possible standard (%.02f)',params.ebbinghaus.small_width_deg,min(params.standard.width_deg));
end


% check for target-inducer overlap in DYNAMIC condition
% (unless we are capping the targ max size on a per-trial basis)
%%if ~params.cap_max_targ_per_trial && max(params.standard.width_deg)/2 > params.ebbinghaus.small_eccen_deg-(params.ebbinghaus.small_width_deg/2)
%%    error('largest standard target (%.02f) will overlap with small inducers (%.02f) at their current eccen (%.02f)',max(params.standard.width_deg),params.ebbinghaus.small_width_deg,params.ebbinghaus.small_eccen_deg)
%%end
%%if ~params.cap_max_targ_per_trial && params.targ.max_dynamic_width_deg/2 > params.ebbinghaus.large_eccen_deg-(params.ebbinghaus.large_width_deg/2)
%%    error('largest target (%.02f) will overlap with large inducers (%.02f) at their current eccen (%.02f) for dynamic conditions',params.targ.max_dynamic_width_deg,params.ebbinghaus.large_width_deg,params.ebbinghaus.large_eccen_deg)
%%end
%%% additional check for target-inducer overlap in STATIC condition
%%if ~params.cap_max_targ_per_trial && params.targ.max_static_width_deg/2 > params.ebbinghaus.small_eccen_deg-(params.ebbinghaus.small_width_deg/2)
%%    error('largest target (%.02f) will overlap with small inducers (%.02f) at their current eccen (%.02f) for static conditions',params.targ.max_dynamic_width_deg,params.ebbinghaus.small_width_deg,params.ebbinghaus.small_eccen_deg)
%%end
%

% verify that we can't exceed screen dimensions
% N.B. we figure out the max extent of stimulus in degrees, and will test
%      against scree dimension once we OpenWindow
tmp.max_radius = params.ebbinghaus.large_eccen_deg + params.ebbinghaus.large_width_deg/2; % max inducer eccen + 1/2 max inducer size
tmp.min_radius = params.ebbinghaus.small_eccen_deg + params.ebbinghaus.small_width_deg/2; % min inducer eccen + 1/2 min inducer size

tmp.hi_start_y = params.start_pos_deg(2) + params.global_pos_jitter_deg; % start y + max jitter
tmp.le_start_x = params.start_pos_deg(1) - params.global_pos_jitter_deg; % start x - max jitter

[tmp.trans_x, tmp.trans_y] = pol2cart(params.trans_ang_rad,params.trans_dist_deg); % translation vector (deg)
tmp.lo_end_y  = params.start_pos_deg(2) - tmp.trans_y - params.global_pos_jitter_deg; % start y - trans y - max jitter

tmp.ri_end_x = params.end_pos_deg(1) + params.global_pos_jitter_deg; % start x - max jitter

params.upper_stimulus_extent = tmp.hi_start_y + tmp.max_radius; % hi_start_y + max_radius (deg)
params.left_stimulus_extent  = tmp.le_start_x - tmp.max_radius; % le_start_y - max_radius (deg)
params.lower_stimulus_extent = tmp.lo_end_y - tmp.max_radius; % low_start_y - max_radius (deg)
params.right_stimulus_extent = tmp.ri_end_x + tmp.max_radius; % ri_start_y - max_radius (deg)


% *** needs adjustment for position of standard and target
% verify that static condition won't ever overlap (given jitter)
params.target_rightmost_edge  = params.start_pos_deg(1) + params.global_pos_jitter_deg + tmp.min_radius; % targ-in-small
params.standard_leftmost_edge = params.end_pos_deg(1) - params.global_pos_jitter_deg - tmp.max_radius; % standard-in-large
if params.target_rightmost_edge > params.standard_leftmost_edge
% ***    error('there is a chance that there will be overlap between the target and standard stimuli in the static condition')
end



%% datafile
% datfiles will be <subjid>_<experimentID>_<datecode>_<fileNumber>
c = clock;
params.time_stamp = sprintf('%02d/%02d/%04d %02d:%02d:%02.0f',c(2),c(3),c(1),c(4),c(5),c(6)); % month/day/year hour:min:sec
params.datecode = datestr(now,'mmddyy');

% get input
params.subjid = 'practice';%input('Enter Subject Code:','s');
params.runid  = 1;%input('Enter Run (reset daily):');

% determine datafile
params.datafile = sprintf('%s_%s_%s_%03d',params.subjid,params.experiment,params.datecode,params.runid);
params.datadir = '../data/';

% % check to see if this file exists
% if exist(fullfile(params.datadir,[params.datafile '.mat']),'file')
%     tmpfile = input('File exists.  Overwrite? y/n:','s');
%     while ~ismember(tmpfile,{'n' 'y'})
%         tmpfile = input('Invalid choice. File exists.  Overwrite? y/n:','s');
%     end
%     if strcmp(tmpfile,'n')
%         display('Bye-bye...');
%         return; % will need to start over for new input
%     end
% end


%% prepare psychtoolbox
% Opens psychtoolbox window
if params.debugging
    Screen('Preference', 'SkipSyncTests', 1);
else
    Screen('Preference', 'SkipSyncTests', 0);
end
[w, rect] = Screen('OpenWindow',params.which_screen,params.bg_color,params.win_size,[],[],[],2); % last argument is multisample - set above zero (~4?) for anti-aliasing using brute-force multisampling method
% turned off for params.debugging
if ~params.debugging
    HideCursor; % Hide mouse cursor from screen
    ListenChar(2); % stops keypresses during experiment
    Priority(MaxPriority(w)); % set to Maximum priority level for the current operating system
end

% measure the frame rate
params.frame_rate = Screen('FrameRate',w); % in seconds.  this is what the operating system is set to, does not work on some systems (e.g., lcd or laptops)
params.flip_interval = Screen('GetFlipInterval',w); % in seconds.  this is an actual measurement, should work on all systems
%params.flip_interval_correction = .80 * params.flip_interval; % this should work even on laptops that don't return a FrameRate value

% Enable alpha blending for contrast manipulations (like images with alpha channels) and anti-aliasing
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% set center of screen
params.x_resolution = rect(3); 
params.x_center = params.x_resolution/2;
xc = params.x_center; % for easy reference
xres = params.x_resolution;
params.y_resolution = rect(4); 
params.y_center = params.y_resolution/2;
yc = params.y_center; % for easy reference
yres = params.y_resolution;

% establish size and ppd, need to change for crt
params.pix_per_deg = (rect(3)/params.mon_width_deg);
ppd = params.pix_per_deg; % for easy reference

% draw to screen so we know things are working...
AddText(w, 'loading and verifying, please wait...', [xc yc], params.text_color, 0, params.instruction_text_size);
Screen(w,'Flip');


% check screen resolution against stimulus extent
err_dir = '';
if (yc - params.upper_stimulus_extent*ppd) < 0
    err_dir = [err_dir ' & upper'];
end
if (xc - params.left_stimulus_extent*ppd) < 0
    err_dir = [err_dir ' & left'];
end
if (yc - params.lower_stimulus_extent*ppd) > yres
    err_dir = [err_dir ' & lower'];
end
if (xc - params.right_stimulus_extent*ppd) > xres
    err_dir = [err_dir ' & right'];
end
% ***&
% if ~strcmp(err_dir,'')
%     KbQueueRelease(params.kb_dev_id); % release keypress queue or process may continue to run in background
%     ShowCursor;
%     ListenChar(0);
%     Priority(0);
%     Screen('CloseAll')
%     error('it appears that the stimulus may extend beyond the %s edge(s) of the screen',err_dir(4:end))
% end


%% convert sizes in deg to pixels (for fixed variables)
%if strcmp(params.targ.min_width_deg,'one_pixel')
%    % now that we know the ppd, we can figure out the minimum size of
%    % target (1 pixel) in degrees.  overwrite params so the datafile is
%    % accurate.
%    params.targ.min_width_deg = 1/ppd;
%end
% and define rects here so we don't have to keep replicating these long bits of code
fix.width_pix     = params.fix.width_deg * ppd;
fix.halfwidth_pix = fix.width_pix / 2;
fix.rect = [xc-fix.halfwidth_pix, yc-fix.halfwidth_pix, xc+fix.halfwidth_pix, yc+fix.halfwidth_pix]; % rect for drawing circular fixation spot
% fix.innerwidth_pix     = params.fix.innerwidth_deg * ppd;
% fix.innerhalfwidth_pix = fix.innerwidth_pix / 2;
% fix.inner_rect = [xc-fix.innerhalfwidth_pix, yc-fix.innerhalfwidth_pix, xc+fix.innerhalfwidth_pix, yc+fix.innerhalfwidth_pix]; % rect for drawing inner portion of fixation spot
% fix.outerwidth_pix = params.fix.outerwidth_deg * ppd;
% fix.outerhalfwidth_pix = fix.outerwidth_pix / 2;
% fix.outer_rect = [xc-fix.outerhalfwidth_pix, yc-fix.outerhalfwidth_pix, xc+fix.outerhalfwidth_pix, yc+fix.outerhalfwidth_pix]; % rect for drawing outer portion of fixation spot


stim.standard.width_pix = params.standard.width_deg * ppd;
stim.standard.halfwidth_pix = stim.standard.width_pix/2;
%%stim.targ.min_width_pix = params.targ.min_width_deg * ppd;
%%stim.targ.min_halfwidth_pix = stim.targ.min_width_pix/2;
%%stim.targ.max_width_pix = params.targ.max_width_deg * ppd;
%%stim.targ.max_halfwidth_pix = stim.targ.max_width_pix/2;
stim.inducers.small_width_pix = params.ebbinghaus.small_width_deg * ppd;
stim.inducers.small_halfwidth_pix = stim.inducers.small_width_pix/2;
stim.inducers.large_width_pix = params.ebbinghaus.large_width_deg * ppd;
stim.inducers.large_halfwidth_pix = stim.inducers.large_width_pix/2;
stim.inducers.small_eccen_pix = params.ebbinghaus.small_eccen_deg * ppd;
stim.inducers.large_eccen_pix = params.ebbinghaus.large_eccen_deg * ppd;

stim.inducers.angs_rad = linspace(0,2*pi,params.ebbinghaus.n+1);
stim.inducers.angs_rad = stim.inducers.angs_rad(1:end-1); % radians

if strcmp(params.corridor.width_deg,'fullscreen')
    stim.corridor.width_pix = params.x_resolution;
else
    stim.corridor.width_pix = params.corridor.width_deg * ppd;
end
stim.corridor.halfwidth_pix = stim.corridor.width_pix/2;
stim.corridor.halfheight_pix = stim.corridor.halfwidth_pix * (900/1440);
tmp.corridor_img = imread(params.corridor.file);
stim.corridor.texid = Screen('MakeTexture',w,tmp.corridor_img);

stim.start_pos_pix = (params.start_pos_deg .* [1 -1]) * ppd; % relative to screen center ([xc yc]) (be sure to flip y dimension for PTB)
[tmp.x,tmp.y] = pol2cart(params.trans_ang_rad,ppd*params.trans_dist_deg);
stim.trans_pos_pix  = [tmp.x tmp.y].*[1 -1]; % relative to starting position (be sure to flip y dimension)
stim.end_pos_pix = (params.end_pos_deg .* [1 -1]) * ppd; % relative to screen center ([xc yc]) (be sure to flip y dimension)



% for a fixation placement in the periphery, we need to define a position
% in space relative to the target/standard positions.  User provides the
% eccentricity from the point half way between target/standard trajectory,
% along a perpendicular path.
tmp.traj_center_point_pos_pix = (stim.start_pos_pix+stim.end_pos_pix)/2; % point half way between target and standard

tmp.periph_eccen_pix = params.fix.placement_args{2} * ppd; % eccentricity of target from the halfway point between target and standard along line perpendicular to trajectory
tmp.periph_ang_rad = params.trans_ang_rad + (pi/2); % angle from horizontal line through start_pos down to line connecting fix point and traj_center_point
tmp.periph_x_from_traj_center_point = tmp.periph_eccen_pix * cos(tmp.periph_ang_rad); % x translation of peripheral fixation from the center point of the trajectory
tmp.periph_y_from_traj_center_point = tmp.periph_eccen_pix * sin(tmp.periph_ang_rad); % y translation of peripheral fixation from the center point of the trajectory
fix.placement.periph_pos_pix = tmp.traj_center_point_pos_pix + [tmp.periph_x_from_traj_center_point -tmp.periph_y_from_traj_center_point]; % relative to screen center ([xc yc]) (be sure to flip y dimension for PTB)



%% conds structure
% conds structure defines various parameters for a given trial-type, or condition.
% much of this structure is hardcoded, or at least the number of conditions
% and how they are defined.
% the trials structure can point to this structure for individual trials
conds.n   = []; % how many total conds for this run?
conds.row = []; % row index for each condition
conds.id  = [101:104,111:114 201:204,211:214 301:304,311:314]'; % number id for each condition (which doesn't have to match rows)
%     101. No Context - Static  - FixTargs
%     102.                      - FixTargs - Jittered
%     103.                      - FixPeriph
%     104.                      - FixPeriph - Jittered
%     111.            - Dynamic - FixTargs
%     112.                      - FixTargs - Jittered
%     113.                      - FixPeriph
%     114.                      - FixPeriph - Jittered
%
%     201. Ebbinghaus - Static  - FixTargs
%     202.                      - FixTargs - Jittered
%     203.                      - FixPeriph
%     204.                      - FixPeriph - Jittered
%     211.            - Dynamic - FixTargs
%     212.                      - FixTargs - Jittered
%     213.                      - FixPeriph
%     214.                      - FixPeriph - Jittered
%
%     301. Corridor   - Static  - FixTargs
%     302.                      - FixTargs - Jittered
%     303.                      - FixPeriph
%     304.                      - FixPeriph - Jittered
%     311.            - Dynamic - FixTargs
%     312.                      - FixTargs - Jittered
%     313.                      - FixPeriph
%     314.                      - FixPeriph - Jittered
conds.label = { % string label for each condition
    'static.normtarg.fixtarg.nocontext'
    'static.jittered.fixtarg.nocontext'
    'static.normtarg.fixperiph.nocontext'
    'static.jittered.fixperiph.nocontext'
    'dynamic.normtarg.fixtarg.nocontext'
    'dynamic.jittered.fixtarg.nocontext'
    'dynamic.normtarg.fixperiph.nocontext'
    'dynamic.jittered.fixperiph.nocontext'
    ...
    'static.normtarg.fixtarg.ebbinghaus'
    'static.jittered.fixtarg.ebbinghaus'
    'static.normtarg.fixperiph.ebbinghaus'
    'static.jittered.fixperiph.ebbinghaus'
    'dynamic.normtarg.fixtarg.ebbinghaus'
    'dynamic.jittered.fixtarg.ebbinghaus'
    'dynamic.normtarg.fixperiph.ebbinghaus'
    'dynamic.jittered.fixperiph.ebbinghaus'
    ...
    'static.normtarg.fixtarg.corridor'
    'static.jittered.fixtarg.corridor'
    'static.normtarg.fixperiph.corridor'
    'static.jittered.fixperiph.corridor'
    'dynamic.normtarg.fixtarg.corridor'
    'dynamic.jittered.fixtarg.corridor'
    'dynamic.normtarg.fixperiph.corridor'
    'dynamic.jittered.fixperiph.corridor'
    };
conds.short_label = conds.label;
conds.short_label = regexprep(conds.short_label,'static','S');
conds.short_label = regexprep(conds.short_label,'dynamic','D');
conds.short_label = regexprep(conds.short_label,'normtarg','NT');
conds.short_label = regexprep(conds.short_label,'jittered','J');
conds.short_label = regexprep(conds.short_label,'fixtarg','FT');
conds.short_label = regexprep(conds.short_label,'tracktarg','TT');
conds.short_label = regexprep(conds.short_label,'fixperiph','FP');
conds.short_label = regexprep(conds.short_label,'nocontext','NC');
conds.short_label = regexprep(conds.short_label,'ebbinghaus','E');
conds.short_label = regexprep(conds.short_label,'corridor','C');

conds.is_dynamic     = repmat([zeros(4,1); ones(4,1)],3,1); % does this condition have dynamic inducers? else static.
conds.is_translating = repmat([repmat(-1,[4,1]); ones(4,1)],3,1); % does this condition translate as a whole? 0=stationary at starting position; 1=translating between starting and ending positions; -1=stationary at end position (for Static conditions)
conds.what_context   = [zeros(8,1); ones(8,1); 2*ones(8,1)]; % what is the context?  0=none, 1=ebbinghaus inducers, 2=corridor

conds.is_jittered    = repmat([0 1]',12,1); % is target phsically jittered?

conds.fix_placement_idx  = repmat([1 1 2 2]',6,1); % 1=target; 2=periphery(empty space); 3=inducer(upperleft); 4=inducer(lowerright); 

% update n and rows
conds.n = length(conds.id);
conds.row = (1:conds.n)'; 

% verify that things are all equal length
for fncell = fieldnames(conds)'
    fn = fncell{1}; % convert to string
    if size(conds.(fn),1) ~= conds.n && ~ismember(fn,{'n'})
        error('it appears that conds.%s is not the appropriate length for the total number of conditions defined',fn)
    end
end


%% staircase structure
% staircase structure contains information to help traack the progress through each of the possible interleaved staircases
% for each 2-dim matrix:
%   rows = condition
%   cols = staircase_idx (includes all anchors and repetitions)
staircase.n = NaN; % how many total staircases will there be?
staircase.idx = []; % id for each staircase (redundant with inherent row organization)

% unique staircases: cross params.included_conds, params.targ_context_sizes, and params.staircase.anchors_idx
[tmp.cond, tmp.targ_context_size, tmp.anchor_idx] = ndgrid(params.included_conds,params.targ_context_sizes,params.staircase.anchors_idx);
staircase.conds_id = tmp.cond(:); % which condition (conds.id) is this trial? (N.B. NOT an index into conds structure)
staircase.conds_row = Remap(staircase.conds_id,conds.id,conds.row); % which condition, as a row index into conds structure
staircase.targ_context_size = tmp.targ_context_size(:);
staircase.anchor_idx  = tmp.anchor_idx(:); % starting target change proportion index for each interleaved staircase
% account for repetitions of each unique staircase (params.staircase.repetitions)
for fncell = fieldnames(staircase)'
    fn = fncell{1}; % convert to string
    if ~ismember(fn,{'n' 'idx'})
        staircase.(fn) = repmat(staircase.(fn),params.staircase.repetitions,1);
    end
end
% update n and id
staircase.n  = length(staircase.conds_id);
staircase.idx = (1:staircase.n)';
% the following will be updated during the experiment
staircase.iteration = zeros(staircase.n,1); % will track how many trials for each interleaved staircase have been run
staircase.cur_targ_change_prop_idx = staircase.anchor_idx; % current target change proportion for each interleaved staircase
staircase.targ_change_prop_idx_seq = NaN(staircase.n,params.staircase.length); % sequence of targ_change_prop_idx for each staircase (init to NaNs, fill as trials come up)


%% trials structure
% initialize trial structure, which will hold all of the trial-specific information in the order that it will be presented
% because of the staircase procedure, much of the trials structure will be filled in as the experiment proceeds
trials.n = staircase.n * params.staircase.length; % total number of trials (we know this upfront from staircase.n and params.staircase.length)
trials.order = (1:trials.n)'; % an ordered list of trial id for the actual presentation order

% psuedo-randomly assign the order of the staircases, with each shown before repetitions
trials.staircase_idx = []; % row index into the corresponding row of the staircase structure
trials.staircase_iter = []; % for the current staircase (trials.staircase_idx), what iteration is this?
for fncell = {'conds_id' 'conds_row' 'targ_context_size'}
  fn = fncell{1}; % convert to string
  trials.(fn) = []; % init
end
for i = 1:params.staircase.length
  this_ref = Shuffle(staircase.idx);
  
  trials.staircase_idx = cat(1,trials.staircase_idx,staircase.idx(this_ref));
  for fncell = {'conds_id' 'conds_row' 'targ_context_size'}
    fn = fncell{1}; % convert to string
    trials.(fn) = cat(1,trials.(fn),staircase.(fn)(this_ref,:));
  end
  trials.staircase_iter = cat(1,trials.staircase_iter,i*ones(staircase.n,1));
end

% randomly filled trial-wise params
trials.standard_width_deg = params.standard.width_deg + params.standard.width_jitter * (2*rand(trials.n,1)-1); % size of standard circle in deg for this trial (accounting for jitter)
trials.standard_halfwidth_pix = (trials.standard_width_deg / 2) * ppd; % size of standard circle (halfwidth) in pix for this trial

trials.global_pos_jitter_eccen_deg    = params.global_pos_jitter_deg*rand([trials.n 1]); % radial jitter of stimulus placement for this trial, in degrees
trials.global_pos_jitter_ang_rad      = 2*pi*rand([trials.n 1]); % angle of jitter of stimulus placement for this trial, in radians
[tmp.x,tmp.y] = pol2cart(trials.global_pos_jitter_ang_rad,trials.global_pos_jitter_eccen_deg); % deg
trials.global_pos_jitter_pix  = ppd * bsxfun(@times,[tmp.x tmp.y],[1 -1]); % x,y jitter of stimulus placement in pixels - flip y when going from degrees to pixels

trials.targ_change_prop = NaN(trials.n,1); % targ_change_prop for each trial. initialize to NaNs - it will be filled in based on the staircase procedure

% verify that things are all equal length
for fncell = fieldnames(trials)'
    fn = fncell{1}; % convert to string
    if size(trials.(fn),1) ~= trials.n && ~ismember(fn,{'n'})
        error('it appears that trials.%s is not the appropriate length for the total number of trials defined',fn)
    end
end
    

%% data structure
% initialize and preallocate memory for data structure, which will hold
% trials-specific data like subjects resp, rts and status (correct/incorrect)
nanfill       = NaN(trials.n,1);
data.resp     = nanfill; % keycode of response button pressed
data.respidx  = nanfill; % index into params.resp_keys of (first) button pressed
data.resptime = nanfill; % response time relative to KbQueue flush (should be cue onset, when stimulus appears and a response would be legitimate for static trials)
data.valid    = nanfill; % valid (1, i.e., completed) or not (0 or NaN)
data.status   = nanfill; % correct (1) or incorrect (0) (but note that this is not a useful performance metric for estimating illusion magnitudes)


%% timing structure
% initialize and preallocate memory for timing structure, which will hold
% timing information for all phases of the run
for scell = {'run'} % one per run
    s = scell{1}; % convert to string
    timing.([s '_start'])  = NaN;
    timing.([s '_stop'])   = NaN;
    timing.([s '_dur'])    = NaN;
end
for scell = {'trial' 'fix' 'cue' 'anim' 'respwin' 'iti'} % one per trial
    s = scell{1}; % convert to string
    timing.([s '_start'])  = nanfill;
    timing.([s '_stop'])   = nanfill;
    timing.([s '_dur'])    = nanfill;
end
breaknanfill = NaN(length(params.break_trials),1);
for scell = {'break'} % one per break period
    s = scell{1}; % convert to string
    timing.([s '_start'])  = breaknanfill;
    timing.([s '_stop'])   = breaknanfill;
    timing.([s '_dur'])    = breaknanfill;
end


 %% setup keypress queue
 % better timing than KbCheck and doesn't require while loops (i.e., works in background)
 allowed_keys = zeros(1,256); % by default, keys are not allowed to influence queue
 allowed_keys([params.resp_keys params.abort_key]) = 1; % only listen for valid response keys and the abort key
 KbQueueCreate(params.kb_dev_id,allowed_keys);
 KbQueueStart(params.kb_dev_id); % start listening for keypresses


%% instructions
% Display instructions
AddText(w, 'On each trial, do your best to indicate', [xc yc], params.text_color, 4.1, params.instruction_text_size);
AddText(w, 'which target circle is larger:', [xc yc], params.text_color, 3, params.instruction_text_size);
AddText(w, 'If top, press ''O''     If bottom, press ''M''', [xc yc], params.text_color, 2, params.instruction_text_size);
AddText(w, 'Please maintain gaze on the green fixation point(s) at all times.', [xc yc], params.text_color, 0, params.instruction_text_size);
AddText(w, 'Press the Return key to begin', [xc yc], params.text_color, -3, round(.80 * params.instruction_text_size)); % slightly smaller...
Screen(w,'Flip');
% wait for subject to press the return key
[~, keycode] = KbPressWait(params.kb_dev_id);
while ~any(keycode(params.return_keys))
    [~, keycode] = KbPressWait(params.kb_dev_id);
end


%% trial loop
timing.run_start = Screen('Flip',w);
abort = 0;
for it = 1:trials.n;
    if abort; break; end; % kick out of trial loop if the abort key was pressed on the last trial
    
    
    % give subject self-timed breaks at certain trials...
    this.break_idx = find(round(params.break_trials * trials.n) == it); % empty if not a break trial
    if this.break_idx
        % display break message
        AddText(w, sprintf('You have completed %d%% of the trials.',round(params.break_trials(this.break_idx)*100)), [xc yc], params.text_color, 1, params.instruction_text_size);
        timing.break_start(this.break_idx,1) = Screen('Flip',w,[],1); % don't clear screen so we can write more
        WaitSecs(1);
        AddText(w, 'Press the Return key when you are ready to continue.', [xc yc], params.text_color, -1, params.instruction_text_size);
        Screen('Flip',w);
        % wait for subject to press the return key
        [~, keycode] = KbPressWait(params.kb_dev_id);
        while ~any(keycode(params.return_keys))
            [~, keycode] = KbPressWait(params.kb_dev_id);
        end
        timing.break_stop(this.break_idx,1) = Screen('Flip',w); % clear screen
    end
    
    
    % mark trial start time
    this.zero_time = GetSecs;
    timing.trial_start(it,1) = this.zero_time;
    
    % extract any trial-specific params into "this" structure
    for fncell = {'conds_id' 'conds_row' 'standard_width_deg' 'standard_halfwidth_pix' 'targ_context_size' 'global_pos_jitter_pix' 'staircase_idx' 'staircase_iter'}
        fn = fncell{1}; % convert to string
        this.(fn) = trials.(fn)(it,:);
    end
    

    % get current trial params based on current staircase
    staircase.iteration(this.staircase_idx,1) = staircase.iteration(this.staircase_idx,1)+1; % update staircase iteration
    if this.staircase_iter ~= staircase.iteration(this.staircase_idx,1)
      KbQueueRelease(params.kb_dev_id); % release keypress queue or process may continue to run in background
      ShowCursor; ListenChar(0); Priority(0); Screen('CloseAll');
      error('hmm, it seems that the trials and staircase structure do not agree on which iteration we are on...')
    end
    this.cur_targ_change_prop_idx = staircase.cur_targ_change_prop_idx(this.staircase_idx);
    this.cur_targ_change_prop = params.staircase.targ_change_props(this.cur_targ_change_prop_idx); % actual target change proportion
    this.targ_halfwidth_pix = this.standard_halfwidth_pix + this.cur_targ_change_prop*this.standard_halfwidth_pix; % size of target (at end pos) in pix
    staircase.targ_change_prop_idx_seq(this.staircase_idx,this.staircase_iter) = this.cur_targ_change_prop_idx; % update staircase's record of target prop change idx

    trials.targ_change_prop(it,1) = this.cur_targ_change_prop; % store this value in the trials structure for easy access during data analysis
    
    % set starting and ending positions accounting for global jitter
    this.stim_start_pos_pix = stim.start_pos_pix + this.global_pos_jitter_pix + [xc yc]; % start pos + jitter (absolute coords)
    this.stim_end_pos_pix   = stim.end_pos_pix + this.global_pos_jitter_pix + [xc yc]; % end pos + jitter (absolute coords)

    % initialize variables needed for draw_fix and draw_stim    
    this.standard_pos_pix = this.stim_start_pos_pix; % starting/standard pos
    if ~conds.is_dynamic(this.conds_row)
      % static condition - target is at end pos
      this.targ_init_pos_pix = this.stim_end_pos_pix;
    else
      % dynamic condition - target is initilized at start pos (and will translate)
      this.targ_init_pos_pix = this.stim_start_pos_pix;
    end   
    
    % add stimulus translation, depending on cycle_prop
    this.half_cycle_prop = 0; % initialize where we are in animation cycle
    if conds.is_translating(this.conds_row) == 1
        % target is translating between starting and ending positions
        this.targ_cur_pos_pix = this.targ_init_pos_pix + this.half_cycle_prop*(stim.trans_pos_pix); % account for translation
    else
        % target is not translating and is stationary, keep initial position
        this.targ_cur_pos_pix = this.targ_init_pos_pix;
    end

    if params.debugging
        fprintf('trial %d\n',it);
        fprintf('\tconds_id = %d (%s)\n',this.conds_id,conds.label{this.conds_row});
        fprintf('\tconds_row = %d (%s)\n',this.conds_row,conds.label{this.conds_row});
        fprintf('\tstandard_width_deg\t%.02f\n',this.standard_width_deg);
        fprintf('\ttarg_context_size\t%d\n',this.targ_context_size);
        fprintf('\n\tcur_targ_change_prop_idx\t%d (%.02f)\n',this.cur_targ_change_prop_idx,this.cur_targ_change_prop);
        fprintf('\tstaircase_idx/iter\t%d/%d\n',this.staircase_idx,this.staircase_iter);
    end
    
    
    % ===== fix period =====
    % fixation only
    % fixation point at start pos ("standard" position)
    draw_fix
    timing.fix_start(it,:) = Screen('Flip',w);
    WaitSecs(params.fix_dur);        

    
    % ===== cue period =====
    % fixation and initial stimulus position
    this.is_cue = 1;
    draw_stim
    draw_fix
    timing.cue_start(it,:) = Screen('Flip',w);
    timing.fix_stop(it,:) = timing.cue_start(it,:);

    if params.allow_resp_during_static_cue && ~conds.is_dynamic(this.conds_row)
        % flush KbQueue
        % N.B. this means that responses prior to this point in a trials will not be recorded
        KbQueueFlush(params.kb_dev_id);
        % reset responese-related variables after flushing KbQueue
        this.resp_keys_time = [];
        this.gotresp    = 0;
        this.firstPress = [];
        this.respidx    = 0; % an invalid index in Matlab
        this.resp       = NaN;
        this.resptime   = NaN;
    end
    
    % tick fixation period duration
    WaitSecs(params.cue_dur);    
    this.is_cue = 0;
    
    if conds.is_dynamic(this.conds_row) || (~params.allow_resp_during_static_cue && ~conds.is_dynamic(this.conds_row))
        % flush KbQueue
        % N.B. this means that responses prior to this point in a trials will not be recorded
        KbQueueFlush(params.kb_dev_id);
        % reset responese-related variables after flushing KbQueue
        this.resp_keys_time = [];
        this.gotresp    = 0;
        this.firstPress = [];
        this.respidx    = 0; % an invalid index in Matlab
        this.resp       = NaN;
        this.resptime   = NaN;
    end
    
    
    % ===== animation period =====
    % dynamic or static stimulus
       
    % animation loop
    tmp.first_frame = 1;
    this.animation_start = GetSecs;
    tmp.last_frame = this.animation_start;
    this.curtime = GetSecs - this.animation_start; % current time of cycle
    while ~this.gotresp && this.curtime < params.animation_dur
        % check to see if we've gotten a response yet
        [this.gotresp, this.firstPress] = KbQueueCheck(params.kb_dev_id); % check response queue
        
        % figure out where we are in the current cycle
        %%this.cycle_prop = mod(this.curtime,params.cycle_dur) / params.cycle_dur; % current proportion of the full cycle
        this.full_cycle_prop = mod(this.curtime,params.full_cycle_dur) / params.full_cycle_dur; % current proportion of the full cycle
        this.half_cycle_prop = 2*(0.5-abs(0.5-this.full_cycle_prop)); % current proportion of the half cycle (works for both increasing and decreasing)

        % add stimulus translation, depending on half_cycle_prop
        if conds.is_translating(this.conds_row) == 1
            % target is translating between starting and ending positions
            this.targ_cur_pos_pix = this.targ_init_pos_pix + this.half_cycle_prop*(stim.trans_pos_pix); % account for translation
        else
            % target is not translating and is stationary, keep initial position
            this.targ_cur_pos_pix = this.targ_init_pos_pix;
        end
        
        % draw the stimulus and flip the screen
        draw_stim
        draw_fix
        tmp.this_frame = Screen('Flip',w);
        if tmp.first_frame
            timing.anim_start(it,:) = tmp.this_frame;
            timing.cue_stop(it,:) = tmp.this_frame;
            tmp.first_frame = 0;
        end
        
        if params.debugging>1
            %fprintf('\t%.02f  %.02f  %.02f\n',this.half_cycle_prop,this.anchor_halfwidth_pix,this.targ_size);
            fprintf('\t\t%.04f ms since last frame\n',tmp.this_frame-tmp.last_frame);
            fprintf('\t\t%.04f ms to response check since last frame\n',this.resptime-tmp.last_frame);
            fprintf('\t\t%.04f ms since response check\n',tmp.this_frame-this.resptime);
            tmp.last_frame = tmp.this_frame;
        end
        
        % update curtime for while check or next loop iteration
        this.curtime = GetSecs - this.animation_start; % current time of animation period
    end
        
    
    % ===== resp period =====
    % N.B. we are using KbQueue, so subject can respond at anytime during trial (after cue, which is when KbQueueFlush occurs).
    %      This is just "extra" post-animation time)
    
    % fixation only (in final stimulus position)
    draw_fix
    timing.respwin_start(it,:) = Screen('Flip',w);
    timing.anim_stop(it,:) = timing.respwin_start(it,:);
        
    this.respwin_start = timing.respwin_start(it,:);
    this.curtime = GetSecs - this.respwin_start; % current time of cycle
    while ~abort && ~this.gotresp && this.curtime < params.resp_dur
        % check to see if we've gotten a response yet
        [this.gotresp, this.firstPress] = KbQueueCheck(params.kb_dev_id); % check response queue
        
        % update curtime for while check or next loop iteration
        this.curtime = GetSecs - this.respwin_start; % current time of cycle
    end
    
    
    % ===== iti =====
    % clear screen
    Screen('FillRect',w,params.bg_color); % draw single circular fix spot
    if params.show_fixspot_during_iti
        draw_fix
    end
    timing.iti_start(it,:) = Screen('Flip',w);
    timing.respwin_stop(it,:) = timing.iti_start(it,:);

    % record data
    if this.gotresp
        % was the abort key pushed?  if so, kick out of this if statement before updating response params
        if this.firstPress(params.abort_key)
            abort = 1;
            continue
        end
        
        % update response-related variables based on the response we obtained
        
        % first detected response key - record and validate depending on current task and trial params
        this.resp_keys_time = this.firstPress(params.resp_keys);
        
        % it is possible that more than one response key was pressed during the trial (if more than one was allowed)
        % in this case, we'll only record the first key that was pressed.
        this.respidx = find(this.resp_keys_time);
        
        % N.B. In the insanely-virtually-impossible case that two keys
        %      were pressed at EXACTLY the same time, we'll arbitrarily
        %      select one as the recorded response.  Not sure if this is the
        %      best thing to do, but (1) it should NEVER happen and (2) it
        %      would screw up the data structure if we try to record multiple
        %      responses.
        if length(this.respidx) > 1
            % then multiple keys have been pressed.  arbitrarily select one as the recorded response.
            this.respidx = Sample(this.respidx);
        end
        
        % this.respidx is the index into params.resp_keys
        this.resp = params.resp_keys(this.respidx); % keycode
        this.resptime = this.resp_keys_time(this.respidx); % response time relative to KbQueue reset (start of cue period)
    end
    data.resp(it,1)      = this.resp;
    data.respidx(it,1)   = this.respidx;
    data.resptime(it,1)  = this.resptime;
    if this.cur_targ_change_prop > 0
        % then either target was LARGER than standard (static trial) or
        % the target was GROWING (dynamic trial).
        % THUS: CIRCLE was BIGGER at END position (bottom of screen)
        % correct response is resp_key idx 2
        this.correct_respidx = 2;
    elseif  this.cur_targ_change_prop < 0
        % then either target was SMALLER than standard (static trial) or
        % the target was SKRINKING (dynamic trial).
        % THUS: CIRCLE was BIGGER at START position (top of screen)
        % correct response is resp_key idx 1
        this.correct_respidx = 1;
    else
        KbQueueRelease(params.kb_dev_id); % release keypress queue or process may continue to run in background
        ShowCursor; ListenChar(0); Priority(0); Screen('CloseAll');
        error('we should not allow the target to be EXACTLY the same size as the standard, else there is no ''correct'' response, and no objective way to update the staircase procedure')
    end
    if ~this.gotresp
         % we did not get a valid response
        data.status(it,1) = -1;
        data.valid(it,1)  = 0;
    elseif this.respidx == this.correct_respidx
        % response was 'correct' and valid
        data.status(it,1) = 1;
        data.valid(it,1)  = 1;
    else
        % response was 'incorrect' and valid
        data.status(it,1) = 0;
        data.valid(it,1)  = 1;
    end

    % update staircase structure based on response
    % figure out how large a step we are going to take based on what iteration we are up to for the current staircase
    this.staircase_step = params.staircase.steps_by_iter.steps(find(this.staircase_iter <= params.staircase.steps_by_iter.iters,1,'first'));

    % N.B. Remember participant is ALWAYS responding whether the circle was
    %      larger at the top or bottom (start or end position).  For both
    %      static and dynamic trials, the STANDARD is at the START position
    %      and the TARGET is at the END position.
    switch this.respidx
        case 1
            % participant reported that the circle at the START POSITION was LARGER
            % circle at the START POSITION is the STANDARD
            % thus participant perceived that TARGER was SMALLER
            % update stairase by making NEXT TARGET LARGER (so INCREASE targ_change_prop)
            staircase.cur_targ_change_prop_idx(this.staircase_idx) = min(this.cur_targ_change_prop_idx + this.staircase_step,length(params.staircase.targ_change_props));
        case 2
            % participant reported that the circle at the END POSITION was LARGER
            % circle at the END POSITION is the TARGET
            % thus participant perceived that TARGER was LARGER
            % update stairase by making NEXT TARGET SMALLER (so DECREASE targ_change_prop)
            staircase.cur_targ_change_prop_idx(this.staircase_idx) = max(this.cur_targ_change_prop_idx - this.staircase_step,1);
        otherwise
            KbQueueRelease(params.kb_dev_id); % release keypress queue or process may continue to run in background
            ShowCursor; ListenChar(0); Priority(0); Screen('CloseAll');
            error('it appears that the participant did not make a valid response.  it is unclear what to do for this trial.  For now, we throw an error.  One possibility would be to not update anyting in the staircase structure and redo the trial later, but this requires altering the trial structure and it counter too') 
    end

    % tick remainig ITI duration
    if params.show_feedback_during_iti
        switch data.status(it,1)
            case -1
                this.feedback_text = 'INVALID RESPONSE';
                this.feedback_color = [0 0 255];
            case 0
                this.feedback_text = 'INCORRECT';
                this.feedback_color = [255 0 0];
            case 1
                this.feedback_text = 'CORRECT';
                this.feedback_color = [0 255 0];
            otherwise
                error('unknown trial status, can''t determine feedback')
        end
        AddText(w, this.feedback_text, [xc yc],this.feedback_color, 1, params.instruction_text_size);
        if params.show_fixspot_during_iti
            draw_fix
        end
        Screen('Flip',w);
    end
    while GetSecs - timing.iti_start(it,:) < params.iti
        WaitSecs('YieldSecs',0.001);
    end

    timing.trial_stop(it,1) = Screen('Flip',w); % clear screen
    timing.iti_start(it,1) = GetSecs;
    timing.trial_stop(it,1) = timing.iti_start(it,1);
end


%% notify subject that run is complete
AddText(w, 'Run Complete!', [xc yc], params.text_color, 1, params.instruction_text_size);
AddText(w, 'Thanks for your participation!', [xc yc], params.text_color, -1, params.instruction_text_size);
timing.run_stop = Screen('Flip',w);

%% post-run calculations and cleanup
% calculate timing durations
for scell = {'run' 'trial' 'fix' 'cue' 'anim' 'respwin' 'iti' 'break'}
    s = scell{1}; % convert to string
    timing.([s '_dur']) = timing.([s '_stop']) - timing.([s '_start']);
end

% save data file
save(fullfile(params.datadir,params.datafile));

% Close window and restore normal functioning
WaitSecs(3); % make sure subject sees the final message
KbQueueRelease(params.kb_dev_id); % release keypress queue or process may continue to run in background
ShowCursor;
ListenChar(0);
Priority(0);
Screen('CloseAll')


%% feedback on console
% ...now that Screen is closed (since PTB writes lots of things to console): some info about last file/run for easy tracking
display(sprintf('\n\n\n*****************************************'));
display(sprintf('   last file: %s', params.datafile));
display(sprintf('   last run : %d', params.runid));
%display(sprintf('   performance = %.2f%%',perf.overall));
display(sprintf('*****************************************\n'));