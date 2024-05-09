% script draw stimuli for ebbcorr_P2
% (to avoid code repetition)

% assumes the following variables are set and up-to-date
%this.half_cycle_prop
%this.targ_cur_pos_pix


% draw target
if conds.is_dynamic(this.conds_row)
    % Dynamic - account for the growth of the target by setting the anchor (i.e.,
    % size of target at the mid-cycle point)
    this.targ_size = this.standard_halfwidth_pix + this.half_cycle_prop*(this.targ_halfwidth_pix-this.standard_halfwidth_pix);
else
    % Static - set the absolute size of the target
    this.targ_size = this.targ_halfwidth_pix;
end
this.trect = [xc-this.targ_size yc-this.targ_size xc+this.targ_size yc+this.targ_size];
this.trect = CenterRectOnPoint(this.trect,this.targ_cur_pos_pix(1),this.targ_cur_pos_pix(2));
this.trect_nojitter = this.trect; % store target position without jitter for fixtarg fixation spot placement
if conds.is_jittered(this.conds_row)
    % add positional jitter
    this.trect = OffsetRect(this.trect,ppd*params.uncertainty.pos_jitter*(2*rand-1),ppd*params.uncertainty.pos_jitter*(2*rand-1)); % cur pos + jitter
end
this.tRGB = params.standard.color;
% N.B. we won't draw the target to the screen until after we deal
%      with the context since the corridor image would cover it


% draw inducers for target
% N.B. we always do the setup for inducers, but only draw them for
%      some conditions, because we may need to know where the
%      inducers would be to draw the fixspot.
if conds.is_dynamic(this.conds_row)
    % Dynamic - inducers grow/shrink (and possibly change eccen) over time
    switch this.targ_context_size
      case 1 % targ-in-small (standard-in-large) [dynamic = large-to-small growth, since adjusted targ is anchor point]
        this.inducer_size  = stim.inducers.large_halfwidth_pix - this.half_cycle_prop*(stim.inducers.large_halfwidth_pix-stim.inducers.small_halfwidth_pix);
        this.inducer_eccen = stim.inducers.large_eccen_pix - this.half_cycle_prop*(stim.inducers.large_eccen_pix-stim.inducers.small_eccen_pix);
      case 2 % targ-in-large (standard-in-small) [dynamic = small-to-large growth, since adjusted targ is anchor point]
        this.inducer_size  = stim.inducers.small_halfwidth_pix + this.half_cycle_prop*(stim.inducers.large_halfwidth_pix-stim.inducers.small_halfwidth_pix);
        this.inducer_eccen = stim.inducers.small_eccen_pix + this.half_cycle_prop*(stim.inducers.large_eccen_pix-stim.inducers.small_eccen_pix);
      otherwise
        error('unrecognized targ_context_size (%d)',this.targ_context_size)
    end
else
    % Static - inducers are static, either large or small
    switch this.targ_context_size
      case 1 % targ-in-small (standard-in-large)
        this.inducer_size  = stim.inducers.small_halfwidth_pix;
        this.inducer_eccen = stim.inducers.small_eccen_pix;
      case 2 % targ-in-large (standard-in-small)
        this.inducer_size  = stim.inducers.large_halfwidth_pix;
        this.inducer_eccen = stim.inducers.large_eccen_pix;
      otherwise
        error('unrecognized targ_context_size (%d)',this.targ_context_size)
    end
end
this.irect_template = [xc-this.inducer_size yc-this.inducer_size xc+this.inducer_size yc+this.inducer_size];
[this.inducer_x,this.inducer_y] = pol2cart(stim.inducers.angs_rad,this.inducer_eccen);
this.irects = OffsetRect(this.irect_template,this.inducer_x',this.inducer_y');
this.irects = OffsetRect(this.irects,this.targ_cur_pos_pix(1)-xc,this.targ_cur_pos_pix(2)-yc); % start pos + global jitter + trans
switch conds.what_context(this.conds_row)
  case 0 % no context
    % do nothing
  case 1 % ebbinghaus
      % actually draw the inducers
    Screen('FillOval',w,params.ebbinghaus.color,this.irects');
  case 2 % corridor
    % draw the corridor background
    this.crect = [0 0 params.x_resolution params.y_resolution]; % init to centered on screen
    %this.crect = [xc-stim.corridor.halfwidth_pix yc-stim.corridor.halfheight_pix xc+stim.corridor.halfwidth_pix yc+stim.corridor.halfheight_pix]; % init to centered on screen (account for requested image size)
    %%this.crect = OffsetRect(this.crect,this.global_pos_jitter_pix(1),this.global_pos_jitter_pix(2)); % account for global offset   
    Screen('DrawTexture',w,stim.corridor.texid,[],this.crect);
  otherwise
    error('unrecognized context %d',conds.what_context(this.conds_row))
end

% draw target (prepared above, drawn to buffer now)
if conds.is_dynamic(this.conds_row) || (~conds.is_dynamic(this.conds_row) && ~this.is_cue) || params.show_targ_during_static_cue
    Screen('FillOval',w,this.tRGB,this.trect);
end

% draw standard (if present)
if ~conds.is_dynamic(this.conds_row)
    % standard is always at the start pos
    this.srect = [xc-this.standard_halfwidth_pix yc-this.standard_halfwidth_pix xc+this.standard_halfwidth_pix yc+this.standard_halfwidth_pix];
    this.srect = CenterRectOnPoint(this.srect,this.stim_start_pos_pix(1),this.stim_start_pos_pix(2)); % start pos + global jitter
    Screen('FillOval',w,params.standard.color,this.srect);    
    
    % add inducers (either large or small, opposite of those surrounding the target)
    switch this.targ_context_size
      case 1 % targ-in-small (standard-in-large)
        this.standard_inducer_size  = stim.inducers.large_halfwidth_pix;
        this.standard_inducer_eccen = stim.inducers.large_eccen_pix;
      case 2 % targ-in-large (standard-in-small)
        this.standard_inducer_size  = stim.inducers.small_halfwidth_pix;
        this.standard_inducer_eccen = stim.inducers.small_eccen_pix;
      otherwise
        error('unrecognized targ_context_size (%d)',this.targ_context_size)
    end
    this.standard_irect_template = [xc-this.standard_inducer_size yc-this.standard_inducer_size xc+this.standard_inducer_size yc+this.standard_inducer_size];
    [this.standard_inducer_x,this.standard_inducer_y] = pol2cart(stim.inducers.angs_rad,this.standard_inducer_eccen);
    this.standard_irects = OffsetRect(this.standard_irect_template,this.standard_inducer_x',this.standard_inducer_y');
    this.standard_irects = OffsetRect(this.standard_irects,this.stim_start_pos_pix(1)-xc,this.stim_start_pos_pix(2)-yc); % start pos + global jitter

    switch conds.what_context(this.conds_row)
      case 0 % no context
        % do nothing
      case 1 % ebbinghaus
        % actually draw the inducers
        Screen('FillOval',w,params.ebbinghaus.color,this.standard_irects');
      case 2 % corridor
        % do nothing, no need to draw a second time
      otherwise
        error('unrecognized context %d',conds.what_context(this.conds_row))
    end
end
