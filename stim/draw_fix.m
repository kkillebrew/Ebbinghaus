% script draw fixation spot for ebbcorr_P2
% (to avoid code repetition)


%******************************************************************
% debugging...(draw lines denoting the trajectory of the dynamic
% stimulus and the line between the trajectory half-way point and
% the peripheral fixation.  angle should be 90 deg)
if params.debugging > 1
    Screen('DrawLine', w,[0 0 255], xc+stim.start_pos_pix(1)+this.global_pos_jitter_pix(1), yc+stim.start_pos_pix(2)+this.global_pos_jitter_pix(2),xc+stim.end_pos_pix(1)+this.global_pos_jitter_pix(1), yc+stim.end_pos_pix(2)+this.global_pos_jitter_pix(2)); % debugging
    Screen('DrawLine', w,[255 0 0], xc+tmp.traj_center_point_pos_pix(1)+this.global_pos_jitter_pix(1), yc+tmp.traj_center_point_pos_pix(2)+this.global_pos_jitter_pix(2), xc+fix.placement.periph_pos_pix(1)+this.global_pos_jitter_pix(1),yc+fix.placement.periph_pos_pix(2)+this.global_pos_jitter_pix(2)); % debugging
    tmp.fix_rect2 = CenterRectOnPoint(fix.rect,tmp.traj_center_point_pos_pix(1),tmp.traj_center_point_pos_pix(2));
    tmp.fix_rect2 = OffsetRect(tmp.fix_rect2,this.global_pos_jitter_pix(1),this.global_pos_jitter_pix(2)); % account for global offset
    Screen('FillOval',w,[0 0 255],tmp.fix_rect2); % draw single circular fix spot
end
%******************************************************************


switch conds.fix_placement_idx(this.conds_row)
    case 0 % = 'none'
        % do nothing (no fixation spot)
    case 1 % 1 = 'target'
        if ~conds.is_dynamic(this.conds_row)
            % STATIC: put a fixation spot on the standard
            tmp.fix_rect = CenterRectOnPoint(fix.rect,this.standard_pos_pix(1),this.standard_pos_pix(2));
            tmp.fix_rect = OffsetRect(tmp.fix_rect,this.global_pos_jitter_pix(1),this.global_pos_jitter_pix(2)); % account for global offset
            Screen('FillOval',w,params.fix.color,tmp.fix_rect); % draw single circular fix spot
        else
            % DYNAMIC: put a fixation spot on the target (will follow target)
            tmp.fix_rect = CenterRectOnPoint(fix.rect,this.targ_cur_pos_pix(1),this.targ_cur_pos_pix(2));
            tmp.fix_rect = OffsetRect(tmp.fix_rect,this.global_pos_jitter_pix(1),this.global_pos_jitter_pix(2)); % account for global offset
            Screen('FillOval',w,params.fix.color,tmp.fix_rect); % draw single circular fix spot
        end            
        
        % the code below will always show a fixation spot on the target,
        %         % and also on the standard for the static condition
        %         % fixation spot centered on target; or in periphery, defined relative to target
        %         tmp.fix_rect = CenterRectOnPoint(fix.rect,this.targ_cur_pos_pix(1),this.targ_cur_pos_pix(2));
        %         tmp.fix_rect = OffsetRect(tmp.fix_rect,this.global_pos_jitter_pix(1),this.global_pos_jitter_pix(2)); % account for global offset
        %         Screen('FillOval',w,params.fix.color,tmp.fix_rect); % draw single circular fix spot
        %
        %         % fixation spot centered on standard (if present)
        %         if ~conds.is_dynamic(this.conds_row)
        %             tmp.fix_rect = CenterRectOnPoint(fix.rect,this.standard_pos_pix(1),this.standard_pos_pix(2));
        %             tmp.fix_rect = OffsetRect(tmp.fix_rect,this.global_pos_jitter_pix(1),this.global_pos_jitter_pix(2)); % account for global offset
        %             Screen('FillOval',w,params.fix.color,tmp.fix_rect); % draw single circular fix spot
        %         end
        
    case 2 % 2 = 'periphery'
        % fixation spot in static periphery, defined relative to the target/standard positions
        % see main code for calculation of fixation spot position.
        tmp.fix_rect = CenterRectOnPoint(fix.rect,xc+fix.placement.periph_pos_pix(1),yc+fix.placement.periph_pos_pix(2));
        tmp.fix_rect = OffsetRect(tmp.fix_rect,this.global_pos_jitter_pix(1),this.global_pos_jitter_pix(2)); % account for global offset
        Screen('FillOval',w,params.fix.color,tmp.fix_rect); % draw single circular fix spot
        

    case {3 4} % ='inducer'
        error('fixation spot placement on inducer is currently disabled.  was used for Frontier''s Dyn Ebb paper. see code to reactivate...')
        
        %         if ~conds.is_dynamic(this.conds_row)
        %             error('fixation spot placement on inducer is currently invalid for static conditions.  see notes at top of code.')
        %         end
        %         % fixation spot centered on an inducer
        %         this.fixinducer_idx = params.fix.placement_args{conds.fix_placement_idx(this.conds_row)};
        %         if this.targ_context_size == 1
        %             % rotate fixinducer_idx by 180 deg to account for
        %             % differenting dynamics for targ-in-small and
        %             % targ-in-large configurations
        %             this.fixinducer_idx = mod(this.fixinducer_idx + params.ebbinghaus.n/2,params.ebbinghaus.n);
        %         end
        %         tmp.fix_rect = AlignRect(fix.rect,this.irects(this.fixinducer_idx,:),'center','center');
        %         Screen('FillOval',w,params.fix.color,tmp.fix_rect); % draw single circular fix spot
    otherwise
        error('fixation spot placement index (%d) undefined.',conds.fix_placement_idx(this.conds_row))
end
