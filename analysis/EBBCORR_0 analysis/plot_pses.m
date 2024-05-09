% quick analysis script for Dynamic Ebbinghaus / Corridor comparison.
% 02.28.15: modified from plot_pses5.m from ebbinghaus_0


% Notes - this needs to be adjusted to account for the change from +/-
% inducers to 3 context, including "no context"



% *********************************************************************
% TO USE WITH GROUP_ANALYSIS SCRIPT, DATA SHOULD BE LOAD BEFORE CALLING
% THIS SCRIPT


%% options
outlier_removal_sd_thresh = 'minmax'; % remove outliers that are mean +/- outlier_removal_sd_thresh*stdev, on a condition-by-condition basis.  set to Inf for no outlier removal.  set to 'minmax' to remove the min and max value for each condition


%% some info on timing
fprintf('== run duration (not including breaks) = %.02f min\n',(timing.run_dur-nansum(timing.break_dur))/60);


%% optional filter to only one targ-context size or one standard size
% params.targ_context_sizes = 2;
% params.standard_width_deg = 2;
% ref = and(ismember(trials.targ_context_size,params.targ_context_sizes),ismember(trials.standard_width_deg,params.standard_width_deg));
% for fncell = {'conds_id' 'conds_row' 'targ_context_size' 'standard_width_deg'}
%     fn = fncell{1}; % convert to string
%     trials.(fn) = trials.(fn)(ref);
% end
% data.PSE_width_prop_standard = data.PSE_width_prop_standard(ref);



%% optional filter to limit total number of trials (cut in half)
% % HACK - lots of inappropriate assumptions
% ref = trials.order <= trials.n/2; % cut data in half
% for fncell = {'conds_id' 'conds_row' 'targ_context_size' 'standard_width_deg'}
%     fn = fncell{1}; % convert to string
%     trials.(fn) = trials.(fn)(ref);
%     trials.(fn) = repmat(trials.(fn),2,1);
% end
% data.PSE_width_prop_standard = data.PSE_width_prop_standard(ref);
% data.PSE_width_prop_standard = repmat(data.PSE_width_prop_standard,2,1);


%% analysis collapsed across standard size and targ-context size
% start with separate data for targ-in-small and targ-in-large
% mean PSE by condition
[n,s] = SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size],'length');
m = SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size],'nanmean');
v = SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size],'nanste');

% remove outliers
vec =  SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size],'vector');
outlier_ref = zeros(size(vec)); % assume no outliers to initialize list
if strcmp(outlier_removal_sd_thresh,'minmax')
    vec2 = vec + .0001*rand(size(vec)); % adding a tiny amount of noise to individual trial values will help avoid selecting more than one trial as the "max" or "min" for that condition
    outlier_ref = or(outlier_ref,bsxfun(@eq,vec2,max(vec2,[],2))); % mark highest value for each condition an outlier
    outlier_ref = or(outlier_ref,bsxfun(@eq,vec2,min(vec2,[],2))); % mark lowest value for each condition an outlier
    % check that this worked and ONLY identified two outliers
    if max(sum(outlier_ref,2)) ~= 2 || min(sum(outlier_ref,2)) ~= 2
        error('outlier detection failed.  see outlier_ref.')
    end
    %sum(outlier_ref,2); % debugging (should be a column vector of 2's)
elseif outlier_removal_sd_thresh < Inf
    sd =  SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size],'std');
    cond_sd_thresh = outlier_removal_sd_thresh * sd;
    sd_from_mean_sorted = sort(bsxfun(@rdivide,bsxfun(@minus,vec,m),sd),2); % how many SD is each value away from the mean, sorted by condiion
    outlier_ref = or(outlier_ref,bsxfun(@gt,vec,m+cond_sd_thresh)); % mark high outliers
    outlier_ref = or(outlier_ref,bsxfun(@lt,vec,m-cond_sd_thresh)); % mark low outliers
end
if sum(outlier_ref(:))
    fprintf('== removing %d outliers (for analysis collapsed across standard size and targ-context size)\n',sum(outlier_ref(:)))
    vec(outlier_ref) = NaN; % replace outliers with NaNs
    m = nanmean(vec,2); % recalculate mean ignoring outliers
    v = nanste(vec,2); % recalculate mean ignoring outliers
end

% convert to percentages, rather than proportions
m = 100*m;
v = 100*v;


% extract with- and without- context conditions, and collapse data for
% without-context conditions
with_context_conds = m(s(:,1) >= 20);
without_context_conds = m(s(:,1) < 20);
without_context_conds = repmat(without_context_conds,2,1); % double since the without_context_conds are used as the baseline for both the Ebbinghaus and Corridoe conditions
% *** we aren't justified in collapsing across the -context targ-in-small
%     and targ-in-large conditions for all manipulations.
% *** however, we ARE justified in collapsing if the experiments only
%     includes one of the targ-in-large or targ-in-small conditions. 
if length(params.targ_context_sizes) > 1
    % collapse the -context condition across the targ-in-small and
    % targ-in-large conditions, since they are not actually different
    without_context_conds = repmat(mean(reshape(without_context_conds,[],2),2),2,1);
end

% subtract -context conditions from +context conditions to normalize for
% response bias
m_norm = with_context_conds - without_context_conds;
s_norm = s(s(:,1) >= 20,:);


dir_labels = {'targ-in-small' 'targ-in-large'};
switch params.experiment
    case {'ebbcorrP1'}
        idx = Remap(s_norm(s_norm(:,2)==min(params.targ_context_sizes),1),conds.id,conds.row); % only the conditions that were actually used in this experiment.  assumes we ALWAYS run the targ-in-large condition
        s_label1 = conds.label(idx);
        s_label1_short = regexprep(s_label1,{'static' 'dynamic' 'normtarg' 'lowcontrast' 'blurred' 'isoluminant' 'jittered' 'fixtarg' 'fixperiph' 'fixcontext' 'tracktarg' 'trackcontext' 'stationary' 'translating' 'nocontext' 'ebbinghaus' 'corridor'},{'S' 'D' 'NT' 'LC' 'B' 'IL' 'J' 'FT' 'FP' 'FI' 'TT' 'TI' 'S' 'T' 'N' 'E' 'C'});
        s_label2 = dir_labels(unique(s(:,2)));
        idx = Remap(s(:,1),conds.id,conds.row);
        allconds = conds.label(idx); % only the conditions that were actually used in this experiment
        allconds_shortlabel = regexprep(allconds,{'static' 'dynamic' 'normtarg' 'lowcontrast' 'blurred' 'isoluminant' 'jittered' 'fixtarg' 'fixperiph' 'fixcontext' 'tracktarg' 'trackcontext' 'stationary' 'translating' 'nocontext' 'ebbinghaus' 'corridor'},{'S' 'D' 'NT' 'LC' 'B' 'IL' 'J' 'FT' 'FP' 'FI' 'TT' 'TI' 'S' 'T' 'N' 'E' 'C'});
        
    otherwise
        error('unrecognized experiment (%s)',params.experiment)
end



% match PSE direction magnitudes for targ-in-small (invert so positive) and targ-in-large (naturally positive)
m_norm_pos_bydir = m_norm .* sign(s_norm(:,2)-1.5);
with_context_conds_matchedPSE = (with_context_conds - 100) .* sign(s_norm(:,2)-1.5) + 100; 

% average over targ-context sizes (after matching PSE direction), if multiple targ-in-X were used
m_norm_pos = mean(reshape(m_norm_pos_bydir,[],length(params.targ_context_sizes)),2);
with_context_conds_collapsed = mean(reshape(with_context_conds_matchedPSE,[],length(params.targ_context_sizes)),2);
without_context_conds_collapsed = mean(reshape(without_context_conds,[],length(params.targ_context_sizes)),2);

% package data for group analysis
subj.adj_pse.lab = s_label1_short';
subj.adj_pse.norm = m_norm_pos';
%subj.adj_pse.norm_bydir  = m_norm_pos_bydir';
subj.adj_pse.context = with_context_conds_collapsed' - 100;  % rescale to same as 'norm_'
%subj.adj_pse.context_bydir = with_context_conds' - 100; % rescale to same as 'norm_'
subj.adj_pse.nocontext = without_context_conds_collapsed' - 100;  % rescale to same as 'norm_'
%subj.adj_pse.nocontext_bydir = without_context_conds' - 100; % rescale to same as 'norm_'



% plot it
figure

% N.B. no variance for a difference metric with a single subject
subplot(2,1,1)
padded_m_norm_pos = cat(1,zeros(3,1),m_norm_pos); % add in filler for the no context conditions
barerr2(1:size(padded_m_norm_pos,1),padded_m_norm_pos,[],'BaseValue',0,'FaceColor','k');
% hold on
% barerr2(1:size(m_norm_pos,1),reshape(m_norm_pos_bydir,[],2),[],'BaseValue',0);
xlim([0 length(padded_m_norm_pos)+1])
set(gca,'XTickLabel',allconds_shortlabel)
ylabel('Illusion Magnitude ({%})')
if length(params.targ_context_sizes) > 1
    title('averaged over standards and equated targ-context sizes')
else
    title('averaged over standards')
end

% show individual trial data (to see spread), by condition (sorted by targ-context-size)
[td td_s] = SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size],'vector'); % trial data sorted by condition
td = 100*td; % convert to percentage, rather than proportion
nconds_per_dir = length(unique(td_s(:,1))); % number of condition per direction (targ-in-small vs. targ-in-large)
subplot(2,1,2)
hold on

td_markers = {'vk' '^k'};
for this_targ_dir = [1 2]; % targ_context_size=1, so targ-in-small, like reverse DISC; targ_context_size=2, so targ-in-large, like DISC
    td_ref = td_s(:,2)==this_targ_dir; % targ-in-small or targ-in-large trials
    if sum(td_ref)
        jittered_x = repmat((1:length(td_s(td_ref,1)))',1,size(td,2)) + (.25*(rand(nconds_per_dir,size(td,2))-.5)); % targ-in-X
        cond_ref = td_s(:,2)==this_targ_dir;
        this_fulldata = td(cond_ref,:); % raw data (conditions are rows)
        this_outliers = outlier_ref(cond_ref,:); % outliers (as defined above)
%         this_mean = mean(this_fulldata,2); % condition means
%         this_std  = std(this_fulldata,[],2);  % condition standard deviations
%         this_outliers = bsxfun(@gt,this_fulldata,this_mean+outlier_removal_sd_thresh*this_std) | bsxfun(@lt,this_fulldata,this_mean-outlier_removal_sd_thresh*this_std);
        plot(jittered_x',this_fulldata',td_markers{this_targ_dir});
        % mark outliers
        if sum(this_outliers(:))
            outlier_jittered_x = jittered_x(this_outliers);
            outlier_this_fulldata = this_fulldata(this_outliers);
            plot(outlier_jittered_x,outlier_this_fulldata,[td_markers{this_targ_dir}(1) 'r'],'MarkerFaceColor','r');
        end
    end
end

xlim([0 nconds_per_dir+1])
set(gca,'XTick',1:nconds_per_dir,'XTickLabel',allconds_shortlabel)
%rotateticklabel(gca,90);
vline(.5:nconds_per_dir+1,'k:');
vline(.5:3:nconds_per_dir+1,'k');
hline(100,':k')
ylabel('PSE ({%})')
title('single-trial data, sorted by condition (targ-context size plotted seperately)')



% % %% analysis collapsed across standard size and targ-context size
% % % start with separate data for targ-in-small and targ-in-large
% % % mean PSE by condition
% % 
% % % only do this if we used categorical levels for standard size (Exp1)
% % if length(unique(trials.standard_width_deg)) == 2
% %     
% %     [n,s] = SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size trials.standard_width_deg],'length');
% %     m = SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size trials.standard_width_deg],'mean');
% %     v = SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size trials.standard_width_deg],'ste');
% %     
% %     % remove outliers
% %     vec =  SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size trials.standard_width_deg],'vector');
% %     outlier_ref = zeros(size(vec)); % assume no outliers to initialize list
% %     if strcmp(outlier_removal_sd_thresh,'minmax')
% %         vec2 = vec + .0001*rand(size(vec)); % adding a tiny amount of noise to individual trial values will help avoid selecting more than one trial as the "max" or "min" for that condition
% %         outlier_ref = or(outlier_ref,bsxfun(@eq,vec2,max(vec2,[],2))); % mark highest value for each condition an outlier
% %         outlier_ref = or(outlier_ref,bsxfun(@eq,vec2,min(vec2,[],2))); % mark lowest value for each condition an outlier
% %         % check that this worked and ONLY identified two outliers
% %         if max(sum(outlier_ref,2)) ~= 2 || min(sum(outlier_ref,2)) ~= 2
% %             error('outlier detection failed.  see outlier_ref.')
% %         end
% %     elseif outlier_removal_sd_thresh < Inf
% %         sd =  SortedFunc(data.PSE_width_prop_standard,[trials.conds_id trials.targ_context_size trials.standard_width_deg],'std');
% %         cond_sd_thresh = outlier_removal_sd_thresh * sd;
% %         sd_from_mean_sorted = sort(bsxfun(@rdivide,bsxfun(@minus,vec,m),sd),2); % how many SD is each value away from the mean, sorted by condiion
% %         outlier_ref = or(outlier_ref,bsxfun(@gt,vec,m+cond_sd_thresh)); % mark high outliers
% %         outlier_ref = or(outlier_ref,bsxfun(@lt,vec,m-cond_sd_thresh)); % mark low outliers
% %     end
% %     if sum(outlier_ref(:))
% %         %    fprintf('== removing %d outliers (for analysis collapsed across standard size and targ-context size)\n',sum(outlier_ref(:)))
% %         vec(outlier_ref) = NaN; % replace outliers with NaNs
% %         m = nanmean(vec,2); % recalculate mean ignoring outliers
% %         v = nanste(vec,2); % recalculate mean ignoring outliers
% %     end
% %     
% %     % convert to percentages, rather than proportions
% %     m = 100*m;
% %     v = 100*v;
% %     
% %     % extract with- (odd conditions) and without- (even conditions) context
% %     % conditions, and collapse data for without-context conditions
% %     with_context_conds = m(1:2:end);
% %     without_context_conds = m(2:2:end);
% %     
% %     % subtract -context conditions from +context conditions to normalize for
% %     % response bias
% %     m_norm = with_context_conds - without_context_conds;
% %     s_norm = s(1:2:end,:);
% %     
% %     m_norm_pos_bydirANDsize = m_norm .* sign(s_norm(:,2)-1.5);
% %     
% %     subj.adj_pse.norm_bydirANDsize = m_norm_pos_bydirANDsize';
% %     subj.adj_pse.context_bydirANDsize = with_context_conds' - 100; % rescale to same as 'norm_'
% %     subj.adj_pse.nocontext_bydirANDsize = without_context_conds' - 100; % rescale to same as 'norm_'
% %     
% % end



MarkPlot(params.datafile);

return
