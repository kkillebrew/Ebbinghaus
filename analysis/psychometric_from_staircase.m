% quick script for plotting a psychometric curve derived from the trials
% of an adaptive staircase procedure.

figure('Position',[715 449 726 335]);

nconds = length(params.included_conds);

nrow = floor(sqrt(nconds)); ncol = ceil(sqrt(nconds)); % number of subplots
if nrow*ncol<nconds
   % not enough cells, add a row (since it was 'floor'ed in the above code)
   nrow = nrow+1;
end

% resorted data
targ_change_props = NaN(params.staircase.n_targ_change_props,nconds); % targ_change_prop x condition
numtrials   = NaN(params.staircase.n_targ_change_props,nconds); % number of trials at each targ_change_prop level
resp_mean   = NaN(params.staircase.n_targ_change_props,nconds); % mean response proportion x condition
resp_sem    = NaN(params.staircase.n_targ_change_props,nconds); % standard error of response proportion x condition
resp_bottom = NaN(params.staircase.n_targ_change_props,nconds); % number of times response was 'bottom is bigger' (i.e., growing for dynamic trials) x condition.  N.B. bottom is actually bigger when targ_change_prop > 0

psy_mag     = 5; % factor to upsample the resolution when drawing fitted psychometric curves (relative to the min staircase step)
psy_targ_change_prop = min(params.staircase.targ_change_props):(min(diff(params.staircase.targ_change_props))/psy_mag):max(params.staircase.targ_change_props); % upsample targ_change_props for psychometric curve fitting
psy_curve   = NaN(length(psy_targ_change_prop),nconds); % psychometric curves estimated at resolution and sampling of targ_change_props
psy_beta    = NaN(2,nconds); % betas for psychometric curve fitting
psy_pse     = NaN(1,nconds); % point of subjetive equality, as extracted from psychometric curve fit crossing 0.5

trim_targ_change_ends = 0; % how many targ_change_props on each side should we ignore when calculating psychometric curves

k = 0;
cond_labels = {};
for cond = params.included_conds
    %fprintf('cond %.01f\n',cond)
    k = k+1;
    cond_labels(k) = conds.short_label(conds.id==cond);
    this_ref = trials.conds_id==cond; % pull out data for this condition only
    
    this_targ_change_props = trials.targ_change_prop(this_ref);
    this_resp = data.respidx(this_ref) - 1; % convert from [1 2] (bigger on [top bottom]) to [0 1] so it represents a proportion of 'bigger on bottom' or 'target larger than standard' responses
    
    % not sure if thie is necessary...needed it during incomplete pilot testing
    valid_ref = ~isnan(this_resp);
    %%%valid_ref = and(valid_ref,(1:length(this_resp))'<=75); % testing how reliable the psychometric curves are with fewer trials (first X trials, set by the selection list provided here)
    this_resp = this_resp(valid_ref);
    this_targ_change_props = this_targ_change_props(valid_ref);
    
    this_used_targ_change_props = unique(this_targ_change_props); % figure out which targ_change_props were actually used in the experiment (because the staircase is allowed to skip some in the early trials)
    used_ref = ismember(params.staircase.targ_change_props',this_used_targ_change_props); % reference list
    
    [resp_mean(used_ref,k) targ_change_props(used_ref,k)] = SortedFunc(this_resp,this_targ_change_props,'mean'); % mean of response for each unique growth value 
    [resp_sem(used_ref,k) tmp] = SortedFunc(this_resp,this_targ_change_props,'ste'); % standard error of response for each unique growth value 
    [numtrials(used_ref,k) tmp] = SortedFunc(this_resp,this_targ_change_props,'length'); % number of trials for each unique growth value 
    [resp_bottom(used_ref,k) tmp] = SortedFunc(this_resp,this_targ_change_props,'sum'); % recall that response has been converted to 0 (shrink) or 1 (grow) at this point.  so sum represents total trials in which subject repsonded 'growing'

    % trim targ_change_prop ends before calculating psychometric curves, by ignoring most extreme values
    used_ref(1:trim_targ_change_ends) = false;
    used_ref(end-trim_targ_change_ends+1:end) = false;

    % plot actual data for this condition
    subplot(nrow,ncol,k)
    errorbar(targ_change_props(used_ref,k),resp_mean(used_ref,k),resp_sem(used_ref,k),'bo'); % used for psychometric
    if sum(~used_ref)
        hold on
        errorbar(targ_change_props(~used_ref,k),resp_mean(~used_ref,k),resp_sem(~used_ref,k),'ko'); % not used for psychometric
    end
    vline(0,'k--');
    hline(0.5,'k--');
    title(sprintf('condition %.01f',cond));
    xlabel('target change proportion (relative to standard)');
    ylabel('proportion reported target is bigger');
    
    % calculate psychometric function
    % adapted from from Gideon's old code (by way of Chris):
    %     PSE(i,j) = -b(i,1,j)/b(i,2,j);
    %     fitdata(j,:,i) = 100* exp(b(i,1,j)+b(i,2,j)*xx_axis')./(1+exp(b(i,1,j)+b(i,2,j)*xx_axis'));
    %     plot(xaxis,100*results(j,:,i)',[color{j} 'o'],'LineWidth',2);
    % N.B. for some cases, the PSE dervied from the beta weights can be slightly biased.
    %      consider the case when responses are perfectly veridical:
    %     n = [3 3 3 3 3 3 3 3 51 51 3 3 3 3 3 3 3 3]';
    %     r = [0 0 0 0 0 0 0 0 0 51 3 3 3 3 3 3 3 3]';
    %     g = linspace(-75,75,18);
    %     b_resp_bottom = glmfit(g,[r n],'binomial','logit');
    %     pse_resp_bottom = -b_resp_bottom(1)/b_resp_bottom(2)
    %     b_resp_shrink = glmfit(g,[n-r n],'binomial','logit');
    %     pse_resp_shrink = -b_resp_shrink(1)/b_resp_shrink(2)
    % pse_resp_bottom and pse_resp_shrink have opposite signs.  In this case, the true PSE should be zero.  this may be related to the "iteration limit reached' warning from glmfit
    b_resp_bottom = glmfit(targ_change_props(used_ref,k),[resp_bottom(used_ref,k) numtrials(used_ref,k)],'binomial', 'logit'); % betas derived from a 'proportion perceived growing' curve
    pse_resp_bottom = -b_resp_bottom(1)/b_resp_bottom(2);
    b_resp_shrink = glmfit(targ_change_props(used_ref,k),[numtrials(used_ref,k)-resp_bottom(used_ref,k) numtrials(used_ref,k)],'binomial', 'logit'); % betas derived from a 'proportion perceived shrinking' curve
    pse_resp_shrink = -b_resp_shrink(1)/b_resp_shrink(2);
    if sign(pse_resp_bottom) ~= sign(pse_resp_shrink)
        fprintf('\n***WARNING: PSE derived from ''prop grow'' and ''prop shrink'' have opposite signs.  Seeting PSE to zero to avoid problems with signrank statistics.\n\n')
        psy_beta(:,k) = [0; b_resp_bottom(2)];
    else 
        psy_beta(:,k) = b_resp_bottom;
    end
    psy_pse(1,k) = -psy_beta(1,k)/psy_beta(2,k);
    
    % the following is the same as:
    % resp_bottom_fit = glmval(b, targ_change_props(:,k), 'logit', 'size', numtrials(:,k));  resp_bottom_fit ./ numtrials(:,k);
    % resp_bottom_fit2 = exp(b(1,k)+b(2,k)*targ_change_props(:,k))./(1+exp(b(1,k)+b(2,k)*targ_change_props(:,k))); % this is taken directly from Gideon's implementation
    psy_curve(:,k) = glmval(psy_beta(:,k),psy_targ_change_prop','logit');
    hold on
    plot(psy_targ_change_prop,psy_curve(:,k), 'r-')
end


% all conditions together
figure('Position',[881 26 560 420]);

% psychometric curves
subplot(2,1,1)
plot(psy_targ_change_prop,psy_curve,'-');
vline(0,'k--');
hline(0.5,'k--');
title('psychometric curves by condition');
xlabel('target change proportion (relative to standard)');
ylabel('proportion reported target is bigger');
%legend(num2str(params.included_conds'),'Location','NorthWest');
legend(cond_labels,'Location','NorthWest');

% PSEs
subplot(2,1,2)
bar(psy_pse)
set(gca,'XTick',1:nconds,'XTickLabel',cond_labels);
title('PSE extracted from psychometric curves')
xlabel('condition')
ylabel('PSE (% rect)')


