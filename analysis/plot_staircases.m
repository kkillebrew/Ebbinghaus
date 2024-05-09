% quick script for plotting staircase results to visually verify that they
% converge and asymtope.

figure('Position',[1 449 759 335]);

nconds = length(params.included_conds);

nrow = floor(sqrt(nconds)); ncol = ceil(sqrt(nconds)); % number of subplots
if nrow*ncol<nconds
   % not enough cells, add a row (since it was 'floor'ed in the above code)
   nrow = nrow+1;
end
 
% resorted data
change_by_cond_sc = NaN(trials_per_stair,nstairs*niter,nconds);
resp_by_cond_sc   = NaN(trials_per_stair,nstairs*niter,nconds);


ignore_last_n  = 0; % same as above, but ignoring trials at end of staircase.  good to test whether PSEs change appreciably if the staircases were shorter (which would be more efficient)
ignore_first_n = 15; % how many trials for each staircase should we ignore when calculating average targ_change_prop for PSE?
avg_pse = NaN(1,nconds); % point of subjective equality: average targ_change_prop rate after the first X trials

rev_n2avg = [5 4]; % average over how many reversals [default min_allowed], counting back from ignore_last_n, when calculating average reversal rate for PSE
rev_pse = NaN(1,nconds); % point of subjective equality: average of reversals after the first X trials


k = 0;
for cond = condition_list
    %fprintf('cond %.01f\n',cond)
    k = k+1;
    this_staircase_ref = staircase.conds_id==cond;
    this_idxseq = staircase.targ_change_prop_idx_seq(this_staircase_ref,:); % pull out all staircases for this condition (index values)
    this_targ_change_prop_seq = params.staircase.targ_change_props(this_idxseq); % actual targ_change_prop values
     
    
    
    % UP TO HERE
    
    
    
     
    j = 0;
    for sc = stairtype_list % start small, start big
        for i = iter_list
            j = j+1;
            ref = and(this_targ_change_prop_seq(:,2)==sc,this_targ_change_prop_seq(:,3)==i);

            % convert targ_change_prop rates from pix/frame to %total change of original bar size
            % ((growth_list*duration)/ppd) / (rmin/ppd)
            % 100 * ((pix/frame*frames) / pix/frame) / (pix/ pix/frame)
            this_growth_pix_per_frame = this_targ_change_prop_seq(ref,4); % pix/frame
            change_by_cond_sc(:,j,k) = 100*((this_growth_pix_per_frame*duration)/ppd) / (rmin/ppd); % %deg of starting bar size (not accounting for bar size jitter...which was not a stored variable)
            
            resp_by_cond_sc(:,j,k)   = this_targ_change_prop_seq(ref,5);
        end
    end
    

    
    subplot(nrow,ncol,k)
    plot(change_by_cond_sc(:,:,k));
    hline(0,'k--');
    vline(ignore_first_n,'k:')
    vline(trials_per_stair-ignore_last_n,'k:')
    title(sprintf('condition %.01f',cond));
    xlabel('trial');
    ylabel('growth rate (% rect)');
end

% calculate PSEs based on reversals
d = diff(resp_by_cond_sc);
d2 = cat(1,zeros(1,nstairs*niter,nconds),d);
rev = ne(d2,0); % reversals occur when the response is different than the previous response for a given staircase

rev_kept = rev(1:end-ignore_last_n,:,:); % ignore some trials at the end if requested
lastn = bsxfun(@gt,cumsum(rev_kept),sum(rev_kept)-rev_n2avg(1)); % reference list for trials we want to average over
rev_kept_lastn = and(rev_kept,lastn); % merge
rev_lastn = cat(1,rev_kept_lastn,false(ignore_last_n,nstairs*niter,nconds)); % pad to be same size as change_by_cond_sc
if any(any(sum(rev_lastn)<rev_n2avg(2)))
    error('at least one condition had too few reversals (less than %d) when ignoring the last %d trials',rev_n2avg(2),ignore_last_n);
end
change_by_cond_sc_rev = change_by_cond_sc;
change_by_cond_sc_rev(~rev_lastn) = NaN;
rev_pse = nanmean(change_by_cond_sc_rev,1); % average over last n reversals
rev_pse = squeeze(mean(rev_pse,2))'; % average over staircases for each condition


%%% the following code is simply averaging all reversals within a given set of trials.  not averaging a defined number of reversals
%%%change_by_cond_sc_rev = change_by_cond_sc;
%%%change_by_cond_sc_rev(~rev) = NaN;
%%%rev_pse = nanmean(change_by_cond_sc_rev(ignore_first_n+1:end-ignore_last_n,:,:),1); % average over trials (ignore some at begining of each staircase)
%%%rev_pse = squeeze(mean(rev_pse,2))'; % average over staircases for each condition

% calculate PSEs based on all staircase steps
avg_pse = nanmean(change_by_cond_sc(ignore_first_n+1:end-ignore_last_n,:,:),1); % average over trials (ignore some at begining of each staircase)
avg_pse = squeeze(mean(avg_pse,2))'; % average over staircases for each condition


% all conditions together
figure('Position',[2 26 560 420]);
mean_growths = squeeze(mean(change_by_cond_sc,2));
% % if k == nrow*ncol
% %     % then there won't be anymore panels to work with
% %     figure
% % else
% %     subplot(nrow,ncol,k+1)
% % end
subplot(2,1,1)
plot(mean_growths);
hline(0,'k--')
vline(ignore_first_n,'k:')
vline(trials_per_stair-ignore_last_n,'k:')
title('mean growth by condition');
xlabel('trial');
ylabel('growth rate (% rect)');
legend(num2str(condition_list'),'Location','Best');

% PSE
subplot(2,1,2)
bar([avg_pse;rev_pse]')
set(gca,'XTick',1:nconds,'XTickLabel',condition_list);
legend({'all steps' 'reversals'},'Location','Best')
title('PSE extracted from staircase')
xlabel('condition')
ylabel('PSE (% rect)')
