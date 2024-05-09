% script to do the paired comparison across same subject datasets for revampHORRM and LC expeirments

% can't clear all if we want to be able to call from group_analysis...
%clear all
%close all


% see group_analysis for more details on the datasets
%
subjlist = {
    'DBD'
    'LXL'
    'TMD' % removed due to high likelihood of response bias
    'UMD'
    'KEN'
    'BKN'
    'REM' % added for VSS presentation
    'DEC' % added for VSS presentation
    };



% save command window output to textfile for quick access to stats
diary_file = './outputs/lc_output_VSS.txt';
if exist(diary_file,'file')
    delete(diary_file); % always start with a blank diary
end
diary(diary_file); 

which_exp = 'lc'

nsubj = length(subjlist)

norm_subj_pses        = 0; % should we normalize subject PSEs to a percentage of the maximum PSE across conditions?
keep_individual_plots = 0; % boolean, should we keep individual figures on screen?
do_stats              = 1; % run various stats on the data and print to screen?
test_for_normality    = 0; % do tests for normality of data in each condition?
do_perm_test          = 1; % run the permutation test?
%extract_params        = 0; % should we extract params?  set to 2 to ONLY extract params, without plotting or stats
%qt_cond_only          = 0; % only extract data for the "QuickTime" condition so we can collapse across all experiments for this condition (always Condition 1)
%update_condition_labels_for_paper = 1; % update condition labels to match the order of presentation in the paper?

which_metrics = {'psy_pse'};%{'rev_pse' 'psy_pse'}; % which PSE metric should we consider? 'grw_pse' 'rev_pse' 'psy_pse'

% init
for metriccell = which_metrics
    g.(metriccell{1}).raw = [];
end

subj_i = 0;
for subjcell = subjlist'
    subj = subjcell{1}; % convert to string
    subj_i = subj_i+1;
    
    fprintf('processing %s...\n',subj)

    % RevampedHORRM - we only care abou conditon 5 here, which will be 1.4 in the paper
    dset = ['../data/' subj '_1_revampHORRM_full.mat'];
    if ~exist(dset,'file')
        dset = ['../data/' subj '_1_revampHOR_full.mat'];
    end
    if ~exist(dset,'file')
        dset = ['../data/' subj '_1_revamp_UPDATED_full.mat'];
    end
    if ~exist(dset,'file') && strcmp(subj,'REM') % added as conditional for VSS presentation
        dset = ['../data/' subj '_1_revamp_full.mat'];
    end
    if ~exist(dset,'file') && strcmp(subj,'DEC') % added as conditional for VSS presentation
        dset = ['../data/' subj '_1_BandE_full.mat'];
    end    
    if ~exist(dset,'file')
        error('unable to locate ''revamp'' dataset for subj %s',subj);
    end
    load(dset);
    
    %plot_staircases
    psychometric_from_staircase
    if ~keep_individual_plots
        close; close; close; close; % two plots from plot_staircases, two from psychometric_from_staircase
    end
    
    % add revamp condition 5 data
    for metriccell = which_metrics
        metric = metriccell{1}; % convert to string
        % for revamp, we only care about condition 1 (full oblique)
        eval(sprintf('g.(metric).raw(subj_i,1) = %s(:,1);',metric));
    end
    
    
    
    % LC (Last Control) - only one condition here, will be 3.1 in the paper
    dset = ['../data/' subj '_1_LC_full.mat'];
    if ~exist(dset,'file')
        error('unable to locate ''LC'' dataset for subj %s',subj);
    end
    load(dset);
    
    %plot_staircases
    psychometric_from_staircase
    if ~keep_individual_plots
        close; close; close; close; % two plots from plot_staircases, two from psychometric_from_staircase
    end
    
    % add revamp condition 5 data
    for metriccell = which_metrics
        metric = metriccell{1}; % convert to string
        % for LC, we only care about condition 1
        eval(sprintf('g.(metric).raw(subj_i,2) = %s(:,1);',metric));
    end
end



% update condition labels to match presentation for the paper
nconds = 2;
condition_list = [1.4 3.1]; % N.B. revampHORRM condition 5 is 1.4 in the paper


if norm_subj_pses
    % normalize all PSEs to the max observed across conditions
    for metriccell = which_metrics
        metric = metriccell{1}; % convert to string
        g.(metric).raw = bsxfun(@rdivide,g.(metric).raw,max(g.(metric).raw,[],2));
    end
end


% plot individual subjct PSEs on one plot
figure
k = 1;
for metriccell = which_metrics
    metric = metriccell{1}; % convert to string
    subplot(length(which_metrics),2,k);
    plot(g.(metric).raw','o-');
    setxrange(gca,0.5,nconds+0.5);
    %setyrange(gca,-0.5,1);
    hline(0,'-k');
    set(gca,'XTick',1:nconds,'XTickLabel',condition_list);
    xlabel('condition')
    ylabel('PSE (% rect)')
    title(regexprep(metric,'_','\\_'));

    
    % also add a plot with subject-means removed (for a visual of between-subj variance removed, see HELP CI)
    subplot(length(which_metrics),2,k+1)
    within_subj_only = bsxfun(@minus,g.(metric).raw,mean(g.(metric).raw,2)) + mean(g.(metric).raw(:)); % data - subj_means + grand mean
    plot(within_subj_only','o-');
    setxrange(gca,0.5,nconds+0.5);
    %setyrange(gca,-0.5,1);
    hline(0,'-k');
    set(gca,'XTick',1:nconds,'XTickLabel',condition_list);
    xlabel('condition')
    ylabel('PSE (% rect)')
    title(sprintf('%s - between subj variance',regexprep(metric,'_','\\_')));
    if k==1
        legend(regexprep(subjlist,'_','\\_'),'Location','Best');
    end
    
    k = k+2;
end


% test normality of the data
if test_for_normality
    fprintf('testing normality of distributions\n')
    for metriccell = which_metrics
        metric = metriccell{1}; % convert to string
        fprintf('\t%s\n',metric)
        for i = 1:nconds
            %%[g.(metric).norm.ks_h(i) g.(metric).norm.ks_p(i)] = kstest(g.(metric).raw(:,i));
            [g.(metric).norm.jb_h(i) g.(metric).norm.jb_p(i)] = jbtest(g.(metric).raw(:,i));
            [g.(metric).norm.lillie_h(i) g.(metric).norm.lillie_p(i)] = lillietest(g.(metric).raw(:,i));
        end
        disp(g.(metric).norm.jb_h)
        disp(g.(metric).norm.lillie_h)
    end
end


% group stats
g.pse_lab  = which_metrics; %{'all steps' 'reversals' 'psychmetric'};
g.pse_mean = [];
g.pse_sem  = [];
g.pse_ci95 = [];
for metriccell = which_metrics
    metric = metriccell{1}; % convert to string
    g.pse_mean = [g.pse_mean; mean(g.(metric).raw)];
    g.pse_sem  = [g.pse_sem; ste(g.(metric).raw)];
    g.pse_ci95 = [g.pse_ci95; ci(g.(metric).raw,95,1,2);];
end
% plot group-averaged PSEs
figure
subplot(2,1,1)
barerr2(condition_list,g.pse_mean',g.pse_ci95'); % grouped by condition
legend(regexprep(g.pse_lab,'_','\\_'),'Location','Best')
title(sprintf('group-average (n=%d)',nsubj))
xlabel('condition')
ylabel('PSE +/- 95%CI (% rect)')
subplot(2,1,2)
barerr2(g.pse_lab,g.pse_mean,g.pse_ci95); % grouped by PSE algorithm
legend(num2str(condition_list'),'Location','Best');
title(sprintf('group-average (n=%d)',nsubj))
xlabel('PSE-algorithm')
ylabel('PSE +/- 95%CI (% rect)')

% nicer figure for presentation/paper?
figure
barerr2(condition_list,g.pse_mean',g.pse_ci95','err_cap_width',0.25); % grouped by condition
legend(regexprep(g.pse_lab,'_','\\_'),'Location','Best')
%title(sprintf('group-average (n=%d)',nsubj))
setxrange(gca,0,6); % match bar width across experiments
setyrange(gca,-10,25); % match bar width across experiments
xlabel('condition')
ylabel('PSE +/- 95%CI')
colormap(repmat(.5,256,3)); % 50% gray bars
box off
print('-dpdf','./figures/lc_group_VSS.pdf')


if do_stats
    mi = 0;
    for metriccell = which_metrics
        metric = metriccell{1}; % convert to string
        mi = mi+1;
        
        % check all conditions against zero (uncorrected)
        [g.(metric).onesamp.ttest_h g.(metric).onesamp.ttest_p] = ttest(g.(metric).raw);
        for i = 1:size(g.(metric).raw,2)
            [g.(metric).onesamp.signrank_p(1,i) g.(metric).onesamp.signrank_h(1,i)] = signrank(g.(metric).raw(:,i));
        end
        % NEED TO ADD ONE-SAMPLE T-TEST.  THIS IS THERE AN STAND-ALONE FUNCITON FOR PERMUTATION TESTING FOR ONE-SAMPLE, PAIRED-SAMPLES, OR INDEPENDENT-SAMPLES (SEP CODE?) WOULD BE REALLY USEFUL!!
        
        
        % ANOVA for Psychometric Curve-Derived PSEs, but N.B.:
        % 1. not actually a repeated measures ANOVA in Matlab
        % 2. using no (LSD) correction for multiple comparisons and non-paired ttests
        % -> so this is all highly inappropriate for our data, but at least I can do it within the script.
        ff = fullfact([nsubj nconds]);
        this_subj = ff(:,1); % arbitrary subject id label (numerals)
        this_cond = condition_list(ff(:,2))'; % condition grouping variable
        this_data = g.(metric).raw(:); % dependent variable - vector
        
        % this isn't really a repeated-measures anova, but it yields the same results as SPSS (assuming sphericity assumption is met)
        [g.(metric).anovan_p g.(metric).anovan_tab g.(metric).anovan_stats] = anovan(this_data,{this_cond this_subj},'varnames',{'cond' 'subj'},'random',2,'display','off');
        %     % this is totally invalid - not based on paired t-tests...
        %     f = figure('MenuBar','none','Color','w');
        %     g.options.posthoc_test = 'hsd';
        %     [g.multcompare.comp g.multcompare.means g.multcompare.h g.multcompare.gnames] = multcompare(g.(metric).anovan_stats,'ctype',g.options.posthoc_test,'dimension',1);
        
        
        % nonparametric Friedman's test (repeated measures ANOVA equivilent)
        % matlab's FRIEDMAN performs a non-parametric two-way ANOVA, adjusting for row effects.
        % so subject should always be the row and ONE other factor will be colums
        if nconds == 1
            % friedman can't execute without multiple conditions
            g.(metric).friedman_p = NaN;
            g.(metric).friedman_tab = '';
            g.(metric).friedman_stats = struct();
        else
            [g.(metric).friedman_p g.(metric).friedman_tab g.(metric).friedman_stats] = friedman(g.(metric).raw,1,'off'); % rows are subj, columns are conditions
            %f = figure('MenuBar','none','Color','w');
            %[g.(metric).multcompare.comp g.(metric).multcompare.means g.(metric).multcompare.h g.(metric).multcompare.gnames] = multcompare(g.(metric).friedman_stats,'ctype','hsd');
        end
        
        
        % paired t-tests (i.e., uncorrected)
        % TODO: this could easily be written as a stand-alone function that takes a variety of arguments for nperms, exhaustive vs. random, metric (or a function handle), tails for p value, etc.  would be very useful in the future...but alas, this MUST wait until after my job apps.
        npairs = nconds*(nconds-1)/2;
        initfill = zeros(1,npairs); % init
        g.(metric).posthoc.t       = initfill;
        g.(metric).posthoc.d       = initfill;
        g.(metric).posthoc.wilcoxon_p = initfill;
        g.(metric).posthoc.wilcoxon_h = initfill;
        g.(metric).posthoc.lsd_h   = initfill;
        g.(metric).posthoc.lsd_p   = initfill;
        g.(metric).posthoc.hsd_h   = initfill;
        g.(metric).posthoc.sidak_h = initfill;
        g.(metric).posthoc.bonf_h  = initfill;
        
        % for paired t-test, only shuffle data within subjects
        if do_perm_test
            max_perms = 100000; % maximum number of permutations to run for every paired comparison.  if total number of possible combinations is less than the max, then we use random sampling
            pos_perms  = 2^nsubj; % total number of possible permutations for paired compairons
            if pos_perms <= max_perms
                g.(metric).posthoc.perm_type = 'exhaustive';
                g.(metric).posthoc.perm_n = pos_perms;
                g.(metric).posthoc.perm_iidx = combn(1:2,nsubj)'; % each column is the index for ivals (below)
            else
                g.(metric).posthoc.perm_type = 'random';
                g.(metric).posthoc.perm_n = max_perms;
                g.(metric).posthoc.perm_iidx = round(rand(nsubj,g.(metric).posthoc.perm_n)+1); % each column is the index for ivals (below)
            end
            g.(metric).posthoc.perm_t  = repmat(initfill,[g.(metric).posthoc.perm_n 1]);
            g.(metric).posthoc.perm_tp  = initfill;
            g.(metric).posthoc.perm_th  = initfill;
            g.(metric).posthoc.perm_d  = repmat(initfill,[g.(metric).posthoc.perm_n 1]);
            g.(metric).posthoc.perm_dp  = initfill;
            g.(metric).posthoc.perm_dh  = initfill;
            %             g.(metric).posthoc.perm_w  = repmat(initfill,[g.(metric).posthoc.perm_n 1]);
            %             g.(metric).posthoc.perm_wp  = initfill;
            %             g.(metric).posthoc.perm_wh  = initfill;
        end
        
        k = 0;
        for i = 1:nconds
            for j = i+1:nconds
                k = k+1;
                ivals = g.(metric).raw(:,i);
                jvals = g.(metric).raw(:,j);
                [g.(metric).posthoc.lsd_h(1,k), g.(metric).posthoc.lsd_p(1,k) ~, s] = ttest(ivals,jvals); % actual data
                g.(metric).posthoc.t(1,k) = s.tstat;
                
                g.(metric).posthoc.d(1,k) = mean(ivals-jvals);
                
                [g.(metric).posthoc.wilcoxon_p(1,k) g.(metric).posthoc.wilcoxon_h(1,k) s] = signrank(ivals,jvals);
                %                 g.(metric).posthoc.w(1,l) = s.signedrank;
                
                % permutations
                if do_perm_test
                    permstart = tic;
                    fprintf('starting permutation test (cond %d vs. %d)...',condition_list(i),condition_list(j))
                    for p = 1:size(g.(metric).posthoc.perm_iidx,2) % loop over total number of perms
                        ijvals = [ivals jvals];
                        pidx = [g.(metric).posthoc.perm_iidx(:,p) 3-g.(metric).posthoc.perm_iidx(:,p)];
                        pvals = NaN(size(ijvals));
                        for ii = 1:length(ivals)
                            pvals(ii,:) = ijvals(ii,pidx(ii,:));
                        end                        
                        [~, ~, ~, ps] = ttest(pvals(:,1),pvals(:,2)); % permuted data
                        g.(metric).posthoc.perm_t(p,k) = ps.tstat;
                        
                        g.(metric).posthoc.perm_d(p,k) = mean(pvals(:,1)-pvals(:,2)); % mean difference
                        
                        %                         [~, ~, ps] = signrank(pvals(:,1),pvals(:,2)); % Wilcoxon
                        %                         g.(metric).posthoc.perm_w(p,k) = ps.signedrank;
                    end
                    fprintf('done in %.02f s\n',toc(permstart))
                end
            end
        end
        
        % tukey hsd
        g.(metric).posthoc.tukey_q = 4.755; % HARDCODED: k=5, df=9; alpha=0.05 <- http://cse.niaes.affrc.go.jp/miwa/probcalc/s-range/srng_tbl.html
        g.(metric).posthoc.pairedq = g.(metric).posthoc.t*sqrt(2);
        g.(metric).posthoc.hsd_h = g.(metric).posthoc.pairedq > g.(metric).posthoc.tukey_q;
        
        % bonferonni
        % index into g.(metric).posthoc. fields for a priori comparisons of interest
        switch nconds
            % a priori, we care about only some of the possible comparisons for each experiment
            case 5 % revampHORRM
                % QuickTime-Horizontal: 1v2; TargetPos: 2v3,2v5,3v5; Pursuit: 2v4
                apriori_idx = [1 5 6 7 9]; 
            case 3 % C3
                % QuickTime-TargetEccen 1v2; QuickTime-StaticBox 1v3
                apriori_idx = [1 2];
            case 2 % BandE
                % QuickTime-BeginEnd 1v2
                apriori_idx = 1;
            case 1 % LC
                apriori_idx = []; % no pairs, would need to compare with revampHORRM
            otherwise
                error('not sure which experiment this is based on nconds (%d).  can''t determine a priori hypotheses',nconds)
        end
        apriori_pairs = length(apriori_idx);
        g.(metric).posthoc.bonf_h = g.(metric).posthoc.lsd_p < 0.05/apriori_pairs;
        
        % dunn-sidak
        g.(metric).posthoc.sidak_h = g.(metric).posthoc.lsd_p < 1-(1-0.05)^(1/apriori_pairs);
        
        % bonferroni-holm stepwise correction...(sort p-vals, check if lowest<alpha/npairs, next<alpha/(npairs-1), ect.)
        
        % permuation test for paired data using t-statistic
        if do_perm_test
             % what is the probability of getting a SOME STAT at least as high as observed, by chance. two-tailed.
             
            % on the t-statistic
            g.(metric).posthoc.perm_tp = mean(bsxfun(@le,abs(g.(metric).posthoc.t),abs(g.(metric).posthoc.perm_t)));
            g.(metric).posthoc.perm_th = g.(metric).posthoc.perm_tp < 0.05;
            % bonferroni and sidak correction on permutation test
            g.(metric).posthoc.perm_th_bonf  = g.(metric).posthoc.perm_tp < 0.05/apriori_pairs;
            g.(metric).posthoc.perm_th_sidak = g.(metric).posthoc.perm_tp < 1-(1-0.05)^(1/apriori_pairs);
            
            % on the mean difference
            g.(metric).posthoc.perm_dp = mean(bsxfun(@le,abs(g.(metric).posthoc.d),abs(g.(metric).posthoc.perm_d)));
            g.(metric).posthoc.perm_dh = g.(metric).posthoc.perm_dp < 0.05;
            % bonferroni and sidak correction on permutation test
            g.(metric).posthoc.perm_dh_bonf  = g.(metric).posthoc.perm_dp < 0.05/apriori_pairs;
            g.(metric).posthoc.perm_dh_sidak = g.(metric).posthoc.perm_dp < 1-(1-0.05)^(1/apriori_pairs);

            %             % on wilcoxon signed rank statistic
            %             g.(metric).posthoc.perm_wp = mean(bsxfun(@le,abs(g.(metric).posthoc.w),abs(g.(metric).posthoc.perm_w))); % what is the probability of getting a t-stat at least as high as observed, by chance.
            %             g.(metric).posthoc.perm_wh = g.(metric).posthoc.perm_wp < 0.05;
            %             % bonferroni and sidak correction on permutation test
            %             g.(metric).posthoc.perm_wh_bonf  = g.(metric).posthoc.perm_wp < 0.05/apriori_pairs;
            %             g.(metric).posthoc.perm_wh_sidak = g.(metric).posthoc.perm_wp < 1-(1-0.05)^(1/apriori_pairs);
        end
        
        
        
        % print some results to command window
        condition_list
        eval(sprintf('%s_mean = g.pse_mean(mi,:)',metric))
        onesamp_ttest = g.(metric).onesamp.ttest_p
        onesamp_signrank = g.(metric).onesamp.signrank_p
        
        rm_anova_p = g.(metric).anovan_p
        g.psy_pse.anovan_tab
        friedman_p = g.(metric).friedman_p
        g.psy_pse.friedman_tab
        if do_perm_test
            %perm_posthoc_t_sig = squareform(g.(metric).posthoc.perm_th)
            perm_posthoc_d_sig = squareform(g.(metric).posthoc.perm_dh)
            perm_posthoc_d_pval = squareform(g.(metric).posthoc.perm_dp)
            %perm_posthoc_w_sig = squareform(g.(metric).posthoc.perm_wh)
            %perm_posthoc_w_pval = squareform(g.(metric).posthoc.perm_wp)
        end
        wilcoxon_posthoc_sig = squareform(g.(metric).posthoc.wilcoxon_h)
        wilcoxon_posthoc_pval = squareform(g.(metric).posthoc.wilcoxon_p)
        lsd_posthoc_sig = squareform(g.(metric).posthoc.lsd_h)
        lsd_posthoc_pval = squareform(g.(metric).posthoc.lsd_p)
        %hsd_posthoc_sig = squareform(g.(metric).posthoc.hsd_h)
        %bonf_posthoc_sig = squareform(g.(metric).posthoc.bonf_h)
        %sidak_posthoc_sig = squareform(g.(metric).posthoc.sidak_h)
    end
end

diary off; % stop storing command window output to text file