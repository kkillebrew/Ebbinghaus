% wrapper for group analysis of the Dynamic Ebbinghaus-Corridor comparison
% expeirment using the adaptive staircase.

clear g
%close all

% Other Early piloting and author's datasets (which we are excluding)
%datasets = {};


% Fall 2015, Pilot (P2)
pilot1_experienced = { % me, chris, etc. EXCLUDE
    'rmru_ebbcorr_P2_110715_001'; % ryan_macair setup, first full run through of P2 (corridor periphfix)
    };
pilot1_sona = { % subjects to include from WSU setup
    };

pilot1_exclude = { % subjects to exclude for some reason
    };
pilot1 = [pilot1_sona;]; % all Experiment 1 (exp1) subjects

% subject info:
% WSU SONA subjs:
%
% OVERALL
%   youngest is ??; oldest is ??
%   NFemale = ??; NMale = ??



% all subjects for a condition 1 grand average (what about subjects in multiple experiments?)
allsubj = [pilot1];


% loop over a set of experiments for batch regeneration of saved figures/stats
all_which_exp = {'pilot1'}; % list of all variables that store datasets to loop over



for which_exp_cell = all_which_exp
    
    % which set of data to analyze
    which_exp = which_exp_cell{1}; % string of variable name that supplies datasets.  this is used so files can be saved with this string in the filename
    eval(sprintf('datasets = %s;',which_exp));
    
    % save command window output to textfile for quick access to stats
    diary_file = sprintf('./outputs/%s_output.txt',which_exp);
    if exist(diary_file,'file')
        delete(diary_file); % always start with a blank diary
    end
    diary(diary_file);
    
    which_exp
    nsubj = length(datasets)
    
    %norm_subj_pses        = 0; % should we normalize subject PSEs to a particular condition?
    keep_individual_plots = 1; % boolean, should we keep individual figures on screen?
    do_stats              = 0; % run various stats on the data and print to screen?
    test_for_normality    = 0; % do tests for normality of data in each condition?
    do_perm_test          = 1; % run the permutation test?
    extract_params        = 0; % should we extract params?  set to 2 to ONLY extract params, without plotting or stats
  
    which_metrics = {'psy_pse'}; % which PSE metric should we consider? 'adj_pse' is only option here
    
    % init
    for metriccell = which_metrics
        g.(metriccell{1}).lab = {};
        g.(metriccell{1}).norm = [];
        g.(metriccell{1}).context = [];
        %g.(metriccell{1}).context_bydir = [];
        %g.(metriccell{1}).context_bydirANDsize = [];
        g.(metriccell{1}).nocontext = [];
        %g.(metriccell{1}).nocontext_bydir = [];
        %g.(metriccell{1}).nocontext_bydirANDsize = [];
    end
    
    
    dset_i = 0;
    for dsetcell = datasets'
        dset = dsetcell{1}; % convert to string
        dset_i = dset_i+1;
        
        fprintf('processing %s...\n',dset)
        
        load(['../data/' dset]);
        
        if extract_params
            % store some useful info about experimental parameters
            for fncell = {'experiment' 'half_cycle_dur' 'standard.widths_deg'}
                fn = fncell{1};
                eval(sprintf('p.%s(dset_i,:) = %s;',fn,fn));
            end
        end
        
        if extract_params~=2
            % plot individual subject data using quickanalysis script
            plot_staircases; % provides avg_pse and rev_pse
            psychometric_from_staircase % provides psy_pse     ?? provides subj structure for subject-specific data
            
            if ~keep_individual_plots
                close; % one plots from plot_pses
            end
            
            % init
            for metriccell = which_metrics
                metric = metriccell{1}; % convert to string
                %                 if qt_cond_only
                %                     % only extract condition 1 (which should be the first column, but we can be extra careful about that...)
                %                     cond_ref = this.cond_lab==1;
                %                 else
                % extract all conditions
                cond_ref = true(size(subj.adj_pse.lab));
                %                 end
                eval(sprintf('g.(metric).lab = cat(1,g.(metric).lab,subj.%s.lab(:,cond_ref));',metric));
                eval(sprintf('g.(metric).norm = cat(1,g.(metric).norm,subj.%s.norm(:,cond_ref));',metric));
                eval(sprintf('g.(metric).context = cat(1,g.(metric).context,subj.%s.context(:,cond_ref));',metric));
                eval(sprintf('g.(metric).nocontext = cat(1,g.(metric).nocontext,subj.%s.nocontext(:,cond_ref));',metric));
                
                %                 cond_ref_for_bydir = repmat(cond_ref,1,length(params.targ_inducer_sizes)); % account for multiple targ-inducer sizes, if necessary
                %                 eval(sprintf('g.(metric).norm_bydir = cat(1,g.(metric).norm_bydir,subj.%s.norm_bydir(:,cond_ref_for_bydir));',metric));
                %                 eval(sprintf('g.(metric).nocontext_bydir = cat(1,g.(metric).nocontext_bydir,subj.%s.nocontext_bydir(:,cond_ref_for_bydir));',metric));
                %                 eval(sprintf('g.(metric).context_bydir = cat(1,g.(metric).context_bydir,subj.%s.context_bydir(:,cond_ref_for_bydir));',metric));
                %
                %                 if isfield(subj.(metric),'norm_bydirANDsize')
                %                     cond_ref_for_bydirANDsize = repmat(cond_ref,1,length(params.targ_inducer_sizes)*length(params.standard.widths_deg)); % account for multiple targ-inducer sizes and standard sizes, if necessary
                %                     eval(sprintf('g.(metric).norm_bydirANDsize = cat(1,g.(metric).norm_bydirANDsize,subj.%s.norm_bydirANDsize(:,cond_ref_for_bydirANDsize));',metric));
                %                     eval(sprintf('g.(metric).nocontext_bydirANDsize = cat(1,g.(metric).nocontext_bydirANDsize,subj.%s.nocontext_bydirANDsize(:,cond_ref_for_bydirANDsize));',metric));
                %                     eval(sprintf('g.(metric).context_bydirANDsize = cat(1,g.(metric).context_bydirANDsize,subj.%s.context_bydirANDsize(:,cond_ref_for_bydirANDsize));',metric));
                %                 end

                fprintf('\n')
            end
        end
    end
    
    
    if extract_params
        % calculate some params and return
        p.init_anchor_jitter = params.init_anchor_jitter;
        p.targ_pos = params.start_pos_deg;
        p.standard_pos = params.standard_pos_deg;
        p.pos_jitter = [params.start_pos_jitter_deg params.standard_pos_jitter_deg];
        p.full_cycle_dur = params.full_cycle_dur;
        p.half_cycle_dur = params.half_cycle_dur;
        p.trans_ang = params.trans_ang_rad * 180/pi;
        p.trans_dist = params.trans_dist_deg;
        p.trans_rate = p.trans_dist / p.half_cycle_dur;
        if int32(100*p.trans_ang)/100 == -60 % avoid rounding errors...
            p.ur_inducer_trans_dist = p.trans_dist; % assumes 60 deg p.trans_ang
        else
            p.ur_inducer_trans_dist = NaN; % geometry is more difficult here
        end
        p.context_small_eccen = params.context.small_eccen_deg;
        p.context_small_width = params.context.small_width_deg;
        p.context_large_eccen = params.context.large_eccen_deg;
        p.context_large_width = params.context.large_width_deg;
        p.min_targ_width = params.targ.min_width_deg;
        p.max_targ_width_small_context = 2*(params.context.small_eccen_deg-params.context.small_width_deg/2);
        p.max_targ_width_large_context = 2*(params.context.large_eccen_deg-params.context.large_width_deg/2);
        p.mouse_targ_width_resolution = (params.targ.max_width_deg - params.targ.min_width_deg)/params.y_resolution;
        p.mouse_targ_rate_resolution = p.mouse_targ_width_resolution / p.half_cycle_dur;
        p.jitter = params.standard.jitter;
        
        if extract_params==2
            return
        end
    end
    
    
    
    
    
    %     if norm_subj_pses
    %         % normalize all PSEs to the max observed across conditions
    %         for metriccell = which_metrics
    %             metric = metriccell{1}; % convert to string
    %             g.(metric).norm = bsxfun(@rdivide,g.(metric).norm,max(g.(metric).norm,[],2));
    %         end
    %     end
    
    
    % plot individual subjct PSEs on one plot
    figure
    k = 1;
    for metriccell = which_metrics
        metric = metriccell{1}; % convert to string
        this.ncond = size(g.(metric).lab,2);
        this.cond_lab = g.(metric).lab(1,:); % assume same for all subjects

        subplot(length(which_metrics),2,k);
        plot(g.(metric).norm','o-');
        xlim([0.5 this.ncond+0.5]);%setxrange(gca,0.5,this.ncond+0.5);
        %setyrange(gca,-0.5,1);
        hline(0,'-k');
        set(gca,'XTick',1:this.ncond,'XTickLabel',this.cond_lab);
        xlabel('condition')
        ylabel('Illusion Magnitude (% standard)')
        title(regexprep(metric,'_','\\_'));
        
        
        % also add a plot with subject-means removed (for a visual of between-subj variance removed, see HELP CI)
        subplot(length(which_metrics),2,k+1)
        within_subj_only = bsxfun(@minus,g.(metric).norm,mean(g.(metric).norm,2)) + mean(g.(metric).norm(:)); % data - subj_means + grand mean
        plot(within_subj_only','o-');
        xlim([0.5 this.ncond+0.5]);%setxrange(gca,0.5,this.ncond+0.5);
        %setyrange(gca,-0.5,1);
        hline(0,'-k');
        set(gca,'XTick',1:this.ncond,'XTickLabel',this.cond_lab);
        xlabel('condition')
        ylabel('Illusion Magnitude (% standard)')
        title(sprintf('%s - between subj variance',regexprep(metric,'_','\\_')));
        if k==1
            legend(regexprep(datasets,'_','\\_'),'Location','Best');
        end
        
        k = k+2;
    end
        
    
    
    % test normality of the data
    % jb and lillie tests both return significant results if data depart
    %    from normality (i.e., Ho = data is normal)
    if test_for_normality
        fprintf('testing normality of distributions\n')
        for metriccell = which_metrics
            metric = metriccell{1}; % convert to string
            this.ncond = size(g.(metric).lab,2);
            this.cond_lab = g.(metric).lab(1,:); % assume same for all subjects

            fprintf('\t%s\n',metric)
            for i = 1:this.ncond
                %%[g.(metric).norm.ks_h(i) g.(metric).norm.ks_p(i)] = kstest(g.(metric).norm(:,i));
                [g.(metric).norm.jb_h(i) g.(metric).norm.jb_p(i)] = jbtest(g.(metric).norm(:,i));
                [g.(metric).norm.lillie_h(i) g.(metric).norm.lillie_p(i)] = lillietest(g.(metric).norm(:,i));
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
    
    g.pse_bydir_mean = [];
    g.pse_bydir_sem  = [];
    g.pse_bydir_ci95 = [];
    
    for metriccell = which_metrics
        metric = metriccell{1}; % convert to string
        g.pse_mean = [g.pse_mean; mean(g.(metric).norm)];
        g.pse_sem  = [g.pse_sem; ste(g.(metric).norm)];
        g.pse_ci95 = [g.pse_ci95; ci(g.(metric).norm,95,1,0);]; % set last argument to 2 for within subject correction
        
        %         g.pse_bydir_mean = [g.pse_bydir_mean; mean(g.(metric).norm_bydir)];
        %         g.pse_bydir_sem  = [g.pse_bydir_sem; ste(g.(metric).norm_bydir)];
        %         g.pse_bydir_ci95 = [g.pse_bydir_ci95; ci(g.(metric).norm_bydir,95,1,0);]; % set last argument to 2 for within subject correction
    end
    % *** assuming all metrics have some conditions
    % *** may not be the case if I use this to separate standard size and
    %     targ-inducer size
    this.ncond = size(g.adj_pse.norm,2);
    this.cond_lab = g.adj_pse.lab(1,:);
    
    % plot group-averaged PSEs
    figure
    %subplot(2,1,1)
    barerr2(this.cond_lab,g.pse_mean',g.pse_ci95','err_cap_width',0,'FaceColor','g'); % grouped by condition
    legend(regexprep(g.pse_lab,'_','\\_'),'Location','Best')
    %title(sprintf('group-average (n=%d)',nsubj))
    xlim([0.5 this.ncond+0.5]);%setxrange(gca,0.5,this.ncond+0.5); % match bar width across experiments
    %setyrange(gca,0,60); % match bar width across experiments
    xlabel('condition')
    ylabel('Illusion Magnitude +/- 95% CI')
    %colormap(repmat(.5,256,3)); % 50% gray bars
    box off
    %saveas(gcf,sprintf('./figures/%s_group.fig',which_exp)); % save figure to a file
    print('-dpdf',sprintf('./figures/%s_group.pdf',which_exp))
    
%     subplot(2,1,2)
%     barerr2(this.cond_lab,reshape(g.pse_bydir_mean,[],length(params.targ_inducer_sizes)),reshape(g.pse_bydir_ci95,[],length(params.targ_inducer_sizes)),'err_cap_width',0); % grouped by condition
%     legend({'targ-in-small' 'targ-in-large'},'Location','Best')
%     %title(sprintf('group-average (n=%d)',nsubj))
%     xlim([0.5 this.ncond+0.5]);%setxrange(gca,0.5,this.ncond+0.5); % match bar width across experiments
%     %setyrange(gca,0,60); % match bar width across experiments
%     xlabel('condition')
%     ylabel('Illusion Magnitude +/- 95% CI')
%     box off
%     %saveas(gcf,sprintf('./figures/%s_group.fig',which_exp)); % save figure to a file
%     print('-dpdf',sprintf('./figures/%s_group.pdf',which_exp))
    
    
    % **** THIS SECTION IS A HACK - VERY-MUCH HARDCODED ****
%     % first, collapse PSEs for noinducer condition across targ-indducer
%     % size (which is meaningless without inducer present).  Note, this is
%     % only necessary (and possible) for Exp 1 (P2)
%     if size(g.adj_pse.lab,2) == 3 % then this is Exp 1 (P2)
%         tmptmp = g.adj_pse.nocontext_bydir; % make a copy
%         g.adj_pse.nocontext_bydir = []; % clear
%         for i = 1:3 % there are 3 conditions in Exp 1 (P2)
%             g.adj_pse.nocontext_bydir(:,i) = mean(tmptmp(:,[i i+3]),2);
%         end
%     end
    
    % plot group-averaged PSEs, but split by with and without context
    figure
    subplot(2,1,1)
    im = mean(g.adj_pse.context);
    ici = ci(g.adj_pse.context,95,1,0); % set last argument to 2 for within subject correction
    COLLAPSE_FOR_NOINDUCER = 1; % should we collapse across the targ-inducer sizes to plot the noinducer data?  we DON'T do this for Illusion Magnitude estimates, but it might be confusing to plot that way (for Exp1)
    if COLLAPSE_FOR_NOINDUCER
        nim = mean(g.adj_pse.nocontext);
        nim = repmat(nim,[1 2]);
        nici = ci(g.adj_pse.nocontext,95,1,0); % set last argument to 2 for within subject correction
        nici = repmat(nici,[1 2]);
    else
        nim = mean(g.adj_pse.nocontext_bydir);
        nici = ci(g.adj_pse.nocontext_bydir,95,1,0); % set last argument to 2 for within subject correction
    end
    

    nim = repmat(nim,[1 2]);
    hold on;
    cs = 'brmygc';
    switch size(g.adj_pse.lab,2)
%         case 3 % Exp 1 (P2)
%             thisn = length(im)/2; % PARSING THIS DEPENDS ON EXP1 VS. EXP2
%             for i = 1:thisn
%                 %plot([1 2]+(i-1)*2,[im(i) nim(i)],[cs(i) '--o']); % targ-in-small
%                 errorbar([1 2]+(i-1)*2,[im(i) nim(i)],[ici(i) nici(i)],[cs(i) '--']); % targ-in-small
%                                 
%                 %plot([1 2]+(i-1)*2,[im(i+thisn) nim(i+thisn)],[cs(i) '-o']); % targ-in-large
%                 errorbar([1 2]+(i-1)*2,[im(i+thisn) nim(i+thisn)],[ici(i+thisn) nici(i+thisn)],[cs(i) '-']); % targ-in-large
%             end
%             legend({'targ-in-small' 'targ-in-large'});
        case 6 % Pilot 1 (P1)
            thisn = length(im); % PARSING THIS DEPENDS ON EXP1 VS. EXP2
            for i = 1:thisn
                %plot([1 2]+(i-1)*2,[im(i) nim(i)],[cs(i) '-o']); % targ-in-large
                errorbar([1 2]+(i-1)*2,[im(i) nim(i)],[ici(i) nici(i)],[cs(i) '-']); % targ-in-large
            end
            legend({'targ-in-large'}); % no targ-in-small for Exp 2 (P4)
        otherwise
            error('can''t figure out which experiment this is based on length of g.adj_pse.lab...')
    end
    xlim([0.5 2*thisn+0.5]);%setxrange(gca,0.5,this.ncond+0.5); % match bar width across experiments
    set(gca,'XTick',1.5:2:2*thisn,'XTickLabel',g.adj_pse.lab(1,:))
    xlabel('condition');
    ylabel('PSE +/- 95% CI');
    vline(2.5:2:2*thisn,':k');
    hline(0,'k-');
    title('group-averaged data, with and without context')
    
    
    % group-averaged PSEs, but split by with and without context, in bar
    % chart form, collapsed across targ-inducer size and standard size
    subplot(2,1,2)
    
    this.data = [mean(g.adj_pse.context); mean(g.adj_pse.nocontext)];
    this.err  = [ci(g.adj_pse.norm,95,1,0); ci(g.adj_pse.nocontext,95,1,0)]; % set last arg of ci() to 2 for within-subject correction.
    barerr2(this.cond_lab,this.data',this.err','err_cap_width',0); % grouped by condition
    hline(0,'k-');
    legend({'with-context' 'without-context'},'Location','Best')
    title('group-averaged data, with and without context')
    xlim([0.5 this.ncond+0.5]); % match bar width across experiments
    %setyrange(gca,0,60); % match bar width across experiments
    xlabel('condition')
    ylabel('PSE +/- 95% CI')
    %colormap(repmat(.5,256,3)); % 50% gray bars
    box off
    print('-dpdf',sprintf('./figures/%s_group_nocontext.pdf',which_exp))
    % **** THIS SECTION IS A HACK - VERY-MUCH HARDCODED ****
    

    
    
    
    
    if do_stats
        mi = 0;
        for metriccell = which_metrics
            metric = metriccell{1}; % convert to string
            mi = mi+1;
            
            % check all norm (with-without context) conditions against zero (uncorrected)
            [g.(metric).onesamp.ttest_h g.(metric).onesamp.ttest_p] = ttest(g.(metric).norm);
            for i = 1:size(g.(metric).norm,2)
                [g.(metric).onesamp.signrank_p(1,i) g.(metric).onesamp.signrank_h(1,i)] = signrank(g.(metric).norm(:,i));
            end
            % NEED TO ADD ONE-SAMPLE PERMTEST.  THIS IS WHERE AN STAND-ALONE FUNCITON FOR PERMUTATION TESTING FOR ONE-SAMPLE, PAIRED-SAMPLES, OR INDEPENDENT-SAMPLES (SEP CODE?) WOULD BE REALLY USEFUL!!

            
            % check all no-inducer conditions against zero (uncorrected)
            [g.(metric).noinducer_onesamp.ttest_h g.(metric).noinducer_onesamp.ttest_p] = ttest(g.(metric).nocontext);
            for i = 1:size(g.(metric).nocontext,2)
                [g.(metric).noinducer_onesamp.signrank_p(1,i) g.(metric).noinducer_onesamp.signrank_h(1,i)] = signrank(g.(metric).nocontext(:,i));
            end
            % NEED TO ADD ONE-SAMPLE PERMTEST.  THIS IS WHERE AN STAND-ALONE FUNCITON FOR PERMUTATION TESTING FOR ONE-SAMPLE, PAIRED-SAMPLES, OR INDEPENDENT-SAMPLES (SEP CODE?) WOULD BE REALLY USEFUL!!

            
            
            % ANOVA for Illusion Magnitudes, but N.B.:
            % 1. not actually a repeated measures ANOVA in Matlab
            % 2. using no (LSD) correction for multiple comparisons and non-paired ttests
            % -> so this is all highly inappropriate for our data, but at least I can do it within the script.
            [ff1 ff2] = ndgrid(1:nsubj,1:this.ncond); ff = [ff1(:) ff2(:)]; % replacement for fullfact (more efficient than fullyfact
            this_subj = ff(:,1); % arbitrary subject id label (numerals)
            this_cond = this.cond_lab(ff(:,2))'; % condition grouping variable
            this.data = g.(metric).norm(:); % dependent variable - vector
            
            % this isn't really a repeated-measures anova, but it yields the same results as SPSS (assuming sphericity assumption is met)
            [g.(metric).anovan_p g.(metric).anovan_tab g.(metric).anovan_stats] = anovan(this.data,{this_cond this_subj},'varnames',{'cond' 'subj'},'random',2,'display','off');
            %     % this is totally invalid - not based on paired t-tests...
            %     f = figure('MenuBar','none','Color','w');
            %     g.options.posthoc_test = 'hsd';
            %     [g.multcompare.comp g.multcompare.means g.multcompare.h g.multcompare.gnames] = multcompare(g.(metric).anovan_stats,'ctype',g.options.posthoc_test,'dimension',1);
            
            
            % nonparametric Friedman's test (repeated measures ANOVA equivilent)
            % matlab's FRIEDMAN performs a non-parametric two-way ANOVA, adjusting for row effects.
            % so subject should always be the row and ONE other factor will be colums
            if this.ncond == 1
                % friedman can't execute without multiple conditions
                g.(metric).friedman_p = NaN;
                g.(metric).friedman_tab = '';
                g.(metric).friedman_stats = struct();
            else
                [g.(metric).friedman_p g.(metric).friedman_tab g.(metric).friedman_stats] = friedman(g.(metric).norm,1,'off'); % rows are subj, columns are conditions
                %f = figure('MenuBar','none','Color','w');
                %[g.(metric).multcompare.comp g.(metric).multcompare.means g.(metric).multcompare.h g.(metric).multcompare.gnames] = multcompare(g.(metric).friedman_stats,'ctype','hsd');
            end
            
            
            % paired t-tests (i.e., uncorrected)
            % TODO: this could easily be written as a stand-alone function that takes a variety of arguments for nperms, exhaustive vs. random, metric (or a function handle), tails for p value, etc.  would be very useful in the future...but alas, this MUST wait until after my job apps.
            npairs = this.ncond*(this.ncond-1)/2;
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
            for i = 1:this.ncond
                for j = i+1:this.ncond
                    k = k+1;
                    ivals = g.(metric).norm(:,i);
                    jvals = g.(metric).norm(:,j);
                    [g.(metric).posthoc.lsd_h(1,k), g.(metric).posthoc.lsd_p(1,k) ~, s] = ttest(ivals,jvals); % actual data
                    g.(metric).posthoc.t(1,k) = s.tstat;
                    
                    g.(metric).posthoc.d(1,k) = mean(ivals-jvals);
                    
                    [g.(metric).posthoc.wilcoxon_p(1,k) g.(metric).posthoc.wilcoxon_h(1,k) s] = signrank(ivals,jvals);
                    %                 g.(metric).posthoc.w(1,l) = s.signedrank;
                    
                    % permutations
                    if do_perm_test
                        permstart = tic;
                        fprintf('starting permutation test (cond %s vs. %s)...',this.cond_lab{i},this.cond_lab{j})
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
            switch this.ncond
                % a priori, we care about only some of the possible comparisons for each experiment
                case 6 % P1
                    % compare all conditions for now
                    apriori_idx = [1 2 3 4 5 6];
                otherwise
                    error('not sure which experiment this is based on this.ncond (%d).  can''t determine a priori hypotheses',this.ncond)
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
            noinducer_mean = mean(g.(metric).nocontext)
            noinducer_std  = std(g.(metric).nocontext)
            noinducer_ste  = ste(g.(metric).nocontext)
            noinducer_onesamp_ttest = g.(metric).noinducer_onesamp.ttest_p
            noinducer_onesamp_signrank = g.(metric).noinducer_onesamp.signrank_p

            this.cond_lab
            %eval(sprintf('%s_mean = g.pse_mean(mi,:)',metric))
            norm_mean = mean(g.(metric).norm)
            norm_std  = std(g.(metric).norm)
            norm_ste  = ste(g.(metric).norm)

            
            onesamp_ttest = g.(metric).onesamp.ttest_p
            onesamp_signrank = g.(metric).onesamp.signrank_p
            
            rm_anova_p = g.(metric).anovan_p
            g.adj_pse.anovan_tab
            friedman_p = g.(metric).friedman_p
            g.adj_pse.friedman_tab
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
    
end % end loop over different experiments (i.e., which_exp variable)