% wrapper for combining two dirty analysis scripts
plot_staircases
psychometric_from_staircase

% plot all PSEs together
figure('Position',[521 27 386 420]);
bar([grw_pse;rev_pse;psy_pse]')
set(gca,'XTick',1:nconds,'XTickLabel',condition_list);
legend({'all steps' 'reversals' 'psychmetric'})
xlabel('condition')
ylabel('PSE')
