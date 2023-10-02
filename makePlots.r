input_directory = 'data/exp1/varyEfficiency'
reduced_data <- read_delim(file.path(input_directory, 'output_summarized.csv'), delim=',')

metrics = c('FP' = 'False Positives',
            'TN' = 'True Negatives',
            'TP' = 'True Postitives',
            'FN' = 'False Negatives',
            'P'  = 'Actual Positive',
            'N'  = 'Actual Negative',
            'PP' = 'Predicted Positive',
            'PN' = 'Predicted Negative',
            'total' = 'Total Population',
            'prevalence' = 'Prevalence of Transmitters',
            'ACC' = 'Accuracy',
            'PPV' = 'Positive Predictive Value',
            'FDR' = 'False Discovery Rate',
            'FOR' = 'False Omission Rate',
            'NPV' = 'Negative Predictive Value',
            'TPR' = 'True Positive Rate',
            'FPR' = 'False Positive Rate',
            'FNR' = 'False Negative Rate',
            'TNR' = 'True Negative Rate',
            'LRP' = 'Positive Likelihood Ratio',
            'LRN' = 'Negative Likelihood Ratio',
            'DOR' = 'Diagnostic Odds Ratio',
            'TS' = 'Threat Score',
            'informedness' = 'Informedness',
            'markedness' = 'Markedness',
            'prevalenceThreshold' = 'Prevalence Threshold',
            'MCC' = 'Matthews Correlation Coefficient',
            'FMI' = 'Fowlkes-Mallows',
            'F1' = 'F1 Score',
            'BA' = 'Balanced Accuracy',
            'RR' = 'Relative Risk',
            'RRall' = 'Relative Risk vs All',
            'contacts_per_index' = 'Contacts Per Index',
            'unique_contacts_per_index' = 'Contacts Per Index (no duplicates)',
            'N_HIV' = 'Number HIV+',
            'N_index' = 'Number of Index cases enrolled',
            'N_untreated' = 'Number Untreated',
            'N_undiagnosed' = 'Number Undiagnosed',
            'contact_mean' = 'Mean Contacts',
            'contact_std' = 'STD of Contacts',
            'percent_HIV' = 'Percent HIV+',
            'NNTTP' = 'NNTT Positive',
            'NNTTT' = 'NNTT Untreated',
            'NNTTD' = 'NNTT Undiagnosed'
           )

titles = c('FP' = 'Number of Contacts Without Trait',
           'TN' = 'Number of Non-Contacts Without Trait',
           'TP' = 'Number of Contacts With Trait',
           'FN' = 'Number of Non-Contacts With Trait',
           'P'  = 'Number of People With Trait',
           'N'  = 'Number of People Without Trait',
           'PP' = 'Number of Contacts of Tracing Group',
           'PN' = 'Number of Non-Contacts of Tracing Group',
           'total' = 'Total Population',
           'prevalence' = 'Prevalence of Trait',
           'ACC' = 'C w/ Trait plus NC w/o Trait divide by Total Population',
           'PPV' = 'Contacts With Trait / All Contacts',
           'FDR' = 'Contacts Without Trait / All Contacts',
           'FOR' = 'Non-Contacts With Trait / All Non-Contacts',
           'NPV' = 'Non-Contacts Without Trait / All Non-Contacts',
           'TPR' = 'Contacts With Trait / # People With Trait',
           'FPR' = 'Contacts Without Trait / # People Without Trait',
           'FNR' = 'Non-Contacts With Trait / # People With Trait',
           'TNR' = 'Non-Contacts Without Trait / # People Without Trait',
           'LRP' = 'True Postive Rate / False Postive Rate',
           'LRN' = 'False Negative Rate / True Negative Rate',
           'DOR' = 'Positive Likelihood Ration / Negative Likelihood Ratio',
           'TS' = 'Contacts With Trait / (All Contacts + Non-Contacts With Trait)',
           'informedness' = 'True Positive Rate + True Negative Rate - 1',
           'markedness' = 'Positive Predictive Value + Negatiev Predictive Value - 1',
           'prevalenceThreshold' = 'Prevalence Threshold',
           'MCC' = '~ Correlation Between Being a Contact and Having Trait',
           'FMI' = 'Fowlkes-Mallows',
           'F1' = 'F1 Score',
           'BA' = 'Mean of True Postive Rate and True Negative Rate',
           'RR' = 'Relative Risk',
           'RRall' = 'Relative Risk vs All',
           'contacts_per_index' = 'Contacts Per Index',
           'unique_contacts_per_index' = 'Contacts Per Index  (Only Contact to One Index)',
           'N_HIV' = 'Number of Contacts HIV+ on Trace',
           'N_index' = 'Number of Index cases enrolled',
           'N_untreated' = 'Number of Contacts Untreated on Trace',
           'N_undiagnosed' = 'Number of Contacts Undiagnosed on Trace',
           'contact_mean' = 'Mean of Contact Distribution',
           'contact_std' = 'STD of Contact Distribution',
           'percent_HIV' = 'Percent of Contacts HIV+ on Trace',
           'NNTTP' = 'Number Needed to Trace to find 1 Positive Contact',
           'NNTTT' = 'Number Needed to Trace to find 1 Untreated Contact',
           'NNTTD' = 'Number Needed to Trace to find 1 Undiagnosed Contact'
)

output_dir = 'figures'

# Plot population stats
for(pop_stat in c('N_index', 'N_HIV', 'N_untreated', 'N_undiagnosed', 'contact_mean',
                  'contact_std', 'percent_HIV', 'NNTTP', 'NNTTT', 'NNTTD', 'total',
                  'P','N','PP','PN')){
  temp = reduced_data %>% filter(trait=="high_risk")
  fig = ggplot(temp, aes(x=sample_rate, y=get(paste0(pop_stat, '_mean')),
                                 color=tracing_rate,
                                 group=tracing_rate)) +
    theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5)) +
    geom_line( ) +
    geom_errorbar( aes(ymin=get(paste0(pop_stat,'_mean'))-get(paste0(pop_stat,'_std')),
                       ymax=get(paste0(pop_stat,'_mean'))+get(paste0(pop_stat,'_std'))),
                   width=3) +
    facet_grid(comparison ~ tracing_group) +
    guides(alpha='none', color=guide_legend(title="% Contacts Traced")) +
    ggtitle(titles[pop_stat]) +
    xlab('Rate of Enrollment') +
    scale_x_continuous(breaks=unique(temp$sample_rate)) +
    ylab(metrics[pop_stat])
  ggsave(file.path(output_dir, paste0(pop_stat, '_singleTrait.png')), plot=fig, height=8, width=12, units='in', dpi=300)
  # Plot single population
  for(population in unique(temp$comparison)){
    fig = ggplot(temp %>% filter(comparison==population),
                 aes(x=sample_rate, y=get(paste0(pop_stat, '_mean')),
                           color=tracing_rate,
                           group=tracing_rate)) +
      theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5)) +
      geom_line( ) +
      geom_errorbar( aes(ymin=get(paste0(pop_stat,'_mean'))-get(paste0(pop_stat,'_std')),
                         ymax=get(paste0(pop_stat,'_mean'))+get(paste0(pop_stat,'_std'))),
                     width=3) +
      facet_wrap(vars(tracing_group)) +
      guides(alpha='none', color=guide_legend(title="% Contacts Traced")) +
      ggtitle(titles[pop_stat]) +
      xlab('Rate of Enrollment') +
      scale_x_continuous(breaks=unique(temp$sample_rate)) +
      ylab(metrics[pop_stat])
    ggsave(file.path(output_dir, paste0(pop_stat, '_singleTrait_', population, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
  }
}

for(i in seq(1,length(metrics))){
  name = metrics[i]
  metric = names(name)
  title = titles[i]
  # fig = ggplot(reduced_data, aes(x=sample_rate, y=get(paste0(metric, '_mean')), fill=tracing_group)) +
  #   theme(text=element_text(size=16), axis.text.x=element_text(angle=-90),
  #         plot.title = element_text(hjust = 0.5)) +
  #   geom_bar(stat='identity',
  #            color="black", position=position_dodge()) +
  #   geom_errorbar( aes(ymin=get(paste0(metric,'_mean'))-get(paste0(metric,'_std')),
  #                      ymax=get(paste0(metric,'_mean'))+get(paste0(metric,'_std'))),
  #                  width=0.4, colour="black", position=position_dodge(0.9)) +
    fig = ggplot(reduced_data, aes(x=sample_rate, y=get(paste0(metric, '_mean')),
                                   color=tracing_group,
                                   alpha=tracing_rate,
                                   group=interaction(tracing_rate, tracing_group))) +
    theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5)) +
    geom_line( aes(size=tracing_rate) ) +
    geom_errorbar( aes(ymin=get(paste0(metric,'_mean'))-get(paste0(metric,'_std')),
                       ymax=get(paste0(metric,'_mean'))+get(paste0(metric,'_std'))),
                   width=3) +
    facet_grid(comparison ~ trait) +
    guides(alpha='none', color=guide_legend(title="Index Group")) +
    ggtitle(title) +
    scale_size("% Contacts Traced", range = c(0.05, 0.7), breaks=unique(reduced_data$tracing_rate)) +
    xlab('Rate of Enrollment') +
    scale_x_continuous(breaks=unique(reduced_data$sample_rate)) +
    ylab(name)
  ggsave(file.path(output_dir, paste0(metric, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
  for(population in unique(reduced_data$comparison)){
    fig = ggplot(reduced_data %>% filter(comparison==population),
                 aes(x=sample_rate, y=get(paste0(metric, '_mean')),
                     color=tracing_rate,
                     group=tracing_rate)) +
      theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5)) +
      geom_line( ) +
      geom_errorbar( aes(ymin=get(paste0(metric,'_mean'))-get(paste0(metric,'_std')),
                         ymax=get(paste0(metric,'_mean'))+get(paste0(metric,'_std'))),
                     width=3) +
      facet_grid(tracing_group ~ trait) +
      guides(alpha='none', color=guide_legend(title="% Contacts Traced")) +
      ggtitle(paste0(title, '(', population, ')')) +
      xlab('Rate of Enrollment') +
      scale_x_continuous(breaks=unique(reduced_data$sample_rate)) +
      ylab(metrics[pop_stat])
    ggsave(file.path(output_dir, paste0(pop_stat, '_', population, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
  }
  for(SR in unique(reduced_data$sample_rate)){
    for(TR in unique(reduced_data$tracing_rate)){
      fig = ggplot(reduced_data %>% filter(sample_rate==SR & tracing_rate==TR),
                   aes(x=trait, y=get(paste0(metric, '_mean')), fill=tracing_group)) +
        theme(text=element_text(size=16), axis.text.x=element_text(angle=-90),
              plot.title = element_text(hjust = 0.5)) +
        geom_bar(stat='identity',
                 color="black", position=position_dodge()) +
        geom_errorbar( aes(ymin=get(paste0(metric,'_mean'))-get(paste0(metric,'_std')),
                           ymax=get(paste0(metric,'_mean'))+get(paste0(metric,'_std'))),
                       width=0.4, colour="black", position=position_dodge(0.9)) +
        facet_wrap(vars(comparison)) +
        guides(fill=guide_legend(title='Index Group')) +
        ggtitle(paste0(title, '\n(Sample Rate: ', SR, ' Tracing Rate: ', TR, ')')) +
        xlab('Trait') +
        ylab(name)
      ggsave(file.path(output_dir, paste0(metric, '_SR', SR, '_TR', TR, '.png')), plot=fig, height=8, width=10, units='in', dpi=300)
    }
  }
}
