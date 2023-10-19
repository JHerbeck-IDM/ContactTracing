library('tidyverse')
input_directory = 'data'
reduced_data <- read_delim(file.path(input_directory, 'output_summarized.csv'), delim=',')

metrics = c('FP' = 'Contacts without Trait',#'False Positives',
            'TN' = 'Non-Contacts without Trait',#'True Negatives',
            'TP' = 'Contacts with Trait',#'True Postitives',
            'FN' = 'Non-Contacts with Trait',#'False Negatives',
            'P'  = 'People with Trait',#'Actual Positive',
            'N'  = 'People without Trait',#'Actual Negative',
            'PP' = 'Sexual Contacts',#'Predicted Positive',
            'PN' = 'Non-Contacts',#'Predicted Negative',
            'total' = 'Total Population',
            'prevalence' = 'Prevalence of Trait',
            'ACC' = 'Accuracy',
            'PPV' = 'Proportion of Contacts with Trait',#'Positive Predictive Value',
            'FDR' = 'False Discovery Rate',
            'FOR' = 'False Omission Rate',
            'NPV' = 'Negative Predictive Value',
            'TPR' = 'Proportion of Trait Amoung Contacts',#'True Positive Rate',
            'FPR' = 'Proportion of Non-Trait Amoung Contacts',#'False Positive Rate',
            'FNR' = 'Proportion of Trait Amoung Non-Contacts',#'False Negative Rate',
            'TNR' = 'Proportion of Non-Trait Amoung Non-Contacts',#'True Negative Rate',
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
            'RRHIV-' = 'Relative Risk vs HIV- Contacts',
            'contacts_per_index' = 'Contacts Per Index',
            'N_HIV' = 'Number HIV+',
            'N_index' = 'Number of Index cases enrolled',
            'N_untreated' = 'Number Untreated',
            'N_undiagnosed' = 'Number Undiagnosed',
            'contact_mean' = 'Number of Contacts (Mean)',
            'contact_std' = 'Number of Contacts (Std Dev)',
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
           'PP' = 'Number of Sexual Contacts of Tracing Group',#'Number of Contacts of Tracing Group',
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
           'RRHIV-' = 'Relative Risk vs HIV- Contacts',
           'contacts_per_index' = 'Contacts Per Index',
           'N_HIV' = 'Number of Contacts HIV+ on Trace',
           'N_index' = 'Number of Index cases enrolled',
           'N_untreated' = 'Number of Contacts Untreated on Trace',
           'N_undiagnosed' = 'Number of Contacts Undiagnosed on Trace',
           'contact_mean' = 'Mean of Contact Distribution',
           'contact_std' = 'Std Dev of Contact Distribution',
           'percent_HIV' = 'Percent of Contacts HIV+ on Trace',
           'NNTTP' = 'Number Needed to Trace to find 1 Positive Contact',
           'NNTTT' = 'Number Needed to Trace to find 1 Untreated Contact',
           'NNTTD' = 'Number Needed to Trace to find 1 Undiagnosed Contact'
)

output_dir = 'figures'

## Make plots with single generation
# Plot population stats
temp = reduced_data %>% filter(trait=="high_risk" & generation==1)
for(pop_stat in c('N_index', 'N_HIV', 'N_untreated', 'N_undiagnosed', 'contact_mean',
                  'contact_std', 'percent_HIV', 'NNTTP', 'NNTTT', 'NNTTD', 'total',
                  'P','N','PP','PN')){
  for(duplicate in c('Duplicates','No_Duplicates')){
    fig = ggplot(temp %>% filter(duplicates==duplicate), aes(x=sample_rate, y=get(paste0(pop_stat, '_mean')),
                           color=tracing_rate,
                           group=tracing_rate)) +
      theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5)) +
      geom_line( ) +
      geom_errorbar( aes(ymin=get(paste0(pop_stat,'_mean'))-get(paste0(pop_stat,'_std')),
                         ymax=get(paste0(pop_stat,'_mean'))+get(paste0(pop_stat,'_std'))),
                     width=3) +
      facet_grid(comparison ~ tracing_group) +
      scale_color_continuous(breaks=unique(temp$tracing_rate)) +
      guides(alpha='none', color=guide_legend(title="% Contacts Traced")) +
      ggtitle(titles[pop_stat]) +
      xlab('Rate of Enrollment') +
      scale_x_continuous(breaks=unique(temp$sample_rate)) +
      ylab(metrics[pop_stat])
    ggsave(file.path(output_dir, paste0(pop_stat, '_singleTrait_gen1_',  duplicate, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
  }
  # Plot single population
  for(population in unique(temp$comparison)){
    fig = ggplot(temp %>% filter(comparison==population),
                 aes(x=tracing_group, y=get(paste0(pop_stat, '_mean')),
                     color=duplicates,
                     group=duplicates)) +
      theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5),
            axis.text.x=element_text(angle=-90)) +
      geom_line( ) +
      geom_errorbar( aes(ymin=get(paste0(pop_stat,'_mean'))-get(paste0(pop_stat,'_std')),
                         ymax=get(paste0(pop_stat,'_mean'))+get(paste0(pop_stat,'_std'))),
                     width=3) +
      # facet_wrap(vars(tracing_group)) +
      facet_grid(sample_rate ~ tracing_rate, labeller = labeller(.rows = label_both, .cols = label_both)) +
      # scale_color_continuous(breaks=unique(temp$tracing_rate)) +
      guides(alpha='none', color=guide_legend(title="Include Duplicates")) +
      ggtitle(titles[pop_stat]) +
      xlab('Rate of Enrollment') +
      # scale_x_continuous(breaks=unique(temp$sample_rate)) +
      ylab(metrics[pop_stat])
    ggsave(file.path(output_dir, paste0(pop_stat, '_singleTrait_', population, '_gen1.png')), plot=fig, height=10, width=12, units='in', dpi=300)
  }
}

# plot metrics
temp = reduced_data %>% filter(generation==1)
for(i in seq(1,length(metrics))){
  name = metrics[i]
  metric = names(name)
  title = titles[i]
  for(duplicate in c('Duplicates','No_Duplicates')){
    fig = ggplot(temp %>% filter(duplicates==duplicate),
                 aes(x=sample_rate, y=get(paste0(metric, '_mean')),
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
      scale_size("% Contacts Traced", range = c(0.05, 0.7), breaks=unique(temp$tracing_rate)) +
      xlab('Rate of Enrollment') +
      scale_x_continuous(breaks=unique(temp$sample_rate)) +
      ylab(name)
    ggsave(file.path(output_dir, paste0(metric, '_gen1_', duplicate, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
    for(population in unique(temp$comparison)){
      fig = ggplot(temp %>% filter(comparison==population & duplicates==duplicate),
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
        scale_x_continuous(breaks=unique(temp$sample_rate)) +
        ylab(name)
      ggsave(file.path(output_dir, paste0(metric, '_', population, '_gen1_', duplicate, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
    }
  }
  SR=100
  TR=100
  fig = ggplot(temp %>% filter(sample_rate==SR & tracing_rate==TR),
               aes(x=trait, y=get(paste0(metric, '_mean')), fill=duplicates)) +
    theme(text=element_text(size=16), axis.text.x=element_text(angle=-90),
          plot.title = element_text(hjust = 0.5)) +
    geom_bar(stat='identity',
             color="black", position=position_dodge()) +
    geom_errorbar( aes(ymin=get(paste0(metric,'_mean'))-get(paste0(metric,'_std')),
                       ymax=get(paste0(metric,'_mean'))+get(paste0(metric,'_std'))),
                   width=0.4, colour="black", position=position_dodge(0.9)) +
    # facet_wrap(vars(comparison), nrow = 2) +
    facet_grid(comparison ~ tracing_group) +
    guides(fill=guide_legend(title='Index Group')) +
    ggtitle(paste0(title, '\n(Sample Rate: ', SR, ' Tracing Rate: ', TR, ')')) +
    xlab('Trait') +
    ylab(name)
  ggsave(file.path(output_dir, paste0(metric, '_SR', SR, '_TR', TR, '_gen1.png')), plot=fig, height=8, width=10, units='in', dpi=300)
  fig = ggplot(temp %>% filter(sample_rate==SR & tracing_rate==TR & tracing_group=='HIV+'),
               aes(x=trait, y=get(paste0(metric, '_mean')), fill=duplicates)) +
    theme(text=element_text(size=16), axis.text.x=element_text(angle=-90),
          plot.title = element_text(hjust = 0.5)) +
    geom_bar(stat='identity',
             color="black", position=position_dodge()) +
    geom_errorbar( aes(ymin=get(paste0(metric,'_mean'))-get(paste0(metric,'_std')),
                       ymax=get(paste0(metric,'_mean'))+get(paste0(metric,'_std'))),
                   width=0.4, colour="black", position=position_dodge(0.9)) +
    facet_wrap(vars(comparison), nrow = 2) +
    guides(fill=guide_legend(title='Population')) +
    ggtitle(paste0(title, '\n(Sample Rate: ', SR, ' Tracing Rate: ', TR, ')')) +
    xlab('Trait') +
    ylab(name)
  ggsave(file.path(output_dir, paste0(metric, '_singleGroup_SR', SR, '_TR', TR, '_gen1.png')), plot=fig, height=8, width=10, units='in', dpi=300)
}

## Plot results from interative tracing
# Plot population stats (trait independent)
temp = reduced_data %>% filter(trait=="high_risk" & comparison=="traceable")
for(pop_stat in c('N_index', 'N_HIV', 'N_untreated', 'N_undiagnosed', 'contact_mean',
                  'contact_std', 'percent_HIV', 'NNTTP', 'NNTTT', 'NNTTD', 'total',
                  'P','N','PP','PN')){
  for(duplicate in c('Duplicates','No_Duplicates')){
    fig = ggplot(temp %>% filter(duplicates==duplicate),
                 aes(x=generation, y=get(paste0(pop_stat, '_mean')),
                     color=tracing_group,
                     group=tracing_group)) +
      theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5),
            axis.text.x=element_text(angle=-90)) +
      geom_line( ) +
      geom_errorbar( aes(ymin=get(paste0(pop_stat,'_mean'))-get(paste0(pop_stat,'_std')),
                         ymax=get(paste0(pop_stat,'_mean'))+get(paste0(pop_stat,'_std'))),
                     width=1) +
      facet_grid(sample_rate ~ tracing_rate, labeller = labeller(.rows = label_both, .cols = label_both)) +
      guides(color=guide_legend(title="Index Group")) +
      ggtitle(paste0(titles[pop_stat], ' (Traceable)')) +
                xlab('Generation of Tracing') +
                # scale_x_continuous(breaks=unique(temp$sample_rate)) +
                ylab(metrics[pop_stat])
    ggsave(file.path(output_dir, paste0(pop_stat, '_singleTrait_traceable_', duplicate, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
  }
}

for(i in seq(1,length(metrics))){
  name = metrics[i]
  metric = names(name)
  title = titles[i]
  for(cur_trait in unique(reduced_data$trait)){
    for(duplicate in c('Duplicates','No_Duplicates')){
      fig = ggplot(reduced_data %>% filter(comparison=='traceable' & trait==cur_trait & duplicates==duplicate),
                   aes(x=generation, y=get(paste0(metric, '_mean')),
                       color=tracing_group,
                       group=tracing_group)) +
        theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5),
              axis.text.x=element_text(angle=-90)) +
        geom_line( ) +
        geom_errorbar( aes(ymin=get(paste0(metric,'_mean'))-get(paste0(metric,'_std')),
                           ymax=get(paste0(metric,'_mean'))+get(paste0(metric,'_std'))),
                       width=1) +
        facet_grid(sample_rate ~ tracing_rate, labeller = labeller(.rows = label_both, .cols = label_both)) +
        guides(color=guide_legend(title="Index Group")) +
        ggtitle(paste0(title, ' (', cur_trait, ', ', 'Traceable)')) +
        xlab('Generation of Tracing') +
        ylab(name)
      max_y = fig$layout$panel_params[[1]]$y.range[2]
      fig + ylim(0, min(50, max_y))
      ggsave(file.path(output_dir, paste0(metric, '_', cur_trait, '_traceable_', duplicate, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
      if(cur_trait %in% c('high_risk', 'trans_source') & metric=='PPV'){
        fig = ggplot(reduced_data %>% filter(comparison=='traceable' & trait==cur_trait & tracing_group %in% c('HIV-', 'HIV+') & duplicates==duplicate),
                     aes(x=generation, y=get(paste0(metric, '_mean')),
                         color=tracing_group,
                         group=tracing_group)) +
          theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5),
                axis.text.x=element_text(angle=-90)) +
          geom_line( ) +
          geom_errorbar( aes(ymin=get(paste0(metric,'_mean'))-get(paste0(metric,'_std')),
                             ymax=get(paste0(metric,'_mean'))+get(paste0(metric,'_std'))),
                         width=1) +
          facet_grid(sample_rate ~ tracing_rate, labeller = labeller(.rows = label_both, .cols = label_both)) +
          guides(color=guide_legend(title="Index Group")) +
          ggtitle(paste0(title, ' (', cur_trait, ', ', 'Traceable)')) +
          xlab('Generation of Tracing') +
          ylab(name)
        max_y = fig$layout$panel_params[[1]]$y.range[2]
        fig + ylim(0, min(50, max_y))
        ggsave(file.path(output_dir, paste0(metric, '_', cur_trait, '_HIV_traceable_', duplicate, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
      }
    }
  }
}

for(duplicate in c('Duplicates','No_Duplicates')){
  fig = ggplot(reduced_data %>% filter(comparison=='traceable' & generation==1 & duplicates==duplicate),
               aes(x=trait, y=`RRHIV-_mean`,
                   fill=tracing_group,
                   group=tracing_group)) +
    theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(angle=-90)) +
    geom_bar(stat='identity',
             color="black", position=position_dodge()) +
    geom_errorbar( aes(ymin=`RRHIV-_mean`-`RRHIV-_std`,
                       ymax=`RRHIV-_mean`+`RRHIV-_std`),
                   width=0.4, colour="black", position=position_dodge(0.9)) +
    facet_grid(sample_rate ~ tracing_rate, labeller = labeller(.rows = label_both, .cols = label_both)) +
    guides(color=guide_legend(title="Index Group")) +
    ggtitle(paste0(titles['RRHIV-'], ' (Traceable)')) +
    xlab('Trait') +
    ylab(metrics['RRHIV-'])
  ggsave(file.path(output_dir, paste0('RRHIV-_traceable_', duplicate, '.png')), plot=fig, height=10, width=12, units='in', dpi=300)
}
