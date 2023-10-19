library('tidyverse')
input_directory = 'data'
reduced_data <- read_delim(file.path(input_directory, 'clusterSummarized.csv'), delim=',')
reduced_data$threshold[reduced_data$threshold==0] = 0.5

metrics = c('percent_clustered' = 'Percent Clustered',
            'num_clusters' = 'Number of Clusters',
            'max_size' = 'Max Cluster Size',
            'median_size' = 'Median Cluster Size',
            'mean_size'  = 'Mean Cluster Size'
)

titles = c('percent_clustered' = 'Percent of Sampled Individuals in Clusters',
           'num_clusters' = 'Number of Clusters',
           'max_size' = 'Maximum Cluster Size',
           'median_size' = 'Median Cluster Size',
           'mean_size'  = 'Mean Cluster Size'
)

output_dir = 'figures'

for(i in seq(1,length(metrics))){
  name = metrics[i]
  metric = names(name)
  title = titles[i]
  
  fig = ggplot(reduced_data,
               aes(x=threshold, y=get(paste0(metric, '_mean')),
                   color=tracing_group,
                   group=tracing_group)) +
    theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(angle=-90)) +
    geom_line( ) +
    geom_errorbar( aes(ymin=get(paste0(metric,'_mean'))-get(paste0(metric,'_std')),
                       ymax=get(paste0(metric,'_mean'))+get(paste0(metric,'_std'))),
                   width=1) +
    facet_grid(n_samples ~ tracing_start_time) +
    guides(color=guide_legend(title="Index Group")) +
    ggtitle(title) +
    xlab('Cluster Threshold (Years)') +
    ylab(name)

  ggsave(file.path(output_dir, paste0(metric, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
}

reduced_data <- read_delim(file.path(input_directory, 'enrichmentSummarized.csv'), delim=',')
reduced_data$threshold[reduced_data$threshold==0] = 0.5

metrics = c('FP' = 'Clustered without Trait',#'False Positives',
            'TN' = 'Non-Clustered without Trait',#'True Negatives',
            'TP' = 'Clustered with Trait',#'True Postitives',
            'FN' = 'Non-Clustered with Trait',#'False Negatives',
            'P'  = 'People with Trait',#'Actual Positive',
            'N'  = 'People without Trait',#'Actual Negative',
            'PP' = 'Clustered',#'Predicted Positive',
            'PN' = 'Non-Clustered',#'Predicted Negative',
            'total' = 'Total Population',
            'prevalence' = 'Prevalence of Trait',
            'ACC' = 'Accuracy',
            'PPV' = 'Proportion of Clustered with Trait',#'Positive Predictive Value',
            'FDR' = 'False Discovery Rate',
            'FOR' = 'False Omission Rate',
            'NPV' = 'Negative Predictive Value',
            'TPR' = 'Proportion of Trait Amoung Clustered',#'True Positive Rate',
            'FPR' = 'Proportion of Non-Trait Amoung Clustered',#'False Positive Rate',
            'FNR' = 'Proportion of Trait Amoung Non-Clustered',#'False Negative Rate',
            'TNR' = 'Proportion of Non-Trait Amoung Non-Clustered',#'True Negative Rate',
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
            'RRrandom' = 'Relative Risk vs Random Sampling'
)

titles = c('FP' = 'Number of Clustered Without Trait',
           'TN' = 'Number of Non-Clustered Without Trait',
           'TP' = 'Number of Clustered With Trait',
           'FN' = 'Number of Non-Clustered With Trait',
           'P'  = 'Number of People With Trait',
           'N'  = 'Number of People Without Trait',
           'PP' = 'Number of Clustered of Tracing Group',#'Number of Clustered of Tracing Group',
           'PN' = 'Number of Non-Clustered of Tracing Group',
           'total' = 'Total Population',
           'prevalence' = 'Prevalence of Trait',
           'ACC' = 'C w/ Trait plus NC w/o Trait divide by Total Population',
           'PPV' = 'Clustered With Trait / All Clustered',
           'FDR' = 'Clustered Without Trait / All Clustered',
           'FOR' = 'Non-Clustered With Trait / All Non-Clustered',
           'NPV' = 'Non-Clustered Without Trait / All Non-Clustered',
           'TPR' = 'Clustered With Trait / # People With Trait',
           'FPR' = 'Clustered Without Trait / # People Without Trait',
           'FNR' = 'Non-Clustered With Trait / # People With Trait',
           'TNR' = 'Non-Clustered Without Trait / # People Without Trait',
           'LRP' = 'True Postive Rate / False Postive Rate',
           'LRN' = 'False Negative Rate / True Negative Rate',
           'DOR' = 'Positive Likelihood Ration / Negative Likelihood Ratio',
           'TS' = 'Clustered With Trait / (All Clustered + Non-Clustered With Trait)',
           'informedness' = 'True Positive Rate + True Negative Rate - 1',
           'markedness' = 'Positive Predictive Value + Negatiev Predictive Value - 1',
           'prevalenceThreshold' = 'Prevalence Threshold',
           'MCC' = '~ Correlation Between Being a Contact and Having Trait',
           'FMI' = 'Fowlkes-Mallows',
           'F1' = 'F1 Score',
           'BA' = 'Mean of True Postive Rate and True Negative Rate',
           'RR' = 'Relative Risk',
           'RRall' = 'Relative Risk vs All',
           'RRrandom' = 'Relative Risk vs Random Sampling'
)

for(i in seq(1,length(metrics))){
  name = metrics[i]
  metric = names(name)
  title = titles[i]
  for(cur_trait in unique(reduced_data$trait)){
    fig = ggplot(reduced_data %>% filter(trait==cur_trait),
                 aes(x=threshold, y=get(paste0(metric, '_mean')),
                     color=tracing_group,
                     group=tracing_group)) +
      theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5),
            axis.text.x=element_text(angle=-90)) +
      geom_line( ) +
      geom_errorbar( aes(ymin=get(paste0(metric,'_mean'))-get(paste0(metric,'_std')),
                         ymax=get(paste0(metric,'_mean'))+get(paste0(metric,'_std'))),
                     width=1) +
      facet_grid(n_samples ~ tracing_start_time) +
      guides(color=guide_legend(title="Index Group")) +
      ggtitle(paste0(title, ' (', cur_trait, ')')) +
      xlab('Cluster Threshold (Years)') +
      ylab(name)
    
    ggsave(file.path(output_dir, paste0(metric, '_', cur_trait, '.png')), plot=fig, height=8, width=12, units='in', dpi=300)
  }
}
