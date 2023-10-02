print('Sourcing tracing functions')
source('HIVContactTracing.r')

directory = readLines("data_path.txt")
tracing_groups=c('HIV+', 'Incident_HIV', 'HIV-', 'Prevalent_HIV',
                 'Undiagnosed_HIV')
comparisons = c('all','traceable','sampleableAtTrace','sampleableEventually')
sample_rates=c(100, 75, 50, 25)
tracing_rates=c(100, 75, 50, 25)
tracing_start_times=2010
tracing_end_times=2012
acute_to_traces=TRUE
look_back_durations=3
tracing_delays=0
traits = c('trans_source', 'yet_to_transmit', 'recently_infected', 'high_risk',
           'not_yet_infected')
unique_individuals = c(FALSE)

agg_data =data.frame()

results = doTracing(directory=directory, tracing_groups,
                    sample_rates, tracing_rates,
                    tracing_start_times, tracing_end_times,
                    acute_to_traces, look_back_durations,
                    tracing_delays, unique_individuals)

# write_csv2(results, 'output.csv')
write.table(results, 'output.csv', sep=',', dec='.', row.names=FALSE)