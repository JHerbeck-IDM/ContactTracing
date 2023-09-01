library('tidyverse')

getGroups <- function(trace_start_time, trace_end_time, look_back_length){
  # Read in file
  eventReport <- read_csv('ReportEventRecorder.csv')
  
  # Set a dataframe with a row per individual
  individuals <- unique(eventReport['Individual_ID'])
  
  # Get the birth time for each individual
  births <- eventReport %>% filter(Event_Name == 'Births') %>% select(Individual_ID, Year)
  individuals <- left_join(individuals, births, by='Individual_ID')
  individuals <- rename(individuals, Birth = Year)
  
  # Get the death time for each individual
  deathEvents = c('DiseaseDeaths',
                  'NonDiseaseDeaths',
                  'OpportunisticInfectionDeath')
  deaths <- eventReport %>% filter(Event_Name %in% deathEvents) %>% select(Individual_ID, Year)
  individuals <- left_join(individuals, deaths, by='Individual_ID')
  individuals <- rename(individuals, Death = Year)
  
  # Get the diagnosis time for each individual
  careCasadeEvents = c('ARTStaging0',
                       'ARTStaging1',
                       'ARTStaging2',
                       'ARTStaging3',
                       'ARTStaging4',
                       'ARTStaging5',
                       'ARTStaging6',
                       'ARTStaging8',
                       'ARTStaging9',
                       'HIV_Positive_at_ANC',
                       'NewInfectionEvent',
                       'OnART0',
                       'OnPreART0',
                       'OnPreART4')
  diagnosis <- eventReport %>% filter(Event_Name %in% careCasadeEvents) %>%
    group_by(Individual_ID) %>%
    summarize(Diagnosis = min(Year[Year > Year[Event_Name=='NewInfectionEvent']]),
              Stage = WHO_Stage[min(which(Year == min(Year[Year > Year[Event_Name=='NewInfectionEvent']])))]) %>%
    mutate(Diagnosis = ifelse(is.infinite(Diagnosis), NA, Diagnosis))
  individuals <- left_join(individuals, diagnosis, by='Individual_ID')
  
  # Reduce to only individuals that are alive during the tracing period
  individuals <- individuals %>% filter(Death < trace_start_time | 
                                          Birth > trace_end_time)
  
  # Create the three subgroups
  individuals['subset'] = 'HIV-'
  individuals$subset[individuals$Diagnosis>=trace_start_time] = 'Incident_HIV'
  individuals$subset[individuals$Diagnosis<trace_start_time] = 'Prevalent_HIV'
  individuals$subset <- factor(individuals$subset, levels=c('Incident_HIV',
                                                            'Prevalent_HIV',
                                                            'HIV-'))
  
  return(individuals)
}