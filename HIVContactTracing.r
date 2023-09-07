library('tidyverse')

getGroups <- function(trace_start_time, trace_end_time, acute_to_trace=TRUE){
  ### This function is focused on getting the population that is alive during
  ### the tracing period as well as information about them such as date infected,
  ### date diagnosed, etc. If these things did not happen for a specific
  ### individual the dates will be NA. Lastly, this fuction flags four groups
  ### within in the population: HIV-, HIV+ but undiagnosed, HIV+ prevalence,
  ### HIV+ incidence.
  ###
  ### Inputs:
  ###    trace_start_time (float): The date when the period for contract tracing
  ###                              started.
  ###    trace_end_time (float)  : The date when the period for contract tracing
  ###                              ended.
  ###    acute_to_trace (bool)   : This flags whether the incidence infections
  ###                              (those whose contacts will be traced) are
  ###                              just those diagnosed during the tracing period
  ###                              or if they also need to be acute at the time
  ###                              of diagnosis.
  ### Return:
  ###    data frame              : Data frame with the aggregated data.
  
  # Read in file
  eventReport <- read_csv('ReportEventRecorder.csv')
  
  # Set a data frame with a row per individual
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
    summarize(Infected = min(Year[Event_Name=='NewInfectionEvent']),
              Diagnosis = min(Year[Year > Year[Event_Name=='NewInfectionEvent']]),
              Stage = WHO_Stage[min(which(Year == min(Year[Year > Year[Event_Name=='NewInfectionEvent']])))]) %>%
    mutate(Diagnosis = ifelse(is.infinite(Diagnosis), NA, Diagnosis),
           Infected = ifelse(is.infinite(Infected), NA, Infected))
  individuals <- left_join(individuals, diagnosis, by='Individual_ID')
  
  # Get time no longer acute
  endAcute <- eventReport %>% filter(Event_Name %in% c('HIVInfectionStageEnteredLatent',
                                                       'HIVInfectionStageEnteredAIDS',
                                                       'HIVInfectionStageEnteredOnART')) %>%
    group_by(Individual_ID) %>%
    summarize(EndAcute = min(Year)) %>%
    mutate(EndAcute = ifelse(is.infinite(EndAcute), NA, EndAcute))
  individuals <- left_join(individuals, endAcute, by='Individual_ID')
  
  # Reduce to only individuals that are alive for any amount of time during the
  # tracing period
  individuals <- individuals %>% filter(Death > trace_start_time & 
                                          Birth < trace_end_time)
  
  # Create the three subgroups
  individuals['subset'] = 'HIV-'
  # Set everyone infected by the end of the trace period to the prevalent group
  individuals$subset[individuals$Infected<trace_end_time] = 'Prevalent_HIV'
  # In the wiki post incident and prevalent are only diagnosed individuals
  # Here moved the undiagnosed to a separate category commenting this out will
  # leave them as prevalent
  individuals$subset[is.na(individuals$Diagnosis) &
                       individuals$Infected<trace_end_time] = 'Undiagnosed_HIV'
  if(acute_to_trace){
    # Only trace (incident) those that are recently diagnosed and acute
    individuals$subset[individuals$Diagnosis>=trace_start_time &
                       individuals$Infected<trace_end_time &
                       individuals$EndAcute > trace_start_time &
                       individuals$Diagnosis <= individuals$EndAcute] = 'Incident_HIV'
  } else{
    # Trace (incident) all recent diagnosed
    individuals$subset[individuals$Diagnosis>=trace_start_time &
                       individuals$Infected<trace_end_time] = 'Incident_HIV'
  }
  individuals$subset <- factor(individuals$subset, levels=c('Incident_HIV',
                                                            'Prevalent_HIV',
                                                            'Undiagnosed_HIV',
                                                            'HIV-'))
  
  return(individuals)
}