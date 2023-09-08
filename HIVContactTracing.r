library('tidyverse')

getGroups <- function(trace_start_time, trace_end_time, acute_to_trace=FALSE){
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
  sim_start_year = min(eventReport['Year'])
  
  # Set a data frame with a row per individual
  individuals <- eventReport %>% select(Individual_ID, Gender) %>%
               distinct(Individual_ID, .keep_all=TRUE)
  
  # Get the birth time for each individual
  births <- eventReport %>% filter(Event_Name == 'Births') %>%
    select(Individual_ID, Year)
  individuals <- left_join(individuals, births, by='Individual_ID')
  individuals <- rename(individuals, Birth = Year)
  rm(births)
  
  # Get the death time for each individual
  deathEvents = c('DiseaseDeaths',
                  'NonDiseaseDeaths',
                  'OpportunisticInfectionDeath')
  deaths <- eventReport %>% filter(Event_Name %in% deathEvents) %>% select(Individual_ID, Year)
  individuals <- left_join(individuals, deaths, by='Individual_ID')
  individuals <- rename(individuals, Death = Year)
  rm(deaths)
  
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
  rm(diagnosis)
  
  # Get time no longer acute
  endAcute <- eventReport %>% filter(Event_Name %in% c('HIVInfectionStageEnteredLatent',
                                                       'HIVInfectionStageEnteredAIDS',
                                                       'HIVInfectionStageEnteredOnART')) %>%
    group_by(Individual_ID) %>%
    summarize(EndAcute = min(Year)) %>%
    mutate(EndAcute = ifelse(is.infinite(EndAcute), NA, EndAcute))
  individuals <- left_join(individuals, endAcute, by='Individual_ID')
  rm(endAcute)
  rm(eventReport)
  
  # Reduce to only individuals that are alive for any amount of time during the
  # tracing period
  individuals <- individuals %>% filter(Death > trace_start_time & 
                                          Birth < trace_end_time)
  
  # Create the three subgroups
  individuals['subset'] = 'HIV-'
  # Set everyone infected by the end of the trace period to the prevalent group
  individuals$subset[individuals$Infected<=trace_end_time] = 'Prevalent_HIV'
  # In the wiki post incident and prevalent are only diagnosed individuals
  # Here moved the undiagnosed to a separate category commenting this out will
  # leave them as prevalent
  individuals$subset[individuals$Infected<=trace_end_time &
                       (is.na(individuals$Diagnosis) |
                          individuals$Diagnosis>=trace_end_time)
                       ] = 'Undiagnosed_HIV'
  if(acute_to_trace){
    # Only trace (incident) those that are recently diagnosed and acute
    individuals$subset[individuals$Diagnosis>=trace_start_time &
                       individuals$Diagnosis<=trace_end_time &
                       individuals$Diagnosis <= individuals$EndAcute] = 'Incident_HIV'
  } else{
    # Trace (incident) all recent diagnosed
    individuals$subset[individuals$Diagnosis>=trace_start_time &
                       individuals$Diagnosis<=trace_end_time] = 'Incident_HIV'
  }
  individuals$subset <- factor(individuals$subset, levels=c('Incident_HIV',
                                                            'Prevalent_HIV',
                                                            'Undiagnosed_HIV',
                                                            'HIV-'))
  
  return(list('df'=individuals, 'start_year'=sim_start_year))
}

getPartners <- function(index_individuals, look_back_duration, sim_start_year=1960.5){
  
  # Read in relationship start file and reduce the amount of data
  relationshipData <- read_csv('RelationshipStart.csv')
  relationshipData <- relationshipData %>% select(Rel_ID,
                                                  Rel_start_time,
                                                  Rel_scheduled_end_time,
                                                  A_ID,
                                                  B_ID,
                                                  A_gender,
                                                  B_gender,
                                                  A_IndividualProperties,
                                                  B_IndividualProperties,
                                                  A_HIV_disease_stage,
                                                  B_HIV_disease_stage)
  
  # Read in relationship end file to get the true end of each relationship
  relationshipEnd <- read_csv('RelationshipEnd.csv')
  relationshipEnd <- relationshipEnd %>% select(Rel_ID,
                                                Rel_actual_end_time)
  relationshipData <- left_join(relationshipData, relationshipEnd, by='Rel_ID')
  rm(relationshipEnd)
  
  # Extract the risk state from the IP column
  relationshipData$A_Risk <- relationshipData$A_IndividualProperties %>%
    str_extract(regex('LOW|MEDIUM|HIGH'))
  relationshipData$B_Risk <- relationshipData$B_IndividualProperties %>%
    str_extract(regex('LOW|MEDIUM|HIGH'))
  relationshipData <- relationshipData %>%
    select(!c(A_IndividualProperties, B_IndividualProperties))
  
  # Get relationships with male index individuals
  maleIndex <- relationshipData %>%
    filter(A_ID %in% index_individuals$Individual_ID) %>% 
    select(A_ID, B_ID, Rel_ID, B_Risk, B_HIV_disease_stage,
           Rel_start_time, Rel_actual_end_time)
  maleIndex <- rename(maleIndex, Individual_ID=A_ID,
                      contact_ID=B_ID,
                      contact_Risk=B_Risk,
                      contact_HIV_stage=B_HIV_disease_stage)
  maleIndex <- left_join(maleIndex,
                         individuals[c('Individual_ID', 'Diagnosis')],
                         by='Individual_ID')
  
  # Get relationships with female index individuals
  femaleIndex <- relationshipData %>%
    filter(B_ID %in% index_individuals$Individual_ID) %>% 
    select(B_ID, A_ID, Rel_ID, A_Risk, A_HIV_disease_stage,
           Rel_start_time, Rel_actual_end_time)
  femaleIndex <- rename(femaleIndex, Individual_ID=B_ID,
                        contact_ID=A_ID,
                        contact_Risk=A_Risk,
                        contact_HIV_stage=A_HIV_disease_stage)
  femaleIndex <- left_join(femaleIndex,
                           individuals[c('Individual_ID', 'Diagnosis', 'subset')],
                           by='Individual_ID')
  rm(relationshipData)
  
  # Combined into a single data frame
  contactData <- bind_rows(maleIndex, femaleIndex)
  rm(maleIndex)
  rm(femaleIndex)
  
  # Reduce to relationships that were active in the look back period
  contactData$Rel_start_time <- contactData$Rel_start_time / 365 + sim_start_year
  contactData$Rel_actual_end_time <- contactData$Rel_actual_end_time / 365 + sim_start_year
  contactData <- contactData %>% filter(Rel_start_time <= Diagnosis &
                                        Rel_actual_end_time >= (Diagnosis - look_back_duration))
  
  return(contactData)
}

getTransData <- function(contactData){
  
  # Read in transmission file and reduce the amount of data
  transData <- read_csv('TransmissionReport.csv')
  transData <- transData %>% select(YEAR,
                                    REL_ID,
                                    SRC_ID,
                                    DEST_ID)
  
  # Get number of transmissions of contacts before and after tracing start
  transData <- left_join(transData, contactData[c('contact_ID',
                                                  'Individual_ID',
                                                  'Diagnosis',
                                                  'contact_Risk',
                                                  'contact_HIV_stage',
                                                  'subset')],
                         by=c('SRC_ID'='contact_ID'), relationship="many-to-many")
  
  transData <- transData %>% group_by(SRC_ID, Individual_ID) %>% 
    summarise(num_trans_before=sum(YEAR<=Diagnosis),
              num_trans_after=sum(YEAR>Diagnosis),
              contact_Risk=contact_Risk[1],
              contact_HIV_stage=contact_HIV_stage[1],
              subset=subset[1]) %>% 
    mutate(num_trans_before = ifelse(is.na(num_trans_before), 0, num_trans_before),
           num_trans_after = ifelse(is.na(num_trans_after), 0, num_trans_after)) %>% 
  
  return(transData)
}
