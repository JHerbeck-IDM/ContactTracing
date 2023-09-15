library('tidyverse')

getDeaths <- function(eventReport){
  ### Get the death time for each individual for the event reporter
  ###
  ### Inputs:
  ###   eventReport (data frame) : The event reporter output
  ### Outputs:
  ###   data frame               : Reduced dataframe containing just individual
  ###                              IDs and the time that they died of any cause
  
  # There are three different death events
  deathEvents = c('DiseaseDeaths',
                  'NonDiseaseDeaths',
                  'OpportunisticInfectionDeath')
  # Filter out the unneeded columns and events
  deaths <- eventReport %>% filter(Event_Name %in% deathEvents) %>% select(Individual_ID, Year)
  # Rename the year column since it is now the time of death
  deaths <- rename(deaths, Death = Year)
  
  return(deaths)
}

getDiagnosis <- function(eventReport){
  ### Get the time of diagnosis for each individual for the event reporter
  ###
  ### Inputs:
  ###   eventReport (data frame) : The event reporter output
  ### Outputs:
  ###   data frame               : Reduced dataframe containing just individual
  ###                              IDs and the time of infection and time of 
  ###                              diagnosis
  
  # The (hopefully) exhaustive list of care cascade events
  careCascadeEvents = c('ARTStaging0',
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
  # Remove events not in the care cascade, group by individual, find the time of
  # infection, and the earliest care event after (or at the same time as) the infection
  diagnosis <- eventReport %>% filter(Event_Name %in% careCascadeEvents) %>%
    group_by(Individual_ID) %>%
    summarize(Infected = min(Year[Event_Name=='NewInfectionEvent']),
              Diagnosis = min(Year[Year >= Year[Event_Name=='NewInfectionEvent'] &
                                     Event_Name!='NewInfectionEvent'])) %>%
    mutate(Diagnosis = ifelse(is.infinite(Diagnosis), NA, Diagnosis),
           Infected = ifelse(is.infinite(Infected), NA, Infected))
  
  return(diagnosis)
}

getGroups <- function(directory, trace_start_time, trace_end_time, acute_to_trace=FALSE){
  ### This function is focused on getting the population that is alive during
  ### the tracing period as well as information about them such as date infected,
  ### date diagnosed, etc. If these things did not happen for a specific
  ### individual the dates will be NA. Lastly, this function flags four groups
  ### within in the population: HIV-, HIV+ but undiagnosed, HIV+ prevalence,
  ### HIV+ incidence.
  ###
  ### Inputs:
  ###    directory (str)         : Path to directory where input files can be
  ###                              found
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
  ###    float                   : The time of the earliest event in the
  ###                              simulation as a proxy for the simulation
  ###                              start time
  
  # Read in file
  eventReport <- read_csv(file.path(directory, 'ReportEventRecorder.csv'))
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
  deaths = getDeaths(eventReport)
  individuals <- left_join(individuals, deaths, by='Individual_ID')
  rm(deaths)
  
  # Get the diagnosis time for each individual
  diagnosis = getDiagnosis(eventReport)
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
  individuals <- individuals %>% filter(Death >= trace_start_time & 
                                          Birth <= trace_end_time)
  
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
  
  # Make sure each person has a time that their contacts are traced
  individuals['begin_tracing'] = trace_start_time
  individuals$begin_tracing[individuals$subset=='Incident_HIV'] = individuals$Diagnosis[individuals$subset=='Incident_HIV']
  
  return(list('df'=individuals, 'start_year'=sim_start_year))
}

getPartners <- function(directory, index_individuals, look_back_duration, sim_start_year=1960.5){
  ### This function is focused on getting the contacts of a given set of index
  ### patients. The index patients could in a data frame with other information
  ### we want to keep with them/their contacts. Due to the size of the
  ### relationship file we will be getting rid of any columns we do not need and
  ### any relationships/individuals that do not effect our time period. 
  ###
  ### Inputs:
  ###    directory (str)               : Path to directory where input files can
  ###                                    be found
  ###    index_individuals (data frame): A data frame with the IDs of the index
  ###                                    individuals and some information about
  ###                                   infection/disease progression
  ###    look_back_duration (float)    : How far into the past individuals are
  ###                                    asked to remember contacts from
  ###    sim_start_time (float)        : Time of the simulations start
  ###
  ### Return:
  ###    data frame                    : Data frame with contacts and their
  ###                                    disease information
  
  # Read in relationship start file and reduce the amount of data
  relationshipData <- read_csv(file.path(directory, 'RelationshipStart.csv'))
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
  relationshipEnd <- read_csv(file.path('RelationshipEnd.csv'))
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
    select(A_ID,
           B_ID,
           Rel_ID,
           B_Risk,
           B_HIV_disease_stage,
           Rel_start_time,
           Rel_actual_end_time)
  maleIndex <- rename(maleIndex,
                      Individual_ID=A_ID,
                      contact_ID=B_ID,
                      contact_Risk=B_Risk,
                      contact_HIV_stage=B_HIV_disease_stage)
  maleIndex <- left_join(maleIndex,
                         index_individuals[c('Individual_ID', 'begin_tracing','subset')],
                         by='Individual_ID')
  
  # Get relationships with female index individuals
  femaleIndex <- relationshipData %>%
    filter(B_ID %in% index_individuals$Individual_ID) %>% 
    select(B_ID,
           A_ID,
           Rel_ID,
           A_Risk,
           A_HIV_disease_stage,
           Rel_start_time,
           Rel_actual_end_time)
  femaleIndex <- rename(femaleIndex,
                        Individual_ID=B_ID,
                        contact_ID=A_ID,
                        contact_Risk=A_Risk,
                        contact_HIV_stage=A_HIV_disease_stage)
  femaleIndex <- left_join(femaleIndex,
                           index_individuals[c('Individual_ID', 'begin_tracing', 'subset')],
                           by='Individual_ID')
  rm(relationshipData)
  
  # Combined into a single data frame
  contactData <- bind_rows(maleIndex, femaleIndex)
  rm(maleIndex)
  rm(femaleIndex)
  
  # Reduce to relationships that were active in the look back period
  contactData$Rel_start_time <- contactData$Rel_start_time / 365 + sim_start_year
  contactData$Rel_actual_end_time <- contactData$Rel_actual_end_time / 365 + sim_start_year
  contactData <- contactData %>% filter(Rel_start_time <= begin_tracing &
                                        Rel_actual_end_time >= (begin_tracing - look_back_duration))
  
  return(contactData)
}

getTransData <- function(directory, contactData, tracing_delay=0){
  ### This function is focused on getting the transmission information for the
  ### contacts of the index individuals. Some transmitters we cannot trace to or
  ### died before we traced to them. These are removed them. Others we can trace
  ### to through multiple index individuals. All these are kept and determining
  ### what path finds them will be left to the sampling function.
  ###
  ### Inputs:
  ###    directory (str)         : Path to directory where input files can be
  ###                             found
  ###    contactData (data frame): Data frame with contacts and their disease
  ###                              information
  ###    tracing_delay (float)   : How long the tracing takes which effects the
  ###                              sample time
  ###
  ### Return:
  ###    data frame              : Data frame with contacts and their
  ###                              transmission and tracing information
  
  # Read in transmission file and reduce the amount of data
  transData <- read_csv(file.path('TransmissionReport.csv'))
  transData <- transData %>% select(YEAR,
                                    REL_ID,
                                    SRC_ID,
                                    DEST_ID)
  
  # Get number of transmissions of contacts before and after tracing start
  transData <- left_join(transData, contactData[c('contact_ID',
                                                  'Individual_ID',
                                                  'begin_tracing',
                                                  'contact_Risk',
                                                  'contact_HIV_stage',
                                                  'subset')],
                         by=c('SRC_ID'='contact_ID'), relationship="many-to-many")
  transData <- transData %>% group_by(SRC_ID, Individual_ID) %>% 
    summarise(num_trans_before = sum(YEAR<=begin_tracing),
              num_trans_after = sum(YEAR>begin_tracing),
              contact_Risk = contact_Risk[1],
              contact_HIV_stage = contact_HIV_stage[1],
              subset = subset[1],
              trace_date = begin_tracing[1] + tracing_delay,
              .groups = 'drop') %>% 
    mutate(num_trans_before = ifelse(is.na(num_trans_before), 0, num_trans_before),
           num_trans_after = ifelse(is.na(num_trans_after), 0, num_trans_after))
  
  
  # Make sure the dead are not found through tracing
  eventReport <- read_csv(file.path('ReportEventRecorder.csv'))
  deaths = getDeaths(eventReport)
  rm(eventReport)
  transData <- left_join(transData, deaths,
                         by=c('SRC_ID'='Individual_ID'))
  transData <- rename(transData, death_date = Death)
  transData <- transData %>% filter(death_date >= trace_date)
  rm(deaths)
  
  # Remove those that we couldn't traced to
  transData <- transData %>% filter(!is.na(Individual_ID))

  return(transData)
}

doSampling <- function(directory='.', tracing_group='Incident_HIV',
                       sample_rate=100, tracing_rate=100,
                       tracing_start_time=2010, tracing_end_time=2012,
                       acute_to_trace=FALSE, look_back_duration=3,
                       tracing_delay=0){
  ### This function is ties every together and prepares the data for phyloModels
  ### to build the tree from.
  ###
  ### Inputs:
  ###    directory (str)           : Path to directory where input files can be
  ###                                found
  ###    tracing_group (str)       : Which group/subset we are using for index
  ###                                individuals
  ###    sample_rate (float)       : Percentage of index individuals who actually
  ###                                respond to tracing efforts
  ###    tracing_rate (float)      : Success rate of tracing contacts gives as a
  ###                                percentage
  ###    trace_start_time (float)  : The date when the period for contract
  ###                                tracing started.
  ###    trace_end_time (float)    : The date when the period for contract
  ###                                tracing ended.
  ###    acute_to_trace (bool)     : This flags whether the incidence infections
  ###                                (those whose contacts will be traced) are
  ###                                just those diagnosed during the tracing
  ###                                period or if they also need to be acute at
  ###                                the time of diagnosis.
  ###    look_back_duration (float): How far into the past individuals are asked
  ###                                to remember contacts from
  ###    tracing_delay (float)     : How long the tracing takes which effects
  ###                                the sample time
  ###
  ### Return:
  ###    data frame                : Transmission line list for tracing
  ###                                information for building the tree
  
  # Do the tracing for everyone
  individuals = getGroups(directory, tracing_start_time, tracing_end_time, acute_to_trace)
  sim_start_year = individuals$start_year
  individuals = individuals$df
  contactData = getPartners(directory, individuals, look_back_duration, sim_start_year)
  rm(individuals)
  transData = getTransData(directory, contactData, tracing_delay)
  rm(contactData)
  
  # Reduce tracing based on desired sampling
  # 1. restrict to the group/subset we are currently looking at
  transData <- transData %>% filter(subset == tracing_group)
  # 2. select those enrolled in tracing
  ids = unique(transData$Individual_ID)
  ids = ids[runif(length(ids), min=0, max=100)<=sample_rate]
  transData <- transData %>% filter(Individual_ID %in% ids)
  # 3. select the contacts that get traced/contacted
  ids = unique(transData$SRC_ID)
  ids = ids[runif(length(ids), min=0, max=100)<=tracing_rate]
  transData <- transData %>% filter(SRC_ID %in% ids)
  
  # Resolve those that are contacts of more then one index person
  sampleData <- transData %>% group_by(SRC_ID) %>% slice_min(order_by=trace_date) %>% 
    filter(row_number() == sample(1:n(),1))
  sampleData <- rename(sampleData, contact_ID=SRC_ID)
  
  # Get the diagnosis time for each individual
  eventReport <- read_csv(file.path(directory, 'ReportEventRecorder.csv'))
  eventReport <- eventReport %>% select(Year, Event_Name, Individual_ID)
  diagnosis <- getDiagnosis(eventReport)
  diagnosis <- rename(diagnosis, sampleTime = Diagnosis)
  rm(eventReport)
  
  # Now that we have the sampled individuals and times add it to the transmission report
  transData <- read_csv(file.path('TransmissionReport.csv'))
  transData <- left_join(transData, diagnosis, by=c('DEST_ID'='Individual_ID'))
  transData <- left_join(transData, sampleData[c('contact_ID',
                                                 'Individual_ID',
                                                 'trace_date',
                                                 'contact_Risk',
                                                 'contact_HIV_stage',
                                                 'subset',
                                                 'num_trans_before',
                                                 'num_trans_after')],
                         by=c('DEST_ID'='contact_ID'))
  
  # Resolve whether tracing or diagnosis happened first and record how was
  # traced into.
  transData['traced_into'] = FALSE
  idx = is.na(transData$sampleTime) & !is.na(transData$trace_date)
  transData$traced_into[idx] = TRUE
  transData$sampleTime[is.na(transData$sampleTime)] = transData$trace_date[is.na(transData$sampleTime)]
  idx = which(transData$trace_date < transData$sampleTime)
  transData$traced_into[idx] = TRUE
  transData$sampleTime[idx] = transData$trace_date[idx]
  
  # Replace character NAs to smooth things over when data is pasted back to R
  transData$contact_Risk[is.na(transData$contact_Risk)] = ''
  
  return(transData)
}