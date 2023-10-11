print('Loading Tidyverse')
library('tidyverse')

directory = 'download/79a4cded-0249-ee11-aa0a-b88303911bc1/'

tracing_groups='Incident_HIV'
comparisons = 'all'
sample_rates=100
tracing_rates=100
tracing_start_times=2010
tracing_end_times=2012
acute_to_traces=TRUE
look_back_durations=3
tracing_delays=0
traits = 'trans_source'

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
    summarize(Infected = suppressWarnings(min(Year[Event_Name=='NewInfectionEvent'])),
              Diagnosis = suppressWarnings(min(Year[Year >= Year[Event_Name=='NewInfectionEvent'] &
                                     Event_Name!='NewInfectionEvent']))) %>%
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
  ###                              just those diagnosed during the tracing
  ###                              period or if they have to have been
  ###                              infected for less than 3 months.
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
  individuals <- rename(individuals, Index_ID = Individual_ID)
  
  # Get the birth time for each individual
  births <- eventReport %>% filter(Event_Name == 'Births') %>%
    select(Individual_ID, Year)
  individuals <- left_join(individuals, births, by=c('Index_ID'='Individual_ID'))
  individuals <- rename(individuals, Index_Born = Year)
  rm(births)
  
  # Get the death time for each individual
  deaths = getDeaths(eventReport)
  individuals <- left_join(individuals, deaths, by=c('Index_ID'='Individual_ID'))
  individuals <- rename(individuals, Index_Dies = Death)
  rm(deaths)
  
  # Get the diagnosis time for each individual
  diagnosis = getDiagnosis(eventReport)
  individuals <- left_join(individuals, diagnosis, by=c('Index_ID'='Individual_ID'))
  individuals <- rename(individuals, Index_Diagnosed=Diagnosis, Index_Infected=Infected)
  rm(diagnosis)
  
  # Get time no longer acute
  endAcute <- eventReport %>% filter(Event_Name %in% c('HIVInfectionStageEnteredLatent',
                                                       'HIVInfectionStageEnteredAIDS',
                                                       'HIVInfectionStageEnteredOnART')) %>%
    group_by(Individual_ID) %>%
    summarize(Index_AcuteEnded = suppressWarnings(min(Year))) %>%
    mutate(Index_AcuteEnded = ifelse(is.infinite(Index_AcuteEnded), NA, Index_AcuteEnded))
  individuals <- left_join(individuals, endAcute, by=c('Index_ID'='Individual_ID'))
  rm(endAcute)
  rm(eventReport)
  
  # Reduce to only individuals that are alive for any amount of time during the
  # tracing period
  individuals <- individuals %>% filter(Index_Dies >= trace_start_time & 
                                        Index_Born <= trace_end_time)
  
  # Create the three subgroups
  individuals['Index_Subset'] = 'HIV-'
  # Set everyone infected by the end of the trace period to the prevalent group
  individuals$Index_Subset[individuals$Index_Infected<=trace_end_time %>%
                             replace_na(FALSE)] = 'Prevalent_HIV'
  # In the wiki post incident and prevalent are only diagnosed individuals
  # Here I move the undiagnosed to a separate category commenting this out will
  # leave them as prevalent
  individuals$Index_Subset[individuals$Index_Infected<=trace_end_time &
                            (is.na(individuals$Index_Diagnosed) |
                             individuals$Index_Diagnosed>trace_end_time) %>%
                             replace_na(FALSE)
                           ] = 'Undiagnosed_HIV'
  if(acute_to_trace){
    # Only trace (incident) those that are recently diagnosed and acute
    individuals$Index_Subset[individuals$Index_Diagnosed>=trace_start_time &
                             individuals$Index_Diagnosed<=trace_end_time &
                             ((individuals$Index_Diagnosed - individuals$Index_Infected) <= 3/12) %>%
                               replace_na(FALSE)] = 'Incident_HIV'
    # Also trace those diagnosed before tracing starts, but still acute when tracing starts
    individuals$Index_Subset[((individuals$Index_Diagnosed - individuals$Index_Infected) <= 3/12) &
                             (trace_start_time <= individuals$Index_Infected + 3/12) &
                               individuals$Index_Diagnosed <= trace_end_time%>%
                               replace_na(FALSE)] = 'Incident_HIV'
  } else{
    # Trace (incident) all recent diagnosed
    individuals$Index_Subset[individuals$Index_Diagnosed>=trace_start_time &
                             individuals$Index_Diagnosed<=trace_end_time %>%
                               replace_na(FALSE)] = 'Incident_HIV'
  }
  individuals$Index_Subset <- factor(individuals$Index_Subset, levels=c('Incident_HIV',
                                                                        'Prevalent_HIV',
                                                                        'Undiagnosed_HIV',
                                                                        'HIV-',
                                                                        'HIV+',
                                                                        'Untraceable'))
  
  # Make sure each person has a time that their contacts are traced
  individuals['begin_tracing'] = trace_start_time
  individuals$begin_tracing[individuals$Index_Subset=='Incident_HIV']   = individuals$Index_Diagnosed[individuals$Index_Subset=='Incident_HIV']
  individuals$begin_tracing[individuals$begin_tracing<trace_start_time] = trace_start_time
  
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
                                                  A_ID,
                                                  B_ID,
                                                  A_IndividualProperties,
                                                  B_IndividualProperties,
                                                  A_HIV_disease_stage,
                                                  B_HIV_disease_stage)
  
  # Read in relationship end file to get the true end of each relationship
  relationshipEnd <- read_csv(file.path(directory, 'RelationshipEnd.csv'))
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
    filter(A_ID %in% index_individuals$Index_ID) %>% 
    select(A_ID,
           B_ID,
           Rel_ID,
           B_Risk,
           B_HIV_disease_stage,
           Rel_start_time,
           Rel_actual_end_time)
  maleIndex <- rename(maleIndex,
                      Index_ID=A_ID,
                      contact_ID=B_ID,
                      contact_Risk=B_Risk,
                      contact_HIV_stage=B_HIV_disease_stage)
  maleIndex <- left_join(maleIndex,
                         index_individuals[c('Index_ID', 'begin_tracing','Index_Subset', 'Index_Diagnosed')],
                         by='Index_ID', relationship="many-to-many")
  
  # Get relationships with female index individuals
  femaleIndex <- relationshipData %>%
    filter(B_ID %in% index_individuals$Index_ID) %>% 
    select(B_ID,
           A_ID,
           Rel_ID,
           A_Risk,
           A_HIV_disease_stage,
           Rel_start_time,
           Rel_actual_end_time)
  femaleIndex <- rename(femaleIndex,
                        Index_ID=B_ID,
                        contact_ID=A_ID,
                        contact_Risk=A_Risk,
                        contact_HIV_stage=A_HIV_disease_stage)
  femaleIndex <- left_join(femaleIndex,
                           index_individuals[c('Index_ID', 'begin_tracing', 'Index_Subset', 'Index_Diagnosed')],
                           by='Index_ID', relationship="many-to-many")
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

getTransData <- function(directory, contactData, tracing_start_time, tracing_delay=0){
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
  transData <- read_csv(file.path(directory, 'TransmissionReport.csv'))
  transData <- transData %>% select(YEAR,
                                    REL_ID,
                                    SRC_ID,
                                    DEST_ID,
                                    SRC_RISK,
                                    SRC_STAGE)
  
  # Get number of transmissions of contacts before and after tracing start
  transData <- full_join(transData, contactData[c('contact_ID',
                                                  'Index_ID',
                                                  'begin_tracing',
                                                  'contact_Risk',
                                                  'contact_HIV_stage',
                                                  'Index_Subset',
                                                  'Index_Diagnosed')],
                         by=c('SRC_ID'='contact_ID'), relationship="many-to-many")
  transData <- transData %>% group_by(SRC_ID, Index_ID) %>% 
    summarise(num_trans_before = sum(YEAR<=begin_tracing + tracing_delay),
              num_trans_after = sum(YEAR>begin_tracing + tracing_delay),
              contact_Risk = contact_Risk[1],
              contact_HIV_stage = contact_HIV_stage[1],
              Index_Subset = Index_Subset[1],
              Index_Diagnosed = Index_Diagnosed[1],
              contact_trace_date = begin_tracing[1] + tracing_delay,
              SRC_RISK = SRC_RISK[1],
              SRC_STAGE = SRC_STAGE[1],
              .groups = 'drop') %>% 
    mutate(num_trans_before = ifelse(is.na(num_trans_before), 0, num_trans_before),
           num_trans_after = ifelse(is.na(num_trans_after), 0, num_trans_after))
 transData['contact_ID'] = transData['SRC_ID']
  
  
  # Make sure the dead are not found through tracing
  eventReport <- read_csv(file.path(directory, 'ReportEventRecorder.csv'))
  deaths = getDeaths(eventReport)
  deaths <- rename(deaths, contact_died = Death)
  transData <- left_join(transData, deaths,
                         by=c('SRC_ID'='Individual_ID'))
transData <- transData %>%
  filter(contact_died >= contact_trace_date %>% replace_na(TRUE)) %>%
  filter(! contact_died < tracing_start_time)
  rm(deaths)
  
  # Get the diagnosis and infection time for each contact
  diagnosis = getDiagnosis(eventReport)
  diagnosis <- rename(diagnosis, Contact_Diagnosed=Diagnosis, Contact_Infected=Infected)
  transData <- left_join(transData, diagnosis, by=c('SRC_ID'='Individual_ID'))
  rm(diagnosis)
  rm(eventReport)
  
  # Make sure we have the risk and stage for the untraceable
  transData$contact_Risk[is.na(transData$contact_Risk)] = transData$SRC_RISK[is.na(transData$contact_Risk)]
  transData$contact_HIV_stage[is.na(transData$contact_HIV_stage)] = transData$SRC_STAGE[is.na(transData$contact_HIV_stage)]
  transData <- transData %>% select(!c('SRC_RISK', 'SRC_STAGE'))
  
  # # Remove those that we couldn't traced to
  # transData <- transData %>% filter(!is.na(Individual_ID))

  return(transData)
}


### This function was written to prepare the data for building a tree.
### Specifically, it does the tracing and assigns sample times to the sampled
### individuals. Will include a random sample for comparison.
doSampling <- function(directory='.', tracing_groups='Incident_HIV',
                       n_samples=5000, tracing_start_times=2010,
                       tracing_duration=2,
                       acute_to_traces=FALSE, look_back_durations=3,
                       tracing_delays=0){
  ### This function is ties every together and prepares the data for phyloModels
  ### to build the tree from.
  ###
  ### Inputs:
  ###    directory (str)            : Path to directory where input files can be
  ###                                 found
  ###    tracing_group (vector)     : Which group(s)/subset(s) we are using for
  ###                                 index individuals. Current options are:
  ###                                 HIV-, HIV+, HIV_Incident, HIV_Prevalent,
  ###                                 HIV_Undiagnosed
  ###    n_samples (vector)         : Number of samples to gather via tracing
  ###    tracing_start_time (vector): The date(s) when the period for contract
  ###                                 tracing starts
  ###    tracing_duration (vector)  : The length of the tracing campaign. All
  ###                                 samples will be gathered during the
  ###                                 duration unless we have to go over to
  ###                                 reach the requested number
  ###    acute_to_trace (vector)    : This flags whether the incidence infections
  ###                                 (those whose contacts will be traced) are
  ###                                 just those diagnosed during the tracing
  ###                                 period or if they have to have been
  ###                                 infected for less than 3 months.
  ###    look_back_duration (vector): How far into the past individuals are
  ###                                 asked to remember contacts from
  ###    tracing_delay (vector)     : How long it takes from a contact being
  ###                                 named to them being contacted and tested
  ###                                 which effects the sample time
  ###
  ### Return:
  ###    data frame                : Transmission line list for tracing
  ###                                information for building the tree
  
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
  transData <- read_csv(file.path(directory, 'TransmissionReport.csv'))
  transData <- left_join(transData, diagnosis, by=c('DEST_ID'='Individual_ID'))
  transData <- left_join(transData, sampleData[c('contact_ID',
                                                 'Individual_ID',
                                                 'trace_date',
                                                 'contact_Risk',
                                                 'contact_HIV_stage',
                                                 'Index_Subset',
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
  
  lineList <- read_csv(file.path(directory, 'TransmissionReport.csv'))
  eventReport <- read_csv(file.path(directory, 'ReportEventRecorder.csv'))
  eventReport <- eventReport %>% select(Year, Event_Name, Individual_ID)
  diagnosis <- getDiagnosis(eventReport)
  diagnosis <- rename(diagnosis, DEST_Diagnosed = Diagnosis, DEST_ID = Individual_ID)
  rm(eventReport)
  lineList <- left_join(lineList, diagnosis, by=c('DEST_ID'='Individual_ID'))
  rm(diagnosis)
  # Iterate over tracing start times
  for(tracing_start_time in tracing_start_times){
    # Iterate over tracing end times
    for(tracing_duration in tracing_durations){
      # Iterate over tracing delays
      for(tracing_delay in tracing_delays){
        # Iterate over look back periods
        for(look_back_duration in look_back_durations){
          # Iterate over definitions of incident definitions
          for(acute_to_trace in acute_to_traces){
            # Do the tracing for everyone using the currently parameters that effect tracing
            # By having this here we can avoid repeating this step in the following loops
            individuals = getGroups(directory, tracing_start_time, tracing_start_time+tracing_duration, acute_to_trace)
            sim_start_year = individuals$start_year
            individuals = individuals$df
            # the index individuals are changing from one generation to the next
            contactData = getPartners(directory, individuals, look_back_duration, sim_start_year)
            rm(individuals)
            transData = getTransData(directory, contactData, tracing_start_time, tracing_delay)
            transData$Index_Subset[is.na(transData$Index_Subset)] = 'Untraceable'
            rm(contactData)
            # Iterate over the number of samples
            for(n_sample in n_samples){
              # Iterate over tracing groups
              for(tracing_group in tracing_groups){
                ## Select samples
                temp = transData
                
                # Restrict to the group/subset we are currently looking at
                if(tracing_group %in% c('HIV+')){
                  temp$Index_Subset[temp$Index_Subset %in% c('Incident_HIV',
                                                             'Prevalent_HIV',
                                                             'Undiagnosed_HIV')] = 'HIV+'
                }
                temp = temp %>% filter(Index_Subset == tracing_group)
                
                # Remove those that are not infected before they are traced into
                temp = temp %>% filter(Contact_Infected <= contact_trace_date %>% 
                                         replace_na(FALSE))
                
                # Get a count of how many people could be sampled
                N_infections = length(unique(c(temp$Index_ID[temp$Index_Diagnosed<=tracing_start_time+tracing_duration],
                                               temp$contact_ID)))
                # Test if we have enough potential samples
                if(N_infections < n_sample){
                  print('WARNING...Not enough infections in tracing window to hit sample request')
                }
                
                # Get an idea of how many index versus contacts to draw to reach the requested number of samples
                N_index = length(unique(temp$Index_ID[temp$Index_Diagnosed<=tracing_start_time+tracing_duration]))
                contacts_per_index = length(unique(temp$contact_ID))/N_index
                N_index_sample = ceiling(n_sample/(contacts_per_index+1))
                
                # Select indexes to sample
                index_sampled = temp %>% filter(temp$Index_Diagnosed<=tracing_start_time+tracing_duration) %>% 
                  select(Index_ID) %>% unique() %>%
                  filter(row_number() %in% sample(1:n(),N_index_sample))
                index_sampled <- left_join(index_sampled, temp %>%
                                             select(Index_ID, Index_Diagnosed) %>% 
                                             distinct(), by='Index_ID')
                samples_remaining = n_sample - length(index_sampled$Index_ID)
                
                # Remove sampled index cases as contacts of other index cases
                temp <- temp %>% filter(! contact_ID %in% index_sampled$Index_ID)
                
                # Sample contacts
                N_contacts = length(unique(temp$contact_ID[temp$Index_ID %in% index_sampled$Index_ID]))
                if(samples_remaining <= N_contacts){
                  contacts_sampled = temp %>% filter(temp$Index_ID %in% index_sampled$Index_ID) %>% 
                    select(contact_ID) %>% unique() %>%
                    filter(row_number() %in% sample(1:n(),samples_remaining))
                  contacts_sampled <- left_join(contacts_sampled, temp %>%
                                                  filter(Index_ID %in% index_sampled$Index_ID) %>% 
                                                  select(contact_ID, contact_trace_date),
                                                by='contact_ID')
                  contacts_sampled <- contacts_sampled %>% group_by(contact_ID) %>%
                    filter(row_number() == sample(1:n(),1))
                } else {
                  
                }
                index_sampled <- rename(index_sampled, DEST_ID = Index_ID, 'Sampled_TG-{tracing_group}_TST-{tracing_start_time}_TL-{tracing_duration}_TD-{tracing_delay}_LBD-{look_back_duration}_NS-{n_sample}' := Index_Diagnosed)
                contacts_sampled <- rename(contacts_sampled, DEST_ID = contact_ID, 'Sampled_TG-{tracing_group}_TST-{tracing_start_time}_TL-{tracing_duration}_TD-{tracing_delay}_LBD-{look_back_duration}_NS-{n_sample}' = contact_trace_date)
                sampledData <- bind_rows(index_sampled, contacts_sampled)
                lineList <- left_join(lineList, sampledData, by='DEST_ID')
              }
            }
          }
        }
      }
    }
  }
  
  return(transData)
}

### This function is for doing the tracing analysis (experiments 1 and 2). It
### takes in the various options for doing tracing, does the tracing, and
### calculates various statistics 
doTracingAnalysis <- function(directory='.', tracing_groups='Incident_HIV',
                              sample_rates=100, tracing_rates=100,
                              tracing_start_times=2010, tracing_end_times=2012,
                              acute_to_traces=FALSE, look_back_durations=3,
                              tracing_delays=0, num_generations=3){
  ### This function does the tracing and calculates various metrics on how the
  ### tracing preforms. All the inputs that are listed as vector can also
  ### receive a single value. TODO: replace tracing_end_times with duration
  ###
  ### Inputs:
  ###    directory (str)            : Path to directory where input files can be
  ###                                 found
  ###    tracing_group (vector)     : Which group(s)/subset(s) we are using for
  ###                                 index individuals. Current options are:
  ###                                 HIV-, HIV+, HIV_Incident, HIV_Prevalent,
  ###                                 HIV_Undiagnosed
  ###    sample_rate (vector)       : Percentage(s) of index individuals who
  ###                                 actually respond to tracing efforts
  ###    tracing_rate (vector)      : Success rate(s) of tracing contacts gives
  ###                                 as a percentage
  ###    trace_start_time (vector)  : The date(s) when the period for contract
  ###                                 tracing starts
  ###    trace_end_time (vector)    : The date(s) when the period for contract
  ###                                 tracing ends
  ###    acute_to_trace (vector)    : This flags whether the incidence infections
  ###                                 (those whose contacts will be traced) are
  ###                                 just those diagnosed during the tracing
  ###                                 period or if they have to have been
  ###                                 infected for less than 3 months.
  ###    look_back_duration (vector): How far into the past individuals are
  ###                                 asked to remember contacts from
  ###    tracing_delay (vector)     : How long it takes from a contact being
  ###                                 named to them being contacted and tested
  ###                                 which effects the sample time
  ###    num_generations (vector)   : Number of generations to trace for
  ###
  ### Return:
  ###    data frame                : A data frame containing the metrics of the
  ###                                tracing analysis
  
  # Create data frame to hold relative risk
  col_names <- c('tracing_group',
                 'comparison',
                 'trait',
                 'sample_rate',
                 'tracing_rate',
                 'tracing_start_time',
                 'tracing_end_time',
                 'campaign_duration',
                 'tracing_delay',
                 'look_back_duration',
                 'acute_to_trace',
                 'generation',
                 'FP',
                 'TN',
                 'TP',
                 'FN',
                 'FPreduced',
                 'TNreduced',
                 'TPreduced',
                 'FNreduced',
                 'N_HIV',
                 'N_untreated',
                 'N_undiagnosed',
                 'N_index',
                 'contact_mean',
                 'contact_std',
                 'contacts_per_index')
  RR = data.frame(matrix(NA,    # Create empty data frame
                         nrow = length(tracing_groups) *
                           length(comparisons) *
                           length(traits) *
                           length(sample_rates) *
                           length(tracing_rates) *
                           length(tracing_start_times) *
                           length(tracing_end_times) *
                           length(tracing_delays) *
                           length(look_back_durations) *
                           length(acute_to_traces) *
                           num_generations,
                         ncol = length(col_names)))
  colnames(RR) <- col_names
  # RR$contact_dist = list()
  # Have a counter for the row
  i = 1
  # Iterate over tracing start times
  for(tracing_start_time in tracing_start_times){
    # Iterate over tracing end times
    for(tracing_end_time in tracing_end_times){
      # Iterate over tracing delays
      for(tracing_delay in tracing_delays){
        # Iterate over look back periods
        for(look_back_duration in look_back_durations){
          # Iterate over definitions of incident definitions
          for(acute_to_trace in acute_to_traces){
            # Do the tracing for everyone using the currently parameters that effect tracing
            # By having this here we can avoid repeating this step in the following loops
            individuals = getGroups(directory, tracing_start_time, tracing_end_time, acute_to_trace)
            sim_start_year = individuals$start_year
            individuals = individuals$df
            # Iteration of the number of generations
            for(generation in seq(num_generations)){
              # For each generation we have to re-generate the contacts list as 
              # the index individuals are changing from one generation to the next
              contactData = getPartners(directory, individuals, look_back_duration, sim_start_year)
              transData = getTransData(directory, contactData, tracing_start_time, tracing_delay)
              rm(contactData)
              # Iterate over the populations.
              for(comparison in comparisons){
                # In this loop members of the population are removed so a copy
                # of the data is made to increase our reuse of previous
                # calculations
                temp = transData
                if(comparison == 'traceable'){
                  if(sum(is.na(temp$Index_ID))==0){
                    print('WARNING: traceable is same as all')
                  }
                    temp = temp %>% filter(!is.na(Index_ID))
                  } else if(comparison == 'sampleableAtTrace'){
                    temp = temp %>% filter(contact_HIV_stage != 0)
                  } else if(comparison == 'sampleableEventually'){
                    temp = temp %>% filter(!is.na(Contact_Infected))
                  }
                # Iterate over tracing groups
                for(tracing_group in tracing_groups){
                  # Iterate over enrollment rates
                  for(sample_rate in sample_rates){
                    # Iterate over tracing rates
                    for(tracing_rate in tracing_rates){
                      # mark the "positive diagnosis" if tracing is the "test"
                      temp['result'] = FALSE
                      if(tracing_group %in% c('Incident_HIV',
                                              'Prevalent_HIV',
                                              'Undiagnosed_HIV',
                                              'HIV-')){
                        temp$result[temp$Index_Subset == tracing_group] = TRUE
                      } else if(tracing_group == 'HIV+'){
                        temp$result[temp$Index_Subset == 'Incident_HIV'] = TRUE
                        temp$result[temp$Index_Subset == 'Prevalent_HIV'] = TRUE
                        temp$result[temp$Index_Subset == 'Undiagnosed_HIV'] = TRUE
                      }
                      
                      # mark the "positive diagnosis" if tracing is the "test"
                      # select those enrolled in tracing
                      ids = unique(temp$Index_ID[temp$result])
                      ids = ids[runif(length(ids), min=0, max=100)<=100*((sample_rate/100)^generation*(tracing_rate/100)^(generation-1))]
                      temp$result[! temp$Index_ID %in% ids] = FALSE
                      N_index = length(unique(temp$Index_ID[temp$result]))
                      # select the contacts that get traced/contacted
                      ids = unique(temp$contact_ID[temp$result])
                      ids = ids[runif(length(ids), min=0, max=100)<=tracing_rate]
                      temp$result[! temp$contact_ID %in% ids] = FALSE
                      N_index_0 = N_index - length(unique(temp$Index_ID[temp$result]))
                      # Iterate over traits to be consider the 'disease'
                      for(trait in traits){
                        # mark the "positive disease" based on trait
                        if(trait == 'trans_source'){
                          temp['condition'] = (temp$num_trans_before + temp$num_trans_after) > 0
                        } else if(trait == 'yet_to_transmit') {
                          temp['condition'] = temp$num_trans_after > 0
                        } else if(trait == 'recently_infected'){
                          temp['condition'] = temp$contact_HIV_stage == 1
                        } else if(trait == 'high_risk'){
                          temp['condition'] = temp$contact_Risk == 'HIGH'
                        } else if(trait == 'not_yet_infected'){
                          temp['condition'] = temp$contact_trace_date < temp$Contact_Infected
                          temp$condition[is.na(temp$condition)] = FALSE
                        }
                        
                        # Save parameters and some population statistics
                        RR$tracing_group[i] = tracing_group
                        RR$comparison[i] = comparison
                        RR$trait[i] = trait
                        RR$sample_rate[i] = sample_rate
                        RR$tracing_rate[i] = tracing_rate
                        RR$tracing_start_time[i] = tracing_start_time
                        RR$tracing_end_time[i] = tracing_end_time
                        RR$campaign_duration[i] = tracing_end_time - tracing_start_time
                        RR$tracing_delay[i] = tracing_delay
                        RR$look_back_duration[i] = look_back_duration
                        RR$acute_to_trace[i] = acute_to_trace
                        RR$generation[i] = generation
                        RR$contacts_per_index[i] = length(temp$Index_ID[temp$result]) / N_index
                        RR$N_HIV[i] = sum(temp$result & (temp$contact_HIV_stage != 0) & (temp$Contact_Infected <= temp$contact_trace_date))
                        RR$N_untreated[i] = sum(temp$result & (temp$contact_HIV_stage != 0) & (temp$contact_HIV_stage != 4))
                        RR$N_undiagnosed[i] = sum(temp$result & (temp$contact_HIV_stage != 0) & (temp$Contact_Diagnosed > temp$contact_trace_date), na.rm=TRUE)
                        RR$N_index[i] = N_index
                        contact_dist = temp %>% filter(result) %>% group_by(Index_ID) %>% summarise(n_contact=n())
                        dist = tabulate(contact_dist$n_contact[!is.na(contact_dist$Index_ID)])
                        # if(any(is.na(contact_dist$Index_ID))){
                        #   dist = c(contact_dist$n_contact[is.na(contact_dist$Index_ID)], dist)  
                        # } else {
                        #   dist = c(0, dist)
                        # }
                        # RR$contact_mean[i] = mean(contact_dist$n_contact[!is.na(contact_dist$Index_ID)])
                        # RR$contact_std[i] = sd(contact_dist$n_contact[!is.na(contact_dist$Index_ID)])
                        dist = c(N_index_0, dist)
                        RR$contact_mean[i] = weighted.mean(seq(0,length(dist)-1), dist, na.rm=TRUE)
                        RR$contact_std[i] = sqrt(weighted.mean( (seq(0,length(dist)-1)-RR$contact_mean[i])^2, dist, na.rm=TRUE ))
                        
                        # set confusion matrix
                        RR$FP[i] = sum( temp$result & !temp$condition)
                        RR$TN[i] = sum(!temp$result & !temp$condition)
                        RR$TP[i] = sum( temp$result &  temp$condition)
                        RR$FN[i] = sum(!temp$result &  temp$condition)
                        
                        # remove duplicates
                        # start by removing non-contacts that also appear as contacts
                        temp2 <- temp %>% filter( !(!result & (contact_ID %in% contact_ID[result])) )
                        # next remove contacts/non-contacts that are also index individuals
                        temp2 <- temp2 %>% filter( !( (contact_ID %in% Index_ID[result])) )
                        # lastly remove any contacts/non-contacts that appear multiple times
                        temp2 <- temp2 %>% group_by(contact_ID) %>% filter(row_number() == sample(1:n(),1))
                        
                        # set reduced confusion matrix
                        RR$FPreduced[i] = sum( temp2$result & !temp2$condition)
                        RR$TNreduced[i] = sum(!temp2$result & !temp2$condition)
                        RR$TPreduced[i] = sum( temp2$result &  temp2$condition)
                        RR$FNreduced[i] = sum(!temp2$result &  temp2$condition)
                        i = i + 1
                      }
                    }
                  }
                }
              }
              # Set up for next generation of tracing
              individuals <- transData %>% select(contact_ID, contact_trace_date, Index_Subset, Contact_Diagnosed)
              individuals <- rename(individuals, c(Index_ID = contact_ID,
                                                   Index_Diagnosed = Contact_Diagnosed,
                                                   begin_tracing = contact_trace_date))
              individuals <- individuals %>% filter(!is.na(Index_ID))
              individuals <- individuals %>% group_by(Index_ID, Index_Subset) %>%
                slice_min(order_by=begin_tracing) %>% filter(row_number() == sample(1:n(),1))
            }
          }
        }
      }
    }
  }
  
  ## Calculate statistics
  
  # Actual Positives
  RR['P'] = RR['TP'] + RR['FN']
  # Actual Negatives
  RR['N'] = RR['FP'] + RR['TN']
  # Predicted Positives
  RR['PP'] = RR['TP'] + RR['FP']
  RR['PPreduced'] = RR['TPreduced'] + RR['FPreduced']
  # Predicted Negatives
  RR['PN'] = RR['FN'] + RR['TN']
  # Total Population
  RR['total'] = RR['P'] + RR['N']
  # Prevalence
  RR['prevalence'] = RR['P'] / RR['total']
  # Percent of contacts HIV+
  RR['percent_HIV'] = RR['N_HIV']/RR['PP']
  # Number needed to trace positive
  RR['NNTTP'] = RR['N_index']/RR['N_HIV']
  # Number needed to trace untreated
  RR['NNTTT'] = RR['N_index']/RR['N_untreated']
  # Number needed to trace undiagnosed
  RR['NNTTD'] = RR['N_index']/RR['N_undiagnosed']
  
  # Accuracy
  RR['ACC'] = (RR['TP'] + RR['TN']) / RR['total']
  # Positive Predictive Value
  RR['PPV'] = RR['TP'] / RR['PP']
  RR['PPVreduced'] = RR['TPreduced'] / RR['PPreduced']
  # False Discovery Rate
  RR['FDR'] = RR['FP'] / RR['PP']
  # False Omission Rate
  RR['FOR'] = RR['FN'] / RR['PN']
  # Negative Predictive Value
  RR['NPV'] = RR['TN'] / RR['PN']
  
  # True Positive Rate
  RR['TPR'] = RR['TP'] / RR['P']
  # False Positive Rate
  RR['FPR'] = RR['FP'] / RR['N']
  # False Negative Rate
  RR['FNR'] = RR['FN'] / RR['P']
  # True Negative Rate
  RR['TNR'] = RR['TN'] / RR['N']
  
  # Positive Likelihood Ratio
  RR['LRP'] = RR['TPR'] / RR['FPR']
  # Negative Likelihood Ratio
  RR['LRN'] = RR['FNR'] / RR['TNR']
  # Diagnostic Odds Ratio
  RR['DOR'] = RR['LRP'] / RR['LRN']
  
  # Threat Score
  RR['TS'] = RR['TP'] / (RR['TP'] + RR['FN'] + RR['FP'])
  # Informedness
  RR['informedness'] = RR['TPR'] + RR['TNR'] - 1
  # Markedness
  RR['markedness'] = RR['PPV'] + RR['NPV'] - 1
  # Prevalence Threshold
  RR['prevalenceThreshold'] = (sqrt(RR['TPR'] * RR['FPR']) - RR['FPR']) / (
    RR['TPR'] - RR['FPR'])
  # Matthews Correlation Coefficient
  RR['MCC'] = sqrt(RR['TPR'] * RR['TNR'] * RR['PPV'] * RR['NPV']) - sqrt(
    RR['FNR'] * RR['FPR'] * RR['FOR'] * RR['FDR'])
  # Fowlkes-Mallows Index
  RR['FMI'] = sqrt(RR['PPV'] * RR['TPR'])
  # F1 Score
  RR['F1'] = 2 * RR['PPV'] * RR['TPR'] / (RR['PPV'] + RR['TPR'])
  # Balanced Accuracy
  RR['BA'] = (RR['TPR'] + RR['TNR']) / 2
  # Relative Risk vs comparison
  RR['RR'] = RR['PPV'] / RR['FOR']
  # Relative Risk vs all
  RR['RRall'] = RR['PPV'] / (RR['P'] / RR['total'])
  
  # Get relative risk compared to the HIV- index group
  RR <- RR %>% group_by(comparison,
                        trait,
                        sample_rate,
                        tracing_rate,
                        tracing_start_time,
                        tracing_end_time,
                        tracing_delay,
                        look_back_duration,
                        acute_to_trace,
                        generation) %>% mutate('RRHIV-' = PPV/PPV[tracing_group=='HIV-'],
                                               'RRHIV-reduced' = PPVreduced/PPVreduced[tracing_group=='HIV-'])
  
  # write.table(RR, file.path(directory, 'RROutput.csv'), sep=',', dec='.', row.names=FALSE)
  return(RR)
}