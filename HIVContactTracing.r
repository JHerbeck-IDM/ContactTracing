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
  individuals$Index_Subset[individuals$Index_Infected<=trace_end_time] = 'Prevalent_HIV'
  # In the wiki post incident and prevalent are only diagnosed individuals
  # Here I move the undiagnosed to a separate category commenting this out will
  # leave them as prevalent
  individuals$Index_Subset[individuals$Index_Infected<=trace_end_time &
                            (is.na(individuals$Index_Diagnosed) |
                             individuals$Index_Diagnosed>=trace_end_time)
                           ] = 'Undiagnosed_HIV'
  if(acute_to_trace){
    # Only trace (incident) those that are recently diagnosed and acute
    individuals$Index_Subset[individuals$Index_Diagnosed>=trace_start_time &
                             individuals$Index_Diagnosed<=trace_end_time &
                             (individuals$Index_Diagnosed - individuals$Index_Infected) <= 3/12] = 'Incident_HIV'
    # Also trace those diagnosed before tracing starts, but still acute when tracing starts
    individuals$Index_Subset[(individuals$Index_Diagnosed - individuals$Index_Infected) <= 3/12 &
                             (trace_start_time - individuals$Index_Infected) <= 3/12] = 'Incident_HIV'
  } else{
    # Trace (incident) all recent diagnosed
    individuals$Index_Subset[individuals$Index_Diagnosed>=trace_start_time &
                             individuals$Index_Diagnosed<=trace_end_time] = 'Incident_HIV'
  }
  individuals$Index_Subset <- factor(individuals$Index_Subset, levels=c('Incident_HIV',
                                                                        'Prevalent_HIV',
                                                                        'Undiagnosed_HIV',
                                                                        'HIV-'))
  
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


### This function was written to prepare the data for building a tree. For now
### it will not be used
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
  
  return(transData)
}


doTracing <- function(directory='.', tracing_groups='Incident_HIV',
                      sample_rates=100, tracing_rates=100,
                      tracing_start_times=2010, tracing_end_times=2012,
                      acute_to_traces=FALSE, look_back_durations=3,
                      tracing_delays=0, unique_individuals=TRUE){
  
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
                 'unique_individuals',
                 'FP',
                 'TN',
                 'TP',
                 'FN',
                 'N_HIV',
                 'N_untreated',
                 'N_undiagnosed',
                 'N_index',
                 'contact_mean',
                 'contact_std',
                 'contacts_per_index',
                 'unique_contacts_per_index')
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
                           length(acute_to_traces),
                         ncol = length(col_names)))
  colnames(RR) <- col_names
  # RR$contact_dist = list()
  i = 1
  for(tracing_start_time in tracing_start_times){
    for(tracing_end_time in tracing_end_times){
      for(tracing_delay in tracing_delays){
        for(look_back_duration in look_back_durations){
          for(acute_to_trace in acute_to_traces){
            # Do the tracing for everyone using the currently parameters that effect tracing
            individuals = getGroups(directory, tracing_start_time, tracing_end_time, acute_to_trace)
            sim_start_year = individuals$start_year
            individuals = individuals$df
            contactData = getPartners(directory, individuals, look_back_duration, sim_start_year)
            rm(individuals)
            transData = getTransData(directory, contactData, tracing_start_time, tracing_delay)
            rm(contactData)
            for(tracing_group in tracing_groups){
              for(comparison in comparisons){
                for(trait in traits){
                  for(sample_rate in sample_rates){
                    for(tracing_rate in tracing_rates){
                      for(unique_individual in unique_individuals){
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
                        RR$unique_individuals[i] = unique_individual
                        if(tracing_group %in% c('Incident_HIV',
                                                'Prevalent_HIV',
                                                'Undiagnosed_HIV',
                                                'HIV-')){
                          
                          # mark the "positive diagnosis" if tracing is the "test"
                          temp = transData
                          temp['result'] = FALSE
                          temp$result[temp$Index_Subset == tracing_group] = TRUE
                        } else if(tracing_group == 'HIV+'){
                          # mark the "positive diagnosis" if tracing is the "test"
                          temp = transData
                          temp['result'] = FALSE
                          temp$result[temp$Index_Subset == 'Incident_HIV'] = TRUE
                          temp$result[temp$Index_Subset == 'Prevalent_HIV'] = TRUE
                          temp$result[temp$Index_Subset == 'Undiagnosed_HIV'] = TRUE
                        }
                        
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
                        
                        # mark the "positive diagnosis" if tracing is the "test"
                        # select those enrolled in tracing
                        ids = unique(temp$Index_ID[temp$result])
                        ids = ids[runif(length(ids), min=0, max=100)<=sample_rate]
                        temp$result[! temp$Index_ID %in% ids] = FALSE
                        # select the contacts that get traced/contacted
                        ids = unique(temp$contact_ID[temp$result])
                        ids = ids[runif(length(ids), min=0, max=100)<=tracing_rate]
                        temp$result[! temp$contact_ID %in% ids] = FALSE
                        RR$contacts_per_index[i] = length(temp$Index_ID[temp$result]) / length(unique(temp$Index_ID[temp$result]))
                        
                        # if(unique_individual){
                        #   # Find which contacts are traced into
                        #   temp <- temp %>% group_by(contact_ID) %>% mutate(traced = any(result))
                        #   # Remove duplicates of traced contacts that are not traced
                        #   temp <- temp %>% filter(!(!result & traced))
                        #   # Resolve those that are contacts of more then one index person
                        #   temp <- temp %>% group_by(contact_ID) %>% slice_min(order_by=contact_trace_date) %>% 
                        #     filter(row_number() == sample(1:n(),1))
                        #   # Make sure people are not indexes and contacts
                        #   ids_both <- unique(temp$Index_ID[!is.na(temp$Index_ID) & temp$Index_ID %in% temp$contact_ID])
                        #   both <- left_join(temp[c('Index_ID', 'Index_Diagnosed')], temp[c('contact_ID', 'contact_trace_date')], by=c('Index_ID'='contact_ID'))
                        #   both <- both %>% filter(Index_ID %in% ids_both)
                        #   # Remove those contacts that are index cases but never tracted into (I think this should remove no one)
                        #   ids_to_remove = temp$contact_ID[(temp$contact_ID %in% ids_both & is.na(temp$contact_trace_date))]
                        #   temp <- temp %>% filter( !(contact_ID %in% ids_to_remove) )
                        #   both <- both %>% filter( !(Index_ID %in% ids_to_remove) )
                        #   ids_both <- ids_both[!(ids_both %in% ids_to_remove)]
                        #   # Remove those contacts that have already been diagnosed
                        #   ids_to_remove = both$Index_ID[Index_Diagnosed <= contact_trace_date]
                        #   temp <- temp %>% filter( !(contact_ID %in% ids_to_remove) )
                        #   both <- both %>% filter( !(Index_ID %in% ids_to_remove) )
                        #   ids_both <- ids_both[!(ids_both %in% ids_to_remove)]
                        # }
                        
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
                        
                        # set confusion matrix
                        RR$FP[i] = sum( temp$result & !temp$condition)
                        RR$TN[i] = sum(!temp$result & !temp$condition)
                        RR$TP[i] = sum( temp$result &  temp$condition)
                        RR$FN[i] = sum(!temp$result &  temp$condition)
                        RR$N_HIV[i] = sum(temp$result & (temp$contact_HIV_stage != 0) & (temp$Contact_Infected <= temp$contact_trace_date))
                        RR$N_untreated[i] = sum(temp$result & (temp$contact_HIV_stage != 0) & (temp$contact_HIV_stage != 4))
                        RR$N_undiagnosed[i] = sum(temp$result & (temp$contact_HIV_stage != 0) & (temp$Contact_Diagnosed > temp$contact_trace_date), na.rm=TRUE)
                        RR$N_index[i] = length(unique(temp$Index_ID[temp$result]))
                        RR$unique_contacts_per_index[i] = length(temp$Index_ID[temp$result]) / RR$N_index[i]
                        contact_dist = temp %>% filter(result) %>% group_by(Index_ID) %>% summarise(n_contact=n())
                        dist = tabulate(contact_dist$n_contact[!is.na(contact_dist$Index_ID)])
                        if(any(is.na(contact_dist$Index_ID))){
                          dist = c(contact_dist$n_contact[is.na(contact_dist$Index_ID)], dist)  
                        } else {
                          dist = c(0, dist)
                        }
                        RR$contact_mean[i] = mean(contact_dist$n_contact[!is.na(contact_dist$Index_ID)])
                        RR$contact_std[i] = sd(contact_dist$n_contact[!is.na(contact_dist$Index_ID)])
                        i = i + 1
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  # Actual Positives
  RR['P'] = RR['TP'] + RR['FN']
  # Actual Negatives
  RR['N'] = RR['FP'] + RR['TN']
  # Predicted Positives
  RR['PP'] = RR['TP'] + RR['FP']
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
  
  # write.table(RR, file.path(directory, 'RROutput.csv'), sep=',', dec='.', row.names=FALSE)
  return(RR)
}

doIteractiveTracing <- function(directory='.', tracing_groups='Incident_HIV',
                                sample_rates=100, tracing_rates=100,
                                tracing_start_times=2010, tracing_end_times=2012,
                                acute_to_traces=FALSE, look_back_durations=3,
                                tracing_delays=0, unique_individuals=TRUE,
                                num_generations=3){
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
                 'unique_individuals',
                 'generation',
                 'FP',
                 'TN',
                 'TP',
                 'FN',
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
  i = 1
  for(tracing_start_time in tracing_start_times){
    for(tracing_end_time in tracing_end_times){
      for(tracing_delay in tracing_delays){
        for(look_back_duration in look_back_durations){
          for(acute_to_trace in acute_to_traces){
            # Do the tracing for everyone using the currently parameters that effect tracing
            individuals = getGroups(directory, tracing_start_time, tracing_end_time, acute_to_trace)
            sim_start_year = individuals$start_year
            individuals = individuals$df
            for(generation in seq(num_generations)){
              contactData = getPartners(directory, individuals, look_back_duration, sim_start_year)
              transData = getTransData(directory, contactData, tracing_start_time, tracing_delay)
              rm(contactData)
              for(comparison in comparisons){
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
                for(tracing_group in tracing_groups){
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
                  for(sample_rate in sample_rates){
                    for(tracing_rate in tracing_rates){
                      # mark the "positive diagnosis" if tracing is the "test"
                      # select those enrolled in tracing
                      ids = unique(temp$Index_ID[temp$result])
                      ids = ids[runif(length(ids), min=0, max=100)<=100*((sample_rate/100)^generation*(tracing_rate/100)^(generation-1))]
                      temp$result[! temp$Index_ID %in% ids] = FALSE
                      # select the contacts that get traced/contacted
                      ids = unique(temp$contact_ID[temp$result])
                      ids = ids[runif(length(ids), min=0, max=100)<=tracing_rate]
                      temp$result[! temp$contact_ID %in% ids] = FALSE
                      for(unique_individual in unique_individuals){
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
                          RR$unique_individuals[i] = unique_individual
                          RR$generation[i] = generation
                          RR$contacts_per_index[i] = length(temp$Index_ID[temp$result]) / length(unique(temp$Index_ID[temp$result]))
                          
                          # set confusion matrix
                          RR$FP[i] = sum( temp$result & !temp$condition)
                          RR$TN[i] = sum(!temp$result & !temp$condition)
                          RR$TP[i] = sum( temp$result &  temp$condition)
                          RR$FN[i] = sum(!temp$result &  temp$condition)
                          RR$N_HIV[i] = sum(temp$result & (temp$contact_HIV_stage != 0) & (temp$Contact_Infected <= temp$contact_trace_date))
                          RR$N_untreated[i] = sum(temp$result & (temp$contact_HIV_stage != 0) & (temp$contact_HIV_stage != 4))
                          RR$N_undiagnosed[i] = sum(temp$result & (temp$contact_HIV_stage != 0) & (temp$Contact_Diagnosed > temp$contact_trace_date), na.rm=TRUE)
                          RR$N_index[i] = length(unique(temp$Index_ID[temp$result]))
                          contact_dist = temp %>% filter(result) %>% group_by(Index_ID) %>% summarise(n_contact=n())
                          dist = tabulate(contact_dist$n_contact[!is.na(contact_dist$Index_ID)])
                          if(any(is.na(contact_dist$Index_ID))){
                            dist = c(contact_dist$n_contact[is.na(contact_dist$Index_ID)], dist)  
                          } else {
                            dist = c(0, dist)
                          }
                          RR$contact_mean[i] = mean(contact_dist$n_contact[!is.na(contact_dist$Index_ID)])
                          RR$contact_std[i] = sd(contact_dist$n_contact[!is.na(contact_dist$Index_ID)])
                          i = i + 1
                        }
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
  
  # Actual Positives
  RR['P'] = RR['TP'] + RR['FN']
  # Actual Negatives
  RR['N'] = RR['FP'] + RR['TN']
  # Predicted Positives
  RR['PP'] = RR['TP'] + RR['FP']
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
  
  RR <- RR %>% group_by(comparison,
                        trait,
                        sample_rate,
                        tracing_rate,
                        tracing_start_time,
                        tracing_end_time,
                        tracing_delay,
                        look_back_duration,
                        acute_to_trace,
                        generation) %>% mutate('RRHIV-' = PPV/PPV[tracing_group=='HIV-'])
  
  # write.table(RR, file.path(directory, 'RROutput.csv'), sep=',', dec='.', row.names=FALSE)
  return(RR)
}