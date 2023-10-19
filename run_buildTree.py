import os
import sys
from COMPS import Client
from COMPS.Data import Experiment, Simulation, Configuration, SimulationFile, QueryCriteria

## All these settings should probably be moved to a yaml file and read in
# The experiment IDs where the trees can be found
data_id  = '9dd0912f-3e69-ee11-aa0b-b88303911bc1'
# The experiment IDs of existing experiments to add simulations too
# exp_id  = ''
# Input parameters
tracing_groups      = ['Incident_HIV', 'HIV+', 'HIV-']
n_samples           = [500, 1000, 1500]
tracing_start_times = [2000, 2005, 2010, 2015]
tracing_durations   = [5]
look_back_windows   = [3]
tracing_delays      = [0]
acute_to_traces     = True
thresholds          = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

# Number cores to use for each simulation
cores = 2

# COMPS settings
compshost = 'https://comps.idmod.org'
compsenv = 'Calculon'
pri = 'Normal'		# Lowest, BelowNormal, Normal, AboveNormal, Highest
ac_id = '26f3cc4f-596c-ee11-aa0b-b88303911bc1'


# Need to make lists strings for passing to the scripts actually doing the work
tracingGroups_str    = ' '.join(str(tracing_group) for tracing_group in tracing_groups)
thresholds_str       = ' '.join(str(threshold) for threshold in thresholds)
nSamples_str         = ' '.join(str(n_sample) for n_sample in n_samples)
traceStartTimes_str  = ' '.join(str(traceStartTime) for traceStartTime in tracing_start_times)
tracingDurations_str = ' '.join(str(tracing_duration) for tracing_duration in tracing_durations)
tracingDelays_str    = ' '.join(str(tracingDelay) for tracingDelay in tracing_delays)
lookBackWindows_str  = ' '.join(str(lookBackWindow) for lookBackWindow in look_back_windows)

exp_name = 'HIV contact tracing clustering  (' + str(tracing_groups) + \
    ', num samples: ' + str(n_samples) + ', tracing start times: ' + \
    str(tracing_start_times) + ', tracing durations: ' + str(tracing_durations) + \
    ', tracing delay: ' + str(tracing_delays) + ', look back window: ' + \
        str(look_back_windows) + ')'

# This is the actual command to be run on COMPS
workorder_str = f"""{{
    "Command": "singularity exec Assets/HIV_contact_tracing.sif python3 buildTrees.py -G {tracingGroups_str} -S {nSamples_str} -C {thresholds_str} -B {traceStartTimes_str} -L {tracingDurations_str} -W {lookBackWindows_str} -D {tracingDelays_str} -A",
    "NumCores": 2,
    "NumNodes": 1,
    "NodeGroupName": "idm_abcd"
}}"""

Client.login(compshost)

print('Getting Data')

# Seting up the path for reading data from COMPS
if 'COMPS_DATA_MAPPING' not in os.environ:
    os.environ['COMPS_DATA_MAPPING'] = '/mnt/idm2;\\\\internal.idm.ctr\\IDM2'

# Get tree data from COMPS experiment
data_exp = Experiment.get(data_id)
data_sims = data_exp.get_simulations(QueryCriteria().select_children(['hpc_jobs']))

# Add to or create the experiment
if 'exp_id' in locals():
    processing_exp = Experiment.get(exp_id)
    # Should add code here to update experiment tags to reflect what is being added
else:
    processing_exp = Experiment(exp_name)
    processing_exp.configuration = Configuration(
        environment_name=compsenv,
        working_directory_root='$COMPS_PATH(USER)/output/',
        priority=pri,
        asset_collection_id=ac_id)
    processing_exp.set_tags({'data_exp_id': data_exp.id, 'traceGroup': tracing_groups, 'nSamples': n_samples, 'threshold': thresholds, 'traceStartTime': tracing_start_times, 'traceDuration': tracing_durations, 'tracingDelay': tracing_delays, 'lookBackWindow': look_back_windows, 'acuteToTrace': acute_to_traces})
    processing_exp.save()

# Create the simulations within the experiment
for i, sim in enumerate(data_sims):
    proc_sim = Simulation('HIV Contact Tracing Clustering')
    proc_sim.experiment_id = processing_exp.id
    proc_sim.set_tags({'buildTree_sim_id': sim.id, 'data_exp_id': data_exp.id, 'traceGroup': tracing_groups, 'nSamples': n_samples, 'threshold': thresholds, 'traceStartTime': tracing_start_times, 'traceDuration': tracing_durations, 'tracingDelay': tracing_delays, 'lookBackWindow': look_back_windows, 'acuteToTrace': acute_to_traces, 'sim_num': i})
    proc_sim.add_file(SimulationFile('WorkOrder.json', 'WorkOrder'), data=bytes(workorder_str, 'utf-8'))
    proc_sim.add_file(SimulationFile('buildTrees.py', 'input'), file_path='buildTrees.py')
    proc_sim.add_file(SimulationFile('HIVContactTracing.r', 'input'), file_path='HIVContactTracing.r')

    working_dir = sim.hpc_jobs[-1].working_directory
    print(working_dir)
    
    # could just save working directory instead of specifically path to .npy if there's other files in that dir you intent to read as well
    proc_sim.add_file(SimulationFile('data_path.txt', 'input'), data=bytes(working_dir, 'utf-8'))

# Send the simulations to COMPS
Simulation.save_all()

# Run the simulations
processing_exp.commission()

# Clean up if adding to an existing experiment
if 'exp_id' in locals():
    del exp_id
