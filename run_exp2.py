import os
import sys
from COMPS import Client
from COMPS.Data import Experiment, Simulation, Configuration, SimulationFile, QueryCriteria

## All these settings should probably be moved to a yaml file and read in
# The experiment IDs where the trees can be found
data_id  = '6bfa88f0-b45c-ee11-aa0a-b88303911bc1'
# The experiment IDs of existing experiments to add simulations too
# exp_id  = ''
# # Input parameters
# traceGroup = 'Incident_HIV'
# sampleRates = [10, 50, 100]
# traceRates = [10, 50, 100]
# traceStartTimes = [2010]
# traceEndTimes = [2012]
# lookBackWindows = [3]
# tracingDelays = [0]
# acuteToTrace = False

# Number cores to use for each simulation
cores = 2

# COMPS settings
compshost = 'https://comps.idmod.org'
compsenv = 'Calculon'
pri = 'Normal'		# Lowest, BelowNormal, Normal, AboveNormal, Highest
ac_id = '7aaf4bfe-0557-ee11-aa0a-b88303911bc1'


# Need to make lists strings for passing to the scripts actually doing the work
tracing_groups='HIV+, Incident_HIV, HIV-, Prevalent_HIV, Undiagnosed_HIV'
comparisons = 'all, traceable, sampleableAtTrace, sampleableEventually'
sample_rates='100, 75, 50, 25'
tracing_rates='100, 75, 50, 25'
tracing_start_times='2010'
tracing_end_times='2012'
acute_to_traces='TRUE'
look_back_durations='3'
tracing_delays='0'
traits = 'trans_source, yet_to_transmit, recently_infected, high_risk, not_yet_infected'
unique_individuals = 'FALSE'
num_generations=10

exp_name = 'HIV contact tracing iteration analysis  (' + tracing_groups + \
     ', sample rates: ' + sample_rates + ', tracing rates: ' + tracing_rates + \
     ', tracing start times: ' + tracing_start_times + ', tracing end times: ' + \
     tracing_end_times + ')'

# This is the actual command to be run on COMPS
workorder_str = f"""{{
    "Command": "singularity exec Assets/tree_stats_analysis.sif Rscript exp2.r",
    "NumCores": 1,
    "NumNodes": 1,
    "NodeGroupName": "idm_abcd"
}}"""
# workorder_str = f"""{{
#     "Command": "singularity exec Assets/tree_stats_analysis.sif Rscript exp2.r",
# }}"""

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
    processing_exp.set_tags({'data_exp_id': data_exp.id, 'traceGroup': tracing_groups, 'sampleRate': sample_rates, 'traceRate': tracing_rates, 'traceStartTime': tracing_start_times, 'traceEndTime': tracing_end_times, 'tracingDelay': tracing_delays, 'lookBackWindow': look_back_durations, 'acuteToTrace': acute_to_traces, 'numGeneration': num_generations})
    processing_exp.save()

# Create the simulations within the experiment
for i, sim in enumerate(data_sims):
    proc_sim = Simulation('HIV Contact Tracing (exp 2)')
    proc_sim.experiment_id = processing_exp.id
    proc_sim.set_tags({'uploadData_sim_id': sim.id, 'data_exp_id': data_exp.id, 'traceGroup': tracing_groups, 'sampleRate': sample_rates, 'traceRate': tracing_rates, 'traceStartTime': tracing_start_times, 'traceEndTime': tracing_end_times, 'tracingDelay': tracing_delays, 'lookBackWindow': look_back_durations, 'acuteToTrace': acute_to_traces, 'numGeneration': num_generations, 'sim_num': i})
    proc_sim.add_file(SimulationFile('WorkOrder.json', 'WorkOrder'), data=bytes(workorder_str, 'utf-8'))
    # proc_sim.add_file(SimulationFile('buildTrees.py', 'input'), file_path='buildTrees.py')
    proc_sim.add_file(SimulationFile('HIVContactTracing.r', 'input'), file_path='HIVContactTracing.r')
    proc_sim.add_file(SimulationFile('exp2.r', 'input'), file_path='exp2.r')

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
