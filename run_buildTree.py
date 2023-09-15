import os
import sys
from COMPS import Client
from COMPS.Data import Experiment, Simulation, Configuration, SimulationFile, QueryCriteria

## All these settings should probably be moved to a yaml file and read in
# Locations/scenarios to run on TODO rename simulation
simulations = ['Atlanta', 'Eswatini']
# The experiment IDs where the trees can be found
data_id  = 'ea624dcb-0249-ee11-aa0a-b88303911bc1'
# The experiment IDs of existing experiments to add simulations too
# exp_id  = ''
# Input parameters
traceGroup = 'Incident_HIV'
sampleRates = [10, 50, 100]
traceRates = [10, 50, 100]
thresholds = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
traceStartTimes = [2010]
traceEndTimes = [2012]
lookBackWindows = [3]
tracingDelays = [0]
acuteToTrace = False

# Number cores to use for each simulation
cores = 2

# COMPS settings
compshost = 'https://comps.idmod.org'
compsenv = 'Calculon'
pri = 'BelowNormal'		# Lowest, BelowNormal, Normal, AboveNormal, Highest
ac_id = '876e396c-5953-ee11-aa0a-b88303911bc1'


# Need to make lists strings for passing to the scripts actually doing the work
thresholds_str      = ' '.join(str(threshold) for threshold in thresholds)
sampleRates_str     = ' '.join(str(sampleRate) for sampleRate in sampleRates)
traceRates_str      = ' '.join(str(traceRate) for traceRate in traceRates)
traceStartTimes_str = ' '.join(str(traceStartTime) for traceStartTime in traceStartTimes)
traceEndTimes_str   = ' '.join(str(traceEndTime) for traceEndTime in traceEndTimes)
tracingDelays_str   = ' '.join(str(tracingDelay) for tracingDelay in tracingDelays)
lookBackWindows_str = ' '.join(str(lookBackWindow) for lookBackWindow in lookBackWindows)

exp_name = 'HIV contact tracing analysis  (' + str(traceGroup) + ', sample rates: ' + \
 str(sampleRates) + ', tracing rates: ' + str(traceRates) + ', tracing start times: ' + \
 str(traceStartTimes) + ', tracing end times: ' + str(traceEndTimes) + ', tracing delay: ' + \
 str(tracingDelays) + ', look back window: ' + str(lookBackWindows) + ')'

# This is the actual command to be run on COMPS
workorder_str = f"""{{
    "Command": "singularity exec Assets/tree_stats_analysis.sif python3 buildTrees.py -G {traceGroup} -S {sampleRates_str} -T {traceRates_str} -C {thresholds_str} -B {traceStartTimes_str} -E {traceEndTimes_str} -W {lookBackWindows_str} -D {tracingDelays_str}",
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
    processing_exp.set_tags({'data_exp_id': data_exp.id, 'traceGroup': traceGroup, 'sampleRate': sampleRates, 'traceRate': traceRates, 'threshold': thresholds, 'traceStartTime': traceStartTimes, 'traceEndTime': traceEndTimes, 'tracingDelay': tracingDelays, 'lookBackWindow': lookBackWindows, 'acuteToTrace': acuteToTrace})
    processing_exp.save()

# Create the simulations within the experiment
for i, sim in enumerate(data_sims):
    proc_sim = Simulation('HIV Contact Tracing')
    proc_sim.experiment_id = processing_exp.id
    proc_sim.set_tags({'buildTree_sim_id': sim.id, 'data_exp_id': data_exp.id, 'traceGroup': traceGroup, 'sampleRate': sampleRates, 'traceRate': traceRates, 'threshold': thresholds, 'traceStartTime': traceStartTimes, 'traceEndTime': traceEndTimes, 'tracingDelay': tracingDelays, 'lookBackWindow': lookBackWindows, 'acuteToTrace': acuteToTrace, 'sim_num': i})
    proc_sim.add_file(SimulationFile('WorkOrder.json', 'WorkOrder'), data=bytes(workorder_str, 'utf-8'))
    proc_sim.add_file(SimulationFile('buildTrees.py', 'input'), file_path='buildTrees.py')
    proc_sim.add_file(SimulationFile('HIVContactTracing.r', 'input'), file_path='HIVContactTracing.r')

    working_dir = sim.hpc_jobs[-1].working_directory
    print(working_dir)
    
    # could just save working directory instead of specifically path to .npy if there's other files in that dir you intent to read as well
    # proc_sim.add_file(SimulationFile('tree_path.txt', 'input'), data=bytes(working_dir, 'utf-8'))

# Send the simulations to COMPS
Simulation.save_all()

# Run the simulations
processing_exp.commission()

# Clean up if adding to an existing experiment
if 'exp_id' in locals():
    del exp_id
