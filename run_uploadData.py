# libraries to do the actually work
import glob, time, os, sys, re
import numpy as np

# libraries for running on COMPS
from COMPS import Client
from COMPS.Data import Suite, Experiment, Simulation, SimulationFile, Configuration, QueryCriteria
from COMPS.Data.Simulation import SimulationState


# Set up the COMPS environment
compshost = 'https://comps.idmod.org'
compsenv = 'Calculon'
pri = 'BelowNormal'		# Lowest, BelowNormal, Normal, AboveNormal, Highest
ac_id = '7aaf4bfe-0557-ee11-aa0a-b88303911bc1'

# Login to COMPS
Client.login(compshost)

# Set up a suite for the HIV contact tracing analysis
suite = Suite('HIV Contact Tracing Analysis')
# Set up the experiment 
experiment = Experiment('Upload Eswatini Data')
experiment.configuration = Configuration(
    environment_name=compsenv,
    working_directory_root='$COMPS_PATH(USER)',
    priority=pri,
    asset_collection_id=ac_id
    )
experiment.save()

workorder_str = """{
    "Command": "echo 'done'",
    "NumCores": 1,
    "NumNodes": 1,
    "NodeGroupName": "idm_abcd"
}"""

files_glob_to_analyze = './download/*/TransmissionReport.csv'
additional_input_filenames = ['RelationshipEnd.csv',
                              'RelationshipStart.csv',
                              'ReportEventRecorder.csv']
    
###############################################################

# Start timing
t1 = time.perf_counter()

print('Searching for files to upload')
files = glob.glob(files_glob_to_analyze,recursive=True)

# Count the number of simulations we create
numsims = 0
# Loop through all realizations of the location/scenario
for f in files:
    d = os.path.dirname(f)
    print(f'\tCreating sim for: {d}')

    s = Simulation('Uploading data')
    s.experiment_id = experiment.id
    s.add_file(SimulationFile('WorkOrder.json', 'WorkOrder'), data=bytes(workorder_str, 'utf-8'))
    s.add_file(SimulationFile(os.path.basename(f), 'input'), file_path=f)
    for af in additional_input_filenames:
        s.add_file(SimulationFile(os.path.basename(af), 'input'), file_path=os.path.join(d, af))

    s.set_tags({'simNum': numsims})

    s.save()

    numsims += 1
    if numsims % 100 == 0:
        print('Number of sims so far: ', numsims)

print(f'Created experiment {experiment.id} with {numsims} simulations')

if numsims == 0:
    print('No sims created.  Exiting.')
else: 
    print('Commissioning...')
    
    experiment.commission()
    
    print('Commissioned')
    
    print("elapsed time: " + str(time.perf_counter() - t1))
