# Example DownloadAnalyzer for EMOD Experiment
# In this example, we will demonstrate how to create an DownloadAnalyzer to download simulation output files locally

# First, import some necessary system and idmtools packages.
from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.analysis.download_analyzer import DownloadAnalyzer
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform
from COMPS.Data import Experiment, WorkItem, WorkItemFile, QueryCriteria

import os
import pandas as pd

if __name__ == '__main__':

    # Set the platform where you want to run your analysis
    exp_id = 'ea624dcb-0249-ee11-aa0a-b88303911bc1'  # comps exp id
    with Platform('IDMcloud') as platform:
        sims = Experiment.get(exp_id).get_simulations(query_criteria=QueryCriteria().select_children(['hpc_jobs']))

        # Arg option for analyzer init are uid, working_dir, data in the method map (aka select_simulation_data),
        # and filenames
        # In this case, we want to provide a filename to analyze
        filenames = ['output/TransmissionReport.csv',
                     'output/ReportEventRecorder.csv',
                     'output/RelationshipStart.csv',
                     'output/RelationshipEnd.csv']
        # Initialize the analyser class with the path of the output files to download
        analyzers = [DownloadAnalyzer(filenames=filenames, output_path='download')]

        directories = next(os.walk('download'))[1]
        
        for sim in sims:
            # Set the experiment you want to analyze
            sim_id = str(sim.id)
    
            if not sim_id in directories:
                # Specify the id Type, in this case an Experiment
                manager = AnalyzeManager(ids=[(sim_id, ItemType.SIMULATION)],
                                         analyzers=analyzers)
                manager.analyze()
                
                # Reduce transmission data
                data = pd.read_csv(os.path.join('download', sim_id, 'TransmissionReport.csv'))
                data['SRC_RISK'] = data['SRC_IP'].str.extract(r'(LOW|MEDIUM|HIGH)')
                data = data[['YEAR', 'REL_ID', 'SRC_ID', 'DEST_ID', 'SRC_STAGE', 'SRC_RISK']]
                data.to_csv(os.path.join('download', sim_id, 'TransmissionReport.csv'),
                                 index=False)
                
                # Reduce relationship start data
                data = pd.read_csv(os.path.join('download', sim_id, 'RelationshipStart.csv'))
                data = data[['Rel_ID', 'Rel_start_time', 'A_ID', 'B_ID',
                            'A_IndividualProperties', 'B_IndividualProperties',
                            'A_HIV_disease_stage', 'B_HIV_disease_stage']]
                data['A_IndividualProperties'] = data['A_IndividualProperties'].str.extract(r'(LOW|MEDIUM|HIGH)')
                data['B_IndividualProperties'] = data['B_IndividualProperties'].str.extract(r'(LOW|MEDIUM|HIGH)')
                data.to_csv(os.path.join('download', sim_id, 'RelationshipStart.csv'),
                                index=False)
                
                # Reduce relationship start data
                data = pd.read_csv(os.path.join('download', sim_id, 'RelationshipEnd.csv'))
                data = data[['Rel_ID', 'Rel_actual_end_time',
                             'A_IndividualProperties', 'B_IndividualProperties',
                             'A_HIV_disease_stage', 'B_HIV_disease_stage']]
                data['A_IndividualProperties'] = data['A_IndividualProperties'].str.extract(r'(LOW|MEDIUM|HIGH)')
                data['B_IndividualProperties'] = data['B_IndividualProperties'].str.extract(r'(LOW|MEDIUM|HIGH)')
                data.to_csv(os.path.join('download', sim_id, 'RelationshipEnd.csv'),
                                index=False)
