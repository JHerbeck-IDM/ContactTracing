from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis
import os
from io import BytesIO
from typing import Dict, Any, Union
from idmtools.entities.ianalyzer import IAnalyzer as BaseAnalyzer

import pandas as pd
from idmtools.entities.iworkflow_item import IWorkflowItem
from idmtools.entities.simulation import Simulation

# Experiment ID for the aggregation
experiment_ID = ['70599d70-e85e-ee11-aa0b-b88303911bc1']
analysis_title = 'Aggergating experiment 1'

# Get aggregate information
agg_columns = {'N_index': ['mean', 'std'],
               'N_HIV': ['mean', 'std'],
               'N_untreated': ['mean', 'std'],
               'N_undiagnosed': ['mean', 'std'],
               'contact_mean': ['mean', 'std'],
               'contact_std': ['mean', 'std'],
               'contacts_per_index': ['mean', 'std'],
               'unique_contacts_per_index': ['mean', 'std'],
               'percent_HIV': ['mean', 'std'],
               'NNTTP': ['mean', 'std'],
               'NNTTT': ['mean', 'std'],
               'NNTTD': ['mean', 'std'],
               'FP': ['mean', 'std'],
               'TN': ['mean', 'std'],
               'TP': ['mean', 'std'],
               'FN': ['mean', 'std'],
               'P':  ['mean', 'std'],
               'N':  ['mean', 'std'], 
               'PP': ['mean', 'std'],
               'PN': ['mean', 'std'],
               'total': ['mean', 'std'],
               'prevalence': ['mean', 'std'],
               'ACC': ['mean', 'std'],
               'PPV': ['mean', 'std'],
               'FDR': ['mean', 'std'],
               'FOR': ['mean', 'std'],
               'NPV': ['mean', 'std'],
               'TPR': ['mean', 'std'],
               'FPR': ['mean', 'std'],
               'FNR': ['mean', 'std'],
               'TNR': ['mean', 'std'],
               'LRP': ['mean', 'std'],
               'LRN': ['mean', 'std'],
               'DOR': ['mean', 'std'],
               'TS':  ['mean', 'std'],
               'informedness': ['mean', 'std'],
               'markedness': ['mean', 'std'],
               'prevalenceThreshold': ['mean', 'std'],
               'MCC': ['mean', 'std'],
               'FMI': ['mean', 'std'],
               'F1':  ['mean', 'std'],
               'BA':  ['mean', 'std'],
               'RR':  ['mean', 'std'],
               'RRall': ['mean', 'std']}

file_name = 'output.csv'
class NNDistAnalyzer(BaseAnalyzer):

    def __init__(self, title='idm'):
        super().__init__(filenames=[file_name], parse=False)
        print(title)

    def initialize(self):
        """
        Initialize our Analyzer. At the moment, this just creates our output folder
        Returns:

        """
        if not os.path.exists(os.path.join(self.working_dir, "output")):
            os.mkdir(os.path.join(self.working_dir, "output"))

    def map(self, data: Dict[str, Any], item: Union[IWorkflowItem, Simulation]) -> Any:
        """
        Extracts the Statistical Population, Data channel from InsetChart.
        Called for Each WorkItem/Simulation.

        Args:
            data: Data mapping str to content of file
            item: Item to Extract Data from (Usually a Simulation)

        Returns:

        """
    
        # Get clustering by attribute information
        dataframe = pd.read_csv(BytesIO(data[file_name]))
        dataframe['sim_num'] = item.tags['sim_num']

        return dataframe

    def reduce(self, all_data: Dict[Union[IWorkflowItem, Simulation], Any]) -> Any:
        """
        Create the Final Population JSON and Plot

        Args:
            all_data: Populate data from all the Simulations

        Returns:
            None
        """

        output_dir = os.path.join(self.working_dir, "output")

        print('Creating output dataframe')
        out_data = pd.DataFrame()
        
        for key, data in all_data.items():
            out_data = pd.concat([out_data, data])
        out_data.to_csv(os.path.join(output_dir, 'aggregateData.csv'), index=False)
        out_data_mean = out_data.groupby(['tracing_group', 
                                          'comparison',
                                          'trait',
                                          'sample_rate',
                                          'tracing_rate',
                                          'tracing_start_time',
                                          'tracing_end_time',
                                          'tracing_delay',
                                          'look_back_duration',
                                          'acute_to_trace',
                                          'unique_individuals'],
                as_index=False).agg(agg_columns)
        out_data_mean.columns = out_data_mean.columns.map("_".join)
        out_data_mean.columns = out_data_mean.columns.str.rstrip('_')
        out_data_mean.to_csv(os.path.join(output_dir, 'output_summarized.csv'), index=False)
        
        return


if __name__ == "__main__":
    docker_image = "docker-staging.packages.idmod.org/greg/idmtools_greg:1.0.0"
    platform = Platform('Calculon', docker_image=docker_image)
    analysis = PlatformAnalysis(
        platform=platform, experiment_ids=experiment_ID,
        analyzers=[NNDistAnalyzer], analyzers_args=[{'title': 'idm'}],
        analysis_name=analysis_title,
        # You can pass any additional arguments needed to AnalyzerManager through the extra_args parameter
        extra_args=dict(max_workers=8)
    )

    analysis.analyze(check_status=True)
    wi = analysis.get_work_item()
    print(wi)
