from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis
import os
from io import BytesIO
from typing import Dict, Any, Union
from idmtools.entities.ianalyzer import IAnalyzer as BaseAnalyzer

import pandas as pd
from numpy import nanmean
from numpy import nanstd
from idmtools.entities.iworkflow_item import IWorkflowItem
from idmtools.entities.simulation import Simulation

# Experiment ID for the aggregation
experiment_ID = ['39705271-4262-ee11-aa0b-b88303911bc1']
analysis_title = 'Aggergating experiment 2'

# Get aggregate information
agg_columns = {'N_index': [nanmean, nanstd],
               'N_HIV': [nanmean, nanstd],
               'N_untreated': [nanmean, nanstd],
               'N_undiagnosed': [nanmean, nanstd],
               'contact_mean': [nanmean, nanstd],
               'contact_std': [nanmean, nanstd],
               'contacts_per_index': [nanmean, nanstd],
               'percent_HIV': [nanmean, nanstd],
               'NNTTP': [nanmean, nanstd],
               'NNTTT': [nanmean, nanstd],
               'NNTTD': [nanmean, nanstd],
               'FP': [nanmean, nanstd],
               'TN': [nanmean, nanstd],
               'TP': [nanmean, nanstd],
               'FN': [nanmean, nanstd],
               'P':  [nanmean, nanstd],
               'N':  [nanmean, nanstd], 
               'PP': [nanmean, nanstd],
               'PN': [nanmean, nanstd],
               'total': [nanmean, nanstd],
               'prevalence': [nanmean, nanstd],
               'ACC': [nanmean, nanstd],
               'PPV': [nanmean, nanstd],
               'FDR': [nanmean, nanstd],
               'FOR': [nanmean, nanstd],
               'NPV': [nanmean, nanstd],
               'TPR': [nanmean, nanstd],
               'FPR': [nanmean, nanstd],
               'FNR': [nanmean, nanstd],
               'TNR': [nanmean, nanstd],
               'LRP': [nanmean, nanstd],
               'LRN': [nanmean, nanstd],
               'DOR': [nanmean, nanstd],
               'TS':  [nanmean, nanstd],
               'informedness': [nanmean, nanstd],
               'markedness': [nanmean, nanstd],
               'prevalenceThreshold': [nanmean, nanstd],
               'MCC': [nanmean, nanstd],
               'FMI': [nanmean, nanstd],
               'F1':  [nanmean, nanstd],
               'BA':  [nanmean, nanstd],
               'RR':  [nanmean, nanstd],
               'RRall': [nanmean, nanstd],
               'RRHIV-': [nanmean, nanstd]}

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
                                          'unique_individuals',
                                          'generation'],
                as_index=False).agg(agg_columns)
        out_data_mean.columns = out_data_mean.columns.map("_".join)
        out_data_mean.columns = out_data_mean.columns.str.rstrip('_')
        out_data_mean.columns = out_data_mean.columns.str.replace('nanmean', 'mean')
        out_data_mean.columns = out_data_mean.columns.str.replace('nanstd', 'std')
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
