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
experiment_ID = ['be786418-a16e-ee11-aa0b-b88303911bc1']
analysis_title = 'Aggergating cluster stats and enrichment'

# Get aggregate information
agg_columns = {'FP': [nanmean, nanstd],
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
               'RRrandom': [nanmean, nanstd]}

file_names = ['enrichmentData.csv', 'clusterData.csv']
class NNDistAnalyzer(BaseAnalyzer):

    def __init__(self, title='idm'):
        super().__init__(filenames=file_names, parse=False)
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
        enrichmentData = pd.read_csv(BytesIO(data[file_names[0]]))
        enrichmentData['sim_num'] = item.tags['sim_num']
        clusterStats = pd.read_csv(BytesIO(data[file_names[1]]))
        clusterStats['sim_num'] = item.tags['sim_num']

        return enrichmentData, clusterStats

    def reduce(self, all_data: Dict[Union[IWorkflowItem, Simulation], Any]) -> Any:
        """
        Create the Final Population JSON and Plot

        Args:
            all_data: Populate data from all the Simulations

        Returns:
            None
        """

        output_dir = os.path.join(self.working_dir, "output")

        print('Creating enrichment dataframe')
        out_data = pd.DataFrame()
        
        for key, data in all_data.items():
            out_data = pd.concat([out_data, data[0]])
        out_data.to_csv(os.path.join(output_dir, 'aggregateEnrichmentData.csv'), index=False)
        out_data_mean = out_data.groupby(['tracing_group',
                                          'threshold',
                                          'n_samples',
                                          'tracing_start_time',
                                          'tracing_duration',
                                          'tracing_delay',
                                          'look_back_window',
                                          'acute_to_trace',
                                          'trait'],
                as_index=False).agg(agg_columns)
        out_data_mean.columns = out_data_mean.columns.map("_".join)
        out_data_mean.columns = out_data_mean.columns.str.rstrip('_')
        out_data_mean.columns = out_data_mean.columns.str.replace('nanmean', 'mean')
        out_data_mean.columns = out_data_mean.columns.str.replace('nanstd', 'std')
        out_data_mean.to_csv(os.path.join(output_dir, 'enrichmentSummarized.csv'), index=False)
        
        print('Creating cluster stats dataframe')
        out_data = pd.DataFrame()
        
        for key, data in all_data.items():
            out_data = pd.concat([out_data, data[1]])
        out_data.to_csv(os.path.join(output_dir, 'aggregateClusterData.csv'), index=False)
        out_data_mean = out_data.groupby(['tracing_group',
                                          'threshold',
                                          'n_samples',
                                          'tracing_start_time',
                                          'tracing_duration',
                                          'tracing_delay',
                                          'look_back_window',
                                          'acute_to_trace'],
                as_index=False).agg({'percent_clustered': [nanmean, nanstd],
                                     'num_clusters': [nanmean, nanstd],
                                     'max_size': [nanmean, nanstd],
                                     'median_size': [nanmean, nanstd],
                                     'mean_size': [nanmean, nanstd]})
        out_data_mean.columns = out_data_mean.columns.map("_".join)
        out_data_mean.columns = out_data_mean.columns.str.rstrip('_')
        out_data_mean.columns = out_data_mean.columns.str.replace('nanmean', 'mean')
        out_data_mean.columns = out_data_mean.columns.str.replace('nanstd', 'std')
        out_data_mean.to_csv(os.path.join(output_dir, 'clusterSummarized.csv'), index=False)
        
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
