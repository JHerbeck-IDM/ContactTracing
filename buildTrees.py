import os
import argparse
import numpy as np
import pandas as pd
from time import perf_counter
from phylomodels.trees import generate_treeFromFile
from phylomodels.trees.transform_joinTrees import transform_joinTrees
from phylomodels.trees.transform_transToPhyloTree import transform_transToPhyloTree
from scipy.sparse.csgraph import connected_components


# R-related packages
import rpy2.robjects as robjects
r = robjects.r

from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# We may need to install some packages
try:
    from rpy2.robjects.packages import importr
    tidyverse = importr('tidyverse')
except RRuntimeError:
    from rpy2.robjects.packages import importr
    utils = importr('utils')
    base = importr('base')
    utils.chooseCRANmirror()
    utils.install_packages('tidyverse')
    
ID = 'DEST_ID'
SRC_ID = 'SRC_ID'

def get_arguments():
    """
    This function parses the input arguments for the rest of the script to use.
    See the information in the add_argument commands or README for more information.
    """

    # Set up the parser and let it know what arguments to expect
    parser = argparse.ArgumentParser(
        prog='Get arguments for contact tracing cluster analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Argument parser for all contact tracing cluster analysis scripts')
    parser.add_argument('--traceGroup', '-G', metavar='TRACEGROUP',
                        help='Choose group (Incident_HIV, Prevalent_HIV, Undiagnosed_HIV, HIV-) to use as index individuals for tracing',
                        default=['Incident_HIV'], nargs='*', type=str)
    parser.add_argument('--numSamples', '-S', metavar='NUMSAMPLES',
                        help='Maximum number of samples to gather',
                        default=[500, 1000], nargs='*', type=float)
    parser.add_argument('--threshold', '-C', metavar='THRESHOLD',
                        help='Threshold to use in defining clusters can be list of values',
                        default=[0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], nargs='*', type=float)
    parser.add_argument('--traceStartTime', '-B', metavar='TRACESTARTTIME',
                        help='The time that the tracing is to start',
                        default=[2000,2005,2010,2015], nargs='*', type=float)
    parser.add_argument('--traceLength', '-L', metavar='TRACELENGTH',
                        help='How long the tracing campaign lasts',
                        default=[5], nargs='*', type=float)
    parser.add_argument('--lookBackWindow', '-W', metavar='LOOKBACKWINDOW',
                        help='The far back from the time of diagnosis the tracing goes',
                        default=[3], nargs='*', type=float)
    parser.add_argument('--tracingDelay', '-D', metavar='TRACINGDELAY',
                        help='How long it takes to trace and test contacts',
                        default=[0], nargs='*', type=float)    
    parser.add_argument('--acuteToTrace', '-A', #metavar='ACUTETOTRACE',
                        help='Only trace those that are acute at diagonsis',
                        action='store_true')#argparse.BooleanOptionalAction)

    input_args = parser.parse_args()

    # Get the arguments as individual variables
    tracing_group      = input_args.traceGroup
    n_samples          = input_args.numSamples
    tracing_start_time = input_args.traceStartTime
    tracing_duration   = input_args.traceLength
    look_back_window   = input_args.lookBackWindow
    tracing_delay      = input_args.tracingDelay
    acute_to_trace     = input_args.acuteToTrace
    thresholds         = input_args.threshold

    # Get the arguements into the right types/formats
    if acute_to_trace:
        acute_to_trace = [True]
    else:
        acute_to_trace = [False]
    n_samples           = list(map(int, n_samples))
    tracing_start_times = list(map(int, tracing_start_time))
    tracing_duration    = list(map(int, tracing_duration))
    look_back_windows   = list(map(int, look_back_window))
    tracing_delays      = list(map(int, tracing_delay))
    thresholds          = list(map(int, thresholds))

    # Print out arguments for reference if needed
    # if any(tracing_group not in ['Incident_HIV', 'Prevalent_HIV', 'Undiagnosed_HIV', 'HIV+', 'HIV-']):
    if 0 > len(set(['Incident_HIV', 'Prevalent_HIV', 'Undiagnosed_HIV', 'HIV+', 'HIV-']).difference(tracing_group)):
        print('ERROR: tracing group "', tracing_group, '" not known.')
    print('Tracing Group:', tracing_group)
    print('Number of Samples:', n_samples)
    print('Tracing Start Time:', tracing_start_time)
    print('Tracing Duration:', tracing_duration)
    print('Look Back Window:', look_back_window)
    print('Tracing Delay:', tracing_delay)
    print('Acute to Trace:', acute_to_trace)
    print('Threshold:', thresholds)


    # Return a dictionary with all the needed variables
    return {'tracing_group': tracing_group,
            'n_samples': n_samples,
            'tracing_start_time': tracing_start_times,
            'tracing_duration': tracing_duration,
            'look_back_window': look_back_windows,
            'tracing_delay': tracing_delays,
            'acute_to_trace': acute_to_trace,
            'thresholds': thresholds}

def walk_down (nodes, curnode, pathlen, topology_only=False):
    """
    Recursive function for walking down the tree and getting the cumulative
    distance from the original nodes to all its descendents.

    Args:
        nodes (list)        : A (pontentially empty) list of tuples where each
                              tuple contains a node name and its distance to
                              the original node.
        curnode (treeNode)  : The node we are currently working on
        pathlen (float)     : The cumulative path length so far (this this
                              node's parent). Because this is suppose to be
                              from teh node's parent and this function adds on
                              the distance from that parent to this current
                              node, if the current node is the node we started
                              with we need to pass in the negative distance
                              from the current node to its parent so that after
                              this function the distance is reported as 0.
        topology_only (bool): Whether to use branch lengths or number of branches
    Returns:
        list                : The updated list of tuples that was passed in
                              (is the nodes input arguement)
    """

    if topology_only:
        # Add the "distance" from the parent to the current node
        pathlen += 1
        # Add the node name/distance pair to the list. If the current node is
        # not a leaf node recurse through the nodes children
        if curnode.is_leaf():
            nodes.append( (curnode.name, pathlen) )
        else:
            nodes.append( (curnode.name, pathlen) )
            for c in curnode.children:
                nodes = walk_down(nodes, c, pathlen, topology_only)
    else:
        # Add the distance from the parent to the current node
        pathlen += curnode.dist
        # Add the node name/distance pair to the list. If the current node is
        # not a leaf node recurse through the nodes children
        if curnode.is_leaf():
            nodes.append( (curnode.name, pathlen) )
        else:
            nodes.append( (curnode.name, pathlen) )
            for c in curnode.children:
                nodes = walk_down(nodes, c, pathlen, topology_only)

    return nodes
    
def find_short_edges(tree, topology_only=False):
    """
    Create the distance matrix for the tree as well as a list of the node names
    in the same order as the matrix. This list of names allows for the selecting
    specific rows and columns out of the matrix later.

    Args:
        tree (ete3 Tree)    : The phylogenetic tree to transform into a
                              distance matrix
        topology_only (bool): Whether to calculate the distance with branch
                              lengths or number of branches
    Returns:
        2darray, list       : Return the distance matrix and the list of node
                              names
    """

    # generate dictionary of child->parent associations in level order to
    # allow for a faster calculation of the distances from nodes on one side
    # of the tree and the other
    parents = {}
    names = []
    for clade in tree.traverse('levelorder'):
        names.append(clade.name)
        for child in clade.children:
            parents.update({child: clade})

    # Get the list of nodes (in levelorder) and preallocate the matrix
    nodes = tree.get_descendants('levelorder')
    num_nodes = len(nodes)
    distance = np.zeros((num_nodes+1, num_nodes+1), dtype=np.float32)
    # Create an dictionary to index between the list of names and actual tree
    index = {tree.name: 0}
    count = 1
    for node in nodes:
        index.update({node.name: count})
        count += 1
        
    # Walk down the tree getting the distance from the root to all the other nodes
    node2 = walk_down([], tree, {True: -1, False: -tree.dist}[topology_only], topology_only)
    
    # Loop through the list of nodes and their distances to the root and add
    # those distances to the matrix in the correct spot
    i = index[tree.name]
    for n2, dist in node2:
        j = index[n2]
        distance[i, j] = dist
    distance[i, i] = 0
    
    # Loop through the rest of the nodes to get their distances to each other
    # Here we are traversing the tree in level order so we can using the
    # distances we already know
    for node in nodes:
        i = index[node.name]
        j = index[parents[node].name]
        dist = {True: 1, False: node.dist}[topology_only]
        # Set the distances for the current node to be the same as the
        # distances of the parent plus the distance between the parent and
        # child. This will get all the correct distances except for the current
        # node to itself and to its descendents.
        distance[i, :] = distance[j, :] + dist
        # Now walk down the tree to get the distances to the descendents correct
        node2 = walk_down([], node, {True: -1, False: -node.dist}[topology_only], topology_only)
        
        # Loop through the list of descents and add their distances to the
        # matrix in the correct spot
        for n2, dist in node2:
            j = index[n2]
            distance[i, j] = dist
        distance[i, i] = 0
            
    return distance, names

def get_cluster_stats(metadata, threshold=5, n_samples=1000,
                      tracing_start_time=2010, tracing_duration=5,
                      tracing_delay=0, tracing_group='Incident_HIV',
                      look_back_window=0, acute_to_trace=False):
    """
    This function calculates basic cluster statistics (number of clusters,
    percent clustered, etc). Many of the inputs are given just to be added to
    the out data. Probably should rework how the data is aggregated.

    Arguments:
        metadata (DataFrame)      : This data frame contains the information on each
                                    individual including what cluster they are in
        threshold (float)         : The threshold used for determining clusters
        tracing_group (str)       : Which group/subset we are using for index
                                    individuals
        N_samples (int)           : Maximum number of samples
        tracing_start_time (float): The date when the period for contract
                                    tracing started.
        tracing_duration (float)  : The date when the period for contract
                                    tracing ended.
        acute_to_trace (bool)     : This flags whether the incidence infections
                                    (those whose contacts will be traced) are
                                    just those diagnosed during the tracing
                                    period or if they also need to be acute at
                                    the time of diagnosis.
        look_back_duration (float): How far into the past individuals are asked
                                    to remember contacts from
        tracing_delay (float)     : How long the tracing takes which effects
                                    the sample time

    Return
        DataFrame                 : The data frame with the results looking at the
                                    full sampling period
    """

    # Pre-allocate the data frame to hold the out put
    data = pd.DataFrame(columns=['threshold', 'n_samples',
                                 'tracing_start_time', 'tracing_duration',
                                 'tracing_delay', 'tracing_group',
                                 'look_back_window', 'acute_to_trace',
                                 'percent_clustered', 'num_clusters',
                                 'max_size', 'median_size', 'mean_size'],
                        index=range(1))
    # Set data types as needed
    data['max_size'] = data['max_size'].astype('float')
    data['mean_size'] = data['mean_size'].astype('float')
    data['median_size'] = data['median_size'].astype('float')
    # Set values passed in
    data['threshold'] = threshold
    data['n_samples'] = n_samples
    data['tracing_start_time'] = tracing_start_time
    data['tracing_duration'] = tracing_duration
    data['tracing_delay'] = tracing_delay
    data['tracing_group'] = tracing_group
    data['look_back_window'] = look_back_window
    data['acute_to_trace'] = acute_to_trace

    # Calculate the percent of individuals in clusters
    metadata['cluster'] = (metadata['cluster_size']>1).astype(bool)
    num_clustered = np.sum(metadata['cluster']==1)
    percent_clustered = num_clustered/metadata.shape[0]
    data['percent_clustered'] = percent_clustered*100

    # Calculate the number of clusters (size >= 2)
    hist, _ = np.histogram(metadata['cluster_size'], [2, 3, 4, 5, 6, 10, 55000], weights = 1/metadata['cluster_size'])
    data['num_clusters'] = np.sum(hist)

    # Calculate the 'middle' of the cluster size distribution
    temp = metadata.loc[metadata['cluster_size']>1,:]
    cluster_sizes = temp.groupby('cluster_name')['cluster_size'].agg('unique').astype(int)
    data['median_size'] = cluster_sizes.median()
    data['mean_size'] = cluster_sizes.mean()
    data['max_size'] = cluster_sizes.max()

    return data

def get_cluster_enrichment(metadata, threshold=5, n_samples=1000,
                           tracing_start_time=2010, tracing_duration=5,
                           tracing_delay=0, tracing_group='Incident_HIV',
                           look_back_window=0, acute_to_trace=False,
                           trans_source=True):
    """
    This function calculates if clusters are 'enriched' with transmitters. To
    do this we calculate all the metrics for evaulating a binary test where we
    consider transmitting to be a true positive and clustering to be a positive
    test result. Many of the inputs are giving just to be added to the out
    data. Probably should rework how the data is aggregated.

    Arguments:
        metadata (DataFrame)  : This data frame contains the information on each
                                individual including what cluster they are in
        threshold (float)         : The threshold used for determining clusters
        tracing_group (str)       : Which group/subset we are using for index
                                    individuals
        N_samples (int)           : Maximum number of samples
        tracing_start_time (float): The date when the period for contract
                                    tracing started.
        tracing_duration (float)  : The date when the period for contract
                                    tracing ended.
        acute_to_trace (bool)     : This flags whether the incidence infections
                                    (those whose contacts will be traced) are
                                    just those diagnosed during the tracing
                                    period or if they also need to be acute at
                                    the time of diagnosis.
        look_back_duration (float): How far into the past individuals are asked
                                    to remember contacts from
        tracing_delay (float)     : How long the tracing takes which effects
                                    the sample time
        trans_source (bool)       : Where the trait is being a transmission
                                    source or high risk

    Return
        DataFrame             : The data frame with the results looking at the
                                full sampling period

    """

    # It was breaking when I started with an empty data frame so here we create
    # a data frame with one row and column
    data = pd.DataFrame({'FP': [0]})

    # Set the values we already know
    data['threshold'] = threshold
    data['n_samples'] = n_samples
    data['tracing_start_time'] = tracing_start_time
    data['tracing_duration'] = tracing_duration
    data['tracing_delay'] = tracing_delay
    data['tracing_group'] = tracing_group
    data['look_back_window'] = look_back_window
    data['acute_to_trace'] = acute_to_trace

    # Here we are treating <trait> as if it is the condition/disease
    # and clusterings is a diagnostic test. We then calculate all the
    # statistics for that test

    if trans_source:
        data['trait'] = 'trans_source'
        # False Positive
        data['FP'] = np.sum(np.logical_and(metadata['cluster_size'] > 1, metadata['num_trans'] < 1))
        # True Negative
        data['TN'] = np.sum(np.logical_and(metadata['cluster_size'] <= 1, metadata['num_trans'] < 1))
        # True Positive
        data['TP'] = np.sum(np.logical_and(metadata['cluster_size'] > 1, metadata['num_trans'] > 0))
        # False Negative
        data['FN'] = np.sum(np.logical_and(metadata['cluster_size'] <= 1, metadata['num_trans'] > 0))
    else:
        data['trait'] = 'high_risk'
        # False Positive
        data['FP'] = np.sum(np.logical_and(metadata['cluster_size'] > 1, metadata['DEST_RISK'] != 'HIGH'))
        # True Negative
        data['TN'] = np.sum(np.logical_and(metadata['cluster_size'] <= 1, metadata['DEST_RISK'] != 'HIGH'))
        # True Positive
        data['TP'] = np.sum(np.logical_and(metadata['cluster_size'] > 1, metadata['DEST_RISK'] == 'HIGH'))
        # False Negative
        data['FN'] = np.sum(np.logical_and(metadata['cluster_size'] <= 1, metadata['DEST_RISK'] == 'HIGH'))

    # Actual Positives
    data['P'] = data['TP'] + data['FN']
    # Actual Negatives
    data['N'] = data['FP'] + data['TN']
    # Predicted Positives
    data['PP'] = data['TP'] + data['FP']
    # Predicted Negatives
    data['PN'] = data['FN'] + data['TN']
    # Total Population
    data['total'] = data['P'] + data['N']
    # Prevalence
    data['prevalence'] = data['P'] / data['total']

    # Accuray
    data['ACC'] = (data['TP'] + data['TN']) / data['total']
    # Positive Predictive Value
    data['PPV'] = data['TP'] / data['PP']
    # False Discovery Rate
    data['FDR'] = data['FP'] / data['PP']
    # False Omission Rate
    data['FOR'] = data['FN'] / data['PN']
    # Negative Predictive Value
    data['NPV'] = data['TN'] / data['PN']

    # True Positive Rate
    data['TPR'] = data['TP'] / data['P']
    # False Positive Rate
    data['FPR'] = data['FP'] / data['N']
    # False Negative Rate
    data['FNR'] = data['FN'] / data['P']
    # True Negative Rate
    data['TNR'] = data['TN'] / data['N']

    # Positive Likelihood Ratio
    data['LRP'] = data['TPR'] / data['FPR']
    # Negative Likelihood Ratio
    data['LRN'] = data['FNR'] / data['TNR']
    # Diagnostic Odds Ratio
    data['DOR'] = data['LRP'] / data['LRN']

    # Threat Score
    data['TS'] = data['TP'] / (data['TP'] + data['FN'] + data['FP'])
    # Informedness
    data['informedness'] = data['TPR'] + data['TNR'] - 1
    # Markedness
    data['markedness'] = data['PPV'] + data['NPV'] - 1
    # Prevalence Threshold
    data['prevalenceThreshold'] = (np.sqrt(data['TPR'] * data['FPR']) - data['FPR']) / (
                data['TPR'] - data['FPR'])
    # Matthews Correlation Coefficient
    data['MCC'] = np.sqrt(data['TPR'] * data['TNR'] * data['PPV'] * data['NPV']) - np.sqrt(
        data['FNR'] * data['FPR'] * data['FOR'] * data['FDR'])
    # Fowlkes-Mallows Index
    data['FMI'] = np.sqrt(data['PPV'] * data['TPR'])
    # F1 Score
    data['F1'] = 2 * data['PPV'] * data['TPR'] / (data['PPV'] + data['TPR'])
    # Balanced Accuracy
    data['BA'] = (data['TPR'] + data['TNR']) / 2
    # Relative Risk vs nonclustered
    data['RR'] = data['PPV'] / data['FOR']
    # Relative Risk vs all
    data['RRall'] = data['PPV'] / (data['P'] / data['total'])

    return data

if __name__ == '__main__':

    input_arguments = get_arguments()
    tracing_groups      = input_arguments['tracing_group']
    n_samples           = input_arguments['n_samples']
    tracing_start_times = input_arguments['tracing_start_time']
    tracing_durations   = input_arguments['tracing_duration']
    look_back_windows   = input_arguments['look_back_window']
    tracing_delays      = input_arguments['tracing_delay']
    acute_to_traces     = input_arguments['acute_to_trace']
    thresholds          = input_arguments['thresholds']
    
    # Read in the directory for the tree
    with open('data_path.txt') as f:
        directory = f.readline()
    # directory = 'download/79a4cded-0249-ee11-aa0a-b88303911bc1/'

    r.source('HIVContactTracing.r')
    print('Starting call to R')
    lineList_r = r.doSampling(directory,
                              tracing_group = tracing_groups,
                              n_samples = n_samples,
                              tracing_start_time = tracing_start_times,
                              tracing_durations=tracing_durations,
                              acute_to_traces=acute_to_traces,
                              look_back_duration=look_back_windows,
                              tracing_delays=tracing_delays)
    print('R is done')
    with localconverter( robjects.default_converter + pandas2ri.converter ):
        lineList  = robjects.conversion.rpy2py( lineList_r  )
    
    # Count the number of times an individual transmits
    temp = pd.DataFrame()
    temp['trans_time'] = lineList.groupby(SRC_ID)['YEAR'].apply(np.array)
    temp['num_trans'] = temp['trans_time'].apply(len)
    temp[ID] = temp.index
    lineList = lineList.merge(temp, on=ID, how='left')
    lineList['num_trans'] = lineList['num_trans'].fillna(0)
    # Make sure ID columns are strings
    lineList[ID] = lineList[ID].astype(int).astype(str)
    lineList[SRC_ID]  = lineList[SRC_ID].astype(int).astype(str)
    lineList.to_csv(os.path.join('lineList.csv'), index=False)
    print('DONE: with data preparation')
    
    # add 'random' as tracing group
    tracing_groups = tracing_groups + ['random']
    # create dataframes to hold output
    cluster_data = []
    enrichment_data = []
    for tracing_start_time in tracing_start_times:
        for tracing_duration in tracing_durations:
            for tracing_delay in tracing_delays:
                for look_back_duration in look_back_windows:
                    for acute_to_trace in acute_to_traces:
                        for n_sample in n_samples:
                            for tracing_group in tracing_groups:
                                # Get the name of the column with sample times
                                sampleTime_name = 'Sampled_TG-' + str(tracing_group) + \
                                    '_TST-' + str(tracing_start_time) + '_TL-' + \
                                    str(tracing_duration) + '_TD-' + str(tracing_delay) + \
                                    '_LBD-' + str(look_back_duration) + '_NS-' + \
                                        str(n_sample)

                                # Generate a transmission tree for each seed infection
                                trees = generate_treeFromFile.read_treeFromLineList(lineList,
                                                                                    ID=ID,
                                                                                    infectorID=SRC_ID,
                                                                                    infectTime='YEAR',
                                                                                    sampleTime=sampleTime_name,
                                                                                    features=['SRC_RISK',
                                                                                              'DEST_RISK',
                                                                                              'num_trans'])

                                # Convert trees into phylo trees
                                for i in np.arange(len(trees)):
                                    trees[i] = transform_transToPhyloTree(trees[i])
                                
                                # Coalesce into a single tree
                                tree = transform_joinTrees(trees)
                                sampleTimes = lineList[[ID, sampleTime_name]]
                                sampleTimes.rename(columns={sampleTime_name: 'SampleTime'}, inplace=True)
    
                                # Start timer
                                tic = perf_counter()
                            
                                # Get distance matrix with every node in the tree
                                distance, names = find_short_edges(tree)
                                # Get leaves
                                leaves = tree.get_leaf_names()
                                del tree
                                names = np.array(names)
                                print(len(names), 'number of nodes')
                                print(len(leaves), 'number of leaves')

                                # Get index for which nodes are leaves
                                idx = np.full((len(names)), False, dtype=bool)
                                for name in leaves:
                                    idx = np.logical_or(idx, names==name)

                                # Reduced the list of names and sample times to just the leaves
                                names = names[idx]
                                names = pd.DataFrame({ID: names})
                                names = names.merge(sampleTimes, on=ID, how='left')

                                # Reduce the distance matrix to just contain leaves
                                distance = distance[idx,:][:,idx]
                                # Fill the diagonal with inf so that the closest node is not itself
                                np.fill_diagonal(distance, float('inf'))

                                # End timer
                                toc = perf_counter()
                                print((toc-tic)/60, " minutes spent on the getting distances")
    
                                names[ID] = names[ID].astype(str)
                                names = names.merge(lineList.loc[lineList[ID].isin(names[ID]),
                                                                         [ID,
                                                                          'SRC_RISK',
                                                                          'DEST_RISK',
                                                                          'YEAR',
                                                                          'num_trans']], on=ID, how='left')
                                
                                for threshold in thresholds:
                                    print('Thershold: ', threshold)

                                    num_clusters, clusters = connected_components(distance < threshold)
                                
                                    # Add basic information: id, time sampled, cluster name, and size of the cluster
                                    names['cluster_name'] = clusters
                                    names['cluster_size'] = names.groupby('cluster_name')['cluster_name'].transform('count')
                                    names['cluster_size'] = names['cluster_size'].astype(int)
                                    
                                    temp = get_cluster_stats(names,
                                                             threshold,
                                                             n_sample,
                                                             tracing_start_time,
                                                             tracing_duration,
                                                             tracing_delay,
                                                             tracing_group,
                                                             look_back_duration,
                                                             acute_to_trace)
                                    cluster_data.append(temp)
                                    temp = get_cluster_enrichment(names,
                                                                  threshold,
                                                                  n_sample,
                                                                  tracing_start_time,
                                                                  tracing_duration,
                                                                  tracing_delay,
                                                                  tracing_group,
                                                                  look_back_duration,
                                                                  acute_to_trace,
                                                                  trans_source=True)
                                    enrichment_data.append(temp)
                                    temp = get_cluster_enrichment(names,
                                                                  threshold,
                                                                  n_sample,
                                                                  tracing_start_time,
                                                                  tracing_duration,
                                                                  tracing_delay,
                                                                  tracing_group,
                                                                  look_back_duration,
                                                                  acute_to_trace,
                                                                  trans_source=False)
                                    enrichment_data.append(temp)
                                
    cluster_data = pd.concat(cluster_data, ignore_index=True)
    enrichment_data = pd.concat(enrichment_data, ignore_index=True)
    enrichment_data['RRrandom'] = enrichment_data.groupby(['threshold',
                                               'n_samples',
                                               'tracing_start_time',
                                               'tracing_duration',
                                               'tracing_delay',
                                               'look_back_window',
                                               'acute_to_trace',
                                               'trait'], as_index=False).apply(lambda df: df['PPV']/df.loc[df['tracing_group']=='random', 'PPV'].values[0]).rename('RRrandom').reset_index(level=0, drop=True)
    cluster_data.to_csv(os.path.join('clusterData.csv'), index=False)
    enrichment_data.to_csv(os.path.join('enrichmentData.csv'), index=False)
    print('DONE')
