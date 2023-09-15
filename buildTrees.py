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
                        default='Incident_HIV')
    parser.add_argument('--sampleRate', '-S', metavar='SAMPLERATE',
                        help='Sample rate (as a percentage) can be list of values',
                        default=[10, 50, 100], nargs='*', type=float)
    parser.add_argument('--traceRate', '-T', metavar='TRACERATE',
                        help='Tracing success rate (as a percentage) can be list of values',
                        default=[10, 50, 100], nargs='*', type=float)
    parser.add_argument('--threshold', '-C', metavar='THRESHOLD',
                        help='Threshold to use in defining clusters can be list of values',
                        default=[0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], nargs='*', type=float)
    parser.add_argument('--traceStartTime', '-B', metavar='TRACESTARTTIME',
                        help='The time that the tracing is to start',
                        default=[2010], nargs='*', type=float)
    parser.add_argument('--traceEndTime', '-E', metavar='TRACEENDTIME',
                        help='The time that the tracing is to end',
                        default=[2012], nargs='*', type=float)
    parser.add_argument('--lookBackWindow', '-W', metavar='LOOKBACKWINDOW',
                        help='The far back from the time of diagnosis the tracing goes',
                        default=[3], nargs='*', type=float)
    parser.add_argument('--tracingDelay', '-D', metavar='TRACINGDELAY',
                        help='How long it takes to trace and test contacts',
                        default=[0], nargs='*', type=float)    
    parser.add_argument('--acuteToTrace', '-A', metavar='ACUTETOTRACE',
                        help='Only trace those that are acute at diagonsis',
                        action=argparse.BooleanOptionalAction)

    input_args = parser.parse_args()

    # Get the arguments as individual variables
    tracing_group      = input_args.traceGroup
    sample_rate        = input_args.sampleRate
    tracing_rate       = input_args.traceRate
    tracing_start_time = input_args.traceStartTime
    tracing_end_time   = input_args.traceEndTime
    look_back_window   = input_args.lookBackWindow
    tracing_delay      = input_args.tracingDelay
    acute_to_trace     = input_args.acuteToTrace
    thresholds         = input_args.threshold

    # Get the arguements into the right types/formats
    if acute_to_trace:
        acute_to_trace = True
    else:
        acute_to_trace = False
    sample_rates        = list(map(float, sample_rate))
    tracing_rates       = list(map(float, tracing_rate))
    tracing_start_times = list(map(float, tracing_start_time))
    tracing_end_times   = list(map(float, tracing_end_time))
    look_back_windows   = list(map(float, look_back_window))
    tracing_delays      = list(map(float, tracing_delay))
    thresholds          = list(map(float, thresholds))

    # Print out arguments for reference if needed
    if tracing_group not in ['Incident_HIV', 'Prevalent_HIV', 'Undiagnosed_HIV', 'HIV']:
        print('ERROR: tracing group "', tracing_group, '" not known.')
    print('Tracing Group:', tracing_group)
    print('Sample Rate:', sample_rate)
    print('Tracing Rate:', tracing_rate)
    print('Tracing Start Time:', tracing_start_time)
    print('Tracing End Time:', tracing_end_time)
    print('Look Back Window:', look_back_window)
    print('Tracing Delay:', tracing_delay)
    print('Acute to Trace:', acute_to_trace)
    print('Threshold:', thresholds)


    # Return a dictionary with all the needed variables
    return {'tracing_group': tracing_group,
            'sample_rate': sample_rates,
            'tracing_rate': tracing_rates,
            'tracing_start_time': tracing_start_times,
            'tracing_end_time': tracing_end_times,
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

def get_cluster_stats(metadata, threshold=5, sample_rate=50, tracing_rate=0,
                      tracing_start_time=0, tracing_end_time=0, tracing_delay=0,
                      tracing_group='Incident_HIV', look_back_window=0,
                      acute_to_trace=False):
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
        sample_rate (float)       : Percentage of index individuals who actually
                                    respond to tracing efforts
        tracing_rate (float)      : Success rate of tracing contacts gives as a
                                    percentage
        trace_start_time (float)  : The date when the period for contract
                                    tracing started.
        trace_end_time (float)    : The date when the period for contract
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
    data = pd.DataFrame(columns=['threshold', 'sample_rate', 'tracing_rate',
                                 'tracing_start_time', 'tracing_end_time',
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
    data['sample_rate'] = sample_rate
    data['tracing_rate'] = tracing_rate
    data['tracing_start_time'] = tracing_start_time
    data['tracing_end_time'] = tracing_end_time
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
    cluster_sizes = temp.groupby('cluster_name')['cluster_size'].agg('unique')
    data['median_size'] = cluster_sizes.median()
    data['mean_size'] = cluster_sizes.mean()
    data['max_size'] = cluster_sizes.max()

    return data

def get_cluster_enrichment(metadata, threshold=0, sample_percent=0, bias=0, burnin=0, start_year=0, end_year=0,
                           window_slide=1, window_len=5, time_dependent=False, trans_yet=0):
    """
    This function calculates if clusters are 'enriched' with transmitters. To
    do this we calculate all the metrics for evaulating a binary test where we
    consider transmitting to be a true positive and clustering to be a positive
    test result. Many of the inputs are giving just to be added to the out
    data. Probably should rework how the data is aggregated.

    Arguments:
        metadata (DataFrame)  : This data frame contains the information on each
                                individual including what cluster they are in
        threshold (float)     : The threshold used for determining clusters
        sample_percent (float): The sampling coverage as a percentage
        bias (float or str)   : The biased used in the sampling
        burnin (float)        : Number of years used for the burn-in
        start_year (float)    : The year in which we began our sampling. This
                                is only needed/used if doing the time dependent
                                analysis
        end_year (float)      : The year in which sampling ended. This is only
                                needed/used if doing the time dependent analysis
        window_slide (float)  : The number of years to move the window for each
                                bin in the time dependent analysis. If this is
                                smaller then the bin width it leds to
                                overlapping bins and gives a smoothing effect
        window_len (float)    : The width of the time bins in years.
        time_dependent (bool) : Whether or not to do the dependent version of
                                the analysis

    Return
        DataFrame             : The data frame with the results looking at the
                                full sampling period
        DataFrame             : The data frame with the results for each time
                                bin. If time_dependent is false then it is an
                                empty data frame

    """

    # It was breaking when I started with an empty data frame us here we create
    # a data frame with one row and column
    data = pd.DataFrame({'FP': [0]})

    # Set the values we already know
    data['bias'] = bias
    data['percent_sampled'] = sample_percent
    data['threshold'] = threshold
    data['burnin'] = burnin

    # Here we are treating transmitter as if it is the condition/disease
    # and clusterings is a diagnostic test. We then calculate all the
    # statistics for that test

    # False Positive
    data['FP'] = np.sum(np.logical_and(metadata['cluster_size'] > 1, metadata['num_trans'] < 1))
    # True Negative
    data['TN'] = np.sum(np.logical_and(metadata['cluster_size'] <= 1, metadata['num_trans'] < 1))
    # True Positive
    data['TP'] = np.sum(np.logical_and(metadata['cluster_size'] > 1, metadata['num_trans'] > 0))
    # False Negative
    data['FN'] = np.sum(np.logical_and(metadata['cluster_size'] <= 1, metadata['num_trans'] > 0))

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
    tracing_group = input_arguments['tracing_group']
    sample_rates = input_arguments['sample_rate']
    tracing_rates = input_arguments['tracing_rate']
    tracing_start_times = input_arguments['tracing_start_time']
    tracing_end_times = input_arguments['tracing_end_time']
    look_back_windows = input_arguments['look_back_window']
    tracing_delays = input_arguments['tracing_delay']
    acute_to_trace = input_arguments['acute_to_trace']
    thresholds = input_arguments['thresholds']

    r.source('HIVContactTracing.r')
    transmissionData_r = r.doSampling(tracing_group=tracing_group,
                                      sample_rate=sample_rates[0],
                                      tracing_rate=tracing_rates[0],
                                      tracing_start_time=tracing_start_times[0],
                                      tracing_end_time=tracing_end_times[0],
                                      acute_to_trace=acute_to_trace,
                                      look_back_duration=look_back_windows[0],
                                      tracing_delay=tracing_delays[0])
    with localconverter( robjects.default_converter + pandas2ri.converter ):
        transmissionData  = robjects.conversion.rpy2py( transmissionData_r  )
    
    # Make sure the tree stops at the end of the campaign
    transmissionData.loc[transmissionData['sampleTime']<=(tracing_end_times[0]+tracing_delays[0]), 'sampleTime'] = np.nan
    # Count the number of times an individual transmits
    temp = pd.DataFrame()
    temp['trans_time'] = transmissionData.groupby(SRC_ID)['YEAR'].apply(np.array)
    temp['num_trans'] = temp['trans_time'].apply(len)
    temp[ID] = temp.index
    transmissionData = transmissionData.merge(temp, on=ID, how='left')
    transmissionData['num_trans'] = transmissionData['num_trans'].fillna(0)
    # Make sure ID columns are strings
    transmissionData[ID] = transmissionData[ID].astype(int).astype(str)
    transmissionData[SRC_ID]  = transmissionData[SRC_ID].astype(int).astype(str)
    # Create Risk column
    transmissionData['SRC_RISK']  = transmissionData['SRC_IP'].str.extract(r'Risk:([A-Za-z]+),')
    transmissionData['DEST_RISK'] = transmissionData['DEST_IP'].str.extract(r'Risk:([A-Za-z]+),')

    # Generate a transmission tree for each seed infection
    trees = generate_treeFromFile.read_treeFromLineList(transmissionData,
                                                        ID=ID,
                                                        infectorID=SRC_ID,
                                                        infectTime='Infected',
                                                        sampleTime='sampleTime',
                                                        features=['Individual_ID',
                                                                  'contact_Risk',
                                                                  'contact_HIV_stage',
                                                                  'subset',
                                                                  'num_trans_before',
                                                                  'num_trans_after',
                                                                  'SRC_RISK',
                                                                  'DEST_RISK',
                                                                  'num_trans',
                                                                  'traced_into'])

    # Convert trees into phylo trees
    for i in np.arange(len(trees)):
        trees[i] = transform_transToPhyloTree(trees[i])
    
    # Coalesce into a single tree
    tree = transform_joinTrees(trees)
    sampleTimes = transmissionData[[ID, 'sampleTime']]
    
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

    # Reduced the list of sames and sample times to just the leaves
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
    names = names.merge(transmissionData.loc[transmissionData[ID].isin(names[ID]),
                                             [ID,
                                              'Individual_ID',
                                              'contact_Risk',
                                              'contact_HIV_stage',
                                              'subset',
                                              'num_trans_before',
                                              'num_trans_after',
                                              'SRC_RISK',
                                              'DEST_RISK',
                                              'Infected',
                                              'num_trans',
                                              'traced_into']], on=ID, how='left')
    del transmissionData

    num_clusters, clusters = connected_components(distance < thresholds[0])
    
    # Create an empty data frame
    meta_data = pd.DataFrame()

    # Add basic information: id, time sampled, cluster name, and size of the cluster
    meta_data[ID] = names[ID]
    meta_data['cluster_name'] = clusters
    meta_data['cluster_size'] = meta_data.groupby('cluster_name')['cluster_name'].transform('count')

    # Add the meta data from the original line list
    meta_data = meta_data.merge(names, on=ID)
    
    print('DONE')
