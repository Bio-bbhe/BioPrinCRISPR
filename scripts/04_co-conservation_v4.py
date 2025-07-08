import networkx as nx
import pandas as pd
import time, os, concurrent.futures,re
from networkx.algorithms.components.connected import connected_components
from itertools import product
from collections import Counter


def tsv2network(tsv_file, mini_node):
    print(f'building {tsv_file} network')
    G = nx.Graph()

    df_net = pd.read_csv(tsv_file, sep='\t', header=None)
    cluster = list(zip(df_net.iloc[:, 0], df_net.iloc[:, 1]))
    print(f'found {len(cluster)} nodes/edges\n')
    # cluster = [('CP012092.1__1835707_1835877', 'CP012092.1__1835707_1835877'),
    # ('DLJL01000280.1__30125_30326', 'DLJL01000280.1__30125_30326')]

    G.add_edges_from(cluster)
    print(f'found {len(list(connected_components(G)))} clusters in {tsv_file}')

    # to_be_removed
    remove_node = [node for cc in connected_components(G) if len(cc) < mini_node for node in cc]
    print(f'removed {len(remove_node)} nodes')
    G.remove_nodes_from(remove_node)
    G.remove_edges_from(nx.selfloop_edges(G))  # remove self-connected nodes
    print(f'found {len(list(connected_components(G)))} clusters in {os.path.basename(tsv_file)}')
    print(f'finished {os.path.basename(tsv_file)} network\n')
    return G


def process_node(sub_cluster):  # sub_cluster is a set = an element from list(connected_components(G))
    df_cluster = pd.DataFrame(list(sub_cluster))
    df_cluster.iloc[:, 0] = df_cluster.iloc[:, 0].str.replace(r'_#\d+', '', regex=True)
    cluster_list = df_cluster.iloc[:, 0].to_list()
    cluster_dict = dict(Counter(cluster_list))
    return set(cluster_list), cluster_dict


def get_intersection(set1, set2):
    intersection = set1.intersection(set2)
    return intersection


def cal_coverage(intersection, processed_prot_dict, prot):
    count = [processed_prot_dict[k] for k in intersection]
    cov = sum(count) / len(prot)
    return round(cov, 2)


def process_df_chunk(chunk_info):  # chunk_infor = [chunk, coverage,array_net,prot_net]
    df_chunk = chunk_info[0]
    cov = chunk_info[1]
    # array_net = chunk_info[2]
    # prot_net = chunk_info[3]
    df_chunk[['processed_prot', 'processed_prot_dict']] = df_chunk.apply(lambda row: process_node(row['prot']),
                                                                         axis=1,
                                                                         result_type='expand'
                                                                         )
    df_chunk['intersection'] = df_chunk.apply(lambda row: get_intersection(row['array'], row['processed_prot']), axis=1)
    df_chunk['cov'] = df_chunk.apply(
        lambda row: cal_coverage(row['intersection'], row['processed_prot_dict'], row['prot']), axis=1)

    # with open('/home/hebeibei/Work/crispr/code/BioPrinCRISPR_pub_data/test_cov_list.txt', 'a') as f:
    #     df_chunk['cov'].to_csv(f, sep='\t', header=False, index=False)
    # only get remaining prot ids
    df_remained = df_chunk[df_chunk['cov'] >= cov]
    keep_prot_ids = set().union(*df_remained['prot'])
    keep_array_ids = set().union(*df_remained['array'])
    # print(df_edges['array_edges'])
    # print(df_edges['prot_edges'])

    return keep_prot_ids, keep_array_ids   # manuscript_v2 using keep_array_ids


def conserved_cluster(array_net, prot_net, coverage,num_cores):
    G_array_updated = nx.Graph()
    G_prot_updated = nx.Graph()

    array_cluster = connected_components(array_net)
    prot_cluster = connected_components(prot_net)

    # version 4, only get remained protein ids  meta: 393 seconds/ NCBI: 18483 seconds
    array_prot_comb = product(list(array_cluster), list(prot_cluster))  # comb[0] = array, comb[1] = prot
    df_comb = pd.DataFrame(list(array_prot_comb), columns=['array', 'prot'])

    print(f'analyzing co-conservation on {num_cores} cores. most time-consuming step')
    chunk_size = 200000 if len(df_comb.index) > 200000 else len(df_comb.index)
    results_prot = []
    results_array = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = [
            executor.submit(process_df_chunk, [
                df_comb.iloc[i:i + chunk_size], coverage
            ]) for i in range(0, len(df_comb.index), chunk_size)
        ]

        for future in concurrent.futures.as_completed(futures):
            keep_prot_ids_result, keep_array_ids_result = future.result()
            results_prot.append(keep_prot_ids_result)
            results_array.append(keep_array_ids_result)

    keep_prot_ids = set().union(*results_prot)
    keep_array_ids = set().union(*results_array)  # manuscript_v2 using keep_array_ids

    keep_edges = list(prot_net.edges(keep_prot_ids))
    G_prot_updated.add_nodes_from(keep_prot_ids)
    G_prot_updated.add_edges_from(keep_edges)

    removed_nodes = set(prot_net.nodes()) - keep_prot_ids
    print(f'removed {len(removed_nodes)} after co-conservation')
    print(f'finished updating co-conservation graphs')
    print(f'writing co-conservation graphs')
    return G_prot_updated, removed_nodes, keep_prot_ids  # G_array_updated, G_prot_updated, removed_nodes,


if __name__ == '__main__':
    start = time.time()
    import argparse

    parser = argparse.ArgumentParser(description='Array-Prot co-conservation analysis')
    parser.add_argument('-a', '--array_cluster', type=str, required=True,
                        help='<-a 0.55_array_out_cluster.tsv> Array cluster')
    parser.add_argument('-p', '--prot_cluster', type=str, required=True,
                        help='<-p 0.35_prot_out_cluster.tsv> Prot cluster')
    parser.add_argument(
        '-c', '--coverage',
        type=lambda x: max(0.0, min(1.0, float(x))),
        required=True,
        default=0.3,
        help='Coverage: [0, 1]. The higher values, the more strict threshold)'
    )
    parser.add_argument('-m', '--mininode_num', type=int, required=False, default=5,
                        help='<-m 5> Number of nodes in a cluster')
    parser.add_argument('-n', '--num_threads', type=int, required=False, default=os.cpu_count(),
                        help='<-n 5> Number of threads (default n = all cpus), all threads  will be used if not given')
    parser.add_argument('-o', '--outdir', type=str, required=True,
                        help='<-o ./output> A path to save output data. e.g., ./output')
    args = parser.parse_args()

    array_tsv = args.array_cluster
    prot_tsv = args.prot_cluster
    coverage = args.coverage
    mininode = args.mininode_num
    num_cores = args.num_threads
    outdir = args.outdir

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # # input for testing
    # prot_tsv = '/home/hebeibei/Work/crispr/code/tmp_test/0.3_prot.tsv_cluster.tsv'
    # array_tsv = '/home/hebeibei/Work/crispr/code/tmp_test/0.5_repeat_cluster.tsv'
    # # output for testing:
    # outdir = '/home/hebeibei/Work/crispr/code/BioPrinCRISPR_pub_data/'
    # coverage = 0
    # mininode = 0  # 5 for < 10k proteins, otherwise 20
    # num_cores = 64

    # # input for working_meta
    # prot_tsv = '/home/hebeibei/Data/minced_output/array_cds_pairs_meta/0.35_prot_cluster.tsv'
    # array_tsv = '/home/hebeibei/Data/minced_output/array_cds_pairs_meta/0.5_array_cluster.tsv'
    # # output for working:
    # outdir = '/home/hebeibei/Data/minced_output/array_cds_pairs_meta'

    # # input for working_NCBI
    # prot_tsv = '/home/hebeibei/Data/minced_output/array_cds_pairs_ncbi/0.35_prot_fasta_cluster.tsv'
    # array_tsv = '/home/hebeibei/Data/minced_output/array_cds_pairs_ncbi/0.5_array_cluster.tsv'
    # # output for working:
    # outdir = '/home/hebeibei/Data/minced_output/array_cds_pairs_ncbi'

    df_prot_id = pd.read_csv(prot_tsv, sep='\t', header=None)
    all_prot = len(df_prot_id.iloc[:, 1].index)
    all_central_prot = set(df_prot_id.iloc[:, 0])
    del df_prot_id

    prot_net = tsv2network(prot_tsv, mininode)
    array_net = tsv2network(array_tsv, mininode)

    G_prot_updated, removed_nodes, keep_prot_ids = conserved_cluster(array_net, prot_net, coverage, num_cores)

    # Function to get the next version number for a given file pattern
    def get_next_version(outdir, filename_prefix, file_extension):
        # Ensure the output directory exists
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # Create a regex pattern to match the existing files
        pattern = re.compile(rf"v(\d+)_{filename_prefix}\.{file_extension}$")
        existing_files = [f for f in os.listdir(outdir) if pattern.match(f)]

        if not existing_files:
            return 1

        versions = []
        for file in existing_files:
            match = pattern.search(file)
            if match:
                versions.append(int(match.group(1)))

        return max(versions) + 1 if versions else 1


    # Get version numbers for each output file
    prot_version = get_next_version(outdir, "prot_cluster", "graphml")
    kept_nodes_version = get_next_version(outdir, "prot_kept","csv")
    kept_array_version = get_next_version(outdir,"array_kept","csv")
    rep_nodes_version = get_next_version(outdir,"prot_representative_nodes","txt")
    cluster_info_version = get_next_version(outdir, "CLUSTER_INFO", "txt")

    nx.write_graphml(G_prot_updated, f'{outdir}/v{prot_version}_prot_cluster.graphml')

    df_kept_nodes = pd.DataFrame(list(keep_prot_ids), columns=[f'remained_nodes'])
    df_kept_nodes.to_csv(f'{outdir}/v{kept_nodes_version}_prot_kept.csv', index=False)

    df_kept_nodes.iloc[:, 0] = df_kept_nodes.iloc[:, 0].str.replace(r'_#\d+', '', regex=True)
    df_kept_array = df_kept_nodes.drop_duplicates(subset=['remained_nodes'], keep='first')
    df_kept_array.to_csv(f'{outdir}/v{kept_array_version}_array_kept.csv', index=False)

    # get central node/representative from updated graph = 1st column in .stv - removed nodes
    prot_representative_nodes = all_central_prot.intersection(keep_prot_ids)
    print(f'prot_representative_nodes:{len(prot_representative_nodes)}')
    with open(f'{outdir}/v{rep_nodes_version}_prot_representative_nodes.txt', 'w') as f:
        f.writelines(f'{node}\n' for node in prot_representative_nodes)

    CLUSTER_INFO = {
        'parameters': '',
        'prot_tsv': prot_tsv,
        'array_tsv': array_tsv,
        'outdir': outdir,
        'df_kept_nodes': f'{outdir}/prot_kept_{len(keep_prot_ids)}_nodes_c_{coverage}_n_{mininode}.csv',
        'prot_representative_nodes': f'{outdir}/prot_representative_nodes.txt',
        'coverage': coverage,
        'mininode': mininode,
        'protein_cluster_info': '',
        'protein_net_clusters': f'{len(list(connected_components(prot_net)))} clusters',
        'number_of_protein_filter_on_mininode': len(prot_net),
        'updated protein_net_clusters': f'{len(list(connected_components(G_prot_updated)))} clusters',
        'removed_protein_number_after_co-conservation': len(removed_nodes),
        'removed_protein_in_total': all_prot - len(G_prot_updated),
        'remaining_proteins': len(G_prot_updated),
        'remaining_arrays': len(df_kept_array['remained_nodes']),
        'prot_representative_nodes:': len(prot_representative_nodes),
        'total_time_used': f'{round(time.time() - start, 2)} seconds'
    }

    with open(f'{outdir}/v{cluster_info_version}_CLUSTER_INFO.txt', 'w') as f:
        for key, value in CLUSTER_INFO.items():
            f.write(f'{key}: {value}\n')

    print(f'Total finished in {round(time.time() - start, 2)} seconds\n')
