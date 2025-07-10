import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import itertools, collections, networkx, json, re,time
from pandas import DataFrame
from Bio import SeqIO
import os


def merge_domain(df):
    # na_rows = df[df['start'].isna()]
    # # print(na_rows)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['c-Evalue'] = df['c-Evalue'].astype(float)
    n2c_sorted = df.sort_values(by='start')
    n2c_sorted = n2c_sorted.reset_index(drop=True)
    while True:
        overlap = False
        rum_index = []
        for i in range(len(n2c_sorted) - 1):
            len_current = n2c_sorted['end'][i] - n2c_sorted['start'][i]
            len_next = n2c_sorted['end'][i + 1] - n2c_sorted['start'][i + 1]
            gap = n2c_sorted['end'][i] - n2c_sorted['start'][i + 1]
            rum_idx = i if n2c_sorted['c-Evalue'][i] >= n2c_sorted['c-Evalue'][i + 1] else i + 1

            if min(len_current, len_next, gap) > 0 and min(gap / len_current, gap / len_next) >= 0.5:
                rum_index.append(rum_idx)
                overlap = True

        if not overlap:
            break

        n2c_sorted = n2c_sorted.drop(rum_index)
        n2c_sorted = n2c_sorted.reset_index(drop=True)

    return n2c_sorted


def parse_pfam(pfam_table, table_from, prot_ids_focus_on=None):
    # pfam_scan.pl -fasta ./repeat_adj_prot.fasta -dir /home/hebeibei/Data/pfam -outfile ./pfamscan_out2.txt  -e_seq
    # 0.01 -e_dom 0.01
    # hmmsearch --domE 0.001 --domtblout /home/hebeibei/Work/crispr/code/hmmsearch_domtblout.txt
    # /home/hebeibei/Data/pfam/Pfam-A.hmm /home/hebeibei/Work/crispr/code/repeat_adj_prot.fasta
    # if using pfam_scan.pl, N-C already sorted in the output, but it's not for hmmsearch need re-sort
    df_pfam = pd.DataFrame()
    if table_from == 'pfamscan':
        df_pfam = pd.read_table(pfam_table,
                                delim_whitespace=True,
                                dtype='str',
                                # comment='#',  # prot_id has #, using 'skiprows' instead
                                skiprows=3,
                                skipfooter=10,
                                engine='python',
                                header=None,
                                usecols=[0, 1, 2, 5, 6, 12],
                                names=['protein_id', 'start', 'end', 'pfam_accession', 'short_name', 'E-value']
                                )
    if table_from == 'hmmsearch':
        df_pfam = pd.read_table(pfam_table,
                                delim_whitespace=True,
                                dtype='str',
                                # comment='#',
                                skiprows=3,
                                skipfooter=10,
                                engine='python',
                                header=None,
                                usecols=[0, 3, 4, 6, 11, 17, 18],
                                names=['protein_id', 'short_name', 'pfam_accession','E-value', 'c-Evalue', 'start', 'end']
                                )
    if prot_ids_focus_on is not None:
        df_prot_focurs_on = pd.read_csv(prot_ids_focus_on, sep='\t')
        df_pfam = df_pfam[df_pfam['protein_id'].isin(df_prot_focurs_on['remained_nodes'])]

    print(df_pfam)
    df_pfam['pfam_accession'] = df_pfam['pfam_accession'].str.replace(r'\..*', '', regex=True)
    df_pfam['E-value'] = df_pfam['E-value'].astype(float)
    df_pfam = df_pfam.loc[df_pfam['E-value'] <= 1e-5]
    # pfam_short_name = df_pfam.groupby(['pfam_accession'])['short_name'].apply(lambda x: '_'.join(set(x))).to_dict()
    df_on_id = df_pfam.groupby('protein_id')
    df_on_pfam = df_pfam.groupby('pfam_accession')

    # get domain_to_protein_list_dict, domain_pair_to_prot_list_dict
    id2pair = collections.defaultdict(list)
    single_nodes = []
    domain_N2C_list = []
    for prot_id, df_gy in df_on_id:
        df_gy = merge_domain(df_gy)  # merge domain if has overlap
        pfam_list = df_gy['pfam_accession'].to_list()

        if len(pfam_list) > 1:
            domain_N2C = [(pfam_list[i], pfam_list[i + 1]) for i in range(len(pfam_list) - 1)]
            domain_N2C_add = [[pfam_list[i], pfam_list[i + 1]] for i in range(len(pfam_list) - 1)]
            domain_N2C_list.extend(domain_N2C_add)
            id2pair[prot_id].extend(domain_N2C)
        if len(pfam_list) == 1:
            single_nodes.append(pfam_list[0])
            # id2pair[prot_id].append((pfam_list[0],pfam_list[0]))

    # pfam_pairs to protein id: for adding edge attributes
    pair2id = collections.defaultdict(list)
    df_id2pair = pd.DataFrame([id2pair])
    df_id2pair = df_id2pair.transpose()
    df_id2pair.columns = ['pfam_pairs']
    df_id2pair.reset_index(inplace=True)
    df_id2pair.rename(columns={'index': 'protein_id'}, inplace=True)
    df_id2pair = df_id2pair.explode('pfam_pairs')
    df_id2pair['pfam_pairs'] = df_id2pair['pfam_pairs'].map(lambda x: ','.join(x))
    for pair_id, df_gy in df_id2pair.groupby('pfam_pairs'):
        pair2id[pair_id].extend(df_gy['protein_id'].to_list())
    # print(pair2id)
    # print(domain_N2C_list)

    return df_pfam, df_on_pfam, id2pair, single_nodes, domain_N2C_list, pair2id


def build_network(pfam_table, pfam_acc, protein_fasta, graph_path, minimal_edges=None, percent=None,reference_pfam_table=None, pfam_needed=None,
                  prot_ids_focus_on=None):
    
    if not os.path.exists(graph_path):
        os.mkdir(os.path.dirname(graph_path))

    df_pfam_acc = pd.read_table(pfam_acc, sep='\t')
    pfam_acc_dict = df_pfam_acc.groupby(['pfam_accession'])['description'].agg(', '.join).to_dict()
    protein_fasta_fict = SeqIO.to_dict(SeqIO.parse(protein_fasta, 'fasta'))

    G = nx.DiGraph()

    df_pfam, df_on_pfam, id2pair, single_nodes, domain_N2C_list, pair2id = parse_pfam(pfam_table, 'hmmsearch',
                                                                                      prot_ids_focus_on)
    pfam_short_name = df_pfam.groupby(['pfam_accession'])['short_name'].apply(lambda x: '_'.join(set(x))).to_dict()

    # get reference info if provided
    ref_pair_tuple = [()]
    ref_pfam_node = set()
    ref_pair2id = {}
    ref_pfam2id = {}
    if reference_pfam_table is not None:
        ref_0, ref_1, ref_2, ref_3, ref_4, ref_5, = parse_pfam(reference_pfam_table, 'hmmsearch')
        # ref_4 = parse_pfam(reference_pfam_table, 'hmmsearch')[4]
        # ref_pair2id = parse_pfam(reference_pfam_table, 'hmmsearch')[5]
        # print(ref_4)
        ref_pair2id = ref_5
        ref_pair_tuple = [tuple(pair) for pair in ref_4]
        ref_pfam_node = set([pfam for pfam, _ in ref_1])
        ref_pfam2id = ref_1['protein_id'].apply(list).to_dict()


    # build graph
    for k, v in id2pair.items():
        G.add_edges_from(v)
    # G.add_nodes_from(set(single_nodes))
    # add node attributes: protein id, num_of_prot, annotation
    node_attr = {}
    pfam2prot_num = {}
    for pfam_id, df_gy in df_on_pfam:
        number_of_proteins = len(df_gy['protein_id'].unique())
        protein_seq_length = [len(protein_fasta_fict[prot_id]) for prot_id in df_gy['protein_id'].to_list()]
        protein_seq_length.sort()
        seq_len_min = protein_seq_length[0]
        seq_len_max = protein_seq_length[-1]
        seq_len_avg = sum(protein_seq_length) / number_of_proteins
        seq_len_mid = protein_seq_length[int(len(protein_seq_length) / 2)]
        node_attr[pfam_id] = {
            'pfam_accession': pfam_id,
            'name': f'{pfam_short_name[pfam_id]}\n{number_of_proteins} proteins\n{seq_len_min}-{seq_len_max}-{int(seq_len_avg)}-{seq_len_mid}',
            'protein_contain_this_domain': '\n'.join(set(df_gy['protein_id'].to_list())),
            'number_of_proteins': number_of_proteins,
            'presence_status': 'New' if pfam_id not in ref_pfam_node else 'Found',
            'present_in': '\n'.join(ref_pfam2id[pfam_id]) if pfam_id in ref_pfam_node else 'Na',
            'present_in_how_many_ref_Cas': len(ref_pfam2id[pfam_id]) if pfam_id in ref_pfam_node else 0,
            'annotation': pfam_acc_dict[pfam_id] if pfam_id in pfam_acc_dict.keys() else None,
            'protein_seq_length': f'min: {seq_len_min}/max: {seq_len_max}/avg: {int(seq_len_avg)}/mid: {seq_len_mid}',
        }

        pfam2prot_num[pfam_id] = number_of_proteins
    nx.set_node_attributes(G, node_attr)

    # attention to specific pfam domain if provided:
    atte_domain = []
    if pfam_needed is not None:
        df_atte = pd.read_table(pfam_needed, sep='\t')
        atte_domain.extend(df_atte['cas_pedia_pfam'].to_list())
        atte_domain.extend(df_atte['endonuclease_nuclease'].to_list())

    # add edge attributes: status(known, or New), num_of_prot
    edge_attr = {}
    missing = [pfam_id for pfam_id in edge if pfam_id not in pfam_acc_dict]
    if missing:
        print(f"Missing pfam_id annotations: {missing}")
    for pfam_pair, prot_id in pair2id.items():
        protein_seq_length = [len(protein_fasta_fict[prot_id]) for prot_id in prot_id]
        protein_seq_length.sort()
        seq_len_min = protein_seq_length[0]
        seq_len_max = protein_seq_length[-1]
        seq_len_avg = sum(protein_seq_length) / len(prot_id)
        seq_len_mid = protein_seq_length[int(len(protein_seq_length) / 2)]
        edge = tuple(pfam_pair.split(','))
        # draw edge on count/percentage
        percentage = round(len(prot_id) / max(pfam2prot_num[edge[0]], pfam2prot_num[edge[1]]), 2)
        include = (edge in ref_pair_tuple or ref_pfam_node.intersection(
            set(edge))) if reference_pfam_table is not None else True
        pfam_atte = set(edge).intersection(set(atte_domain)) if pfam_needed is not None else True
        if len(prot_id) >= minimal_edges and percentage >= percent:  # or include or pfam_atte
            edge_attr[edge] = {
                'pfam_accession': ','.join(list(edge)),
                #'annotation': ','.join([pfam_acc_dict[pfam_id] if pfam_id in pfam_acc_dict.keys() else None for pfam_id in edge]),
                #'annotation': ','.join([pfam_acc_dict[pfam_id] for pfam_id in edge if pfam_id in pfam_acc_dict]),
                'annotation': ','.join([pfam_acc_dict[pfam_id] if pfam_id in pfam_acc_dict else 'Unknown' for pfam_id in edge]),
                'name': f'C:{len(prot_id)}/P:{percentage}',
                'protein_contain_this_domain_pair': '\n'.join(prot_id),
                'number_of_proteins': len(prot_id),
                'presence_status': 'New' if edge not in ref_pair_tuple else 'Found',
                'present_in': '\n'.join(ref_pair2id[pfam_pair]) if edge in ref_pair_tuple else 'Na',
                'present_in_how_many_ref_Cas': len(ref_pair2id[pfam_pair]) if edge in ref_pair_tuple else 0,
                'protein_seq_length': f'min: {seq_len_min}/max: {seq_len_max}/avg: {int(seq_len_avg)}/mid: {seq_len_mid}',
            }
        else:
            G.remove_edge(*edge)
    nx.set_edge_attributes(G, edge_attr)
    removed_nodes = list(nx.isolates(G))
    G.remove_nodes_from(removed_nodes)

    # G.remove_edges_from(nx.selfloop_edges(G))  # remove self-connected nodes
    nx.write_graphml(G, graph_path)

    return node_attr, id2pair


if __name__ == '__main__':
    start = time.time()
    import argparse

    parser = argparse.ArgumentParser(description='Domain co-occurrence analysis of protein sequences.')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='<-i hmmout.txt> Domain output generated by "hmmsearch -domtblout "')
    parser.add_argument('-f', '--prot_fasta', type=str, required=True,
                        help='<-f protein.fasta> Protein fasta file.')
    parser.add_argument('-r', '--reference', type=str, required=False,
                        help='<-r ref_prot_hmmout.txt> Optional, hmmout txt of reference protein fasta.')
    parser.add_argument('-a', '--prot_id_focused', type=str, required=False, default=None,
                        help='<-a focused_prot_id.txt> Optional, protein ids that pay attention to,if not given, all proteins '
                             'all be analyzed')
    parser.add_argument(
        '-p', '--percent',
        type=lambda x: max(0.0, min(1.0, float(x))),
        required=False,
        default=0,
        help='Coverage: [0, 1]. Default=0, the higher values, the more strict threshold)'
    )
    parser.add_argument('-m', '--mini_edge_num', type=int, required=False, default=0,
                        help='<-m 0> Default=0, number of edges in a cluster')
    parser.add_argument('-o', '--outdir', type=str, required=True,
                        help='<-o ./output> A path to save output data. e.g., ./output')
    args = parser.parse_args()

    # input constant:
    # pfam_acc = '/home/hebeibei/Data/pfam/pfam_acc2des.txt'
    # pfam_acc = '/media/Data/qichen/pfam_acc2des.txt'
    #pfam_acc = '/media/Data/qichen/BioPrin/pfam_acc2des.txt'
    pfam_acc = '/media/Data/qichen/pfam_acc2des.37.1.txt'
    # # testing input if using pfamscan:
    # pfam_table = '/home/hebeibei/Work/crispr/code/pfamscan_out2.txt'
    # known_pfam_table = '/home/hebeibei/Work/crispr/code/pfamscan_out_known.txt'
    # # or using hmm if using hmmsearch:
    # hmm_table = '/home/hebeibei/Work/crispr/code/hmmsearch_domtblout.txt'
    # known_hmm_table = '/home/hebeibei/Data/Database_cas_db/caspdb_and_caspedia_hmmout.txt'
    # prot_fasta = '/home/hebeibei/Work/crispr/code/repeat_adj_prot.fasta'
    # graph_path = '/home/hebeibei/Work/crispr/code/tmp_test/domain_filter_on_edge.graphml'

    # # ----------------------NCBI input and output---------------------------- #
    # hmm_table = '/home/hebeibei/Data/minced_output/array_cds_pairs_ncbi/prot_fasta_hmmout.txt'
    # prot_fasta = '/home/hebeibei/Data/minced_output/array_cds_pairs_ncbi/prot_fasta_unique.csv'
    # graph_path = '/home/hebeibei/Data/minced_output/array_cds_pairs_ncbi/NCBI_domain_filter_edge_20_0.02_prot_kept_13597804_v4.graphml'
    # # optional - NCBI
    # known_hmm_table = '/home/hebeibei/Data/Database_cas_db/caspdb_and_caspedia_hmmout.txt'
    # pfam_needed = '/home/hebeibei/Data/cas_db/cas_pedia_endo-_nuclease.txt'  # currently not used
    # prot_ids_focus_on = '/home/hebeibei/Data/minced_output/array_cds_pairs_ncbi/prot_kept_13597804_nodes_c_0.5_n_20_v4.csv'
    # percent = 0.02
    # minimal_edges = 20

    # # ----------------------meta input and output---------------------------- #
    # hmm_table = '/home/hebeibei/Work/crispr/code/meer/miced_out/meer_hmm_out.txt'
    # prot_fasta = '/home/hebeibei/Work/crispr/code/meer/miced_out/prot_fasta.csv'
    # graph_path = '/home/hebeibei/Work/crispr/code/meer/miced_out/meer.graphml'
    # # optional - meta
    # known_hmm_table = '/home/hebeibei/Data/Database_cas_db/caspdb_and_caspedia_hmmout.txt'
    # pfam_needed = '/home/hebeibei/Data/cas_db/cas_pedia_endo-_nuclease.txt'  # currently not used
    # prot_ids_focus_on = '/home/hebeibei/Work/crispr/code/meer/miced_out/kept_prot.txt'
    # percent = 0
    # minimal_edges = 0
    hmm_table = args.input
    prot_fasta = args.prot_fasta
    known_hmm_table = args.reference
    pfam_needed = None
    prot_ids_focus_on = args.prot_id_focused
    percent = args.percent
    minimal_edges = args.mini_edge_num
    graph_path = f'{args.outdir}/domain_p-{percent}_m-{minimal_edges}.graphml'

    build_network(pfam_table=hmm_table,
                  pfam_acc=pfam_acc,
                  protein_fasta=prot_fasta,
                  graph_path=graph_path,
                  minimal_edges=minimal_edges,
                  percent=percent,
                  reference_pfam_table=known_hmm_table,
                  pfam_needed=pfam_needed,
                  prot_ids_focus_on=prot_ids_focus_on
                  )
    print(f'Total finished in {round(time.time() - start, 2)} seconds\n')
    # print(build_network(pfam_table))
