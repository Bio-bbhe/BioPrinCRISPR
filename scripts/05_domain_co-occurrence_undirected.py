import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import itertools, collections, json, re, time
from pandas import DataFrame
from Bio import SeqIO
import os

def merge_domain(df):
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['c-Evalue'] = df['c-Evalue'].astype(float)
    n2c_sorted = df.sort_values(by='start').reset_index(drop=True)
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
        n2c_sorted = n2c_sorted.drop(rum_index).reset_index(drop=True)
    return n2c_sorted


def parse_pfam(pfam_table, table_from, prot_ids_focus_on=None):
    if table_from == 'pfamscan':
        df_pfam = pd.read_table(pfam_table, delim_whitespace=True, dtype='str',
                                skiprows=3, skipfooter=10, engine='python',
                                header=None, usecols=[0, 1, 2, 5, 6, 12],
                                names=['protein_id', 'start', 'end', 'pfam_accession', 'short_name', 'E-value'])
    elif table_from == 'hmmsearch':
        df_pfam = pd.read_table(pfam_table, delim_whitespace=True, dtype='str',
                                skiprows=3, skipfooter=10, engine='python',
                                header=None, usecols=[0, 3, 4, 6, 11, 17, 18],
                                names=['protein_id', 'short_name', 'pfam_accession', 'E-value', 'c-Evalue', 'start', 'end'])

    if prot_ids_focus_on is not None:
        df_focus = pd.read_csv(prot_ids_focus_on, sep='\t')
        df_pfam = df_pfam[df_pfam['protein_id'].isin(df_focus['remained_nodes'])]

    df_pfam['pfam_accession'] = df_pfam['pfam_accession'].str.replace(r'\..*', '', regex=True)
    df_pfam['E-value'] = df_pfam['E-value'].astype(float)
    ### Evalue threshold
    df_pfam = df_pfam.loc[df_pfam['E-value'] <= 1e-5]
    df_on_id = df_pfam.groupby('protein_id')
    df_on_pfam = df_pfam.groupby('pfam_accession')

    id2pair = collections.defaultdict(list)
    domain_pairs_all = []
    single_nodes = []

    for prot_id, df_gy in df_on_id:
        df_gy = merge_domain(df_gy)
        pfam_list = df_gy['pfam_accession'].tolist()

        if len(pfam_list) > 1:
            all_combinations = list(itertools.combinations(sorted(set(pfam_list)), 2))
            id2pair[prot_id].extend(all_combinations)
            domain_pairs_all.extend(all_combinations)
        elif len(pfam_list) == 1:
            single_nodes.append(pfam_list[0])

    pair2id = collections.defaultdict(list)
    for prot_id, pfam_pairs in id2pair.items():
        for pair in pfam_pairs:
            pair2id[','.join(pair)].append(prot_id)

    return df_pfam, df_on_pfam, id2pair, single_nodes, domain_pairs_all, pair2id


def build_network(pfam_table, pfam_acc, protein_fasta, graph_path, minimal_edges=0, percent=0,
                  reference_pfam_table=None, pfam_needed=None, prot_ids_focus_on=None):
    df_pfam_acc = pd.read_table(pfam_acc, sep='\t')
    pfam_acc_dict = df_pfam_acc.groupby('pfam_accession')['description'].agg(', '.join).to_dict()
    protein_fasta_dict = SeqIO.to_dict(SeqIO.parse(protein_fasta, 'fasta'))

    G = nx.Graph()  # Use an undirected graph here

    df_pfam, df_on_pfam, id2pair, single_nodes, domain_pairs_all, pair2id = parse_pfam(pfam_table, 'hmmsearch', prot_ids_focus_on)
    pfam_short_name = df_pfam.groupby('pfam_accession')['short_name'].apply(lambda x: '_'.join(set(x))).to_dict()

    ref_pair_tuple = set()
    ref_pfam_node = set()
    ref_pair2id = {}
    ref_pfam2id = {}

    if reference_pfam_table is not None:
        ref_df, ref_on_pfam, _, _, ref_pairs, ref_pair2id = parse_pfam(reference_pfam_table, 'hmmsearch')
        ref_pair_tuple = set(ref_pairs)
        ref_pfam_node = set(ref_df['pfam_accession'])
        ref_pfam2id = ref_df.groupby('pfam_accession')['protein_id'].apply(list).to_dict()

    for prot_id, pairs in id2pair.items():
        for pair in pairs:
            G.add_edge(pair[0], pair[1])  # Add edges between domains that co-occur

    node_attr = {}
    pfam2prot_num = {}
    for pfam_id, df_gy in df_on_pfam:
        prot_ids = df_gy['protein_id'].unique()
        lengths = [len(protein_fasta_dict[p]) for p in prot_ids if p in protein_fasta_dict]
        if not lengths:
            continue
        lengths.sort()
        stats = {
            'min': lengths[0],
            'max': lengths[-1],
            'avg': int(sum(lengths) / len(lengths)),
            'mid': lengths[len(lengths) // 2],
        }

        node_attr[pfam_id] = {
            'pfam_accession': pfam_id,
            'name': f'{pfam_short_name.get(pfam_id, pfam_id)}\n{len(prot_ids)} proteins\n{stats["min"]}-{stats["max"]}-{stats["avg"]}-{stats["mid"]}',
            'protein_contain_this_domain': '\n'.join(prot_ids),
            'number_of_proteins': len(prot_ids),
            'presence_status': 'New' if pfam_id not in ref_pfam_node else 'Found',
            'present_in': '\n'.join(ref_pfam2id.get(pfam_id, [])),
            'present_in_how_many_ref_Cas': len(ref_pfam2id.get(pfam_id, [])),
            'annotation': pfam_acc_dict.get(pfam_id, ''),
            'protein_seq_length': f'min: {stats["min"]}/max: {stats["max"]}/avg: {stats["avg"]}/mid: {stats["mid"]}',
        }
        pfam2prot_num[pfam_id] = len(prot_ids)

    nx.set_node_attributes(G, node_attr)

    # read target pfam if any
    atte_domain = []
    if pfam_needed is not None:
        df_atte = pd.read_table(pfam_needed, sep='\t')
        atte_domain.extend(df_atte['cas_pedia_pfam'].dropna().tolist())
        atte_domain.extend(df_atte['endonuclease_nuclease'].dropna().tolist())

    edge_attr = {}
    for pair_str, prot_ids in pair2id.items():
        pfam1, pfam2 = pair_str.split(',')
        edge = tuple(sorted((pfam1, pfam2)))
        lengths = [len(protein_fasta_dict[p]) for p in prot_ids if p in protein_fasta_dict]
        if not lengths:
            continue
        lengths.sort()
        stats = {
            'min': lengths[0],
            'max': lengths[-1],
            'avg': int(sum(lengths) / len(lengths)),
            'mid': lengths[len(lengths) // 2],
        }

        percentage = round(len(prot_ids) / max(pfam2prot_num.get(pfam1, 1), pfam2prot_num.get(pfam2, 1)), 2)
        include = (edge in ref_pair_tuple or ref_pfam_node.intersection(edge)) if reference_pfam_table else True
        atte = set(edge).intersection(atte_domain) if pfam_needed else True

        if len(prot_ids) >= minimal_edges and percentage >= percent:
            edge_attr[edge] = {
                'pfam_accession': ','.join(edge),
                'annotation': ','.join([pfam_acc_dict.get(p, '') for p in edge]),
                'name': f'C:{len(prot_ids)}/P:{percentage}',
                'protein_contain_this_domain_pair': '\n'.join(prot_ids),
                'number_of_proteins': len(prot_ids),
                'presence_status': 'New' if edge not in ref_pair_tuple else 'Found',
                'present_in': '\n'.join(ref_pair2id.get(','.join(edge), [])),
                'present_in_how_many_ref_Cas': len(ref_pair2id.get(','.join(edge), [])),
                'protein_seq_length': f'min: {stats["min"]}/max: {stats["max"]}/avg: {stats["avg"]}/mid: {stats["mid"]}',
            }
        else:
            G.remove_edge(*edge)

    nx.set_edge_attributes(G, edge_attr)
    G.remove_nodes_from(list(nx.isolates(G)))
    nx.write_graphml(G, graph_path)

    return node_attr, id2pair

if __name__ == '__main__':
    import argparse
    start = time.time()
    parser = argparse.ArgumentParser(description='Domain co-occurrence analysis (undirected).')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='<-i hmmout.txt>')
    parser.add_argument('-f', '--prot_fasta', type=str, required=True,
                        help='<-f protein.fasta>')
    parser.add_argument('-r', '--reference', type=str, required=False,
                        help='<-r ref_prot_hmmout.txt>')
    parser.add_argument('-a', '--prot_id_focused', type=str, required=False, default=None,
                        help='<-a focused_prot_id.txt>')
    parser.add_argument('-p', '--percent', type=lambda x: max(0.0, min(1.0, float(x))), default=0,
                        help='Percentage threshold [0-1], default=0')
    parser.add_argument('-m', '--mini_edge_num', type=int, default=0,
                        help='Minimal number of proteins per domain-pair')
    parser.add_argument('-o', '--outdir', type=str, required=True,
                        help='<-o ./output>')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    #pfam_acc = '/media/Data/qichen/BioPrin/pfam_acc2des.txt'
    pfam_acc = '/media/Data/qichen/pfam_acc2des.37.1.txt'
    graph_path = f'{args.outdir}/domain_undirected_p-{args.percent}_m-{args.mini_edge_num}.graphml'

    build_network(
        pfam_table=args.input,
        pfam_acc=pfam_acc,
        protein_fasta=args.prot_fasta,
        graph_path=graph_path,
        minimal_edges=args.mini_edge_num,
        percent=args.percent,
        reference_pfam_table=args.reference,
        pfam_needed=None,
        prot_ids_focus_on=args.prot_id_focused
    )
    print(f'Total finished in {round(time.time() - start, 2)} seconds\n')

