import os.path, sys
import pandas as pd
import time, subprocess, ast


def prot_fasta(column):
    prot = ast.literal_eval(column['cds_list'])
    if len(prot) > 0:
        prot_id_prefix = column['faa_file_name']
        prot_fasta_list = ([f">{prot_id_prefix}_#{n}\n{seq.replace('*', '')}" for n, seq in enumerate(prot) if seq])

        return prot_fasta_list, prot


def pairs_to_fasta(pairs_file,repeat_number=3, repeat_len=20):
    df = pd.read_csv(pairs_file, sep='\t')
    df_pairs = df.loc[df['cds_list'] != '[]']
    df_pairs = df_pairs.dropna(subset=['cds_list'])
    df_pairs['array_fasta'] = df_pairs.apply(
        lambda row: f">{row['faa_file_name']}\n{row['repeat_region']}", axis=1
    )
    df_pairs.drop_duplicates(['repeat_region'], keep='first', inplace=True)

    df_pairs = df_pairs[(df_pairs['rpt_unit_seq'].str.len() <= repeat_len) &
                        (df_pairs['repeat_number'] >= repeat_number)
                        ]

    df_pairs['prot_fasta'] = df_pairs.apply(lambda column: prot_fasta(column)[0], axis=1)
    df_pairs['prot_list'] = df_pairs.apply(lambda column: prot_fasta(column)[1], axis=1)
    df_pairs['cds2rep_dist'] = df_pairs['cds2rep_dist'].apply(ast.literal_eval)
    df_pairs['num_th_orf'] = df_pairs['num_th_orf'].apply(ast.literal_eval)

    df_pairs_explored = df_pairs.explode(['prot_fasta', 'prot_list', 'cds2rep_dist', 'num_th_orf'])
    df_pairs_explored.reset_index(drop=True, inplace=True)

    # dedup on concatenated array_prot_concat
    df_pairs_explored['array_prot_concat'] = df_pairs_explored['repeat_region'] + df_pairs_explored['prot_list']
    df_pairs_explored.drop_duplicates(['array_prot_concat'], keep='first', inplace=True)

    return df_pairs_explored


if __name__ == '__main__':
    start = time.perf_counter()
    import argparse

    def validate_arguments(args):
        # Validate ranges based on mode
        if args.mode == 'bp':
            if not (1000 <= args.base_pair_num <= 20000):
                print("Error: --base_pair_num must be between 1000 and 20000 when --mode is 'bp'")
                sys.exit(1)
        elif args.mode == 'orf':
            if not (1 <= args.orf_num <= 50):
                print("Error: --orf_num must be between 1 and 50 when --mode is 'orf'")
                sys.exit(1)


    parser = argparse.ArgumentParser(description='Extract array and protein fasta and deduplicate on array-protein '
                                                 'concatenation')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), required=True,
                        help='<-i array2prot_pairs.csv> The array2prot_pairs.csv obtained from 02_ap_pairs_v2.py')
    parser.add_argument('-mode', choices=['bp', 'orf'], required=True,
                        help="Specify the mode of operation: 'bp' or 'orf'")
    parser.add_argument('-rN', '--repeat_number', type=int, default=3, required=False,
                        help='<-rN 3> Lower limit of repeat number, default=3')
    parser.add_argument('-rL', '--repeat_length', type=int, default=20, required=False,
                        help='<-rL 20> Lower limit of repeat length, default=20')
    parser.add_argument('-bp_num', '--base_pair_num', type=int, default=8000, required=False,
                        help='<-bp_num 8000> Range [1000, 20000], default=8000, fetch orfs on distance')
    parser.add_argument('-orf_num', '--orf_num', type=int, default=5, required=False,
                        help='<-orf_num 5> Range [1, 50], default=5, fetch orf_num orfs up-/downstream of array')
    parser.add_argument('-o', '--outdir', type=str, required=True,
                        help='<-o /test/output> A path to save output data. e.g., ./output')
    args = parser.parse_args()
    validate_arguments(args)

    pairs_file = args.input
    fasta_dir = args.outdir
    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)

    df_pairs_final = pairs_to_fasta(pairs_file,args.repeat_number,args.repeat_length)

    # print(df_pairs_final)
    # df_pairs_final['array_prot_concat'].to_csv(f'{fasta_dir}/array_prot_concat.csv', sep='\t', index=False)
    with open(f'{fasta_dir}/array_fasta_all.csv', 'w') as f_array, open(f'{fasta_dir}/prot_fasta_all.csv', 'w') as f_prot:
        for _, row in df_pairs_final.iterrows():
            f_array.write(f"{row['array_fasta']}\n")
            f_prot.write(f"{row['prot_fasta']}\n")
    f_prot.close()
    f_array.close()

    input_type, value = '', ''
    if args.mode == 'orf':
        df_pairs_final = df_pairs_final[df_pairs_final['num_th_orf'] <= args.orf_num]
        input_type = 'orf'
        value = args.orf_num
    if args.mode == 'bp':
        df_pairs_final = df_pairs_final[df_pairs_final['cds2rep_dist'] <= args.base_pair_num]
        input_type = 'bp'
        value = args.base_pair_num

    filtered_array = f'{fasta_dir}/array_fasta_filter_by_{value}_{input_type}.csv'
    filtered_prot = f'{fasta_dir}/prot_fasta_filter_by_{value}_{input_type}.csv'
    with open(filtered_array, 'w') as f_array, open(filtered_prot, 'w') as f_prot:
        for _, row in df_pairs_final.iterrows():
            f_array.write(f"{row['array_fasta']}\n")
            f_prot.write(f"{row['prot_fasta']}\n")
    f_prot.close()
    f_array.close()

    seqkit = ['seqkit', 'rmdup', '-w', '0', '-n', f'{fasta_dir}/array_fasta_all.csv', '-o',
              f'{fasta_dir}/array_fasta_unique.csv']
    subprocess.run(seqkit)

    seqkit = ['seqkit', 'rmdup', '-w', '0', '-n', f'{filtered_array}', '-o',
              f'{filtered_array.replace(".csv", "_unique.csv")}']
    subprocess.run(seqkit)

    finish = time.perf_counter()

    print(f'finished in {round(finish - start, 3)} seconds')

