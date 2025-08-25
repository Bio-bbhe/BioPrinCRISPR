import os.path, itertools,re
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import time, subprocess
from pathlib import Path
from Bio import SeqIO


def ensure_paths(paths):
    [Path(path).mkdir(parents=True, exist_ok=True) for path in paths]


def minced_to_fna(minced_out_path, save_dir=None):
    df_minced_path = pd.read_csv(minced_out_path)
    df_minced_path['fna_path'] = df_minced_path['file_path'].map(
        lambda file_path: file_path.replace('minced_output/metadatabase', 'metadata_from_QiChen/metadata_fna').replace(
            '_minced_out.txt', ''))
    df_minced_path = df_minced_path.dropna()
    if save_dir is not None:
        df_minced_path.to_csv(f'{save_dir}minced_to_fna.csv', sep='\t', index=False)

    return df_minced_path


def get_sequences(column, prodigal_faa_path, window):
    # input infor
    start = column['start'] - 1
    end = column['end']
    sequence = column['sequence']
    repeat_sequence = sequence.seq[start: end]  # sequence object
    faa_file_name = column['faa_file_name']

    cds_start = max(0, start - window)
    cds_end = min(end + window, int(column['sequence_length']))
    cds_masked = sequence.seq[cds_start: start] + (end - start) * 'N' + sequence.seq[end:cds_end]

    # new start and end of repeat region
    new_positions = [start - cds_start, end - cds_start]

    # run prodigal
    seq = f'>seq\n{cds_masked}'
    faa_file_path = f"{prodigal_faa_path}/{faa_file_name}.faa"
    command = ['prodigal', '-a', faa_file_path, '-q', '-p', 'meta', '-f', 'gff', '-m']
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    input_data = seq.encode('utf-8')
    output_data, _ = process.communicate(input_data)

    # get cds list and distance to array
    cds_dict = SeqIO.to_dict(SeqIO.parse(faa_file_path, 'fasta'))
    cds_list = [str(cds.seq) for cds in cds_dict.values()]

    if cds_list:
        cds_start_end = [[cds.description.split(' # ')[1],cds.description.split(' # ')[2]] for cds in cds_dict.values()]
        cds2rep_dist_abs = [
            min([abs(a - int(b)) for a, b in itertools.product(new_positions, cds_start_end[i])])
            for i in range(len(cds_start_end))
        ]

        # get cds label, e.g., -4,-3,-2, -1, [repeat], 1, 2, 3, 4......
        cds_start_end = [[int(x), int(y)] for x, y in cds_start_end]
        cds_start_end.append(new_positions)
        cds_start_end.sort(key=lambda x: x[0])
        index_of_new_positions = cds_start_end.index(new_positions)
        elements_before = [i for i in range(index_of_new_positions, 0, -1)]
        elements_after = len(cds_start_end) - index_of_new_positions - 1
        elements_after = [i for i in range(elements_after + 1) if i > 0]
        elements_before.extend(elements_after)  # number of orf relative to array

        return repeat_sequence, sequence.seq[cds_start:cds_end], cds_masked, cds_list, cds2rep_dist_abs, elements_before,new_positions


def x(minced_path, fna_path, faa_path, minced2fna_csv_path, window):
    try:
        # Parse minced_output
        df_minced_out = pd.read_csv(minced_path, sep='\t', comment='#', header=None,
                                    names=['scaffold', 'source', 'repeat_region', 'start', 'end',
                                           'repeat_number', 'unk_1', 'unk_2', 'attributes']
                                    )
        # fna fasta to dict
        record_dict = SeqIO.to_dict(SeqIO.parse(open(fna_path, 'r', encoding='ISO-8859-1'), 'fasta'))
        df_minced_out['repeat_unit'] = df_minced_out['attributes'].apply(lambda x: x.split('rpt_unit_seq=')[-1])
        df_minced_out['sequence'] = df_minced_out['scaffold'].map(record_dict)
        df_minced_out['sequence_length'] = df_minced_out['sequence'].map(lambda sequence: len(sequence))
        df_minced_out['faa_file_name'] = df_minced_out.apply(
            lambda row: f"{row['scaffold']}__{row['start']}_{row['end']}", axis=1
        )

        # Process sequences and extract information
        seq_infor = df_minced_out.apply(lambda column: pd.Series(get_sequences(column, faa_path, window)), axis=1)
        if not seq_infor.empty:
            seq_infor.columns = ['repeat_sequence', 'cds_region', 'cds_masked', 'cds_list', 'cds2rep_dist',
                                 'num_th_orf', 'repeat_start_end']
            df_minced_out = pd.concat([df_minced_out, seq_infor], axis=1)
            minced_out_to_fna_name = os.path.basename(fna_path)
            df_minced_out.dropna()
            df_minced_out.to_csv(f'{minced2fna_csv_path}/{minced_out_to_fna_name}.csv', sep='\t', index=False)
            return df_minced_out[['faa_file_name', 'cds_region', 'repeat_sequence', 'cds_list', 'cds2rep_dist',
                                  'num_th_orf', 'repeat_number', 'repeat_unit', 'repeat_start_end']]
    except Exception as e:
        return fna_path, e  # Return the file path and the error to be logged


if __name__ == '__main__':
    start = time.perf_counter()
    import argparse

    parser = argparse.ArgumentParser(description='Define windows size and run prodigal')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='<-i ./minced_out2fna.tsv> The minced_out2fna.tsv obtained from 01_array.py')
    parser.add_argument('-w', '--window', type=int, required=True, default=10000,
                        help='<-w 10000> Base pairs to the array')
    parser.add_argument('-n', '--num_threads', type=int, required=False, default=os.cpu_count(),
                        help='<-n 5> Number of threads (default n=all cpus), all threads  will be used if not given')
    parser.add_argument('-o', '--outdir', type=str, required=True,
                        help='<-o /test/output> A path to save output data. e.g., ./output')
    args = parser.parse_args()

    minced_to_fna = args.input
    df_minced2fna = pd.read_csv(minced_to_fna, sep='\t')
    window = args.window
    faa_path = f'{args.outdir}/prodigal_faa/'
    single_csv_faa = f'{args.outdir}/single_csv_and_faa/'
    array2prot_pairs_path = f'{args.outdir}/'
    error_log_path = f'{args.outdir}/error_log.txt'
    paths = [faa_path, single_csv_faa, array2prot_pairs_path]
    ensure_paths(paths)

    batch_size = 10000
    with open(f'{array2prot_pairs_path}/array2prot_pairs.csv', 'w') as f, open(error_log_path, 'w') as error_log:
        f.write('faa_file_name\tcds_region\trepeat_region\tcds_list\tcds2rep_dist\tnum_th_orf\trepeat_number'
                '\trpt_unit_seq\trepeat_start_end\n')
        with ProcessPoolExecutor(max_workers=args.num_threads) as executor:
            for i in range(0, len(df_minced2fna), batch_size):
                batch = df_minced2fna.iloc[i:i + batch_size]
                futures = [
                    executor.submit(
                        x, row['file_path'], row['fna_path'], faa_path, single_csv_faa, window
                    ) for _, row in batch.iterrows()
                ]
                for future in as_completed(futures):
                    result = future.result()
                    if isinstance(result, tuple) and isinstance(result[1], Exception):  # Check if result is an error
                        fna_path, error = result
                        error_log.write(f"Error in file {fna_path}: {error}\n")  # Write error and file path to log
                    elif result is not None:
                        try:
                            result.to_csv(f, header=False, index=False, sep='\t')
                        except Exception as e:
                            error_log.write(f"Error while saving data: {e}\n")
    f.close()
    error_log.close()
    finish = time.perf_counter()
    print(f'finished in {round(finish - start, 3)} seconds')

