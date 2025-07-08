import pandas as pd
import os, subprocess, concurrent.futures, time, logging


def locate_array_kmers(num_of_repeats, rep_min, rep_max, spa_min, spa_max, input_file, sequence_out, logger):
    try:
        minced_run = [
            'minced',
            '-minNR', str(num_of_repeats),
            '-minRL', str(rep_min),
            '-maxRL', str(rep_max),
            '-minSL', str(spa_min),
            '-maxSL', str(spa_max),
            '-gff',
            input_file,
            sequence_out
        ]
        subprocess.run(minced_run)
        if os.path.getsize(sequence_out) > 0:
            return input_file, sequence_out
        else:
            os.remove(sequence_out)
            return input_file, None
    except Exception as e:
        logger.exception(f"An error occurred while processing the {input_file}: {e}")
        return input_file, None  # Return None for sequence_out if an error occurs


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Use minced to find array in each fasta fna file.')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='<-i /home/sequence_folder/file_path.txt> A file that stores the DNA sequence file path '
                             'with "file_path" in the first row')
    parser.add_argument('-rN', '--repeat_number', type=int, required=True,
                        help='<-rN 3> The minimal number of repeat')
    parser.add_argument('-rMin', '--repeat_minimal_length', type=int, required=False, default=11,
                        help='<-rMin 20> The minimal length of the repeat, default=11')
    parser.add_argument('-rMax', '--repeat_maximal_length', type=int, required=False, default=80,
                        help='<-rMax 80> The maximal length of the repeat, default=80')
    parser.add_argument('-sMin', '--spacer_minimal_length', type=int, required=False, default=11,
                        help='<-sMin 20> The minimal length of the spacer, default=11')
    parser.add_argument('-sMax', '--spacer_maximal_length', type=int, required=False, default=80,
                        help='<-sMax 80> The maximal length of the spacer, default=80')
    parser.add_argument('-n', '--num_threads', type=int, required=False, default=os.cpu_count(),
                        help='<-n 5> Number of threads (default n = all cpus), all threads  will be used if not given')
    parser.add_argument('-o', '--outdir', type=str, required=True,
                        help='<-o /test/output> An absolute directory to save output files. e.g., user/output')
    # parser.add_argument('-tsv', '--output_tsv', type=str, required=True,
    #                     help='Path to the TSV file for saving input-output paths')

    args = parser.parse_args()
    if not os.path.isabs(args.outdir):
        parser.error("The output directory must be an absolute path.")
        sys.exit(1)
        
    if not os.path.exists(f'{args.outdir}/repeat_info/'):
        os.makedirs(f'{args.outdir}/repeat_info/')

    logging.basicConfig(filename=f'{args.outdir}/output.log', level=logging.INFO)
    logger = logging.getLogger(__name__)

    df = pd.read_csv(args.input, sep='\t')
    df_file_list = df['file_path'].to_list()

    start = time.time()
    batch_size = 10000

    # List to store input-output file path pairs
    file_pairs = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.num_threads) as executor:
        for i in range(0, len(df_file_list), batch_size):
            batch = df_file_list[i:i + batch_size]
            try:
                futures = [
                    executor.submit(
                        locate_array_kmers,
                        args.repeat_number,
                        args.repeat_minimal_length,
                        args.repeat_maximal_length,
                        args.spacer_minimal_length,
                        args.spacer_maximal_length,
                        fna,
                        f'{args.outdir}/repeat_info/{fna.split("/")[-1]}_minced_out.txt',
                        logger
                    ) for fna in batch
                ]

                for future in concurrent.futures.as_completed(futures):
                    input_file, output_file = future.result()
                    if output_file is not None:  # Only append if output_file is not None
                        file_pairs.append((input_file, output_file))

            except Exception as e:
                logger.exception(f"An error occurred while processing the batch: {e}")

    # Save the input-output file paths to a TSV file
    output_df = pd.DataFrame(file_pairs, columns=['fna_path', 'file_path'])
    output_df.to_csv(f'{args.outdir}/minced_out2fna.tsv', sep='\t', index=False)

    print(f'Total finished in {round(time.time() - start, 2)} seconds\n')
    # 90 Gb fna file -> 4090.61 s