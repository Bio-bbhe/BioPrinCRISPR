import torch
import math, time, os
import pandas as pd
from Bio import SeqIO
from transformers import EsmModel, EsmTokenizer
from itertools import islice
import argparse as ap

def get_embeddings(sequences, esm_model, tokenizer,emb_path,cds_dic, batch_size=32):  # Reduced batch size to reduce GPU memory usage
    with open(emb_path, 'w') as f:
        f.write('Sequence_ID\t' + '\t'.join([f'Embedding_{i}' for i in range(esm_model.config.hidden_size)]) + '\n')  # Write header
        # f.write('Sequence_ID\t' + 'embedding' + '\n')
        for i in range(0, len(sequences), batch_size):
            batch_seqs = sequences[i:i + batch_size]
            inputs = tokenizer(batch_seqs, return_tensors="pt", padding=True, truncation=True)
            input_ids = inputs["input_ids"].to(device)
            attention_mask = inputs["attention_mask"].to(device)

            with torch.no_grad():
                torch.cuda.empty_cache()  # Clear cache before model inference
                outputs = esm_model(input_ids=input_ids, attention_mask=attention_mask)
                seq_embedding = outputs.last_hidden_state.mean(dim=1).to(dtype=torch.float32)

            # Append embeddings to CSV file
            for idx, embedding in enumerate(seq_embedding.cpu().numpy()):
                seq_id = list(cds_dic.keys())[i + idx]
                f.write(seq_id + '\t' + '\t'.join(map(str, embedding)) + '\n')
                # f.write(seq_id + '\t' + str(list(embedding)) + '\n')

            del input_ids, attention_mask, outputs, seq_embedding  # Free GPU memory
            torch.cuda.empty_cache()  # Release unused memory


def compute_distances_from_csv(query_csv, ref_csv, batch_size=1024):
    embeddings_query = pd.read_csv(query_csv, sep='\t', index_col=0)
    embeddings_ref = pd.read_csv(ref_csv, sep='\t', index_col=0)

    query_tensors = torch.tensor(embeddings_query.values.astype(float), dtype=torch.float32)
    ref_tensors = torch.tensor(embeddings_ref.values.astype(float), dtype=torch.float32)

    distances = []
    for i in range(0, query_tensors.size(0), batch_size):
        batch_query = query_tensors[i:i + batch_size].to(device)
        with torch.no_grad():
            batch_distances = torch.cdist(batch_query, ref_tensors.to(device), p=2)
        distances.append(batch_distances.cpu())  # Move distances to CPU
        del batch_query, batch_distances  # Free GPU memory
        torch.cuda.empty_cache()  # Release unused memory

    return torch.cat(distances)


if __name__ == '__main__':
    # ----------------------------- input ----------------------------- #

    start = time.time()
    
    #model_path = '/home/hebeibei/Data/model/esm2_t33_650M_UR50D'
    #faa_query_dirs = '/home/hebeibei/Work/dl/code/seq2seq/inference/e.coli_PPI'

    opts = ap.ArgumentParser()
    opts.add_argument("-m", "--model_path", type=str, help="Path to the ESM model directory", required=True)
    opts.add_argument("-q", "--faa_query_dirs", type=str, help="Path to the query FASTA directory", required=True)
    opts.add_argument("-e", "--emb_path", type=str, help="Path to the output embedding file", required=True)

    ap = opts.parse_args()

    model_path = ap.model_path
    faa_query_dirs = ap.faa_query_dirs
    emb_path = ap.emb_path

    if not os.path.exists(emb_path):
        os.makedirs(emb_path)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    esm_model = EsmModel.from_pretrained(model_path)
    tokenizer = EsmTokenizer.from_pretrained(model_path)
    esm_model = esm_model.to(device)

    for faa in os.listdir(faa_query_dirs):
        print(faa)
        faa_file_query = os.path.join(faa_query_dirs, faa)
        outcsv = os.path.join(emb_path, faa.replace('.fasta', '_emb.csv'))
        cds_dict_query = SeqIO.to_dict(SeqIO.parse(faa_file_query, 'fasta'))
        # cds_dict_query = dict(islice(cds_dict_query.items(), 128))  # for testing, uncomment for full seq
        sequences_query = [str(v.seq).upper() for v in cds_dict_query.values()]
        get_embeddings(sequences_query, esm_model, tokenizer,outcsv,cds_dict_query)

    print(f'Total finished in {round(time.time() - start, 2)} seconds\n')
