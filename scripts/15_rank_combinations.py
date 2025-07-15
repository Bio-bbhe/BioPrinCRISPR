import networkx as nx
import csv
import argparse
import re

def parse_length_summary(length_str):
    match = re.search(r"min:\s*(\d+)/max:\s*(\d+)/avg:\s*(\d+)", length_str)
    if match:
        return int(match.group(1)), int(match.group(2)), int(match.group(3))
    else:
        return 0, 0, 0

def load_pfam_descriptions(desc_file):
    pfam_desc = {}
    with open(desc_file, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split(None, 1)
            if len(parts) == 2:
                pfam_id, desc = parts
                pfam_desc[pfam_id] = desc
    return pfam_desc

def extract_domain_combinations(graphml_path, desc_file, output_csv):
    pfam_desc = load_pfam_descriptions(desc_file)
    G = nx.read_graphml(graphml_path)
    combo_stats = []

    for u, v, data in G.edges(data=True):
        domain_pair = data.get("pfam_accession", f"{u},{v}")
        pf1, pf2 = domain_pair.split(",")
        desc1 = pfam_desc.get(pf1, pf1)
        desc2 = pfam_desc.get(pf2, pf2)
        domain_combination = f"{pf1} - {pf2}"
        domain_description = f"{desc1} - {desc2}"

        protein_count = int(data.get("number_of_proteins", 0))
        length_summary = data.get("protein_seq_length", "")
        min_len, max_len, avg_len = parse_length_summary(length_summary)

        present_in = data.get("present_in_how_many_ref_Cas", "")
        presence_status = data.get("presence_status", "")

        combo_stats.append({
            'domain_combination': domain_combination,
            'domain_description': domain_description,
            'protein_count': protein_count,
            'avg_protein_length': avg_len,
            'min_length': min_len,
            'max_length': max_len,
            'present_in_how_many_ref_Cas': present_in,
            'presence_status': presence_status
        })

    combo_stats.sort(key=lambda x: x['protein_count'], reverse=True)
    rank = 1
    prev_count = None
    for i, row in enumerate(combo_stats):
        if row['protein_count'] != prev_count:
            rank = i + 1
            prev_count = row['protein_count']
        row['rank'] = rank

    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'rank', 'domain_combination', 'domain_description', 'protein_count',
            'avg_protein_length', 'min_length', 'max_length',
            'present_in_how_many_ref_Cas', 'presence_status'
        ])
        writer.writeheader()
        writer.writerows(combo_stats)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract domain co-occurrence combinations from GraphML and map to descriptions.")
    parser.add_argument("-i", "--input", required=True, help="Path to input .graphml file")
    parser.add_argument("-d", "--desc", required=True, help="Path to pfam_accession-to-description file")
    parser.add_argument("-o", "--output", required=True, help="Path to output .csv file")
    args = parser.parse_args()

    extract_domain_combinations(args.input, args.desc, args.output)

