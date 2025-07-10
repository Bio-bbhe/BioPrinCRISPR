import torch
import pandas as pd
import numpy as np
import torch.nn.functional as F
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
import os
import argparse as ap
import umap.umap_ as umap
import seaborn as sns
from adjustText import adjust_text


def perform_kmeans_clustering(embeddings_df: pd.DataFrame,
                              n_clusters: int = None,
                              run_k_optimization: bool = False,
                              k_range: range = None,
                              random_state: int = 42,
                              n_init: int = 10,
                              output_dir: str = '.',
                              plot_filename: str = 'kmeans_k_optimization.png'):
    if not isinstance(embeddings_df, pd.DataFrame):
        raise TypeError("embeddings_df must be a pandas DataFrame.")
    if embeddings_df.empty:
        raise ValueError("embeddings_df cannot be empty.")
    if run_k_optimization and k_range is None:
        k_range = range(2, min(11, len(embeddings_df) // 2))  # Default k_range up to 10 or half the data if smaller
        print(f"No k_range specified for optimization, defaulting to {list(k_range)}")
    if not run_k_optimization and n_clusters is None:
        raise ValueError("Either n_clusters must be specified or run_k_optimization must be True.")

    embeddings_np = embeddings_df.values.astype(float)

    print("Ensuring L2 normalization of embeddings...")
    embeddings_np = normalize(embeddings_np, axis=1, norm='l2')
    print("Normalization complete.")

    if run_k_optimization:
        print(f"\n--- Running K-Means optimization for k in {list(k_range)} ---")
        inertia_values = []
        silhouette_scores = []

        # Ensure k_range is valid for silhouette_score
        valid_k_range = [k for k in k_range if k > 1 and k < len(embeddings_np)]

        if not valid_k_range:
            print(
                "Warning: No valid k values in range for optimization (k must be > 1 and < num_samples). Skipping optimization.")
            return None, None, None

        for k in valid_k_range:
            print(f"  Clustering with k={k}...")
            kmeans = KMeans(n_clusters=k, random_state=random_state, n_init=n_init, verbose=False)
            kmeans.fit(embeddings_np)
            inertia_values.append(kmeans.inertia_)

            # Calculate silhouette score
            # Only if k is valid for silhouette (k > 1 and k < num_samples)
            if k > 1:
                current_labels = kmeans.labels_
                # Avoid calculating silhouette for trivial cases or if all points are in one cluster
                if len(np.unique(current_labels)) > 1 and len(np.unique(current_labels)) < len(embeddings_np):
                    silhouette_avg = silhouette_score(embeddings_np, current_labels, metric='euclidean')
                    silhouette_scores.append(silhouette_avg)
                else:
                    silhouette_scores.append(np.nan)  # Mark as Not a Number
            else:
                silhouette_scores.append(np.nan)

        metrics_df = pd.DataFrame({
            'k': valid_k_range,
            'Inertia': inertia_values,
            'Silhouette_Score': silhouette_scores
        }).set_index('k')

        # Plotting the results
        os.makedirs(output_dir, exist_ok=True)
        plot_path = os.path.join(output_dir, plot_filename)

        fig, ax1 = plt.subplots(figsize=(10, 6))

        ax1.plot(metrics_df.index, metrics_df['Inertia'], marker='o', linestyle='-', color='b')
        ax1.set_xlabel('Number of Clusters (k)')
        ax1.set_ylabel('Inertia (Elbow Method)', color='b')
        ax1.tick_params(axis='y', labelcolor='b')
        ax1.set_title('K-Means Optimization: Elbow Method and Silhouette Score')
        ax1.grid(True)

        ax2 = ax1.twinx()
        ax2.plot(metrics_df.index, metrics_df['Silhouette_Score'], marker='x', linestyle='--', color='r')
        ax2.set_ylabel('Silhouette Score', color='r')
        ax2.tick_params(axis='y', labelcolor='r')
        ax2.set_ylim([-0.1, 1.0])  # Silhouette scores range from -1 to 1

        fig.tight_layout()
        plt.savefig(plot_path)
        plt.close(fig)
        print(f"Optimization plot saved to {plot_path}")
        print("\n--- K-Means Optimization Results ---")
        print(metrics_df)
        print("\nExamine the plot and table to choose an optimal 'k'.")
        print("Look for an 'elbow' in the Inertia curve and the highest Silhouette Score.")

        if n_clusters is None:
            return None, None, metrics_df  # Return optimization results, no final labels/model yet.

    # Perform final clustering with the chosen n_clusters (either from input or after optimization)
    print(f"\n--- Performing final K-Means clustering with k={n_clusters} ---")
    kmeans_final = KMeans(n_clusters=n_clusters, random_state=random_state,
                          n_init=n_init, verbose=True)

    labels = kmeans_final.fit_predict(embeddings_np)
    print("Final K-Means clustering complete.")

    # Calculate final silhouette score
    final_silhouette_score = np.nan
    if len(np.unique(labels)) > 1 and len(np.unique(labels)) < len(embeddings_np):
        final_silhouette_score = silhouette_score(embeddings_np, labels, metric='euclidean')
        print(f"Final Silhouette Score: {final_silhouette_score}")
    else:
        print("Cannot calculate Silhouette Score (too few or too many clusters for final run).")

    metrics_df_final = pd.DataFrame({
        'k': [n_clusters],
        'Inertia': [kmeans_final.inertia_],
        'Silhouette_Score': [final_silhouette_score]
    }).set_index('k')

    return labels, kmeans_final, metrics_df_final

def run_visualization(embedding_csv_ref: str,
                      embedding_csv_query: str,
                      clusters_filepath: str,
                      visualization_png: str,
                      no_clustered_png: str):
    print("Loading embeddings...")
    embeddings_query_df = pd.read_csv(embedding_csv_query, sep='\t', index_col=0)
    embeddings_ref_df = pd.read_csv(embedding_csv_ref, sep='\t', index_col=0)

    combined_embeddings_df = pd.concat([embeddings_query_df, embeddings_ref_df])

    if not combined_embeddings_df.index.is_unique:
        print(
            "Warning: Duplicate Sequence IDs found after combining. Removing duplicate index entries (keeping first occurrence).")
        combined_embeddings_df = combined_embeddings_df[~combined_embeddings_df.index.duplicated(keep='first')]

    print(f"Combined embeddings shape: {combined_embeddings_df.shape}")

    print(f"Loading cluster assignments from: {clusters_filepath}")
    cluster_assignments_df = pd.read_csv(clusters_filepath, sep='\t', index_col='Sequence_ID')

    merged_data_df = combined_embeddings_df.join(cluster_assignments_df, how='inner')

    if merged_data_df.empty:
        raise ValueError(
            "No common Sequence_IDs found between combined embeddings and cluster assignments. Check your data files and indices.")

    print(f"Shape after merging embeddings and cluster assignments: {merged_data_df.shape}")

    embeddings_for_umap_df = merged_data_df.drop(columns=['Cluster_ID'])
    cluster_labels_aligned = merged_data_df['Cluster_ID'].astype(int) 

    embeddings_for_umap_np = normalize(embeddings_for_umap_df.values.astype(float), axis=1, norm='l2')

    print(f"Normalized embedding shape for UMAP: {embeddings_for_umap_np.shape}")

    print("Starting UMAP dimensionality reduction (this might take a few minutes)...")
    reducer = umap.UMAP(
        n_components=2,  # We want 2D for standard scatter plot
        n_neighbors=15,  # Common default, adjust if needed (e.g., 5-50)
        min_dist=0.1,  # Common default, adjust if needed (e.0.0-0.99)
        metric='euclidean',  # Use Euclidean distance as embeddings are L2-normalized
        random_state=42,  # For reproducibility
        verbose=True  # To see progress
    )
    umap_embeddings_2d = reducer.fit_transform(embeddings_for_umap_np)
    print("UMAP reduction complete.")
    print(f"UMAP 2D embedding shape: {umap_embeddings_2d.shape}")

    umap_plot_df = pd.DataFrame(umap_embeddings_2d,
                                columns=['UMAP_X', 'UMAP_Y'],
                                index=embeddings_for_umap_df.index)

    umap_plot_df['Cluster_ID'] = cluster_labels_aligned
    umap_plot_df['Sequence_ID'] = umap_plot_df.index.to_list() 
    if umap_plot_df['Cluster_ID'].isnull().any():
        print(
            f"ERROR: {umap_plot_df['Cluster_ID'].isnull().sum()} NaNs found in 'Cluster_ID' column after merge. Plotting might be incorrect.")

    num_unique_clusters = len(umap_plot_df['Cluster_ID'].unique())

    plt.figure(figsize=(12, 10))

    highlight_ids = umap_plot_df.index[
        umap_plot_df.index.str.startswith('Cas13an') |
        (umap_plot_df.index == 'MEQ1724827.1') |
        (umap_plot_df.index == 'MEK6541296.1')
    ].tolist()

    umap_plot_df['Highlight'] = umap_plot_df.index.isin(highlight_ids)

    palette = sns.color_palette("tab20", n_colors=num_unique_clusters)

    sns.scatterplot(
        x='UMAP_X',
        y='UMAP_Y',
        hue='Cluster_ID',
        palette=palette,
        data=umap_plot_df[~umap_plot_df['Highlight']],
        legend=False,
        alpha=0.6,
        s=10,
        linewidth=0
    )

    sns.scatterplot(
        x='UMAP_X',
        y='UMAP_Y',
        hue='Cluster_ID',
        data=umap_plot_df[umap_plot_df['Highlight']],
        legend=False,
        palette=palette,
        edgecolor='black',
        linewidth=1,
        alpha=1,
        s=40
    )

    texts = []
    for _, row in umap_plot_df[umap_plot_df['Highlight']].iterrows():
        texts.append(
            plt.text(row['UMAP_X'], row['UMAP_Y'], row['Sequence_ID'],
                    fontsize=8, color='black')
        )
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

    plt.title(f'UMAP Projection of Protein Embeddings with K-Means Clusters (k={num_unique_clusters})')
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    plt.grid(True)

    output_path = os.path.join(output_dir, 'umap_clusters_2d_highlighted.png')
    plt.savefig(output_path)
    plt.show()
    print(f"2D UMAP cluster plot with highlights saved to {output_path}")

    

if __name__ == '__main__':
    #embedding_csv_ref = '/media/Data/qichen/project/11.forEmb/caspedia_emb/caspedia_type2_5_and_6_emb.csv'
    #embedding_csv_query = '/media/Data/qichen/project/11.forEmb/ReadyQuery_emb/ReadyQuery_emb.csv'
    #output_dir = '/home/hebeibei/Work/dl/code/crispr_part/data/kmeans_results'

    ap = ap.ArgumentParser(description="K-Means Clustering for Protein Embeddings")
    ap.add_argument('-q','--embedding_csv_query', type=str, 
                    default="/media/Data/qichen/project/11.forEmb/caspedia_emb/caspedia_type2_5_and_6_emb.csv",
                    help='Path to the query embeddings CSV file.')
    ap.add_argument('-r','--embedding_csv_ref', type=str, 
                    default="/media/Data/qichen/project/11.forEmb/ReadyQuery_emb/ReadyQuery_emb.csv",
                    help='Path to the reference embeddings CSV file.')
    ap.add_argument('-o','--output_dir', type=str, required=True,
                    help='Directory to save output files (default: current directory).')
    ap.add_argument('-c','--n_clusters', type=int, required=True,
                    help='Number of clusters for K-Means.')
    args = ap.parse_args()

    embedding_csv_query = args.embedding_csv_query
    embedding_csv_ref = args.embedding_csv_ref
    output_dir = args.output_dir
    n_clusters = args.n_clusters

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    embeddings_query = pd.read_csv(embedding_csv_query, sep='\t', index_col=0)
    embeddings_ref = pd.read_csv(embedding_csv_ref, sep='\t', index_col=0)
    combined_embeddings_df = pd.concat([embeddings_query, embeddings_ref])

    cluster_labels, kmeans_final_model, final_metrics = perform_kmeans_clustering(
        embeddings_df=combined_embeddings_df,
        n_clusters=n_clusters,
        output_dir=output_dir
    )

    if cluster_labels is not None:
        # Assign clusters back to original protein IDs
        cluster_assignments = pd.DataFrame({
            'Sequence_ID': combined_embeddings_df.index,
            'Cluster_ID': cluster_labels
        })

        # Save assignments
        assignments_path = os.path.join(output_dir, 'protein_cluster_assignments.tsv')
        cluster_assignments.to_csv(assignments_path, sep='\t', index=False)
        print(f"Cluster assignments saved to {assignments_path}")

        # Inspect cluster sizes
        print("\nCluster sizes:")
        print(cluster_assignments['Cluster_ID'].value_counts().sort_index())

    ### Visualization
    clustered_png = os.path.join(output_dir, 'umap_clusters_2d_final.png')
    no_clustered_png = os.path.join(output_dir, 'umap_no_clusters_2d_final.png')
    run_visualization(embedding_csv_ref=embedding_csv_ref,
                      embedding_csv_query=embedding_csv_query,
                      clusters_filepath=assignments_path,
                      visualization_png=clustered_png,
                      no_clustered_png=no_clustered_png
                      )
