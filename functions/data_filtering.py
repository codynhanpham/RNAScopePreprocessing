import os, json, math
import typing
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot
import matplotlib.figure

from functions import helpers


def adata_metadata(adata: anndata.AnnData) -> str:
    """
    Return the basic cells and transcripts metadata of the AnnData object.

    Parameters
    ----------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix.

    Returns
    -------
    - `metadata` : **str**
        Metadata of the AnnData object.
    """
    
    metadata = f"AnnData object with {adata.shape[0]} cells and {adata.shape[1]} genes (with Blanks if exist).\n"

    gene_count = adata.var["Genes"].sum()
    blank_count = adata.var["Blanks"].sum()
    metadata += f"-> Genes: {gene_count} | Blanks: {blank_count}\n"

    # Cell, gene, and transcript metadata:
    cell_count = adata.shape[0]
    cell_volume = adata.obs["volume"]

    metadata += f"\nCells Volume Info: ({cell_count})\n"
    metadata += f"  Min: {round(cell_volume.min(), 2)} | Max: {round(cell_volume.max(), 2)}\n  Mean: {round(np.mean(cell_volume), 2)} | Std: {round(np.std(cell_volume), 2)}\n  Median: {round(np.median(cell_volume), 2)} | SEM: {round(np.std(cell_volume) / np.sqrt(cell_count), 2)}\n"
    metadata += f"    1st  - 99th Percentile: {round(np.percentile(cell_volume, 1), 2)} - {round(np.percentile(cell_volume, 99), 2)}\n"
    metadata += f"    3rd  - 97th Percentile: {round(np.percentile(cell_volume, 3), 2)} - {round(np.percentile(cell_volume, 97), 2)}\n"
    metadata += f"    5th  - 95th Percentile: {round(np.percentile(cell_volume, 5), 2)} - {round(np.percentile(cell_volume, 95), 2)}\n"
    metadata += f"    10th - 90th Percentile: {round(np.percentile(cell_volume, 10), 2)} - {round(np.percentile(cell_volume, 90), 2)}\n"
    metadata += f"    25th - 75th Percentile: {round(np.percentile(cell_volume, 25), 2)} - {round(np.percentile(cell_volume, 75), 2)}\n"

    del cell_count, cell_volume

    transcripts = adata.X.sum(axis=1)
    transcripts = transcripts[transcripts > 0]
    total_transcripts = transcripts.sum()
    metadata += f"\nTranscripts Info: ({total_transcripts} non-zero)\n"
    metadata += f"  Min: {round(transcripts.min(), 2)} | Max: {round(transcripts.max(), 2)}\n  Mean: {round(np.mean(transcripts))} | Std: {round(np.std(transcripts), 2)}\n  Median: {round(np.median(transcripts))} | SEM: {round(np.std(transcripts) / np.sqrt(len(transcripts)), 2)}\n"
    metadata += f"    1st  - 99th Percentile: {round(np.percentile(transcripts, 1), 2)} - {round(np.percentile(transcripts, 99), 2)}\n"
    metadata += f"    3rd  - 97th Percentile: {round(np.percentile(transcripts, 3), 2)} - {round(np.percentile(transcripts, 97), 2)}\n"
    metadata += f"    5th  - 95th Percentile: {round(np.percentile(transcripts, 5), 2)} - {round(np.percentile(transcripts, 95), 2)}\n"
    metadata += f"    10th - 90th Percentile: {round(np.percentile(transcripts, 10), 2)} - {round(np.percentile(transcripts, 90), 2)}\n"
    metadata += f"    25th - 75th Percentile: {round(np.percentile(transcripts, 25), 2)} - {round(np.percentile(transcripts, 75), 2)}\n"


    del transcripts, total_transcripts

    gene_per_cell = (adata.X > 0).sum(axis=1)
    # Subset adata to only include non Blanks genes
    adata = adata[:, adata.var["Genes"]]
    metadata += f"\nGenes per Cell Info: ({len(adata.var['Genes'])})\n"
    metadata += f"  Min: {gene_per_cell.min()} | Max: {gene_per_cell.max()}\n  Mean: {round(np.mean(gene_per_cell))} | Std: {round(np.std(gene_per_cell), 2)}\n  Median: {round(np.median(gene_per_cell))} | SEM: {round(np.std(gene_per_cell) / np.sqrt(len(gene_per_cell)), 2)}\n"
    metadata += f"    1st  - 99th Percentile: {round(np.percentile(gene_per_cell, 1), 2)} - {round(np.percentile(gene_per_cell, 99), 2)}\n"
    metadata += f"    3rd  - 97th Percentile: {round(np.percentile(gene_per_cell, 3), 2)} - {round(np.percentile(gene_per_cell, 97), 2)}\n"
    metadata += f"    5th  - 95th Percentile: {round(np.percentile(gene_per_cell, 5), 2)} - {round(np.percentile(gene_per_cell, 95), 2)}\n"
    metadata += f"    10th  - 90th Percentile: {round(np.percentile(gene_per_cell, 10), 2)} - {round(np.percentile(gene_per_cell, 90), 2)}\n"
    metadata += f"    25th - 75th Percentile: {round(np.percentile(gene_per_cell, 25), 2)} - {round(np.percentile(gene_per_cell, 75), 2)}\n"


    del gene_per_cell

    return metadata

    

def plot_cell_metadata(adata: anndata.AnnData, metadata_str: typing.Union[str, None], title="Cell Metadata") -> matplotlib.figure.Figure:
    import matplotlib.pyplot as plt
    plt = helpers.applyBasicPlotStyle(plt)
    
    mainfig = plt.figure(layout='constrained', figsize=(26, 14))
    subfigs = mainfig.subfigures(1, 2, wspace=0.07, width_ratios=[1.8, 1])

    axs1 = subfigs[0].subplots(2, 2, sharex=False, sharey=False)
    axs2 = subfigs[1].subplots(1, 1, sharex=False, sharey=False)

    # Plot the histogram of transcript counts per cell
    # Need to calculate the sum of transcripts from X
    transcript_counts = adata.X.sum(axis=1)
    axs1[0, 0].hist(transcript_counts, bins=100, color='#997543', alpha=1)
    axs1[0, 0].set_title("Transcripts per Cell")
    axs1[0, 0].set_xlabel("Transcript count")
    axs1[0, 0].set_ylabel("Number of cells")


    # Remove cells with zero transcript counts
    transcript_counts_nonzero = transcript_counts[transcript_counts > 0]
    axs1[0, 1].hist(transcript_counts_nonzero, bins=100, color='#997543', alpha=1)
    axs1[0, 1].set_title("Transcripts per Cell (Non-zero)")
    axs1[0, 1].set_xlabel("Transcript count")
    axs1[0, 1].set_ylabel("Number of cells")

    # Plot the histogram of cell volumes
    cell_volumes = adata.obs["volume"]
    axs1[1, 0].hist(cell_volumes, bins=100, color='#9F343D', alpha=1)
    axs1[1, 0].set_title("Cell Volumes Distribution")
    axs1[1, 0].set_xlabel("Cell volume (μm³)")
    axs1[1, 0].set_ylabel("Number of cells")

    # Make the last plot a scatter plot of cell volumes vs transcript counts
    axs1[1, 1].scatter(cell_volumes, transcript_counts, color='#0065FF', alpha=0.4, s=4)
    axs1[1, 1].set_title("Cell volumes vs. transcript counts")
    axs1[1, 1].set_xlabel("Cell volume (μm³)")
    axs1[1, 1].set_ylabel("Transcript counts")

    if metadata_str is not None:
        axs2.text(0, 0.5, metadata_str, fontsize=14, ha='left', va='center', wrap=True, fontweight='medium', fontfamily='monospace')
        axs2.axis('off')

    mainfig.suptitle(title, fontsize=22, fontweight='bold')

    return mainfig
    

def filter_cell_volume(adata: anndata.AnnData, min_volume: int, max_volume: int, format: typing.Literal["percentile", "literal"]) -> anndata.AnnData:
    """
    Filter the cells based on their volume.

    Parameters
    ----------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix.
    - `min_volume` : **int**
        Minimum volume of the cell.
    - `max_volume` : **int**
        Maximum volume of the cell.

    Returns
    -------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix after filtering the cells based on their volume.
    """
    adata_filtered = adata.copy()

    if format == "percentile":
        min_volume = np.percentile(adata_filtered.obs["volume"], min_volume)
        max_volume = np.percentile(adata_filtered.obs["volume"], max_volume)

    adata_filtered = adata_filtered[(adata_filtered.obs["volume"] >= min_volume) & (adata_filtered.obs["volume"] <= max_volume)]

    return adata_filtered


def filter_blank_thresholding(adata: anndata.AnnData, set_threshold: typing.Literal["blank_max", "blank_min", "blank_mean", "blank_median"], filter_action: typing.Literal["remove", "subtract"]) -> anndata.AnnData:
    """
    Filter the cells based on the blank thresholding.

    Parameters
    ----------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix.
    - `set_threshold` : **str**
        Method to determine the threshold: `blank_max`, `blank_min`, `blank_mean`, `blank_median`
    - `filter_action` : **str**
        Action to take when the threshold is exceeded: `remove` (set the transcript count to 0), `subtract` (subtract the current transcript count by the threshold)
    """
    adata_filtered = adata.copy()

    blanks = adata_filtered.var["Blanks"]
    adata_blanks_only = adata_filtered[:, blanks == True]

    # If no blanks are present, return the original AnnData object
    if adata_blanks_only.shape[1] == 0:
        return adata_filtered
    
    # Make an array of the blank thresholds for each cell, based on the set_threshold and the blanks in the AnnData object
    if set_threshold == "blank_max":
        blanks = adata_blanks_only.X.max(axis=1)
    elif set_threshold == "blank_min":
        blanks = adata_blanks_only.X.min(axis=1)
    elif set_threshold == "blank_mean":
        blanks = adata_blanks_only.X.mean(axis=1)
        blanks = blanks.astype(int)
    elif set_threshold == "blank_median":
        blanks = adata_blanks_only.X.median(axis=1)
        blanks = blanks.astype(int)

    # Repeat the blank threshold for every genes in cell (to match the shape of the AnnData object)
    blanks = np.tile(blanks, (adata_filtered.shape[1], 1))
    blanks = blanks.T
    blanks = blanks.astype(int)

    # Take action based on the filter_action
    if filter_action == "remove":
        # If gene transcript count is less than the threshold, set the transcript count to 0
        adata_filtered.X[adata_filtered.X < blanks] = 0
    elif filter_action == "subtract":
        # If gene transcript count is less than the threshold, subtract the current transcript count by the blank threshold
        adata_filtered.X = adata_filtered.X - blanks

    return adata_filtered


def filter_transcript_count(adata: anndata.AnnData, min_transcript: int, max_transcript: int, format: typing.Literal["percentile", "literal"]) -> anndata.AnnData:
    """
    Filter the cells based on the number of transcripts.

    Parameters
    ----------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix.
    - `min_transcript` : **int**
        Minimum number of transcripts.
    - `max_transcript` : **int**
        Maximum number of transcripts.

    Returns
    -------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix after filtering the cells based on the number of transcripts.
    """
    adata_filtered = adata.copy()

    if format == "percentile":
        transcripts = adata_filtered.X.sum(axis=1)
        transcripts = transcripts[transcripts > 0]
        min_transcript = np.percentile(transcripts, min_transcript)
        max_transcript = np.percentile(transcripts, max_transcript)

    adata_filtered = adata_filtered[(adata_filtered.X.sum(axis=1) >= min_transcript) & (adata_filtered.X.sum(axis=1) <= max_transcript)]

    return adata_filtered


def filter_gene_per_cell(adata: anndata.AnnData, min_genes: int, max_genes: int, format: typing.Literal["percentile", "literal"]) -> anndata.AnnData:
    """
    Filter the cells based on the number of genes expressed.

    Parameters
    ----------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix.
    - `min_genes` : **int**
        Minimum number of genes expressed.
    - `max_genes` : **int**
        Maximum number of genes expressed.

    Returns
    -------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix after filtering the cells based on the number of genes expressed.
    """
    adata_filtered = adata.copy()

    if format == "percentile":
        min_genes = np.percentile((adata_filtered.X > 0).sum(axis=1), min_genes)
        max_genes = np.percentile((adata_filtered.X > 0).sum(axis=1), max_genes)

    # Subset to only include non Blanks genes
    genes = adata_filtered.var["Genes"]

    min_genes = max(0, min_genes)

    # sum of unique, non zero genes per each cell, of genes that are not Blanks
    adata_filtered = adata_filtered[((adata_filtered.X[:, genes] > 0).sum(axis=1) >= min_genes) & ((adata_filtered.X[:, genes] > 0).sum(axis=1) <= max_genes)]

    return adata_filtered


def normalize_transcript_by_volume(adata: anndata.AnnData) -> anndata.AnnData:
    """
    Normalize the transcript count by the largest cell volume.

    Parameters
    ----------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix.

    Returns
    -------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix after normalizing the transcript count by the cell volume.
    """
    adata_normalized = adata.copy()

    max_volume = adata_normalized.obs["volume"].max()

    # Find the scale factor and scale up the transcript count
    scale_factor = max_volume / adata_normalized.obs["volume"].to_numpy()

    # format the scale_factor to var length (number of genes) dimemsion
    scale_factor = scale_factor[:, np.newaxis]
    transcript_matrix = adata_normalized.X
    # Scale up the transcript count
    transcript_matrix = transcript_matrix * scale_factor
    # transcript_matrix should be int
    transcript_matrix = transcript_matrix.astype(int)

    # Overwrite the transcript count matrix
    adata_normalized.X = transcript_matrix

    return adata_normalized




def apply_filters(adata: anndata.AnnData, filtering_procedure: dict, retain_raw_transcript_count: bool) -> dict:
    """
    Read through the filtering procedure and apply the filters to the AnnData object.

    Parameters
    ----------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix.
    - `filtering_procedure` : **dict**
        Dictionary containing the filtering procedure. Provided through the pipeline configuration.

    Returns
    -------
    - `dict` : Output dictionary with the following fields:
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix after applying the filters.
    - `stats` : **list**
        List of dictionaries containing the statistics and figures of the filtered cells of every filtering step.
        - `plt` : **matplotlib.figure.Figure**
            Matplotlib plot of the cell metadata.
        - `metadata` : **str**
            Metadata of the AnnData object.
    """

    adata_filtered = adata.copy()
    if retain_raw_transcript_count:
        # Save the raw data in .raw; when anndata is sliced, the raw data is also sliced https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.raw.html
        adata_filtered.raw = adata.copy()

    adata_rows_len = adata_filtered.shape[0]

    stats = []

    filtering_params_str = f"Filtering Procedure:\n{json.dumps(filtering_procedure, indent=None)}\n\n\n"
    metadata_str = filtering_params_str + adata_metadata(adata_filtered)
    plt = plot_cell_metadata(adata_filtered, metadata_str, title=f"[{adata_filtered.uns['info']['experiment_name']}] Raw Cell Metadata")
    stats.append(dict(plt=plt, metadata=metadata_str))


    for i, filter_step in enumerate(filtering_procedure):
        filter_name = list(filter_step.keys())[0]
        filter_params = list(filter_step.values())[0]
    
        print(f"\t[{i+1}/{len(filtering_procedure)}] Applying filter: '{filter_name}' with parameters: {filter_params}")

        if filter_name == "cell_volume":
            adata_filtered = filter_cell_volume(adata_filtered, filter_params["min"], filter_params["max"], filter_params["format"])
        elif filter_name == "blank_thresholding":
            adata_filtered = filter_blank_thresholding(adata_filtered, filter_params["set_threshold"], filter_params["filter_action"])
        elif filter_name == "transcript_count":
            adata_filtered = filter_transcript_count(adata_filtered, filter_params["min"], filter_params["max"], filter_params["format"])
        elif filter_name == "gene_per_cell":
            adata_filtered = filter_gene_per_cell(adata_filtered, filter_params["min"], filter_params["max"], filter_params["format"])
        elif filter_name == "normalize_transcript_by_volume":
            if filter_params:
                adata_filtered = normalize_transcript_by_volume(adata_filtered)

        percentage = round((adata_filtered.shape[0] / adata_rows_len) * 100, 2)
        print(f"\tFilter '{filter_name}' applied successfully. Current size: {adata_filtered.shape[0]} cells ({percentage}% of pre-filtered cells)")

        current_filtering_params = f"[{i+1}/{len(filtering_procedure)}] ({filter_name}) Parameters:\n{filter_params}\n\n\n"
        metadata_str = current_filtering_params + adata_metadata(adata_filtered)
        plt = plot_cell_metadata(adata_filtered, metadata_str, title=f"[{adata_filtered.uns['info']['experiment_name']}] Cell Metadata after '{filter_name}' filter ({i+1}/{len(filtering_procedure)})")
        stats.append(dict(plt=plt, metadata=metadata_str))


    if retain_raw_transcript_count and any(filter_step.get("normalize_transcript_by_volume") for filter_step in filtering_procedure):
        temp_adata = adata_filtered.copy()
        temp_adata.X = adata_filtered.raw.X
        # Metadata for the filtered cells with raw transcript counts
        metadata_str = adata_metadata(temp_adata)
        plt = plot_cell_metadata(temp_adata, metadata_str, title=f"[{adata_filtered.uns['info']['experiment_name']}] Filtered Cells Metadata with Raw (Unnormalized) Transcript Counts")
        stats.append(dict(plt=plt, metadata=metadata_str))

    adata_filtered = adata_filtered.copy()
    adata_filtered.uns["filtered"] = True

    return dict(adata=adata_filtered, stats=stats)



def plotSpatialCellFilteringStats(pre_filter_adata: anndata.AnnData, post_filter_adata: anndata.AnnData, rotation_deg: typing.Union[int, float], flip_x: bool) -> matplotlib.figure.Figure:
    """
    Plot the spatial statistics of the cells before and after filtering.

    Parameters
    ----------
    - `pre_filter_adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix before filtering.
    - `post_filter_adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix after filtering.
    - `rotation_deg` : **int** or **float**
        The degree of rotation to apply to the spatial coordinates

    Returns
    -------
    - `plt` : **matplotlib.figure.Figure**
        Matplotlib plot of the spatial statistics of the cells before and after filtering.
    """
    import matplotlib.pyplot as plt

    plt = helpers.applyBasicPlotStyle(plt)

    mainfig = plt.figure(layout='tight', figsize=(26, 10))

    # 3 axes for the 3 plots
    axs = mainfig.subplots(1, 3, sharex=False, sharey=False)

    # Plot of pre-filter (raw), post-filter (filtered), and the difference (got filtered out) cells
    pre_filter_coords = pre_filter_adata.obsm["spatial"]
    pre_filter_coords = helpers.rotateCoordinates(pre_filter_coords, rotation_deg, flip_x)

    # ax1: burlywood color for the pre-filtered cells
    # ax2: mediumaquamarine color, alpha 1, for the passed filtered + linen color, alpha 0.5, for the filtered out cells
    # ax3: lightcoral color for the filtered out cells, alpha 1, and the passed filtered cells are linen color, alpha 0.5

    axs[0].scatter(pre_filter_coords[:, 0], pre_filter_coords[:, 1], color='burlywood', alpha=0.9, s=0.01)
    axs[0].set_title("Pre-filtered Cells")
    axs[0].axis('off')
    axs[0].set_aspect('equal', 'box')

    post_filter_coords = post_filter_adata.obsm["spatial"]
    # Find the cells that got filtered out
    
    post_filter_coords = helpers.rotateCoordinates(post_filter_coords, rotation_deg, flip_x)
    filtered_out_coords = pre_filter_coords[~np.isin(pre_filter_coords, post_filter_coords).all(1)]

    axs[1].scatter(filtered_out_coords[:, 0], filtered_out_coords[:, 1], color='linen', alpha=0.3, s=0.01)
    axs[1].scatter(post_filter_coords[:, 0], post_filter_coords[:, 1], color='mediumaquamarine', alpha=0.9, s=0.01)
    axs[1].set_title("Post-filtered Cells")
    axs[1].axis('off')
    axs[1].set_aspect('equal', 'box')

    axs[2].scatter(post_filter_coords[:, 0], post_filter_coords[:, 1], color='linen', alpha=0.3, s=0.01)
    axs[2].scatter(filtered_out_coords[:, 0], filtered_out_coords[:, 1], color='lightcoral', alpha=1, s=0.01)
    axs[2].set_title("Removed Cells")
    axs[2].axis('off')
    axs[2].set_aspect('equal', 'box')

    mainfig.suptitle(f"[{post_filter_adata.uns['info']['experiment_name']}] Spatial Cell Metadata", fontsize=22, fontweight='bold')

    return mainfig