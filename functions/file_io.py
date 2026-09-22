import os, json
import typing
import anndata
import matplotlib.figure
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages

import shelve




def saveGenericFigures(output_dir: str, figures: typing.List[matplotlib.figure.Figure], experiment_name: str, filename_prefix: str, filename_suffix: str) -> str:
    """
    Save the a list of matplotlib figures to a PDF file. The filename is formatted as `{filename_prefix}_{experiment_name}_{filename_suffix}`.pdf.

    Parameters
    ----------
    - `output_dir` : **str**
        Path to the output directory.
    - `figures` : list of **matplotlib.figure.Figure**
        The `matplotlib.Figure.figure` figure(s) to save.
    - `experiment_name` : **str**
        Name of the experiment.
    - `filename_prefix` : **str**
        Prefix to the filename.
    - `filename_suffix` : **str**
        Suffix to the filename.

    Returns
    ----------
    - `pdf_path` : **str**
        Path to the saved PDF report file.
    """

    # filename = experiment_name + "_" + filename_prefix + "_" + filename_suffix + ".pdf"
    if filename_prefix != "":
        filename_prefix = filename_prefix + "_"
    if filename_suffix != "":
        filename_suffix = "_" + filename_suffix
    filename = filename_prefix + experiment_name + filename_suffix + ".pdf"
    pdf_file = os.path.join(output_dir, filename)
    pdf_file = os.path.normpath(pdf_file)

    with PdfPages(pdf_file) as pdf:
        for fig in figures:
            pdf.savefig(fig)

    return pdf_file




def loadDataTables(cell_metadata_csv: str, cell_by_gene_csv: str, prioritized_hdf5: typing.Union[str, None], experiment_name: str, additional_metadata: dict, drop_zero_transcript_cells: bool = False) -> anndata.AnnData:
    """
    Load the cell metadata and cell-by-gene data matrix. Return an AnnData object.
    
    Parameters
    ----------
    - `cell_metadata_csv` : **str**
        Path to the cell metadata CSV file.
    - `cell_by_gene_csv` : **str**
        Path to the cell-by-gene data matrix file.
    - `prioritized_hdf5`: **str**
        Path to the HDF5 to load instead of the CSV files. Only works if the HDF5 file is valid
    - `experiment_name` : **str**
        Name of the experiment.
    - `additional_metadata`: **dict**
        Additional metadata to add to anndata.uns["info"]["metadata"]
    
    Returns
    -------
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix. Additional **`info`** attribute is added to the AnnData object to store experiment specific information.
    """

    info = dict(
        experiment_name=experiment_name,
        cell_metadata_csv=cell_metadata_csv,
        cell_by_gene_csv=cell_by_gene_csv,
        metadata=additional_metadata,
    )


    # Check if prioritized_hdf5 is provided
    if prioritized_hdf5 is not None:
        if os.path.exists(prioritized_hdf5):
            print(f"Loading HDF5 file: {prioritized_hdf5}")
            adata = sc.read_h5ad(prioritized_hdf5)
            return adata
        else:
            print(f"HDF5 file not found: {prioritized_hdf5}")
            print("Loading data from CSV files instead...")


    cell_metadata_df = pd.read_csv(
        cell_metadata_csv,
        sep=",",
        header=0,
        index_col=0,
    )
    cell_by_gene_df = pd.read_csv(
        cell_by_gene_csv,
        sep=",",
        header=0,
        index_col=0,
    )

    # Validation: the cell_by_gene matrix must contain at least one data column
    # (genes, blanks, or intensity values). An index-only file cannot form a valid
    # AnnData matrix and would crash later with an opaque pandas/anndata error.
    if len(cell_by_gene_df.columns) == 0:
        raise ValueError(
            f"cell_by_gene file '{cell_by_gene_csv}' contains no data columns "
            "(only the index). The pipeline requires at least one gene, blank, or intensity column."
        )

    cell_metadata_df.index = [str(x) for x in cell_metadata_df.index]
    cell_by_gene_df.index = [str(x) for x in cell_by_gene_df.index]
    info["detected_cells"] = len(cell_metadata_df)

    # Identify intensity columns (raw histology intensity values, NOT gene transcript counts)
    # A column whose name contains the substring 'intensity_' is an intensity value.
    intensity_columns = [col for col in cell_by_gene_df.columns if "intensity_" in col]
    gene_blank_columns = [col for col in cell_by_gene_df.columns if "intensity_" not in col]
    has_transcript_data = len(gene_blank_columns) > 0

    # Optionally filter for cells with at least one transcript.
    # drop_zero_transcript_cells=True preserves the legacy behavior of dropping
    # zero-transcript cells before any yaml filtering runs. When False (default),
    # zero-transcript cells are kept and filtering is left to the yaml
    # filtering_procedure. Intensity-only datasets always keep all cells.
    if has_transcript_data and drop_zero_transcript_cells:
        cell_metadata_df = cell_metadata_df[cell_metadata_df["transcript_count"] > 0]
        cell_by_gene_df["sum"] = cell_by_gene_df[gene_blank_columns].sum(axis=1)
        cell_by_gene_df = cell_by_gene_df[cell_by_gene_df["sum"] > 0]
        cell_by_gene_df.drop(columns=["sum"], inplace=True)

    # Only keep the cells that are present in both the cell metadata and cell-by-gene data matrix
    cellOverlap = list(set(cell_metadata_df.index) & set(cell_by_gene_df.index))
    cellOverlap.sort()
    cell_metadata_df = cell_metadata_df.loc[cellOverlap]
    cell_by_gene_df = cell_by_gene_df.loc[cellOverlap]

    if has_transcript_data:
        # Transcript counts are computed over gene/blank columns only (intensity excluded)
        cell_metadata_df["transcript_count"] = cell_by_gene_df[gene_blank_columns].sum(axis=1)
        cell_metadata_df["genes_count"] = (cell_by_gene_df[gene_blank_columns] > 0).sum(axis=1)
        total_transcripts = cell_by_gene_df[gene_blank_columns].sum().sum()
    else:
        # Intensity-only dataset: no transcripts exist, keep schema consistency with zeros
        cell_metadata_df["transcript_count"] = 0
        cell_metadata_df["genes_count"] = 0
        total_transcripts = 0

    transcripts_coords = cell_metadata_df[["center_x", "center_y"]]

    info["detected_cells_with_transcripts"] = len(cell_metadata_df)
    info["total_intracellular_transcripts"] = total_transcripts
    info["intensity_only"] = not has_transcript_data

    adata = anndata.AnnData(
        X=cell_by_gene_df,
        obs=cell_metadata_df,
        obsm={"spatial": transcripts_coords.to_numpy()},
    )
    adata.var["Blanks"] = ["Blank" in x for x in adata.var.index]
    adata.var["Intensity"] = ["intensity_" in x for x in adata.var.index]
    adata.var["Genes"] = [not x and not y for x, y in zip(adata.var["Blanks"], adata.var["Intensity"])]

    # Per-cell intensity channel info (additive obs columns; null value = channel absent for that cell)
    if len(intensity_columns) > 0:
        intensity_channels = sorted({col.split("__")[1] for col in intensity_columns if "__" in col})
        info["intensity_channels"] = intensity_channels

        intensity_values = cell_by_gene_df[intensity_columns]
        # A channel is present for a cell if it has at least one non-null, non-NaN
        # intensity value across all of its stat columns
        channel_present = pd.DataFrame(
            {ch: intensity_values[[c for c in intensity_columns if c.split("__")[1] == ch]].notna().any(axis=1)
             for ch in intensity_channels},
            index=cell_by_gene_df.index,
        )
        cell_metadata_df["intensity_channel_count"] = channel_present.sum(axis=1).astype(int)
        cell_metadata_df["intensity_channel"] = [
            ",".join(intensity_channels[i] for i in np.flatnonzero(row)) if row.any() else None
            for row in channel_present.to_numpy()
        ]
        adata.obs["intensity_channel_count"] = cell_metadata_df["intensity_channel_count"]
        adata.obs["intensity_channel"] = cell_metadata_df["intensity_channel"]
    else:
        info["intensity_channels"] = []

    adata.uns["info"] = info
    adata.uns["spatial"] = { experiment_name: {} }
    adata.uns["spatial"][experiment_name]["images"] = {}

    return adata



def saveDeltaUsableTranscriptsReport(output_dir: str, diff_output: dict) -> typing.Tuple[str, str]:
    """
    Save the difference in usable transcripts report to a PDF file.

    Parameters
    ----------
    - `output_dir` : **str**
        Path to the output directory. The names are generated based on the experiment name.
    - `diff_output` : **dict**
        The output from the `deltaUsableTranscripts()` function.

    Returns
    ----------
    - `pdf_path` : **str**
        Path to the saved PDF report file.
    - `json_path` : **str**
        Path to the saved JSON report file.
    """

    filename = diff_output["experiment_name"] + "_transcripts_metadata_QC.pdf"
    pdf_file = os.path.join(output_dir, filename)
    pdf_file = os.path.normpath(pdf_file)

    with PdfPages(pdf_file) as pdf:
        pdf.savefig(diff_output["axes"][0].figure)

    # Save the difference data to a JSON file
    filename = diff_output["experiment_name"] + "_transcripts_metadata_QC.json"
    json_file = os.path.join(output_dir, filename)
    json_file = os.path.normpath(json_file)

    with open(json_file, "w") as f:
        json.dump(diff_output["diff"], f, indent=4)

    return pdf_file, json_file


def saveFilteringCellMetadataReport(output_dir: str, cell_metadata_stats: list, experiment_name: typing.Union[str, None]) -> str:
    """
    Save the cell filtering report to a PDF file.

    Parameters
    ----------
    - `output_dir` : **str**
        Path to the output directory. The names are generated based on the experiment name.
    - `cell_metadata_stats` : **matplotlib.pyplot**
        The `stats` field from the output of `apply_filters()` function.

    Returns
    ----------
    - `pdf_path` : **str**
        Path to the saved PDF report file.
    """

    if experiment_name is None:
        experiment_name = "exp"

    filename = experiment_name + "_cell_metadata.pdf"
    pdf_file = os.path.join(output_dir, filename)
    pdf_file = os.path.normpath(pdf_file)

    with PdfPages(pdf_file) as pdf:
        for stat in cell_metadata_stats:
            pdf.savefig(stat['plt'].figure)

    return pdf_file



def exportMapMyCellsInput(output_dir: str, adata: anndata.AnnData, gene_panel: dict) -> str:
    """
    Export the H5AD file for MapMyCells processing.

    Parameters
    ----------
    - `output_dir` : **str**
        Path to the output directory.
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix.
    - `gene_panel` : **dict**
        The dictionary provided in the pipeline configuration file.

    Returns
    ----------
    - `h5ad_path` : **str**
        Path to the saved H5AD file.
    """

    gene_panel_df = pd.read_csv(
        gene_panel["gene_panel_csv"],
        sep=",",
        header=0,
        index_col=0,
        dtype=str,
        usecols=[gene_panel["gene_panel_vizgen_gene_header"], gene_panel["gene_panel_ensemble_id_header"]],
    )

    adata_formatted = adata.copy()

    # No genes to export (e.g. intensity-only dataset): skip MapMyCells export
    if adata_formatted.var["Genes"].sum() == 0:
        print("No gene columns available for MapMyCells export. Skipping.")
        return ""

    adata_formatted = adata_formatted[:, adata_formatted.var["Genes"]]
    adata_formatted.obs.drop(columns=adata_formatted.obs.columns, inplace=True)

    # Rename the genes from Vizgen to Ensembl
    adata_formatted = anndata.AnnData(
        X=adata_formatted.X,
        obs=adata_formatted.obs,
        var=adata_formatted.var["Genes"].to_frame(),
    )

    # Map the gene names to Ensembl IDs
    adata_formatted.var.index = gene_panel_df.loc[adata_formatted.var.index, gene_panel["gene_panel_ensemble_id_header"]].values

    del adata_formatted.var["Genes"]

    h5ad_path = os.path.join(output_dir, adata.uns["info"]["experiment_name"] + "_MapMyCells_input.h5ad")
    h5ad_path = os.path.normpath(h5ad_path)

    adata_formatted.write(h5ad_path, compression="gzip")

    return h5ad_path



def exportAnnData(output_dir: str, adata: anndata.AnnData, experiment_name: str, filename_suffix: str) -> str:
    """
    Export the AnnData object to a HDF5 file.

    Parameters
    ----------
    - `output_dir` : **str**
        Path to the output directory.
    - `adata` : **anndata.AnnData**
        An AnnData object containing the cell metadata and cell-by-gene data matrix.
    - `experiment_name` : **str**
        Name of the experiment.

    Returns
    ----------
    - `hdf5_path` : **str**
        Path to the saved HDF5 file.
    """

    hdf5_path = os.path.join(output_dir, experiment_name + "_" + filename_suffix +".hdf5")
    hdf5_path = os.path.normpath(hdf5_path)

    adata.write(hdf5_path)

    return hdf5_path


def exportCellClusteringAnnotationsSpatialReport(output_dir: str, figures: typing.List[matplotlib.figure.Figure], experiment_name: str) -> str:
    """
    Export the spatial MapMyCells annotations report to a PDF file.

    Parameters
    ----------
    - `output_dir` : **str**
        Path to the output directory.
    - `figures` : list of **matplotlib.figure.Figure**
        The figure(s) generated by `plotSpatialMapMyCellsHierarchy()` and/or `plotSpatialLeidenClusters()`.
    - `experiment_name` : **str**
        Name of the experiment.

    Returns
    ----------
    - `pdf_path` : **str**
        Path to the saved PDF report file.
    """

    filename = experiment_name + "_clustered_celltyped_spatial.pdf"
    pdf_file = os.path.join(output_dir, filename)
    pdf_file = os.path.normpath(pdf_file)

    with PdfPages(pdf_file) as pdf:
        for fig in figures:
            pdf.savefig(fig)

    return pdf_file