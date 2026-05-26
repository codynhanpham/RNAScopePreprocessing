import os, json
import typing
import anndata
import matplotlib.patches
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot
import matplotlib.figure
from mpl_toolkits.axes_grid1 import make_axes_locatable


from functions import helpers


def leidenClustering(adata: anndata.AnnData, clustering_params: dict) -> anndata.AnnData:
    """
    Perform Leiden clustering on the AnnData object.

    Parameters
    ----------
    - `adata` : **anndata.AnnData**
        The AnnData object to perform clustering on.
    - `clustering_params` : **dict**
        The parameters for the Leiden clustering specified in the pipeline configuration.

    Returns
    ----------
    - `adata` : **anndata.AnnData**
        The AnnData object with the Leiden clustering results added.
    """
    
    original_adata = adata.copy()
    adata_clustered = adata.copy()

    # Leiden Clustering
    sc.pp.normalize_total(adata_clustered)
    sc.pp.log1p(adata_clustered)
    sc.pp.scale(adata_clustered, max_value=10)
    sc.tl.pca(adata_clustered, svd_solver='arpack')
    sc.pp.neighbors(adata_clustered, n_neighbors=clustering_params["n_neighbors"], n_pcs=clustering_params["n_pcs"], method="umap")
    sc.tl.leiden(adata_clustered, resolution=clustering_params["resolution"], flavor="leidenalg", n_iterations=3, directed=False, key_added="leiden_clusters")
    sc.tl.umap(adata_clustered, min_dist=clustering_params["min_dist"], spread=clustering_params["spread"])

    # Save new data to the original adata object (since scanpy modifies the object in place)
    original_adata.obsm["X_umap"] = adata_clustered.obsm["X_umap"].copy()
    original_adata.obs["leiden_clusters"] = adata_clustered.obs["leiden_clusters"].copy()
    original_adata.obs["leiden_clusters"] = original_adata.obs["leiden_clusters"].astype("category")

    return original_adata



def mirrorLeidenClusteringResults(raw_adata: anndata.AnnData, adata_clustered: anndata.AnnData) -> anndata.AnnData:
    """
    Transfer the clustering results from the clustered AnnData object to the raw AnnData object.

    Parameters
    ----------
    - `raw_adata` : **anndata.AnnData**
        The raw AnnData object.
    - `adata_clustered` : **anndata.AnnData**
        The AnnData object with the clustering results.

    Returns
    ----------
    - `raw_adata` : **anndata.AnnData**
        The raw AnnData object with the clustering results transferred. Cells that were not clustered are marked as 'Unassigned (filtered)'.
    """

    # Transfer the clustering results to the original AnnData object
    # Also transfer the UMAP coordinates, match the shape of the raw_adata object by adding NaNs
    raw_adata.obs["leiden_clusters"] = "Unassigned (filtered)"
    raw_adata.obsm["X_umap"] = np.full((len(raw_adata.obs), 2), np.nan)

    for cell_id in adata_clustered.obs.index:
        if cell_id in raw_adata.obs.index:
            raw_adata.obs.loc[cell_id, "leiden_clusters"] = adata_clustered.obs.loc[cell_id, "leiden_clusters"]
            raw_adata.obsm["X_umap"][raw_adata.obs.index.get_loc(cell_id)] = adata_clustered.obsm["X_umap"][adata_clustered.obs.index.get_loc(cell_id)]

    return raw_adata



def plotSpatialLeidenClusters(adata: anndata.AnnData, rotation_deg: typing.Union[int, float], flip_x: bool) -> matplotlib.figure.Figure:
    """
    Plot the Leiden clusters on the spatial coordinates and UMAP.

    Parameters
    ----------
    - `adata` : **anndata.AnnData**
        The AnnData object with the Leiden clustering results.
    - `rotation_deg` : **int**
        The degree to rotate the spatial coordinates.
    - `flip_x` : **bool**
        Whether to flip the x-axis.

    Returns
    ----------
    - `plt` : **matplotlib.figure.Figure**
    """
    import matplotlib.pyplot as plt
    plt = helpers.applyBasicPlotStyle(plt)

    # 40 for left, 60 for right
    fig, axs = plt.subplots(1, 2, figsize=(20, 10), gridspec_kw={"width_ratios": [4, 6]})

    # Plot the spatial coordinates on the left
    ax = axs[0]
    coords = adata.obsm["spatial"]
    coords = helpers.rotateCoordinates(coords, rotation_deg, flip_x)
    leiden_colors = matplotlib.cm.get_cmap("tab20", len(adata.obs["leiden_clusters"].cat.categories))

    # legend needs artist and label
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=[leiden_colors.colors[i] for i in adata.obs["leiden_clusters"].cat.codes],
        s=0.1,
        alpha=0.8,
    )

    ax.set_title("Leiden Clusters (Spatial)", fontsize=16)
    ax.axis("off")
    ax.set_aspect("equal", "box")


    # Plot the UMAP on the right
    ax = axs[1]
    ax.scatter(
        adata.obsm["X_umap"][:, 0],
        adata.obsm["X_umap"][:, 1],
        c=[leiden_colors.colors[i] for i in adata.obs["leiden_clusters"].cat.codes],
        s=0.1,
        alpha=0.8,
        label="Leiden Clusters",
    )

    # Calculate and plot cluster centroids for UMAP plot
    for cluster in adata.obs["leiden_clusters"].cat.categories:
        cluster_coords = adata.obsm["X_umap"][adata.obs["leiden_clusters"] == cluster]
        centroid = cluster_coords.mean(axis=0)
        ax.text(centroid[0], centroid[1], str(cluster), fontsize=11, ha='center', va='center', color='black', path_effects=[matplotlib.patheffects.withStroke(linewidth=3, foreground='white')])

    ax.set_title("Leiden Clusters UMAP", fontsize=16)
    ax.axis("off")
    ax.set_aspect("equal", "box")
    
    # Legend noting the cluster ID and color
    legend_elements = [
        matplotlib.patches.Patch(color=leiden_colors.colors[i], label=f"Cluster {i}") for i in range(len(adata.obs["leiden_clusters"].cat.categories))
    ]
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12, ncol=2)


    fig.suptitle(f"[{adata.uns['info']['experiment_name']}] Leiden Clusters", fontsize=18, fontweight='bold')
    fig.tight_layout()

    return fig



def reformatMapMyCellsResults(MapMyCells_JSON_path: str, MapMyCells_CSV_path: str) -> dict:
    """
    Reformat the JSON and CSV data from MapMyCells to a more usable format.

    Parameters
    ----------
    - `MapMyCells_JSON_path` : **str**
        Path to the JSON file from MapMyCells.
    - `MapMyCells_CSV_path` : **str**
        Path to the CSV file from MapMyCells.

    Returns
    ----------
    - `MapMyCells_annotations` : **dict**
        A dictionary containing the reformatted JSON data from MapMyCells.
    """

    # File path validation, just in case...
    if not os.path.exists(MapMyCells_JSON_path):
        print("File not found!")
        exit(1)
    if not MapMyCells_JSON_path.endswith(".json"):
        print("File is not JSON")
        exit(1)

    MapMyCells_JSON_path = os.path.abspath(MapMyCells_JSON_path)

    if not os.path.exists(MapMyCells_CSV_path):
        print("File not found!")
        exit(1)
    if not MapMyCells_CSV_path.endswith(".csv"):
        print("File is not CSV")
        exit(1)

    MapMyCells_CSV_path = os.path.abspath(MapMyCells_CSV_path)


    # Load the JSON file
    with open(MapMyCells_JSON_path, "r") as json_file:
        MapMyCells_annotations = json.load(json_file)

    if not all(key in MapMyCells_annotations for key in ['results', 'marker_genes', 'taxonomy_tree', 'config']):
            print("JSON data is invalid. The JSON file from MapMyCells should contain at least these keys: ['results', 'marker_genes', 'taxonomy_tree', 'config']")
            exit(1)

    # Load the CSV file
    # only headers: cell_id,class_name,subclass_name,supertype_name,cluster_name,cluster_alias
    MapMyCells_df = pd.read_csv(MapMyCells_CSV_path, sep=",", header=0, index_col=0, comment="#", dtype=str, usecols=["cell_id", "class_name", "subclass_name", "supertype_name", "cluster_name", "cluster_alias"])

    # MAKE SURE THE INDEX IS A STRING!!!
    MapMyCells_df.index = MapMyCells_df.index.astype(str)

    # Reformat MapMyCells.results data to a dictionary with cell_id as key (instead of an array) for faster lookup
    MapMyCells_results = {}
    for cell in MapMyCells_annotations['results']:
        cell_id = cell['cell_id']
        MapMyCells_results[cell_id] = cell

        # Create a new 'hierarchy' key for each cell_id, 'hierarchy' will be a dictionary with the keys ('class', 'subclass', 'superclass', 'cluster',...)
        # This will be used for looking up the taxonomy name from the taxonomy_tree later on
        hierarchy = {}
        for key in cell:
            # Probably better to do match-case here instead of if-elif, but it's Python 3.9 and no match-case yet

            if key.endswith("_CLAS"):
                hierarchy['class'] = {
                    "name": key, # for example, CCN20230722_CLAS
                    "id": cell[key]['assignment'] or "Unassigned" # for example, CS20230722_CLAS_01
                }
            elif key.endswith("_SUBC"):
                hierarchy['subclass'] = {
                    "name": key,
                    "id": cell[key]['assignment'] or "Unassigned"
                }
            elif key.endswith("_SUPT"):
                hierarchy['supertype'] = {
                    "name": key,
                    "id": cell[key]['assignment'] or "Unassigned"
                }
            elif key.endswith("_CLUS"):
                hierarchy['cluster'] = {
                    "name": key,
                    "id": cell[key]['assignment'] or "Unassigned"
                }
            elif key.endswith("_NEUR"):
                hierarchy['neurotransmitter'] = {
                    "name": key,
                    "id": cell[key]['assignment'] or "Unassigned"
                }

        # Add the hierarchy to the cell_id
        MapMyCells_results[cell_id]['hierarchy'] = hierarchy

        # Add the cluster alias to the cell_id based on the CSV data
        if str(cell_id) in MapMyCells_df.index:
            MapMyCells_results[cell_id]['cluster_alias'] = MapMyCells_df.loc[str(cell_id), 'cluster_alias']
        else:
            MapMyCells_results[cell_id]['cluster_alias'] = "Unassigned"


    MapMyCells_annotations['results'] = MapMyCells_results # Replace the original results data with the reformatted one

    return MapMyCells_annotations


def addMapMyCellsAnnotationsToAnnData(adata: anndata.AnnData, MapMyCells_annotations: dict, bootstrap_probability_threshold: float) -> anndata.AnnData:
    """
    Add the annotations from MapMyCells to the AnnData object.

    Parameters
    ----------
    - `adata` : **anndata.AnnData**
        The AnnData object to add the annotations to.
    - `MapMyCells_annotations` : **dict**
        The annotations from MapMyCells.
    - `bootstrap_probability_threshold` : **float**
        The threshold for the bootstrap probability to consider a cell type as valid.

    Returns
    ----------
    - `adata` : **anndata.AnnData**
        The AnnData object with the annotations from MapMyCells added to `obs`.
    """

    MapMyCells_results = MapMyCells_annotations['results']
    adata_added = adata.copy()

    default_colors = {
        "Unassigned": "#D3D3D3",
        "Unassigned (filtered)": "#D3D3D3",
        "Low Bootstrapping Probability (filtered)": "#FFFF00",
        "Low Transcripts (filtered)": "#FF0000",
    }

    print("\tNumber of cells in provided AnnData: ", len(adata_added.obs))
    print("\tNumber of cells in the MapMyCells results: ", len(MapMyCells_results))

    print("\n\tAdding cell type annotations to the AnnData object...")

    for key in ['class', 'subclass', 'supertype', 'cluster', 'neurotransmitter']:
        column = []

        cluster_aliases_list = [""] * len(adata_added.obs.index)
        for i, cell_id in enumerate(adata_added.obs.index):
            if cell_id in MapMyCells_results:
                hierarchy = MapMyCells_results[cell_id]['hierarchy']
                if key in hierarchy:
                    name = MapMyCells_annotations['taxonomy_tree']['name_mapper']\
                        [hierarchy[key]['name']]\
                            [hierarchy[key]['id']]['name']
                    # Check if the bootstrapping probability is above the threshold
                    if MapMyCells_results[cell_id][hierarchy[key]['name']]\
                        ['bootstrapping_probability'] >= bootstrap_probability_threshold:
                        column.append(name)
                        cluster_aliases_list[i] = str(MapMyCells_results[cell_id]['cluster_alias'])

                    else:
                        # If the bootstrapping probability is below the threshold, the classification is unreliable
                        column.append("Low Bootstrapping Probability (filtered)")
                        cluster_aliases_list[i] = "Low Bootstrapping Probability (filtered)"
                else:
                    column.append("Unassigned")
                    cluster_aliases_list[i] = "Unassigned"
            else:
                # These cells were filtered out before MapMyCells analysis
                column.append("Low Transcripts (filtered)")
                cluster_aliases_list[i] = "Low Transcripts (filtered)"

        adata_added.obs['MapMyCells_' + key] = column
        adata_added.obs['MapMyCells_' + key] = adata_added.obs['MapMyCells_' + key].astype('category')

        # Add colors to the uns metadata for plotting using Yao's color scheme
        if key == 'neurotransmitter':
            continue
        
        color_map = helpers.YAO_MAPMYCELLS_COLORS # dataframe after read_csv
        colors = []
        cm = color_map[f"{key}_color"]
        cm.index = cm.index.astype(str)

        for i, _ in enumerate(adata_added.obs.index):
            if cluster_aliases_list[i] in cm.index:
                colors.append(cm[cluster_aliases_list[i]])
            else:
                colors.append(default_colors[cluster_aliases_list[i]])

        adata_added.uns['MapMyCells_{key}_colors'.format(key=key)] = colors

    print(f"\n\tCell type annotations added with a bootstrapping probability threshold of {bootstrap_probability_threshold}.")

    # Number of reliably mapped cells for each key
    print("\n\tNumber of cells with reliable classification for each key:")
    faulty_values = ['Low Bootstrapping Probability (filtered)', 'Unassigned', 'Low Transcripts (filtered)']
    for key in ['class', 'subclass', 'supertype', 'cluster', 'neurotransmitter']:
        reliable_cells = adata_added.obs[adata_added.obs['MapMyCells_{key}'.format(key=key)].isin(faulty_values) == False]
        keep_percentage_total = len(reliable_cells) / len(adata_added.obs) * 100
        keep_percentage_filtered = len(reliable_cells) / len(MapMyCells_results) * 100
        print(f"\t\t{key}: {len(reliable_cells)} ({keep_percentage_filtered:.2f}% of filtered cells | {keep_percentage_total:.2f}% of total cells)")

    print()

    return adata_added



def plotSpatialMapMyCellsHierarchy(adata_with_mapmycells: anndata.AnnData, rotation_deg: typing.Union[int, float], flip_x: bool) -> matplotlib.figure.Figure:
    """
    Plot the MapMyCells `class`, `subclass`, `supertype`, and `cluster` annotations.

    Parameters
    ----------
    - `adata_with_mapmycells` : **anndata.AnnData**
        The AnnData object with the MapMyCells annotations in `obs`.
    - `output_dir` : **str**
        The output directory to save the plots.
    - `experiment_name` : **str**
        The name of the experiment.

    Returns
    ----------
    - `plt` : **matplotlib.figure.Figure**
    """
    import matplotlib.pyplot as plt
    plt = helpers.applyBasicPlotStyle(plt)

    # 2x2 grid of plots for class, subclass, supertype, and cluster
    fig, axs = plt.subplots(2, 2, figsize=(20, 20))

    for i, key in enumerate(['class', 'subclass', 'supertype', 'cluster']):
        ax = axs[i // 2, i % 2]
        
        this_adata = adata_with_mapmycells.copy()
        colors = this_adata.uns['MapMyCells_{key}_colors'.format(key=key)]
        
        # Filter for cells with reliable classification (not 'Low Bootstrapping Probability (filtered)', 'Unassigned', 'Low Transcripts (filtered)')
        this_adata = adata_with_mapmycells[adata_with_mapmycells.obs['MapMyCells_{key}'.format(key=key)].isin(['Low Bootstrapping Probability (filtered)', 'Unassigned', 'Low Transcripts (filtered)']) == False]

        # Also filter the colors
        colors = [colors[i] for i in range(len(colors)) if adata_with_mapmycells.obs.index[i] in this_adata.obs.index]


        coords = this_adata.obsm["spatial"]
        coords = helpers.rotateCoordinates(coords, rotation_deg, flip_x)

        # Plot the spatial coordinates
        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=colors,
            s=0.1,
            alpha=0.8,
        )

        ax.set_title(f"MapMyCells {key}", fontsize=16)
        ax.axis("off")
        ax.set_aspect("equal", "box")

    fig.suptitle(f"[{adata_with_mapmycells.uns['info']['experiment_name']}] MapMyCells Annotations", fontsize=18, fontweight='bold')
    fig.tight_layout()
    
    return fig



def plotNeurotransmitterSpatialDistribution(adata_with_mapmycells: anndata.AnnData, rotation_deg: typing.Union[int, float], flip_x: bool) -> matplotlib.figure.Figure:
    """
    Plot the MapMyCells types my its main neurotransmitter (Glutamatergic, GABAergic, or Non-Neuronal) based on the `class` level annotations.

    Parameters
    ----------
    - `adata_with_mapmycells` : **anndata.AnnData**
        The AnnData object with the MapMyCells annotations in `obs`.
    - `output_dir` : **str**
        The output directory to save the plots.
    - `experiment_name` : **str**
        The name of the experiment.

    Returns
    ----------
    - `plt` : **matplotlib.figure.Figure**
    """
    import matplotlib.pyplot as plt
    plt = helpers.applyBasicPlotStyle(plt)

    # Filter for cells with reliable classification (not 'Low Bootstrapping Probability (filtered)', 'Unassigned', 'Low Transcripts (filtered)')
    adata = adata_with_mapmycells.copy()
    adata = adata[adata.obs['MapMyCells_class'].isin(['Low Bootstrapping Probability (filtered)', 'Unassigned', 'Low Transcripts (filtered)']) == False]


    # 2x2 grid of plots for Overlayed, Glutamatergic, GABAergic, and Non-Neuronal
    fig, axs = plt.subplots(2, 2, figsize=(20, 20))


    class_names = adata.obs['MapMyCells_class'].astype('category').cat.categories.tolist()
    class_names = [name.lower() for name in class_names]

    for i, key in enumerate(['Overlayed', 'Glutamatergic', 'GABAergic', 'Non-Neuronal']):
        ax = axs[i // 2, i % 2]

        # Colors: if overlayed, use indianred for GABA, mediumseagreen for Glut, and royalblue for Non-Neuronal
        # For individual neurotransmitters, use the above colors, and lightgray for the rest
        colors = []
        for _, class_name in enumerate(class_names):
            if key == 'Overlayed':
                if class_name.endswith('gaba'):
                    colors.append('indianred')
                elif class_name.endswith('glut'):
                    colors.append('mediumseagreen')
                else:
                    colors.append('royalblue')
            else:
                if class_name.endswith('gaba') and key == 'GABAergic':
                    colors.append('indianred')
                elif class_name.endswith('glut') and key == 'Glutamatergic':
                    colors.append('mediumseagreen')
                elif key == 'Non-Neuronal' and not (class_name.endswith('gaba') or class_name.endswith('glut')):
                    colors.append('royalblue')
                else:
                    colors.append('lightgray')

        # Since colors only map to unique class names, we need to map the colors to the cells
        full_colors = []
        for cell_class in adata.obs['MapMyCells_class']:
            full_colors.append(colors[class_names.index(cell_class.lower())])

        coords = adata.obsm["spatial"]
        coords = helpers.rotateCoordinates(coords, rotation_deg, flip_x)

        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=full_colors,
            s=0.1,
            alpha=0.8,
        )

        ax.set_title(key, fontsize=16)
        ax.axis("off")
        ax.set_aspect("equal", "box")

    fig.suptitle(f"[{adata.uns['info']['experiment_name']}] Neurotransmitter Spatial Distribution", fontsize=18, fontweight='bold')
    fig.tight_layout()

    return fig



def roiDifferentialGeneByCellTypeBubbles(adata: anndata.AnnData, comparison_options: dict, roi: typing.List[dict], whitelisted_types: typing.Union[typing.List[str], None]) -> matplotlib.figure.Figure:
    """
    Generate a bubble plot for the differential gene expression between cell types in a pair of ROIs. The plot will be cell type (rows) by gene (columns) with the bubble size representing the log2 total expression and color representing the percent different between the two ROIs (range: -1 to 1; 0 will be white, -1 will be blue, 1 will be red).

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object containing the gene expression data and cell type annotations.
    comparison_options : dict
        The comparison options provided in one of `roi_differential_gene_by_celltype` in the pipeline configuration.
    roi : list of dict
        The list of ROI sets to compare, set in the pipeline configuration.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object containing the bubble plot.
    """
    import matplotlib.pyplot as plt
    plt = helpers.applyBasicPlotStyle(plt)

    roi_set_names = comparison_options["roi"]
    sort_types = comparison_options["sort_types"]
    sort_genes = comparison_options["sort_genes"]

    # Default sort_types to True if they are None
    if sort_types is None:
        sort_types = True
    if sort_genes is None:
        sort_genes = True

    roi_data = dict()
    total_cells_in_comparison = 0

    # Filter for the two ROIs:
    for roi_set_name in roi_set_names:
        # Get the ROI set values
        # find the ROI set in the pipeline config
        roi_set_values = [item["roi_set_values"] for item in roi if item["roi_set_name"] == roi_set_name][0]

        # Get the ROI set obs names
        roi_obs_names = [f"ROI__{roi_name}" for roi_name in roi_set_values]

        roi_cells = adata.obs[roi_obs_names].any(axis=1)
        roi_data[roi_set_name] = adata[roi_cells]
        # number of cells in this ROI set (roi_cells != 0)
        total_cells_in_comparison += len(roi_cells[roi_cells == True])


    print("\t\t[1/3] Generating the cell type by gene expression matrices...")

    # Cell type by gene expression matrix for each ROI:
    roi_celltype_by_gene = dict()

    # Get the cell types and Genes
    cell_types = adata.obs["MapMyCells_" + comparison_options["celltype_hierarchy_level"]].cat.categories
    cell_types = cell_types[~cell_types.isin(["Unassigned", "Unassigned (filtered)", "Low Bootstrapping Probability (filtered)", "Low Transcripts (filtered)"])]

    genes = adata.var_names[adata.var["Genes"] == True]

    # Fill the matrix for each ROI
    for roi_name in roi_set_names:
        roi_celltype_by_gene[roi_name] = pd.DataFrame(index=cell_types, columns=genes)
        # Fill the matrix by summing the expression for each gene of all cells in each cell type
        for cell_type in cell_types:
            this_cell_by_gene_expression = roi_data[roi_name][roi_data[roi_name].obs["MapMyCells_" + comparison_options["celltype_hierarchy_level"]] == cell_type].X[:, adata.var["Genes"] == True].sum(axis=0)

            roi_celltype_by_gene[roi_name].loc[cell_type] = this_cell_by_gene_expression


    print("\t\t[2/3] Calculating the differential gene expression between the two ROIs...")
    # The expression matrix
    expression_df = pd.DataFrame(index=cell_types, columns=genes)
    for cell_type in cell_types:
        for gene in genes:
            # Get the expression in each ROI
            expression_roi1 = roi_celltype_by_gene[roi_set_names[0]].loc[cell_type, gene]
            expression_roi2 = roi_celltype_by_gene[roi_set_names[1]].loc[cell_type, gene]

            # log of the total expression, or 0 if both are 0
            if expression_roi1 == 0 and expression_roi2 == 0:
                expression = 0.0
            else:
                expression:float = np.log(expression_roi1 + expression_roi2)

            expression_df.loc[cell_type, gene] = expression

    # The percent difference matrix
    percent_diff_df = pd.DataFrame(index=cell_types, columns=genes)
    for cell_type in cell_types:
        for gene in genes:
            # Get the expression in each ROI
            expression_roi1 = roi_celltype_by_gene[roi_set_names[0]].loc[cell_type, gene]
            expression_roi2 = roi_celltype_by_gene[roi_set_names[1]].loc[cell_type, gene]

            # percent difference
            if expression_roi1 == 0 and expression_roi2 == 0:
                percent_diff = 0
            else:
                percent_diff = (expression_roi2 - expression_roi1) / ((expression_roi1 + expression_roi2) / 2)

            percent_diff_df.loc[cell_type, gene] = percent_diff


    # Sort both expression df by value and percent_diff df by absolute value, get the top K overlapping genes
    total_genes = len(genes)

    expression_df = expression_df.loc[:, expression_df.mean().sort_values(ascending=False).index]
    percent_diff_df = percent_diff_df.loc[:, percent_diff_df.abs().mean().sort_values(ascending=False).index]
    top_k_genes = comparison_options["top_k_genes"] or 40
    if top_k_genes >= len(genes):
        top_k_genes = len(genes)

    # Get the top k genes pair by pair in the two dataframes: get unique until top_k_genes reached
    genes = set()
    for i in range(expression_df.shape[1]):
        genes.add(expression_df.columns[i])
        genes.add(percent_diff_df.columns[i])
        if len(genes) >= top_k_genes:
            break

    genes = list(genes)
    genes = sorted(genes)

    # Remove last genes if exceeds the top_k_genes
    while len(genes) > top_k_genes:
        genes.pop()

    # Subset the dataframes to the top k genes
    expression_df = expression_df.loc[:, genes]
    percent_diff_df = percent_diff_df.loc[:, genes]


    # If whitelisted_types is not None, subset the cell types to only those in the whitelisted_types
    if whitelisted_types is not None:
        cell_types = cell_types[cell_types.isin(whitelisted_types)]
        expression_df = expression_df.loc[cell_types]
        percent_diff_df = percent_diff_df.loc[cell_types]

    # Default sort_types by the cell types names
    cell_types = cell_types.sort_values(ascending=False)

    # Sort the cell types and genes by the abs of percent difference
    if sort_types:
        cell_types = percent_diff_df.abs().sum(axis=1).sort_values(ascending=True).index
    if sort_genes:
        genes = percent_diff_df.abs().sum(axis=0).sort_values(ascending=False).index

    # Sort the cell types and genes by the total of expression
    if sort_types:
        cell_types = expression_df.sum(axis=1).sort_values(ascending=True).index
    if sort_genes:
        genes = expression_df.sum(axis=0).sort_values(ascending=False).index

    print("\t\t[3/3] Plotting the bubble plot...")

    # Bubble plot
    x, y = np.meshgrid(np.arange(len(genes)), np.arange(len(cell_types)))

    plt.subplots_adjust(bottom=0.1)
    fig, ax = plt.subplots(figsize=(0.27 * len(genes) + 9, 0.27 * len(cell_types) + 7.5))

    sizes = np.zeros((len(cell_types), len(genes)))
    colors = np.zeros((len(cell_types), len(genes)))

    for i, cell_type in enumerate(cell_types):
        for j, gene in enumerate(genes):
            sizes[i, j] = expression_df.loc[cell_type, gene] * 28
            colors[i, j] = percent_diff_df.loc[cell_type, gene]

    # Plot the bubbles
    sc = ax.scatter(x, y, s=sizes, c=colors, cmap="coolwarm", alpha=1)


    # Set the ticks and labels
    ax.set_xticks(np.arange(len(genes)))
    ax.xaxis.set_tick_params(labeltop=True, labelbottom=False, top=True, bottom=False)
    ax.set_xticklabels(genes, fontstyle="italic", rotation=60, ha="left")
    ax.set_yticks(np.arange(len(cell_types)))
    ax.set_yticklabels(cell_types)
    ax.autoscale(enable=True, axis='both', tight=True)


    ax.margins(0.7/len(genes), 0.7/len(cell_types))
    ax.set_title(f"[{adata.uns['info']['experiment_name']}] Differential Gene Expression between {roi_set_names[0]} and {roi_set_names[1]} ({total_cells_in_comparison} cells | {total_genes} genes)")
    ax.set_aspect("equal", "box")


    # Add the colorbar (horizontal, bottom of plot). Limit the size of the colorbar to the axes width, height of 5% the axes height
    divider = make_axes_locatable(ax)
    # barsize, min of (0.4 and 5% of axes height)
    barsize = min(0.4, ax.get_window_extent().height * 0.05)
    cax = divider.append_axes("bottom", size=barsize, pad=0.1)

    cbar = fig.colorbar(sc, cax=cax, orientation="horizontal")
    # label to the left
    min_diff = percent_diff_df.min().min()
    max_diff = percent_diff_df.max().max()
    cbar.set_label("Percent Difference")
    cbar.set_ticks([-2, 0, 2])
    cbar.set_ticklabels([f"Only in {roi_set_names[0]}", "No Difference", f"Only in {roi_set_names[1]}"])

    return fig