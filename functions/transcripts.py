import sys, os, json, shutil, subprocess
sys.path.append('../')

import anndata
import numpy as np
import matplotlib
import matplotlib.figure
import scanpy as sc
import typing

from functions import helpers


def getDetectedTranscriptsStats(detected_transcripts_file: str, experiment_name: str, chunksize_KiB: int = 2048, n_threads=8) -> typing.Union[dict, None]:
    """
    Load the detected transcripts data matrix. Return an AnnData object.

    Required the binary `detected_transcripts_metadata` file to be in the same directory as the `main.py` script.
    
    Parameters
    ----------
    - `detected_transcripts_file` : **str**
        Path to the detected transcripts data matrix file. Can be a .csv or a pre-processed .json file.
    - `experiment_name` : **str**
        Name of the experiment.
    
    Returns
    -------
    - `accumulated_data` : **dict**
        A dictionary containing the accumulated data:
        - `cell_count`: **int**
            Number of unique cells with transcripts.
        - `transcripts_per_cell`: **array**
            Number of transcripts per cell.
        - `total_transcripts` : **int**
            Total number of transcripts.
        - `transcripts_within_cells` : **int**
            Number of transcripts within cells. (cell_id is not null or -1)
        - `layers_transcripts_count` : **dict**
            Number of transcripts per layer. (key: z index value, value: frequency of that z index)
        - `layers_transcripts_within_cells_count` : **dict**
            Number of transcripts per layer that belong to a cell. (key: z index value, value: frequency of that z index)
        - `genes_frequency` : **dict**
            Frequency of genes. (key: gene name, value: frequency of that gene)
        - `genes_frequency_per_layer` : **dict**
            Frequency of genes per layer. (key: z index value, value: dictionary of `gene frequency`)
        - `transcripts_mask_area`: **float**
    """
    

    # Check if the file exists
    if not os.path.exists(detected_transcripts_file):
        raise FileNotFoundError(f"The detected transcripts data matrix file '{detected_transcripts_file}' does not exist.")
    
    # If JSON file, load the data and return
    if detected_transcripts_file.endswith(".json"):
        try:
            # Load the json file if error occurs, return None
            json_file = open(detected_transcripts_file, "r")
            data = json.load(json_file)
            json_file.close()

            # Add the experiment name to the data
            data["experiment_name"] = experiment_name

            # Sort the z-layers and gene names dicts
            data["layers_transcripts_count"] = dict(sorted(data["layers_transcripts_count"].items()))
            data["layers_transcripts_within_cells_count"] = dict(sorted(data["layers_transcripts_within_cells_count"].items()))
            data["genes_frequency"] = dict(sorted(data["genes_frequency"].items()))
            data["genes_frequency_per_layer"] = {k: dict(sorted(v.items())) for k, v in data["genes_frequency_per_layer"].items()}

            return data
        except Exception as e:
            print(f"Error: {e}")
            return None


    # Extend the path to help find the binary file
    project_root = helpers.getProjectRoot()

    os.environ["PATH"] += os.pathsep + os.path.join(project_root)
    os.environ["PATH"] += os.pathsep + os.path.join(project_root, "bin")
    os.environ["PATH"] += os.pathsep + os.path.join(project_root, "functions")
    os.environ["PATH"] += os.pathsep + os.path.join(project_root, "functions", "bin")

    # Check the binary file exist before calling with subprocess
    try:
        binpath = shutil.which("detected_transcripts_metadata")
        if binpath is None:
            raise FileNotFoundError
        
    except: 
        print("No reachable executable for 'detected_transcripts_metadata', which is required for efficient transcripts metadata processing.")
        print("Contact nhanp@wustl.edu for instructions.")
        return None

    # Call the binary to process the detected transcripts data file
    try:
        print(f"Processing detected transcripts data...")
        # Use the resolved absolute executable path; with shell=False, each arg is passed safely.
        cmd = [binpath, "-i", detected_transcripts_file, "-q", "-o", "-", "-c", f"{chunksize_KiB}", "-p", f"{n_threads}"]
        json_string = subprocess.run(cmd, shell=False, capture_output=True)
    except Exception as e:
        print(f"Error: {e}")

        return None

    # Load the json string to a dictionary
    data = json.loads(json_string.stdout)

    # Sort the z-layers and gene names dicts
    data["layers_transcripts_count"] = dict(sorted(data["layers_transcripts_count"].items()))
    data["layers_transcripts_within_cells_count"] = dict(sorted(data["layers_transcripts_within_cells_count"].items()))
    data["genes_frequency"] = dict(sorted(data["genes_frequency"].items()))
    data["genes_frequency_per_layer"] = {k: dict(sorted(v.items())) for k, v in data["genes_frequency_per_layer"].items()}


    # Get the csv file name without extension, save data to {csv_file_name}_metadata.json 
    csv_file_name = os.path.basename(detected_transcripts_file).split(".")[0]
    metadata_file = os.path.join(os.path.dirname(detected_transcripts_file), f"{csv_file_name}_metadata.json")

    with open(metadata_file, "w") as json_file:
        json.dump(data, json_file, separators=(",", ":"))

    # Add the experiment name to the data
    data["experiment_name"] = experiment_name

    print(f"Metadata saved to '{metadata_file}'.")

    return data


# Load detected transcripts from pipeline
def pipelineLoadDetectedTranscriptsFiles(detected_transcripts: dict, experiment_name: str) -> typing.Union[dict, None]:
    """
    Load the detected transcripts data matrix for all the detected transcripts files listed in the pipeline configuration.

    Return an array of `accumulated_data` : **dict**

    Parameters
    ----------
    

    """



def deltaUsableTranscripts(transcripts_metadata_1: dict, transcripts_metadata_2: dict, experiment_name: str, plot_options: dict = dict()) -> typing.Union[dict, None]:
    """
    Transcripts are only usable if they are within a cell!

    Calculate and plot the difference in the number of usable transcripts between two different segmentation result. The total transcripts count must be the same as this only evaluate the segmentation quality.

    The report includes:
    - Info:
        - The number of cell count in set 1
        - The number of cell count in set 2
        - Percentage difference in cell count between the two metadata sets, as (2) compared to (1)
        - Number of usable transcripts in set 1
        - Number of usable transcripts in set 2
        - Percentage difference in usable transcripts between the two metadata sets, as (2) compared to (1)
        - The same info above for every layer

    - Plots:
        - Overlayed line chart of the number of usable transcripts per layer for both metadata sets. The difference is filled with `graph_options.line_diff_fill[0]` color when (2) > (1) and `graph_options.line_diff_fill[1]` color when (2) < (1).
        - Bar chart of the number of usable transcripts in each metadata set.
        - Stacked bubble chart to see intracellular transcripts of each metadata as part of the total transcripts.

    Parameters
    ----------
    - `transcripts_metadata_1` : **dict**
        The 1st detected transcripts metadata.
    - `transcripts_metadata_2` : **dict**
        The 2nd detected transcripts metadata.
    - `experiment_name` : **str**
        Name of the experiment.
    - `plot_options` : **dict**
        - `line_total_color`: **str**
            Color for the total line on the line chart. Default: `"black"`
        - `line_set_color`: **list**
            Colors for the line chart, in the order of (set_1, set_2). Default: `["chocolate", "teal"]`
        - `line_diff_fill`: **list**
            Colors for the fill between the two lines, in the order of (set_2 > set_1, set_2 < set_1). Default: `["green", "red"]`
        - `bar_color`: **list**
            Colors for the bar chart, in the order of (set_1, set_2). Default: `["chocolate", "teal"]`
        - `bubble_color`: **list**
            Colors for the Venn diagram, in the order of (set_1, set_2). Default: `["chocolate", "teal"]`
        - `bubble_radius`: **int**
            Radius of the Venn diagram bubbles. Default: `50`
        - `bubble_align_degree`: **int**
            Degree to align the stacked bubbles. Default: `-90` (0 to align the stacked bubbles at the right side of the total circle, 90 at the top, etc.)
        - `bubble_text_offset_deg`: **int**
            Degree to offset the text from the bubble. Default: `-90-60`

    Returns
    -------
    - `output` : **dict**
        A dictionary containing the difference and plot data for saving to a report.
        - `experiment_name` : **str**
            Name of the experiment.
        - `diff` : **dict**
            The difference in cell count and usable transcripts.
        - `axes` : **list**
            The matplotlib axes for the plots.
    """
    import matplotlib.pyplot as plt

    required_keys = ["total_transcripts", "transcripts_within_cells", "layers_transcripts_count", "layers_transcripts_within_cells_count"]

    # Check if the metadata is valid and that the two metadata have the same keys and total transcripts count
    if not all(key in transcripts_metadata_1 for key in required_keys) or not all(key in transcripts_metadata_2 for key in required_keys):
        print("Invalid metadata content.")
        return

    # if transcripts_metadata_1.keys() != transcripts_metadata_2.keys():
    #     print("Metadata keys mismatch.")
    #     return
    
    if transcripts_metadata_1["total_transcripts"] != transcripts_metadata_2["total_transcripts"]:
        print("Total transcripts count mismatch.")
        return None
    

    # Experiment names is string
    if not isinstance(experiment_name, str):
        print("Experiment name must be a string.")
        return None


    # Default metadata set names
    if "experiment_name" not in transcripts_metadata_1 or transcripts_metadata_1["experiment_name"] == "" or transcripts_metadata_1["experiment_name"] is None:
        name_1 = "Set 1"
    else:
        name_1 = transcripts_metadata_1["experiment_name"]
    
    if "experiment_name" not in transcripts_metadata_2 or transcripts_metadata_2["experiment_name"] == "" or transcripts_metadata_2["experiment_name"] is None:
        name_2 = "Set 2"
    else:
        name_2 = transcripts_metadata_2["experiment_name"]

    
    diff = dict()
    
    # Calculate the difference in cell count and usable transcripts
    diff["cell_count_diff"] = transcripts_metadata_2["cell_count"] - transcripts_metadata_1["cell_count"]
    diff["cell_count_diff_percent"] = diff["cell_count_diff"] / transcripts_metadata_1["cell_count"]

    diff["usable_transcripts_diff"] = transcripts_metadata_2["transcripts_within_cells"] - transcripts_metadata_1["transcripts_within_cells"]
    diff["usable_transcripts_diff_percent"] = diff["usable_transcripts_diff"] / transcripts_metadata_1["transcripts_within_cells"]

    # Calculate the difference in usable transcripts per layer
    diff["layers_usable_transcripts_diff"] = {
        z: transcripts_metadata_2["layers_transcripts_within_cells_count"][z] - transcripts_metadata_1["layers_transcripts_within_cells_count"][z]
        for z in transcripts_metadata_1["layers_transcripts_within_cells_count"]
    }
    diff["layers_usable_transcripts_diff_percent"] = {
        z: diff["layers_usable_transcripts_diff"][z] / transcripts_metadata_1["layers_transcripts_within_cells_count"][z]
        for z in transcripts_metadata_1["layers_transcripts_within_cells_count"]
    }


    # Plotting

    # Default plot options
    plot_options_defaults = dict(
        line_total_color="black",
        line_set_color=["chocolate", "teal"],
        line_diff_fill=["green", "red"],
        bar_color=["chocolate", "teal"],
        bubble_radius=50,
        bubble_color=["chocolate", "teal"],
        bubble_align_degree=-90, # 0 to align the stacked bubbles at the right side of the total circle, 90 at the top, etc.
        bubble_text_offset_deg=-90-60
    )

    # Update the plot options with the defaults if any of the keys are missing
    plot_options = {**plot_options_defaults, **plot_options}
    
    # Default plot styling
    plt = helpers.applyBasicPlotStyle(plt)
    
    fig, axs = plt.subplots(3, 1, figsize=(14, 24), gridspec_kw={"hspace": 0.3})

    # Line chart of the difference in usable transcripts per layer
    
    x = list(transcripts_metadata_1["layers_transcripts_count"].keys())
    y_total = list(transcripts_metadata_1["layers_transcripts_count"].values())
    y_usable_1 = list(transcripts_metadata_1["layers_transcripts_within_cells_count"].values())
    y_usable_2 = list(transcripts_metadata_2["layers_transcripts_within_cells_count"].values())

    axs[0].plot(x, y_total, label="All Transcripts", color=plot_options["line_total_color"], linestyle="--")
    axs[0].plot(x, y_usable_1, label=f"Intracellular: {name_1}", color=plot_options["line_set_color"][0])
    axs[0].plot(x, y_usable_2, label=f"Intracellular: {name_2}", color=plot_options["line_set_color"][1])

    # Dots for every [x, y] point
    axs[0].scatter(x, y_total, color=plot_options["line_total_color"], s=30)
    axs[0].scatter(x, y_usable_1, color=plot_options["line_set_color"][0], s=30)
    axs[0].scatter(x, y_usable_2, color=plot_options["line_set_color"][1], s=30)

    axs[0].fill_between(x, y_usable_1, y_usable_2, where=[y_usable_2[i] > y_usable_1[i] for i in range(len(x))], color=plot_options["line_diff_fill"][0], alpha=0.1)
    axs[0].fill_between(x, y_usable_1, y_usable_2, where=[y_usable_2[i] < y_usable_1[i] for i in range(len(x))], color=plot_options["line_diff_fill"][1], alpha=0.1)

    axs[0].set_title("Intracellular Transcripts per Layer", pad=10)
    axs[0].set_xlabel("Z Layer")
    axs[0].set_ylabel("Transcripts")

    # Y range [0, max(y_total) * 1.15] (for extra padding)
    axs[0].set_ylim(bottom=0)
    axs[0].set_ylim(top=max(max(y_total), max(y_usable_1), max(y_usable_2)) * 1.15)

    axs[0].grid(True, axis="y", alpha=0.4)
    axs[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)
    axs[0].yaxis.major.formatter._useMathText = True
    
    axs[0].legend()
    


    # Bar chart of total usable transcripts for both metadata sets

    x = [name_1, name_2]
    y = [transcripts_metadata_1["transcripts_within_cells"], transcripts_metadata_2["transcripts_within_cells"]]
    colors = [plot_options["bar_color"][0], plot_options["bar_color"][1]]

    axs[1].bar(x, y, color=colors)
    axs[1].set_title("Total Intracellular Transcripts", pad=10)
    axs[1].set_ylabel("Transcripts")
    axs[1].grid(True, axis="y", alpha=0.4)
    axs[1].ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)
    axs[1].yaxis.major.formatter._useMathText = True


    # Stacked bubble chart to show the intracellular transcripts as part of the total transcripts

    bubble_radius = plot_options["bubble_radius"]
    total_bubble_area = np.pi * bubble_radius ** 2

    # Ratio of intracellular transcripts to total transcripts for both metadata sets
    ratio_1 = transcripts_metadata_1["transcripts_within_cells"] / transcripts_metadata_1["total_transcripts"]
    ratio_2 = transcripts_metadata_2["transcripts_within_cells"] / transcripts_metadata_1["total_transcripts"]
    bubble_1_area = total_bubble_area * ratio_1
    bubble_2_area = total_bubble_area * ratio_2
    bubble_1_radius = np.sqrt(bubble_1_area / np.pi)
    bubble_2_radius = np.sqrt(bubble_2_area / np.pi)

    # Find the xy of the aligned location on the total circle
    aligned_x = np.cos(np.radians(plot_options["bubble_align_degree"])) * bubble_radius
    aligned_y = np.sin(np.radians(plot_options["bubble_align_degree"])) * bubble_radius

    # Find the xy of the bubble 1 and bubble 2 so that after drawing, the edge of the bubble is aligned with the total circle
    bubble_1_x = aligned_x - bubble_1_radius * np.cos(np.radians(plot_options["bubble_align_degree"]))
    bubble_1_y = aligned_y - bubble_1_radius * np.sin(np.radians(plot_options["bubble_align_degree"]))
    bubble_2_x = aligned_x - bubble_2_radius * np.cos(np.radians(plot_options["bubble_align_degree"]))
    bubble_2_y = aligned_y - bubble_2_radius * np.sin(np.radians(plot_options["bubble_align_degree"]))

    # Draw the total bubble
    total_circle = plt.Circle((0, 0), bubble_radius, color=plot_options["line_total_color"], fill=False, linewidth=2)
    axs[2].add_artist(total_circle)

    # Draw the intracellular bubbles
    bubble_colors = plot_options["bubble_color"]
    bubble_1 = plt.Circle((bubble_1_x, bubble_1_y), bubble_1_radius, color=bubble_colors[0], fill=False, linewidth=2)
    bubble_2 = plt.Circle((bubble_2_x, bubble_2_y), bubble_2_radius, color=bubble_colors[1], fill=False, linewidth=2)

    axs[2].add_artist(bubble_1)
    axs[2].add_artist(bubble_2)

    N = 800
    text_rot_deg = abs(plot_options["bubble_text_offset_deg"])
    text_offset = min(bubble_radius, bubble_1_radius, bubble_2_radius) * 0.05
    bubble_0_points = [
        [0 + (bubble_radius - text_offset) * np.cos(np.radians(text_rot_deg + i * 360 / N)) for i in range(N)],
        [0 + (bubble_radius - text_offset) * np.sin(np.radians(text_rot_deg + i * 360 / N)) for i in range(N)]
    ]
    bubble_1_points = [
        [bubble_1_x + (bubble_1_radius - text_offset) * np.cos(np.radians(text_rot_deg + i * 360 / N)) for i in range(N)],
        [bubble_1_y + (bubble_1_radius - text_offset) * np.sin(np.radians(text_rot_deg + i * 360 / N)) for i in range(N)]
    ]
    bubble_2_points = [
        [bubble_2_x + (bubble_2_radius - text_offset) * np.cos(np.radians(text_rot_deg + i * 360 / N)) for i in range(N)],
        [bubble_2_y + (bubble_2_radius - text_offset) * np.sin(np.radians(text_rot_deg + i * 360 / N)) for i in range(N)]
    ]

    # Flip the plotting direction so the text is upright
    bubble_0_points = [bubble_0_points[0][::-1], bubble_0_points[1][::-1]]
    bubble_1_points = [bubble_1_points[0][::-1], bubble_1_points[1][::-1]]
    bubble_2_points = [bubble_2_points[0][::-1], bubble_2_points[1][::-1]]

    helpers.CurvedText(
        x=bubble_0_points[0],
        y=bubble_0_points[1],
        text="All Transcripts",
        color=plot_options["line_total_color"],
        axes=axs[2],
        va="top",
    )
    helpers.CurvedText(
        x=bubble_1_points[0],
        y=bubble_1_points[1],
        text=name_1,
        color=bubble_colors[0],
        axes=axs[2],
        va="top",
    )
    helpers.CurvedText(
        x=bubble_2_points[0],
        y=bubble_2_points[1],
        text=name_2,
        color=bubble_colors[1],
        axes=axs[2],
        va="top",
    
    )
    
    axs[2].set_title("Intracellular Transcripts Proportion", pad=-10)

    axs[2].set_aspect("equal")
    axs[2].set_xlim(-bubble_radius * 1.15, bubble_radius * 1.15)
    axs[2].set_ylim(-bubble_radius * 1.15, bubble_radius * 1.15)
    axs[2].axis("off")


    fig.align_ylabels()
    
    fig.suptitle(f"[{experiment_name}] Detected Transcripts Comparison", fontsize=18, y=0.95)

    # Return the diff and plot data for saving to a report
    return dict(
        experiment_name = experiment_name,
        diff = diff,
        axes = axs
    )




def datasetsTranscriptsPerLayerNormalized(datasets_transcripts_metadata: typing.List[dict], colors: typing.List[str]) -> typing.List[dict]:
    # Line chart of the normalized number of usable transcripts per layer for all datasets

    import matplotlib.pyplot as plt

    plt = helpers.applyBasicPlotStyle(plt)

    data_metadata = datasets_transcripts_metadata.copy()

    # Each data will have the maximum number of total transcripts per layer normalize 0-1
    # The number of usable transcripts will be normalized to the same scale
    for data in data_metadata:
        max_total_transcripts = max(data["layers_transcripts_count"].values())
        data["layers_transcripts_count"] = {k: v / max_total_transcripts for k, v in data["layers_transcripts_count"].items()}
        data["layers_transcripts_within_cells_count"] = {k: v / max_total_transcripts for k, v in data["layers_transcripts_within_cells_count"].items()}


    fig, ax = plt.subplots(1, 1, figsize=(14, 8))

    if colors is None:
        colors = []

    # Extend/append the colors to match the length of the data_metadata
    colormap = matplotlib.cm.get_cmap("tab10")
    while len(colors) < len(data_metadata):
        colors.append(colormap(len(colors) % 10))
        

    for i, data in enumerate(data_metadata):
        x = list(data["layers_transcripts_count"].keys())
        y_total = list(data["layers_transcripts_count"].values())
        y_usable = list(data["layers_transcripts_within_cells_count"].values())

        ax.plot(x, y_total, label=f"All Transcripts: {data['experiment_name']}", color=colors[i], linestyle="--")
        ax.plot(x, y_usable, label=f"Intracellular: {data['experiment_name']}", color=colors[i])

        # Dots for every [x, y] point
        ax.scatter(x, y_total, color=colors[i], s=30)
        ax.scatter(x, y_usable, color=colors[i], s=30)




    ax.set_title("Normalized Intracellular Transcripts per Layer", pad=10)
    ax.set_xlabel("Z Layer")
    ax.set_ylabel("Transcripts")

    # Y range [0, 1.15] (for extra padding)
    ax.set_ylim(bottom=0)
    ax.set_ylim(top=1.15)

    ax.grid(True, axis="y", alpha=0.4)
    ax.legend()
    plt.show()


    # Return normalized metadata
    return data_metadata


def qc_exp_set_expression(exp_set_anndata: anndata.AnnData, sort_max_expression_first=False, fig_title: str= 'Normalized Gene Expression', output_pdf_path:typing.Optional[str] = None) -> typing.Optional[matplotlib.figure.Figure]:
    """
        Given an AnnData object with the un-normalized expression data, return a plot to visualize the expression levels of the genes across the experiments/runs.

        Parameters:
            exp_set_anndata (AnnData): AnnData object with un-normalized transcript counts in X, columns are genes and rows are cells, with an obs called `exp_name` that contains the experiment name.
            sort_max_expression_first (bool): Whether to sort the cells with maximum expression level first within each experiment. Default is `False`, which put the cell with highest expression in the center of each experiment group and fade to the min expression cells on the edges to make the experiment groups blend together better.
            fig_title (str): Title of the figure. Default is `Normalized Gene Expression`.
            output_pdf_path (str): Path to save the PDF output. Default is `None`, show without saving the plot.
    """

    import distinctipy
    import matplotlib.pyplot as plt
    import matplotlib.patches

    if exp_set_anndata.obs is None or 'exp_name' not in exp_set_anndata.obs.columns or exp_set_anndata.obs['exp_name'].isnull().all():
        raise ValueError("AnnData object must have an `exp_name` column in obs.")
    if output_pdf_path is not None and (not output_pdf_path.endswith(".pdf") or not isinstance(output_pdf_path, str)):
        raise ValueError("output_pdf_path, when specfied, must be a string and end with .pdf")

    # exp_names = data_summary['exp_name'].to_list()
    exp_names = exp_set_anndata.obs['exp_name'].unique().tolist()

    count_data_norm = sc.pp.normalize_total(exp_set_anndata, target_sum=1e6, copy=True)

    # Create a blank dataframe and update it with the per experiment data to make sure all cells in same exp are together
    ordered_data_X = np.zeros((count_data_norm.X.shape))
    exp_cell_counts = np.zeros((len(exp_names), 1))

    current_start_index = 0
    for i, exp_name in enumerate(exp_names):
        this_count_data = count_data_norm[count_data_norm.obs['exp_name'] == exp_name]
        print(this_count_data.obs['exp_name'].unique(), end=': ')

        exp_X = this_count_data.X
        # sort X by std of rows
        exp_X = exp_X[np.argsort(np.std(exp_X, axis=1))]

        # OPTIONAL, but is default: Make the lowest std row the center of the data and the min rows the edges
        if not sort_max_expression_first:
            exp_X_odd = exp_X[::2, :]
            exp_X_even = exp_X[1::2, :]
            exp_X_odd = exp_X_odd[::-1, :]
            exp_X = np.concatenate((exp_X_odd, exp_X_even), axis=0)
        else:
            exp_X = exp_X[::-1, :]

        # Update the aggregated data
        ordered_data_X[current_start_index:(current_start_index + exp_X.shape[0]), :] = exp_X
        current_start_index += exp_X.shape[0]
        exp_cell_counts[i] = exp_X.shape[0]
        print(f"{exp_X.shape[0]} samples", end='\n')

    ordered_data_X = np.log(ordered_data_X + 1.0)
    gene_list = np.array(list(count_data_norm.var_names))
    # sort the gene list by alphabetical order, then sort the data using the sorted gene list
    gene_list_sort_idx = np.argsort(gene_list)
    ordered_data_X = ordered_data_X[:, gene_list_sort_idx]



    fig, ax = plt.subplots(figsize=(24, 30), nrows=1, ncols=2, width_ratios=[0.97, 0.03], layout="constrained", sharey=True)
    ax[0].imshow(ordered_data_X, cmap='viridis', aspect='auto')

    # Get the first letters of the gene names
    gene_first_letters = np.array([gene_list[i][0] for i in gene_list_sort_idx])
    # Get the index of the first new letter in the gene_first_letters array
    gene_first_letters, gene_first_letters_idx = np.unique(gene_first_letters, return_index=True)
    ax[0].set_xticks(gene_first_letters_idx+0.5)
    xticks_labels = np.array([
        f"\n{1}" if i == 0 else 
        f"\n{ordered_data_X.shape[1]} " if i == (ordered_data_X.shape[1] - 1) else 
        "" for i in range(ordered_data_X.shape[1])])
    xticks_labels[gene_first_letters_idx] = np.char.add(gene_first_letters, xticks_labels[gene_first_letters_idx])
    xticks_labels = xticks_labels[xticks_labels != ""]
    ax[0].set_xticklabels(xticks_labels.tolist(), fontsize=14)
    ax[0].set_xlabel('Genes', fontsize=20)

    ax[0].set_yticks([])
    ax[0].set_yticklabels([])

    # Draw the bar outside of plot to indicate the location/span of each experiment in the plot
    current_start_index = 0
    colors = distinctipy.get_colors(len(exp_cell_counts))
    for i, exp_name in enumerate(exp_names):
        exp_y_pos = [current_start_index, (current_start_index + exp_cell_counts[i])[0]]
        ax[1].add_patch(
            matplotlib.patches.Rectangle(
                (0, exp_y_pos[0]),
                1, # width
                exp_y_pos[1] - exp_y_pos[0],
                facecolor=colors[i],
                edgecolor=None,
                linewidth=0,
            )
        )
        ax[1].text(
            2, # this is the x-position, give it some gap from the Rectangle width
            np.mean(exp_y_pos),
            exp_name,
            rotation=-90,
            fontsize=18,
            ha='center',
            va='center',
        )
        current_start_index += exp_cell_counts[i][0]
    ax[1].set_ylim([0, ordered_data_X.shape[0]])
    ax[1].set_xlim([0, 5])
    ax[1].axis('off')

    fig.suptitle(fig_title, fontsize=24, y=1.02)
    # plt.show()


    if output_pdf_path is not None and isinstance(output_pdf_path, str):
        try:
            fig.savefig(output_pdf_path, bbox_inches='tight', pad_inches=0.2)
        except Exception as e:
            print("Error saving figure as PDF")
            raise e
        
    plt.close()
    return fig