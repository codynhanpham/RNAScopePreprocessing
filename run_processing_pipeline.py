from functions import helpers, file_io, data_hashing, transcripts, data_filtering, cell_types, yaml_instructions
import yaml

import os
import sys
from typing import cast


def main():
    # CLI arguments for the pipeline yaml file
    if len(sys.argv) > 1:
        yaml_file = sys.argv[1]
    else:
        print("Please provide a yaml file for the pipeline configuration.")
        return

    print()

    # Load the pipeline configuration
    print(f"Loading pipeline configuration from: \x1b[35m{yaml_file}\x1b[0m")
    (PIPELINE_CONFIG, config_same_as_last_run) = yaml_instructions.load_pipeline_yaml(yaml_file)
    print("Pipeline configuration loaded successfully.\n")
    hash_file_path = data_hashing.getHashFilePath(PIPELINE_CONFIG['data_outputs']['output_dir'], PIPELINE_CONFIG['data_inputs']['experiment_name'])
    if not config_same_as_last_run:
        print("\x1b[33mThe pipeline configuration has changed since the last run.\x1b[0m\n")
        data_hashing.updateHashFileLog(hash_file_path, "PIPELINE_CONFIG", data_hashing.hashVar(PIPELINE_CONFIG), "vars")

    # Make the output directory
    os.makedirs(PIPELINE_CONFIG["data_outputs"]["output_dir"], exist_ok=True)

    print("Start the Processing Pipeline on Experiment: \x1b[32m", PIPELINE_CONFIG["data_inputs"]["experiment_name"], "\x1b[0m\n")


    # QC of Transcripts:
    # Load the transcripts data
    print("Loading transcripts data...")
    chunk_multiplier = 4
    transcripts_csv_len = len(PIPELINE_CONFIG["data_inputs"]["detected_transcripts"]["detected_transcripts_csv"])
    metadata: list[dict | None] = [None] * transcripts_csv_len
    for i, detected_transcripts_csv in enumerate(PIPELINE_CONFIG["data_inputs"]["detected_transcripts"]["detected_transcripts_csv"]):
        csv = detected_transcripts_csv["csv"]
        json = detected_transcripts_csv["json"]
        name = detected_transcripts_csv["name"]
        prioritize_json = detected_transcripts_csv["prioritize_json"]

        if (json is None or os.path.exists(json) is False):
            print(f"\t[{i+1}/{transcripts_csv_len}] No pre-processed JSON file for {name}. Processing the CSV file...")
            metadata[i] = transcripts.getDetectedTranscriptsStats(csv, name, chunk_multiplier * 10**4)
        else:
            if prioritize_json:
                print(f"\t[{i+1}/{transcripts_csv_len}] Loading pre-processed JSON file for {name}...")
                metadata[i] = transcripts.getDetectedTranscriptsStats(json, name, chunk_multiplier * 10**4)
            else:
                print(f"\t[{i+1}/{transcripts_csv_len}] Processing the CSV file for {name}...")
                metadata[i] = transcripts.getDetectedTranscriptsStats(csv, name, chunk_multiplier * 10**4)

    assert all(m is not None for m in metadata), "Not all metadata entries were populated by the loop."
    metadata_list = cast(list[dict], metadata)
    del metadata

    print("Transcripts data loaded successfully.\n")

    # Generate the QC report
    print("Generating the QC report...")
    if len(metadata_list) >= 2:
        delta_usable_transcripts = transcripts.deltaUsableTranscripts(metadata_list[-2], metadata_list[-1], PIPELINE_CONFIG["data_inputs"]["experiment_name"], dict(bubble_align_degree=-120, bubble_text_offset_deg=90))
        assert delta_usable_transcripts is not None, "deltaUsableTranscripts returned None."
        file_output = file_io.saveDeltaUsableTranscriptsReport(PIPELINE_CONFIG["data_outputs"]["output_dir"], delta_usable_transcripts)
        print(f"QC report generated successfully.\nPDF saved to:\n\t{file_output[0]}\n\t{file_output[1]}\n")
        del delta_usable_transcripts

    # Only keep the metadata set for PIPELINE_CONFIG["data_inputs"]["detected_transcripts"]["use_for_analysis"] == PIPELINE_CONFIG["data_inputs"]["detected_transcripts"]["detected_transcripts_csv"][i]["name"]
    # Find the index of the metadata set to use for analysis
    csv_names = [detected_transcripts_csv["name"] for detected_transcripts_csv in PIPELINE_CONFIG["data_inputs"]["detected_transcripts"]["detected_transcripts_csv"]]
    use_for_analysis_index = csv_names.index(PIPELINE_CONFIG["data_inputs"]["detected_transcripts"]["use_for_analysis"])
    metadata_single = metadata_list[use_for_analysis_index]

    if "transcripts_per_cell" in metadata_single:
        del metadata_single["transcripts_per_cell"]
    if "transcripts_histogram2d_data" in metadata_single:
        del metadata_single["transcripts_histogram2d_data"]

    # Extend PIPELINE_CONFIG["data_inputs"]["metadata"] with the metadata in a "transcripts_metadata" key
    if "transcripts_metadata" not in PIPELINE_CONFIG["data_inputs"]["metadata"]:
        PIPELINE_CONFIG["data_inputs"]["metadata"]["transcripts_metadata"] = metadata_single

    del(transcripts_csv_len, chunk_multiplier, use_for_analysis_index, csv_names, metadata_list, metadata_single)



    # Load the cell metadata and gene expression data:

    print("Loading cell metadata, gene expression data, and ROI definitions...")
    cell_anndata = file_io.loadDataTables(PIPELINE_CONFIG["data_inputs"]["cell_metadata_csv"], PIPELINE_CONFIG["data_inputs"]["cell_by_gene_csv"], PIPELINE_CONFIG["data_inputs"]["prioritized_hdf5"], PIPELINE_CONFIG["data_inputs"]["experiment_name"], PIPELINE_CONFIG["data_inputs"]["metadata"])
    print("Cell metadata and gene expression data loaded successfully.\n")

    cell_anndata_filtered = True if cell_anndata.uns.get("filtered", False) == True else False


    # Cell Filtering and Clustering:

    # Filter the cells according to pipeline configuration
    print("Filtering the cells...")
    # if PIPELINE_CONFIG["data_filters"] and PIPELINE_CONFIG["data_filters"]["filtering_procedure"] and len(PIPELINE_CONFIG["data_filters"]["filtering_procedure"]) > 0 and (PIPELINE_CONFIG["data_inputs"]["prioritized_hdf5"] is None or PIPELINE_CONFIG["data_filters"]["prefiltered_hdf5"] != True):
    if not cell_anndata_filtered and PIPELINE_CONFIG["data_filters"] and PIPELINE_CONFIG["data_filters"]["filtering_procedure"] and len(PIPELINE_CONFIG["data_filters"]["filtering_procedure"]) > 0:
        filtering_result = data_filtering.apply_filters(cell_anndata, PIPELINE_CONFIG["data_filters"]["filtering_procedure"], PIPELINE_CONFIG["data_filters"]["retain_raw_transcript_count"])
        stats = filtering_result["stats"]
        cell_anndata_filtered = filtering_result["adata"]
        print("Cells filtered successfully.\n")

        cell_filtering_spatial = data_filtering.plotSpatialCellFilteringStats(cell_anndata, cell_anndata_filtered, PIPELINE_CONFIG["data_inputs"]["spatial_rotation_degrees"], PIPELINE_CONFIG["data_inputs"]["spatial_flip_x"])

        # Save the filtering report
        print("Saving the cell filter report...")
        stats.append(dict(plt=cell_filtering_spatial))
        file_output = file_io.saveFilteringCellMetadataReport(
            PIPELINE_CONFIG["data_outputs"]["output_dir"],
            stats,
            PIPELINE_CONFIG["data_inputs"]["experiment_name"]
        )
        print(f"Filtering report PDF saved to {file_output}.\n")

        del(stats, cell_filtering_spatial)
    else:
        if cell_anndata_filtered:
            print("Input HDF5 file is marked as filtered. Skipping the filtering step.\n")
        else:
            print("No filtering procedure provided in the pipeline configuration. Skipping the filtering step.\n")
        cell_anndata_filtered = cell_anndata.copy()

    file_output = file_io.exportAnnData(PIPELINE_CONFIG["data_outputs"]["output_dir"], cell_anndata_filtered, PIPELINE_CONFIG["data_inputs"]["experiment_name"],"filtered")
    print(f"\tFiltered AnnData saved to {file_output}.\n")




if __name__ == '__main__':
    main()