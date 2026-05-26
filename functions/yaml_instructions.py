from schema import Schema, SchemaError, And, Or, Optional
import yaml

import os
import typing

from functions import data_hashing

pipeline_schema = Schema({
    "data_inputs": {
        "experiment_name": str,

        Optional("metadata"): Or(None, dict),  # type: ignore[arg-type]

        "detected_transcripts": {
            "detected_transcripts_csv": [{
                "csv": And(str, lambda x: x.endswith(".csv") and os.path.exists(x), error="'detected_transcripts_csv' -> 'csv' path must end in '.csv' and lead to a valid existing file"),
                "json": Or(None, And(str, lambda x: x.endswith(".json") and os.path.exists(x), error="'detected_transcripts_csv' -> 'json' path must end in '.json' and lead to a valid existing file")),  # type: ignore[arg-type]
                "name": str,
                "prioritize_json": bool,
            }],
            "use_for_analysis": str,
        },
        "cell_metadata_csv": And(str, lambda x: x.endswith(".csv") and os.path.exists(x), error="'cell_metadata_csv' path must end in '.csv' and lead to a valid existing file"),
        "cell_by_gene_csv": And(str, lambda x: x.endswith(".csv") and os.path.exists(x), error="'cell_by_gene_csv' path must end in '.csv' and lead to a valid existing file"),
        Optional("prioritized_hdf5"): Or(None, And(str, lambda x: x.endswith(".hdf5") and os.path.exists(x), error="'prioritized_hdf5' path must end in '.hdf5' and lead to a valid existing file")),  # type: ignore[arg-type]

        "spatial_rotation_degrees": Or(float, int),
        "spatial_flip_x": bool,

        Optional("roi"): Or(None, [{  # type: ignore[arg-type]
            "roi_set_name": str,
            "roi_set_values": [str],
        }]),
    },

    "data_outputs": {
        "output_dir": And(str, lambda x: (os.makedirs(x, exist_ok=True) or True), error="'output_dir' must be a valid directory path"),
        Optional("normalize_mean_transcript_per_cell"): {
            "enable": bool,
            "target_value": Or(float, int),
        },
    },

    "data_filters": {
        Optional("filtering_procedure"): Or(None, [  # type: ignore[arg-type]
            {
                Optional("cell_volume"): {
                    "format": And(str, lambda x: x in ("percentile", "literal"), error="'cell_volume' -> 'format' must be one of ['percentile', 'literal']"),
                    "min": Or(float, int),
                    "max": Or(float, int),
                },
            },
            {
                Optional("blank_thresholding"): {
                    "set_threshold": And(str, lambda x: x in ("blank_max", "blank_min", "blank_mean", "blank_median"), error="'blank_thresholding' -> 'set_threshold' must be one of ['blank_max', 'blank_min', 'blank_mean', 'blank_median']"),
                    "filter_action": And(str, lambda x: x in ("remove", "subtract"), error="'blank_thresholding' -> 'filter_action' must be one of ['remove', 'subtract']"),
                },
            },
            {
                Optional("transcript_count"): {
                    "format": And(str, lambda x: x in ("percentile", "literal"), error="'transcript_count' -> 'format' must be one of ['percentile', 'literal']"),
                    "min": Or(float, int),
                    "max": Or(float, int),
                },
            },
            {
                Optional("normalize_transcript_by_volume"): bool,
            },
            {
                Optional("gene_per_cell"): {
                    "format": And(str, lambda x: x in ("percentile", "literal"), error="'gene_per_cell' -> 'format' must be one of ['percentile', 'literal']"),
                    "min": Or(float, int),
                    "max": Or(float, int),
                },
            },
        ]),

        Optional("retain_raw_transcript_count"): bool,
    },

}, ignore_extra_keys=True)


def validate_pipeline_yaml(loaded_yaml: dict) -> dict:
    try:
        # Basic type validation
        pipeline_schema.validate(loaded_yaml)

        # --- Defaults for optional keys ---

        if "prioritized_hdf5" not in loaded_yaml["data_inputs"] or loaded_yaml["data_inputs"]["prioritized_hdf5"] == "":
            loaded_yaml["data_inputs"]["prioritized_hdf5"] = None

        if "roi" not in loaded_yaml["data_inputs"]:
            loaded_yaml["data_inputs"]["roi"] = None

        if "normalize_mean_transcript_per_cell" not in loaded_yaml["data_outputs"]:
            loaded_yaml["data_outputs"]["normalize_mean_transcript_per_cell"] = {
                "enable": False,
                "target_value": -1,  # Does not matter, since enable is False
            }

        if "filtering_procedure" not in loaded_yaml["data_filters"]:
            loaded_yaml["data_filters"]["filtering_procedure"] = None

        if "retain_raw_transcript_count" not in loaded_yaml["data_filters"]:
            loaded_yaml["data_filters"]["retain_raw_transcript_count"] = False

        # --- Additional cross-field validation ---

        # detected_transcripts names must be unique
        transcript_names = [item["name"] for item in loaded_yaml["data_inputs"]["detected_transcripts"]["detected_transcripts_csv"]]
        if len(transcript_names) != len(set(transcript_names)):
            raise SchemaError(f"'name' in 'detected_transcripts_csv' must be unique. Duplicate names found: {transcript_names}")

        # use_for_analysis must refer to one of the declared names
        if loaded_yaml["data_inputs"]["detected_transcripts"]["use_for_analysis"] not in transcript_names:
            raise SchemaError(f"'use_for_analysis' must be one of the names in 'detected_transcripts_csv'. Valid names are: {transcript_names}")

        # roi_set_name values must be unique
        if loaded_yaml["data_inputs"]["roi"] is not None:
            roi_names = [item["roi_set_name"] for item in loaded_yaml["data_inputs"]["roi"]]
            if len(roi_names) != len(set(roi_names)):
                raise SchemaError(f"'roi_set_name' in 'roi' must be unique. Duplicate names found: {roi_names}")

        # Percentile bounds and min < max checks on filtering steps
        if loaded_yaml["data_filters"]["filtering_procedure"] is not None:
            for step in loaded_yaml["data_filters"]["filtering_procedure"]:
                # Each step is {filter_type: {format, min, max}} — unwrap the inner dict
                inner = next(iter(step.values())) if step else None
                if not isinstance(inner, dict):
                    continue
                if inner.get("format") == "percentile":
                    for bound in ("min", "max"):
                        if bound in inner and not (0 <= inner[bound] <= 100):
                            raise SchemaError(f"'{bound}' in 'filtering_procedure' must be between 0 and 100 when format is 'percentile'.")
                if "min" in inner and "max" in inner and inner["min"] > inner["max"]:
                    raise SchemaError(f"'min' in 'filtering_procedure' must be less than 'max'.")

        return loaded_yaml

    except SchemaError as e:
        print(f"\nError in the pipeline configuration file:\n{e}\n")
        if "'experiment_name'" in str(e):
            print("Please name this experiment and ensure that 'experiment_name' is a string.\n")
        exit(1)



def load_pipeline_yaml(yaml_path: str) -> typing.Tuple[dict, bool]:
    """
    Load and validate the pipeline configuration YAML file.

    Parameters
    ----------
    yaml_path : str
        Path to the pipeline configuration YAML file.

    Returns
    ----------
    PIPELINE_CONFIG: **dict**
        The validated pipeline configuration dictionary.
    same_as_last_run: **bool**
        Whether the pipeline configuration is the same as the last run, checked by the hash logged in the hash file.
    """

    try:
        with open(yaml_path, 'r') as yaml_file:
            loaded_yaml = yaml.safe_load(yaml_file)
            validated_yaml = validate_pipeline_yaml(loaded_yaml)
            config_hash = data_hashing.hashVar(validated_yaml)
            same_as_last_run = data_hashing.checkHashExists(data_hashing.getHashFilePath(validated_yaml['data_outputs']['output_dir'], validated_yaml['data_inputs']['experiment_name']), "PIPELINE_CONFIG", config_hash, "vars")
            return validated_yaml, same_as_last_run
    except FileNotFoundError:
        print(f"\nError: The pipeline configuration file was not found at: {yaml_path}\n")
        exit(1)
    except yaml.YAMLError:
        print(f"\nError: The pipeline configuration file is not a valid YAML file: {yaml_path}\n")
        exit(1)
    except Exception as e:
        print(f"\n{e}\nError while parsing the pipeline configuration file: {yaml_path}\nPlease make sure it is a valid YAML file.\n")
        exit(1)