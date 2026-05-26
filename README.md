# RNAScope Preprocessing

> Work-in-progress preprocessing pipeline for RNAScope data
>
> For now, except for RNA-Scope specific steps (see [notebook](./notebooks/scratch_imapose2vizgen.ipynb)), the main processing pipeline is adapted from https://github.com/codynhanpham/MERFISH_Analysis_Pipeline

## Installation
0. Install [uv](https://docs.astral.sh/uv/getting-started/installation/) for Python environment and dependency management.

1. Clone the repository
    ```bash
    git clone
    cd RNAScopePreprocessing
    ```
2. Install dependencies
    ```bash
    uv sync
    ```


## Usage
### Imaris spots + Cellpose segmentation to Vizgen-compatible formats
See [notebook](./notebooks/scratch_imapose2vizgen.ipynb) for converting Imaris puncta detection output and Cellpose segmentation into Vizgen-compatible formats for downstream analysis.

### Setup experiment config
Copy and edit [`pipeline_yaml/template_pipeline.yaml`](./pipeline_yaml/template_pipeline.yaml) to create a config file for a specific experiment.

### Run the pipeline
```bash
uv run ./run_processing_pipeline.py </path/to/config.yaml>

# For example
uv run ./run_processing_pipeline.py pipeline_yaml/your_pipeline.yaml
```