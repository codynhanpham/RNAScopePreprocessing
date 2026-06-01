# RNAScope Preprocessing

> [!WARNING]
> Work-in-progress preprocessing pipeline for RNAScope data
>
> For now, except for RNA-Scope specific steps (see [notebook](./notebooks/imapose2vizgen/scratch_imapose2vizgen.ipynb)), the main processing pipeline is adapted from https://github.com/codynhanpham/MERFISH_Analysis_Pipeline

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
In addition to this repository, you will also need:
1. [ZEISS ZEN Lite](https://www.zeiss.com/microscopy/us/products/software/zeiss-zen-lite.html) to export the raw CZI image file to a single-plane, downsampled TIFF file as input for Cellpose 

2. The source [Cellpose](https://github.com/mouseland/cellpose) package, or alternatively, the preconfigured and one-command install [cellpose-local](https://github.com/codynhanpham/cellpose-local) repository, set up elsewhere on your system to run Cellpose segmentation

> [!IMPORTANT]
> It is recommended that you have Cellpose installed separately from this repository's environment to avoid dependency conflicts. Cellpose relies on older versions of many packages and is strongly recommended to be run in its own isolated environment. Use the [cellpose-local](https://github.com/codynhanpham/cellpose-local) repository to skip the hassle and set up Cellpose with a single `uv sync` command.

See [notebook](./notebooks/imapose2vizgen/scratch_imapose2vizgen.ipynb) for converting Imaris puncta detection output and Cellpose segmentation into Vizgen-compatible formats for downstream analysis.


### Setup experiment config
Copy and edit [`pipeline_yaml/template_pipeline.yaml`](./pipeline_yaml/template_pipeline.yaml) to create a config file for a specific experiment.

### Run the pipeline
```bash
uv run ./run_processing_pipeline.py </path/to/config.yaml>

# For example
uv run ./run_processing_pipeline.py pipeline_yaml/your_pipeline.yaml
```