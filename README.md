# Patient Deduplicator

A Python package for identifying potential duplicate patients across different clinical studies based on both clinical characteristics and genetic mutation profiles.

## Features

- Load and preprocess patient clinical data and mutation data
- Create patient profiles based on demographic and clinical information
- Create mutation profiles from genetic data
- Compute similarity between patients using various methods:
  - Clinical similarity based on demographic matching
  - Mutation similarity based on common genetic alterations
  - Combined similarity using both clinical and mutation data
- Find potential duplicate patients across different studies
- Export results to CSV for further analysis
- Generate visualizations of potential duplicate pairs

## Installation

### Option 1: Install from source

```bash
# Clone the repository
git clone https://github.com/Almogkm/patient_deduplicator.git
cd patient_deduplicator

# Install in development mode
pip install -e .
```

### Option 2: Create a conda environment (recommended)

```bash
# Create conda environment
conda create -n patient_dedup python=3.9 pandas numpy matplotlib pyarrow
conda activate patient_dedup

# Install the package
cd patient_deduplicator
pip install -e .
```

## Usage

### Python API

```python
from patient_deduplicator import PatientDeduplicator

# Initialize the deduplicator
deduplicator = PatientDeduplicator(
    mutations_file="path/to/mutations.parquet",
    patients_file="path/to/patients.parquet",
    output_dir="results"
)

# Configure parameters
deduplicator.set_configuration(
    similarity_threshold=0.7,
    clinical_threshold=0.8,
    min_matches=2,
    method="combined"
)

# Run the complete pipeline
results = deduplicator.run_pipeline(
    sample_size=None,  # Use entire dataset
    max_comparisons=None,  # No limit on comparisons
    visualize=True
)

# Access results
potential_duplicates = results['potential_duplicates']
csv_path = results['csv_output']
visualization_paths = results['visualizations']

print(f"Found {len(potential_duplicates)} potential duplicate pairs")
print(f"Results saved to {csv_path}")
```

### Command Line

The package includes a command-line interface for running the deduplication pipeline:

```bash
python -m patient_deduplicator.run \
    --mutations path/to/mutations.parquet \
    --patients path/to/patients.parquet \
    --output results \
    --method combined \
    --min-matches 2 \
    --similarity 0.7 \
    --clinical 0.8
```

Run with `--help` to see all available options:

```bash
python -m patient_deduplicator.run --help
```

### SLURM Cluster

For running on a SLURM-based computational cluster, use the provided `run.sh` script:

```bash
sbatch run.sh \
    --mutations /path/to/mutations.parquet \
    --patients /path/to/patients.parquet \
    --output results \
    --method combined
```

You can customize SLURM parameters (memory, CPUs, etc.) by editing the `#SBATCH` directives at the top of the `run.sh` file.

## Pipeline Steps

1. **Load Data**: Load patient clinical data and mutation data from parquet files
2. **Preprocess Data**: Handle missing values, standardize formats, and prepare data for analysis
3. **Create Patient Profiles**: Build patient profiles based on clinical data
4. **Create Mutation Profiles**: Generate mutation profiles from genetic data
5. **Find Potential Duplicates**: Compare patients across studies to identify potential duplicates
6. **Export Results**: Save potential duplicate pairs to CSV
7. **Visualize Results**: Generate visualizations of potential duplicate pairs

## Requirements

- Python 3.8 or higher
- pandas
- numpy
- matplotlib
- pyarrow (for parquet support)

## License

[MIT License](LICENSE)