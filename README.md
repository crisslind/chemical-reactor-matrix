# High-Performance Chemical Reaction Enumerator

A multi-threaded chemical library generator built on RDKit. This tool performs virtual library matrix enumeration by processing combinations of core scaffolds and building blocks across diverse reaction SMARTS, offering real-time progress bars, internal product deduplication, and chemical property filtering.

## Features
* **High Performance:** Dynamic multiprocessing worker pooling scales to large building-block catalogs.
* **Smart Deduplication:** Identical product compounds generated via symmetrical sites or multiple pathways are merged into a single row, cleanly concatenating input IDs and reagent SMILES.
* **Real-time Monitoring:** Nested terminal progress tracking windows showing core scaffold progress alongside live building block parsing.
* **Property Filters:** Filter out unviable products on-the-fly using Molecular Weight (MW) and LogP thresholds.

## Options and Parameters:
-s, --scaffolds: Path to your target scaffolds file (.csv, .tsv, .sdf) or a single raw SMILES string.

-b, --building_blocks: Path to your building block library metadata matrix (.csv, .tsv, .json).

-o, --output: Name/destination path of the output CSV file (Default: library_output.csv).

-w, --workers: Explicitly override default CPU core allocation.

--max_mw: Cutoff ceiling for maximum product Molecular Weight (Default: 2000.0).

--max_logp: Cutoff ceiling for maximum product LogP (Default: 20.0).

## 📋 Prerequisites & Installation

Ensure you have a modern Python environment installed. Then, clone the repository and install the dependencies:

```bash
git clone [https://github.com/crisslind/chemical-reactor-matrix.git](https://github.com/crisslind/chemical-reactor-matrix.git)
cd chemical-reactor-matrix
pip install -r requirements.txt


