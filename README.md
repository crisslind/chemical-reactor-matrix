# High-Performance Chemical Reaction Enumerator

A multi-threaded chemical library generator built on RDKit. This tool performs virtual library matrix enumeration by processing combinations of core scaffolds and building blocks across diverse reaction SMARTS, offering real-time progress bars, internal product deduplication, and chemical property filtering.

## Features
* **High Performance:** Dynamic multiprocessing worker pooling scales to large building-block catalogs.
* **Smart Deduplication:** Identical product compounds generated via symmetrical sites or multiple pathways are merged into a single row, cleanly concatenating input IDs and reagent SMILES.
* **Real-time Monitoring:** Nested terminal progress tracking windows showing core scaffold progress alongside live building block parsing.
* **Property Filters:** Filter out unviable products on-the-fly using Molecular Weight (MW) and LogP thresholds (rounded to 1 decimal place).

## Prerequisites & Installation

Ensure you have a modern Python environment installed. Clone the repository and install the dependencies:

```bash
git clone [https://github.com/crisslind/chemical-reactor-matrix.git](https://github.com/crisslind/chemical-reactor-matrix.git)
cd chemical-reactor-matrix
pip install -r requirements.txt
