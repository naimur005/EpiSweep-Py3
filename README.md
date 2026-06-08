# EpiSweep-Py3
EpiSweep-Py3: Computationally-driven reengineering of therapeutic proteins to  reduce immunogenicity while maintaining function. Modernized Python 3 codebase with rotamer-based protein design and Docker support.

## Overview
EpiSweep-Py3 is a Python 3 modernization of the original EpiSweep framework that computationally designs therapeutic protein variants with reduced T cell epitopes. The toolkit integrates epitope prediction with structure-guided protein design to generate Pareto-optimal protein candidates that balance immunogenicity reduction with functional preservation.

## Key Features
```
Python 3.8+ - Modern codebase (upgraded from Python 2)
Rotamer-Based Design - Structure-guided optimization using OSPREY 3
Pareto Optimization - Multi-objective design for immunogenicity vs. function trade-offs
Docker Ready - Containerized environment for reproducibility
Sequence & Structure-Based Design - Flexible design methodologies
```
## How to run
Before you run this project make sure you have Docker and necessary tools inside it.

Step 1. Cloning the repository.
```
git clone git clone https://github.com/naimur005/EpiSweep-Py3.git
```
Step 2. Create a Docker environment.
```
docker compose build && docker compose up -d
```
To work inside EpiSweep-Py3
```
docker exec -it episweep /bin/bash
```
To restart Docker
```
docker restart episweep
```
Step 3. Install the requirements
```
pip install -r requirements.txt
```
## Dependencies
Python Packages
```
numpy, scipy, pandas, matplotlib, seaborn
biopython, PyYAML, scikit-learn
```
External Tools
```
IBM ILOG CPLEX (required for optimization)
OSPREY 3.x (for rotamer-based design)
Java JDK 17+ (for OSPREY)
R 3.x+ (for visualization)
```
Optional: 
```
CLEVER 1.0, Rosetta, PyMOL
```

## Design Methods
```
1. Sequence-Based Design: Uses MSA-derived amino acid frequencies. Fast and suitable for initial exploration.
2. Rotamer-Based Design : Incorporates rotamer conformations from OSPREY. More accurate structural predictions.
3. Structure-Based Design: CLEVER-derived energetic potentials combining rotamer and sequence information.
4. Combinatorial Library Design: Generates enriched variant libraries using "tube" potentials and degenerate codons.
```
## Usage
Basic Workflow
```
from episweep import Design

#Initialize design
design = Design()
design.load_spec('design_spec.csv')

#Filter mutations from MSA
design.initialize_filter()
design.msa_filter(low_thresh=0.30, high_thresh=0.95)
design.background_frequency_filter()
design.epitope_score_filter(min_epi_del=1)

#Generate potentials and optimize
design.generate_seq_potential()
```
## Run EpiSweep Optimizer
```
python bin/sweep.py design_spec.csv
```
## Analyze Results
```
#Plot Pareto frontier
Rscript bin/pareto_curve.R -i design_spec.csv -o pareto.pdf

#Predict and visualize epitopes
python bin/epitopes.py -i design_spec.csv -n 0 -o epitopes.csv
Rscript bin/epitope_map.R -i epitopes.csv -o epitope_map.pdf
```
## Citation
If you use EpiSweep-Py3 in your research, please cite:

> Choi, Y., Verma, D., Griswold, K. E., & Bailey-Kellogg, C. (2017).
> *EpiSweep: Computationally-driven reengineering of therapeutic proteins to reduce immunogenicity while maintaining function.*
> **Methods in Molecular Biology**, 1529, 375–398.
>
> Journal Paper: https://doi.org/10.1007/978-1-4939-6637-0_20