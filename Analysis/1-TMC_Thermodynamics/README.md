# 1 - TMC Thermodynamics

This analysis layer describes the isolated-molecule thermodynamics of the SCO
ATLAS and supports the accompanying manuscript. It reads SCOPE System objects
from the external ATLAS data store and writes all derived products within this
directory.

## Contents

- `Notebooks/`: numbered production notebooks, diagnostic `X-` notebooks, and
  `atlas_functions.py`.
- `Outputs/`: derived CSV and LaTeX tables used by the analysis and manuscript.
- `Figures/`: locally generated static PNG/PDF figures.
- `Figures/HTML/`: versioned interactive ligand maps and their CSV data tables.
- `Bibliography/`: Crossref/OpenAlex utilities and journal-name resources.

## Workflow

Run the production notebooks in numerical order:

```text
0-Prepare_Systems.ipynb
1-Crystal_Breakdown.ipynb
2-Bibliographic_Data.ipynb
3-Thermodynamics.ipynb
4-Ligand_Composition.ipynb
5-Structural_Correlations.ipynb
6-VNM_Analysis.ipynb
7-Finalize_ATLAS.ipynb
```

The configuration cells expect:

```python
data_path = "/path/to/ATLAS/Data/"
analysis_path = "/path/to/SCO-Atlas/Analysis/1-TMC_Thermodynamics/"
env_path = "/path/to/ATLAS/Data/scope_env_atlas.npy"
```

The notebooks target SCOPE 0.9.6. Files beginning with `X-` are diagnostic or
exploratory and are not required for the standard workflow.
