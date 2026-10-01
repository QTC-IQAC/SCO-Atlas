# 1 - TMC Thermodynamics

This analysis layer describes the gas-phase thermodynamics of the SCO
ATLAS and supports the accompanying manuscript. It reads SCOPE `System` objects
from the external ATLAS data, and executes the analyses discussed in the manuscript

## Requirements

The notebooks require:

- a Python environment containing SCOPE 0.9.6, including the `scope` and `scope_sco` packages;
- the complete SCO ATLAS `Data` directory, including `1-Sources/`, `2-Systems/`, and `3-Computations/`; and
- a local `scope_env_atlas.npy` environment generated as described in the ATLAS `README.md`.

The released notebooks target SCOPE 0.9.6. The repository version may be
updated when necessary to remain compatible with later SCOPE releases.

Each notebook begins with a configuration cell that collects its external data
paths and the additional files it reads or writes. Adapt the following root
paths to the local installation before running the workflow:

```python
data_path     = "/path/to/ATLAS/Data/"
env_path      = "/path/to/ATLAS/Data/scope_env_atlas.npy"
analysis_path = "/path/to/SCO-Atlas/Analysis/1-TMC_Thermodynamics/"
figures_path  = "/path/to/SCO-Atlas/Analysis/1-TMC_Thermodynamics/Figures/"
```

## Contents

- `Notebooks/`: production and diagnostic notebooks.
- `Notebooks/local_functions.py`: helper functions developed specifically for
  this analysis layer. They support the notebooks but are not part of SCOPE
- `Outputs/`: output CSV and LaTeX tables used by the analysis and manuscript.
- `Figures/`: folder to hosts generated PNG/PDF figures locally by the user.
- `Figures/HTML/`: hosts interactive ligand maps and their CSV data tables.
- `Bibliography/`: Crossref/OpenAlex utilities and journal-name resources.

## Notebooks

The documented workflow contains two notebook categories. Numbered production
notebooks form the reproducible analysis sequence, whereas notebooks prefixed
with `D-` are focused diagnostic examples used to inspect or validate individual
procedures and are not required for production.

Run the production notebooks in numerical order:

```text
0-Prepare_Systems.ipynb
1-Crystal_Breakdown.ipynb
2-Bibliographic_Data.ipynb
3-Thermodynamics.ipynb
4-Ligand_Composition.ipynb
5-Structural_Correlations.ipynb
6-VNM_Analysis.ipynb
```

The CSV and LaTeX files in `Outputs/` are the primary records of this analysis
layer. The workflow does not add a parallel set of analysis-specific
`atlas_*` results to the SCOPE `System` objects; underlying thermodynamic values
remain accessible from their original SCOPE `Collection` objects.

The current correspondence with the main manuscript is:

| Notebook | Main role | Main-text section(s) |
| --- | --- | --- |
| `0-Prepare_Systems.ipynb` | System preparation, completion checks, structural filtering, and representative selection | 2.1, Dataset Construction; 3.1, Dataset Analysis |
| `1-Crystal_Breakdown.ipynb` | Crystal and calculation-status composition | 3.1, Dataset Analysis |
| `2-Bibliographic_Data.ipynb` | Bibliographic provenance and publication coverage | 3.1, Dataset Analysis |
| `3-Thermodynamics.ipynb` | Thermodynamic terms, transition temperatures, and frequency-resolved contributions | 2.2-2.3, Thermodynamic Framework and Computational Details; 3.2, Thermodynamics of the Spin Transition |
| `4-Ligand_Composition.ipynb` | Ligand descriptors, interactive family maps, and family-resolved statistics | 2.5, Ligand-composition Split and Mapping; 3.4, Ligand Families Guide Thermodynamics |
| `5-Structural_Correlations.ipynb` | Structure-property correlations and ligand-family structural ranges | 3.5, Structure-Property Correlations |
| `6-VNM_Analysis.ipynb` | Normal-mode mapping and group-resolved vibrational analysis | 2.4, Vibrational Normal-Mode Mapping and Analysis; 3.3, Vibrational Origin of the Denticity Dependence |

### Diagnostic Notebooks

The diagnostic notebooks document the ligand-core assignment and the
normal-mode mapping for individual systems:

```text
D-Core_Ring_Assignment.ipynb
D-VNM_Mapping.ipynb
```
