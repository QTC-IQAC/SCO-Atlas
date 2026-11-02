# SCO ATLAS Analysis 

This repository contains reproducible analyses built on the SCO ATLAS, 
a chemically aware data infrastructure represented and navigated with SCOPE. 
The underlying data objects and calculations are archived on Zenodo, 
while this repository contains the notebooks, derived outputs and 
interactive resources used to interrogate these data.

The first analysis layer accompanies *Mapping How Ligands Shape Spin-State
Thermodynamics in Fe(II)N₆ Complexes: Insights from the SCO ATLAS*.

The ATLAS data itself, including the SCOPE System objects and computation files,
is archived on Zenodo at [https://doi.org/10.5281/zenodo.23083798](https://doi.org/10.5281/zenodo.23083798).

## Repository Structure

```text
SCO-Atlas/                 Repository root (local folder: ATLAS/Analyses)
├── 1-TMC_Thermodynamics/
│   ├── Bibliography/      Bibliographic lookup and retrieval utilities
│   ├── Notebooks/         Production and diagnostic notebooks, including helpers
│   ├── Outputs/           Derived CSV tables
│   ├── Figures/
│   │   └── HTML/          Versioned interactive maps and associated CSV tables
│   └── README.md          Layer-specific workflow documentation
├── docs/                  Analysis-layer and publication workflow guidance
├── AGENTS.md              Project context for assistants
└── README.md
```

Each scientific layer (e.g. 1-TMC_Thermodynamics) has its own numbered
directory directly at the repository root. Static PNG (and PDF) figures are generated locally by the notebooks
and are not versioned. `1-TMC_Thermodynamics/Figures/HTML/` is the
exception: it contains the configurable ligand maps and their associated CSV
tables, which are included in the repository for direct exploration.

See [Analysis layers](docs/analysis_layers.md) and
[Publication workflow](docs/publication_workflow.md) for adding curated analyses.

## ATLAS Data Organization

The external data is expected to have the following structure:

```text
ATLAS/
└── Data/
    ├── 1-Sources/          Experimental structures and SCOPE Source objects
    ├── 2-Systems/          Serialized SCOPE System objects
    ├── 3-Computations/
    │   └── <system>/
    │       ├── isolated/   Gas-phase molecular calculations
```

Each file under `Data/2-Systems/<system>/` links the corresponding sources,
branches, workflows, calculations, and derived results. These serialized System
objects are the principal entry points for navigating the ATLAS with SCOPE.

Instructions for extracting, verifying, and configuring the archived data are
provided in the README distributed with the Zenodo record.

## Requirements

- Python 3.12.
- SCOPE 0.9.6, including the core and SCO packages importable as `scope` and `scope_sco`.
- The complete [SCO ATLAS data archive](https://doi.org/10.5281/zenodo.23083798).

If you're unsure, you can check the installed SCOPE version with:

```bash
python -c "import scope; print(scope.__version__)"
```

The notebooks in this release target SCOPE 0.9.6. The repository may be updated
as SCOPE evolves, so users of later releases should consult the current notebook
version.

## Notebook Configuration

All notebooks have configuration cells that defines the paths

The main notebooks are:

```text
0-Prepare_Systems.ipynb
1-Crystal_Breakdown.ipynb
2-Bibliographic_Data.ipynb
3-Thermodynamics.ipynb
4-Ligand_Composition.ipynb
5-Structural_Correlations.ipynb
6-VNM_Analysis.ipynb
```

Notebooks beginning with `D-` are focused diagnostic examples and are not
required for the standard analysis sequence. Layer-specific notebook helpers
are collected in `1-TMC_Thermodynamics/Notebooks/local_functions.py`.

## Trust and Distribution

SCOPE System and Environment files are serialized Python objects and should only
be loaded from a trusted copy of the dataset. Before public distribution, confirm
that all crystallographic source material is shared in accordance with the
applicable CSD/CCDC licensing conditions.
