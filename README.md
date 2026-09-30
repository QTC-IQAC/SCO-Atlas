# SCO ATLAS

This repository contains the analysis notebooks, derived tabular data, figures,
and manuscript associated with the **SCO ATLAS**. The ATLAS is a collection of
experimentally grounded Fe(II)N6 transition-metal complexes represented and
analysed with SCOPE.

The large source, system, and computation files are stored separately from this
Git repository. In the current local layout, the two locations are:

```text
/Users/sergivela/Documents/SCOPE/Sub-Projects/1-ATLAS/  # This repository
/Users/sergivela/Documents/SCOPE/ATLAS/Data/            # External data store
```

The manuscript accompanying the present release analyses the `isolated` branch,
which contains gas-phase HS and LS calculations. The data store also contains a
`solid` branch that will be analysed separately.

## Repository Structure

```text
1-ATLAS/
├── Analysis/
│   ├── README.md       Convention for numbered analysis layers
│   └── 1-TMC_Thermodynamics/
│       ├── Bibliography/  Bibliographic lookup and retrieval utilities
│       ├── Notebooks/     Numbered notebooks and shared functions
│       ├── Outputs/       Derived CSV and LaTeX data products
│       ├── Figures/       Locally generated static figures
│       │   └── HTML/      Versioned interactive maps and CSV tables
│       └── README.md      Layer-specific workflow documentation
├── Manuscript/         Main text, ESI, bibliography, and LaTeX build files
├── Notes/              Scientific and development notes
├── AGENTS.md            Project context for coding assistants
└── README.md
```

Development material under `Old/`, notebook checkpoints, caches, and local
backups are not part of the released analysis.

Each scientific layer receives a chronological, descriptive directory under
`Analysis/`. Static PNG and PDF figures are generated locally by the notebooks
and are not versioned. `Analysis/1-TMC_Thermodynamics/Figures/HTML/` is the
exception: it contains the configurable ligand maps and their associated CSV
tables, which are included in the repository for direct exploration.

## Data Store

The external data store is expected to have the following structure:

```text
ATLAS/
└── Data/
    ├── 1-Sources/          Experimental structures and SCOPE Source objects
    ├── 2-Systems/          Serialized SCOPE System objects
    ├── 3-Computations/
    │   └── <system>/
    │       ├── isolated/   Gas-phase molecular calculations
    │       └── solid/      Periodic solid-state calculations, when available
    ├── scope_env_atlas.npy Serialized SCOPE Environment object
    ├── *.scope             SCOPE workflow definitions
    └── *.py                Workflow execution scripts
```

Each file under `Data/2-Systems/<system>/` links the corresponding sources,
branches, workflows, calculations, and derived results. These serialized System
objects are the principal entry points for navigating the ATLAS with SCOPE.

Instructions for exporting the shared data and computation branches as separate
archives are provided in `/Users/sergivela/Documents/SCOPE/ATLAS/README.md`.

## Requirements

- Python 3.12.
- SCOPE 0.9.6, including the core and SCO packages importable as `scope` and
  `scope_sco`.
- The complete shared ATLAS data and the computation branch required by the
  analysis being run.
- The supplied `Data/scope_env_atlas.npy` environment file.

Check the installed SCOPE version with:

```bash
python -c "import scope; print(scope.__version__)"
```

The notebooks in this release target SCOPE 0.9.6. The repository may be updated
as SCOPE evolves, so users of later releases should consult the current notebook
version.

## Notebook Configuration

Notebook `Analysis/1-TMC_Thermodynamics/Notebooks/0-Prepare_Systems.ipynb` is the analysis entry point.
Its configuration cell defines the external data store and this repository's
analysis directory:

```python
data_path     = "/path/to/ATLAS/Data/"
analysis_path = "/path/to/1-ATLAS/Analysis/1-TMC_Thermodynamics/"
env_path      = "/path/to/ATLAS/Data/scope_env_atlas.npy"
```

Notebook 0 checks the SCOPE version and required data directories before loading
the serialized environment. It then creates the CSV indexes used by subsequent
notebooks, including `finished_systems.csv`, `not_finished_systems.csv`,
`unique_systems.csv`, and `filtered_structure_systems.csv` in
`Analysis/1-TMC_Thermodynamics/Outputs/`.

System paths in these indexes are relative to the external `Data/` directory:

```text
2-Systems/ABITEM/ABITEM.npy
```

This keeps the indexes relocatable. After loading a System, the notebooks call
`system.set_paths_from_environment(env)` to remap its source and computation
paths to the configured local data store.

The main notebooks are:

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

Notebooks beginning with `X-` are diagnostic or exploratory and are not required
for the standard analysis sequence. Shared notebook helpers are collected in
`Analysis/1-TMC_Thermodynamics/Notebooks/atlas_functions.py`.

## Manuscript

The manuscript uses XeLaTeX, Biber, and `biblatex`. Build the ESI followed by the
main document from the repository root with:

```bash
make -C Manuscript
```

The current PDFs are written to:

```text
Manuscript/build/main.pdf
Manuscript/build/esi.pdf
```

See `Manuscript/README.md` for the LaTeX Workshop configuration and troubleshooting
instructions.

## Trust and Distribution

SCOPE System and Environment files are serialized Python objects and should only
be loaded from a trusted copy of the dataset. Before public distribution, confirm
that all crystallographic source material is shared in accordance with the
applicable CSD/CCDC licensing conditions.
