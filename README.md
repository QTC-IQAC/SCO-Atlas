# SCO ATLAS Data

This directory contains the data layer of **SCO ATLAS**, an experimentally grounded collection of spin-crossover (SCO) systems represented and analysed with SCOPE. The analysis accompanying the present manuscript focuses on the isolated, gas-phase high-spin (HS) and low-spin (LS) states.

## Requirements

- A Python environment with **SCOPE 0.9.6** installed.
- The SCOPE core and SCO packages must both be importable as `scope` and `scope_sco`.
- The complete `Data` directory, including `scope_env_atlas.npy`.

The installed version can be checked with:

```bash
python -c "import scope; print(scope.__version__)"
```

The analysis notebooks accompanying this dataset are intended to be run with SCOPE 0.9.6, and Notebook 0 validates this version before loading the data. The notebooks will also be published in the project repository and updated, when necessary, as SCOPE evolves; users of later SCOPE releases should therefore consult the current repository version of the analysis.

## Contents

```text
Data/
├── 1-Sources/          Experimental structures and SCOPE source objects
├── 2-Systems/          Serialized SCOPE System objects
├── 3-Computations/     Quantum-chemistry calculation files
├── scope_env_atlas.npy Serialized SCOPE Environment object
├── *.scope             SCOPE workflow definitions
└── *.py                Scripts used to run the workflows
```

Each directory below `2-Systems/` contains the binary `System` object that connects the corresponding sources, HS and LS states, computations, workflows, and derived results. These System files are the principal entry points for navigating the ATLAS.

The development-only `old/` directory is not part of the distributed dataset.

## Local Configuration

The serialized environment records the paths used when the ATLAS was created. Notebook `0-Create_System_Files.ipynb` therefore remaps it to the local data copy. Set the following two paths in its configuration cell:

```python
from pathlib import Path

data_path = Path("/path/to/Database_SCO/Data").expanduser().resolve()
env_path = Path("/path/to/Database_SCO/Data/scope_env_atlas.npy").expanduser().resolve()
```

The notebook checks the SCOPE version and the presence of `1-Sources`, `2-Systems`, `3-Computations`, and the environment file before continuing.

## Analysis Entry Points

The notebooks under `Analysis/1-Analyse_Database/` operate on selected collections of System objects rather than duplicating those binaries. Notebook 0 prepares CSV indexes such as `finished_systems.csv`, `not_finished_systems.csv`, `unique_systems.csv`, and `filtered_structure_systems.csv` in the project-level `Outputs/` directory.

The `path` fields in newly generated indexes are stored relative to the project root, for example:

```text
Data/2-Systems/ABITEM/ABITEM.npy
```

This keeps the data relocatable while preserving the real on-disk names of systems whose CSD refcodes contain suffixes such as `01` or `02`. Subsequent notebooks load the indexed System objects and call `system.set_paths_from_environment(env)` so their internal source and computation paths refer to the configured local data copy.

SCOPE binary files should only be loaded from a trusted copy of the dataset.
