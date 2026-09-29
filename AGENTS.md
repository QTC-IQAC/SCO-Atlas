# SCO ATLAS Notes

## Project Identity

- Project name: `SCO ATLAS`.
- Repository root: `/Users/sergivela/Documents/SCOPE/Sub-Projects/1-ATLAS`.
- External data root: `/Users/sergivela/Documents/SCOPE/ATLAS/Data`.
- Purpose: construct, document, and analyse a chemically aware dataset of
  Fe(II)N6 spin-crossover complexes.
- The current manuscript analyses the gas-phase `isolated` branch. Periodic
  `solid` computations belong to a later analysis layer.

## Relationship to SCOPE

- SCOPE source repository: `/Users/sergivela/Documents/SCOPE/Program/Scope`.
- The analysis notebooks target SCOPE 0.9.6 and require both `scope` and
  `scope_sco`.
- Most production calculations were orchestrated with earlier frozen SCOPE
  releases; do not rewrite or regenerate them unless explicitly requested.
- Prefer analysis changes in this repository. Changes to SCOPE should be made
  only when the task explicitly concerns a reusable SCOPE bug or feature.

## Repository Structure

- Analysis notebooks: `Analysis/Notebooks/`.
- Shared notebook helpers: `Analysis/Notebooks/atlas_functions.py`.
- Derived tables and CSV files: `Analysis/Outputs/`.
- Bibliographic tools and journal data: `Analysis/Bibliography/`.
- Figures and interactive maps: `Figures/`.
- Scientific and development notes: `Notes/`.
- Main manuscript and ESI: `Manuscript/main.tex` and `Manuscript/esi.tex`.
- LaTeX PDFs: `Manuscript/build/main.pdf` and `Manuscript/build/esi.pdf`.

The standard analysis sequence is:

```text
0-Prepare_Systems.ipynb
1-Crystal_Breakdown.ipynb
2-Bibliographic_Data.ipynb
3-Thermodynamics.ipynb
4-Ligand_Composition.ipynb
5-Structural_Correlations.ipynb
```

Files beginning with `X-` are diagnostic or exploratory notebooks. Some of them
are computationally expensive; inspect them before executing cells. In
particular, rebuilding the full HO/QRRHO comparison in `X-Test_HO.ipynb` takes
approximately 42 minutes.

## External Data Structure

The large ATLAS data are deliberately kept outside the Git repository:

```text
/Users/sergivela/Documents/SCOPE/ATLAS/Data/
├── 1-Sources/
├── 2-Systems/
├── 3-Computations/
│   └── <system>/
│       ├── isolated/
│       └── solid/
├── scope_env_atlas.npy
├── *.scope
└── *.py
```

System indexes in `Analysis/Outputs/` store paths relative to the external
`Data/` directory, for example `2-Systems/ABITEM/ABITEM.npy`. Notebooks load the
serialized environment and call `system.set_paths_from_environment(env)` to
restore local source and computation paths.

The data-store release and branch-compression instructions live in
`/Users/sergivela/Documents/SCOPE/ATLAS/README.md`.

## Working Conventions

- Preserve the numbered, readable notebook workflow. Put reusable or complex
  helpers in `atlas_functions.py`, but keep short analysis loops visible in the
  notebooks.
- Avoid pandas and `pathlib.Path` in the analysis notebooks unless a task
  specifically requires them; prefer `csv`, ordinary dictionaries, and string
  paths.
- Do not execute long-running notebooks or regenerate expensive data unless the
  user explicitly asks for it.
- Do not overwrite raw `Sources`, `Systems`, or `Computations` while performing
  analysis.
- Preserve six-letter refcodes in analysis tables while retaining physical file
  paths containing derived suffixes such as `01` or `02`.
- Treat serialized SCOPE objects as trusted-data-only files.

## Manuscript

- Build the ESI before the main manuscript with `make -C Manuscript`.
- Figures are read through `Manuscript/figures`, a link to the repository-level
  `Figures/` directory.
- Do not edit prose unless requested. When discussing manuscript changes, propose
  LaTeX text in chat unless the user explicitly asks for direct modification.
- Keep gas-phase molecular results distinct from experimental crystal or
  solution observables. The computed molecular crossing temperature is an
  analogue of experimental `T1/2`, not a direct prediction of it.
