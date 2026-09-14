# SCO ATLAS Notes

## Project Identity

- Workspace name: `DATABASE_SCO`
- Project name: `SCO ATLAS`
- Purpose: analyze results of quantum chemistry computations for a dataset of spin crossover (SCO) molecules and materials

## Relationship To SCOPE

- Related software repo: `/Users/sergivela/Documents/SCOPE/Program/Scope`
- `SCOPE` is under review and should stay as frozen as possible
- Prefer doing analysis work in this repository/workspace instead of modifying `SCOPE`
- If a future task appears to require changes in `SCOPE`, pause and reassess before touching that codebase

## Main Analysis Area

- Analysis notebooks live in:
  - `/Users/sergivela/Documents/SCOPE/Database_SCO/Analysis/1-Analyse_Database`
- The intended structure is numerated notebooks that are being polished progressively
- Consolidated analysis outputs now live in:
  - `/Users/sergivela/Documents/SCOPE/Database_SCO/Outputs`
- `Analysis/1-Analyse_Database/Outputs` is a compatibility link to the workspace-level output directory for existing notebooks

## Main Data Sources

- Raw/project data lives in:
  - `/Users/sergivela/Documents/SCOPE/Database_SCO/Data`
- This includes:
  - computations
  - `SCOPE` system objects derived from parsing those computations

## Current Analysis Entry Points

- Consolidated system collections for notebook analysis:
  - `/Users/sergivela/Documents/SCOPE/Database_SCO/Analysis/1-Analyse_Database/finished_systems.npy`
  - `/Users/sergivela/Documents/SCOPE/Database_SCO/Analysis/1-Analyse_Database/not_finished_systems.npy`
  - `/Users/sergivela/Documents/SCOPE/Database_SCO/Analysis/1-Analyse_Database/unique_systems.npy`
    - currently being created

## Manuscript

- LaTeX sources live in:
  - `/Users/sergivela/Documents/SCOPE/Database_SCO/Manuscript`
- Main documents:
  - `Manuscript/main.tex`
  - `Manuscript/esi.tex`
  - `Manuscript/references.bib`
- The manuscript uses XeLaTeX, Biber, and `biblatex` with the `chem-acs` style
- Build both documents with `make -C Manuscript`; PDFs are written to `Manuscript/build/`
- VS Code LaTeX Workshop is configured to build the ESI before the main document so external references resolve
- Manuscript figures are read from `Manuscript/figures/`; generated LaTeX tables are read directly from the workspace-level `Outputs/` directory
- The funded SCOPE project proposal is an important source of scientific context and literature:
  - `/Users/sergivela/ownCloud - Sergi Vela Llausi (svlqtc)@nube.cid.csic.es/Archius_Sync/Feina/My_Current_Projects/SCOPE/Memoria/Memoria-RS.pdf`
  - Use its scientific motivation, hypotheses, descriptor framework, and references; ignore implementation, staffing, scheduling, infrastructure, and budget material
  - A scientific-only digest is stored in `Notes/Memoria_RS_Scientific_Context.md`

## Working Assumptions

- Prefer reading pre-gathered system arrays when building notebook analyses
- Keep notebook organization clear, incremental, and numbered
- Treat this file as a lightweight memory aid for future assistance in this workspace
