# SCO ATLAS Publication Repository

## Identity and boundaries

- Repository root: `/Users/sergivela/Documents/SCOPE/ATLAS/Analyses`.
- Remote: `https://github.com/QTC-IQAC/SCO-Atlas.git`.
- Purpose: curate reproducible, publication-worthy ATLAS analyses under
  `<analysis>/`.
- Development projects live separately under
  `/Users/sergivela/Documents/SCOPE/Sub-Projects/`, with independent private
  version tracking. Do not turn this public checkout into a mirror of them.
- The parent ATLAS store and its `Data/` directory are not part of this Git
  repository. Read the store's guidance before accessing its data.

Read `/Users/sergivela/Documents/SCOPE/Docs/agent_guidance.md` and applicable
analysis-specific guidance. Default to the local SCOPE development checkout
unless the user specifies a version-pinned task; published layers document the
software versions needed to reproduce their original results.

## Publication workflow

Select the necessary code, notebooks, portable configuration, documentation and
final derived results from the private source project. Record its source
revision and software/data dependencies. Copy only that reviewed selection to
the appropriate analysis directory and inspect the diff and ignore rules.
Validate the exported workflow independently of private machine paths before
publication. Do not copy private `.git` metadata or merge private histories
into this repository. See `docs/publication_workflow.md`.

## Existing first analysis

`1-TMC_Thermodynamics/` contains the first gas-phase study. Its
development source remains in `Sub-Projects/1-ATLAS`, whose existing Git setup
is unchanged pending a separate migration. Its published notebook sequence is
numbered 0 through 6, with focused `D-` diagnostics and existing helpers in
`Notebooks/local_functions.py`. Preserve its established layout. Its original
release targets SCOPE 0.9.6 with both `scope` and `scope_sco`. Do not execute
expensive workflows or regenerate outputs without task scope.

Preserve six-letter refcodes and physical paths with derived suffixes. Keep
gas-phase molecular crossing temperatures, periodic results and experimental
observables distinct. Never overwrite raw ATLAS sources, computations, Systems
or published data releases while preparing an analysis publication.
