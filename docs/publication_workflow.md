# Curated analysis publication

The repository at `/Users/sergivela/Documents/SCOPE/ATLAS/Analyses` tracks
`https://github.com/QTC-IQAC/SCO-Atlas.git`. It is the central local checkout for
the shared publication repository. The parent ATLAS directory is a data store,
not a Git repository.

Private development projects live under `SCOPE/Sub-Projects/`. Each is intended
to have its own history and private remote. The first study's existing Git setup
is migrated separately; its local manuscript, notes, figures and other
development files must be retained. The solid-state project has a fresh
independent history.

Publishing an analysis means selecting reproducible code/notebooks, necessary
portable settings, documentation and final results from its private source.
Copy the selection into `<analysis>/`; include the source revision and
software/data provenance. Adapt local paths and helper dependencies for the
public environment. Review the diff and validate the exported material before
committing or pushing it here. Export updates in the same way.

Keep private history and draft material in the source project. Do not mirror
entire development folders, copy their `.git` directories, delete destination
files automatically, or merge unrelated private histories. Publication copies
are intentional curated snapshots; editing private work does not update them.
