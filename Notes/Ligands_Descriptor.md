### Distance-Weighted Core--Substituent Ligand Descriptor

Ligand similarity was encoded separately for each denticity using coordination-core and substituent fingerprints. For homoleptic complexes, one representative ligand was used.

To emphasize the ligand region closest to the metal center, Morgan-type atom environments were weighted according to the topological distance between each ligand atom and the metal atom in the molecular graph. For a ligand atom at graph distance \(d_i\) from the metal, an exponential weight was assigned as:

\[
w_i = \exp[-\alpha(d_i - 1)]
\]

where atoms directly bound to the metal have \(w_i = 1\), and more distant atoms contribute progressively less.

For each ligand atom, Morgan environments up to a fixed radius were generated using the atom as the environment center. Instead of constructing a binary fingerprint, each activated fingerprint bit was incremented by the weight of the corresponding center atom, yielding a distance-weighted count fingerprint.

The ligand was divided into a coordination core and its substituents. For ring-based ligands, the core contains the metal-connected ring systems and its fingerprint is computed exclusively on the isolated core fragment. For ligands without a metal-connected ring, an acyclic core is defined from the shortest graph paths connecting the donor atoms; its fingerprint is evaluated in the complete ligand so radius-2 environments can retain information about the immediate substituent environment. The remaining atoms define the substituent fingerprint.

Families are assigned from a ring-composition label or, for acyclic cores, from the donor-backbone topology. Formal charges are ignored. Bond orders within an acyclic backbone are also ignored except for bonds immediately adjacent to the metal-bound donor atoms, thereby distinguishing donor environments such as N--C and N=C without allowing remote unsaturation to split a family. Each family is represented by the mean core fingerprint of its members. Cosine distances between these vectors are projected into two dimensions by metric multidimensional scaling, after which family centres are arranged on a sparse grid. The substituent fingerprints determine the local arrangement of compounds within each family.
