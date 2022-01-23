#pragma once
#include "atom.h"
#include "lattice.h"
#include "molecule.h"

/* Append cd2 mol name to lattice name
 * Returns malloced string
 */
char *append_cd2_mol_name(Lattice *lat, Adsorbate *ads, int c1, int c2);
// Align the adsorbate with the carbon chain
void ads_align_with_carbon_chain(Lattice *lat, Adsorbate *ads);
// Move the adsorbate to the target position
Lattice *Add_cd2_mol_to_carbon_chain(Lattice *, Adsorbate *, int, int);
