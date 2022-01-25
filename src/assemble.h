#pragma once
#include "atom.h"
#include "lattice.h"
#include "molecule.h"

/* Append cd2 mol name to lattice name
 * Returns malloced string
 */
char *append_cd2_mol_name(Lattice *lat, Adsorbate *ads, int c1, int c2);

/* Align the adsorbate with the carbon chain */
void ads_align_with_carbon_chain(Lattice *lat, Adsorbate *ads);

/* Align the adsorbate to carbon-metal direction.
 * Allowed carbonId: 41 or 42
 */
void ads_align_with_carbon_metal(Lattice *lat, Adsorbate *ads, int carbonId);

/* Move the adsorbate to the target position (carbon chain ver)
 * The carbonId 1 and 2 must be the adjacent carbon atoms
 */
Lattice *Add_cd2_mol_to_carbon_chain(Lattice *, Adsorbate *, int, int);
/* Move the adsorbate to the carbon metal position */
Lattice *Add_cd2_mol_to_carbon_metal(Lattice *, Adsorbate *, int);
